# Changelog

All notable changes to ElectromagneticFields.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Releases before v0.6.0 are
not covered here; see the git history for those.


## [0.7.1] - 2026-08-10

### Added

- `code` emits a nullary `orientation()` into the generated module, so a loaded equilibrium can
  recover the handedness of its own chart and reconstruct `det DF` as `orientation() * J`. It is the
  one generated function that takes no arguments — the sign depends on neither `t` nor `ξ` — and the
  value is interpolated as a literal `-1` or `1` at generation time. It carries its own `export`, so
  both `Mod.@code` and `load_equilibrium` work. `ElectromagneticFields.orientation` itself stays
  unexported, as `J` is, to avoid a collision in `Main`.

- **A do-block form of `load_equilibrium`**, for callers that are not at top level. The methods
  `load_equilibrium` evaluates are defined in a world age newer than the frame that called it, so
  that frame cannot call them — from inside a function they are "too new", which Julia reports as a
  `MethodError` on a name that plainly exists. At top level the world age advances between
  statements, so the plain form is fine there; everywhere else,

  ```julia
  load_equilibrium(equ; target_module = MyModule) do mod
      mod.orientation()
      mod.DF(t, ξ)
  end
  ```

  passes the target module to `f` through `Base.invokelatest`, so the body runs in the world age the
  generated methods live in, and returns `f`'s value. Both the name lookup and the call have to
  happen inside `f`. The `@code` macros are unaffected either way — they splice the definitions in
  at expansion time.

- Test coverage for `load_equilibrium`, which previously had none: one left-handed and one
  right-handed equilibrium are loaded into throwaway modules and checked for the sign, agreement
  with the trait, `det(DF) ≈ orientation() * J`, and that `orientation` is exported. The calls sit
  inside the `@testset` in do-block form, so the world-age case is what the test exercises.

### Changed

- `load_equilibrium` returns the target module, rather than whatever `Core.eval` gave back for the
  last generated statement (the `orientation` method object, which nothing used).
- `code` no longer assigns its parameter list to a module-level global on every call. Nothing read
  it, and on Julia 1.12 creating and reassigning a global binding is a world-age event.
- `code(...; output = 1)` now reports `orientation` among the generated functions.
- The `orientation`, `code` and `load_equilibrium` docstrings describe the generated `orientation()`
  and the world-age constraint.


## [0.7.0] - 2026-08-07

**This release is breaking.** The magnetic field changes sign in every left-handed chart. Downstream
results computed in `AxisymmetricTokamakCylindrical`, `AxisymmetricTokamakToroidal`,
`AxisymmetricTokamakToroidalRegularization` and every `Solovev*` equilibrium except
`SolovevSymmetric` will change accordingly; the cartesian charts are unaffected.

### Fixed

- **`B¹ = ⋆dA` used the unsigned volume element.** `hodge²¹` was handed `J = √|g| = |det DF|`, but
  the Hodge star is orientation-dependent, and four of the charts here are left-handed: `(R, Z, ϕ)`
  and `(r, θ, ϕ)` both have `det DF < 0`, the right-handed orderings being `(R, ϕ, Z)` and
  `(r, ϕ, θ)`. In every one of them the magnetic field came out reversed — the same vector potential
  produced a `B` antiparallel to the one the cartesian chart gives at the same physical point. The
  Hodge star and the cross product, the only two orientation-dependent operations, now use the
  signed determinant `orientation(equ) * J(x, equ)`.
- **`AxisymmetricTokamakToroidalRegularization` had the wrong sign on `A₂`**, making its toroidal
  field antiparallel to the other three axisymmetric tokamak charts. The chart was internally
  consistent — it satisfied `B = ∇×A` for its own `A` — so it was wrong rather than inconsistent,
  which is why the orientation work did not surface it. All four axisymmetric tokamak charts now
  agree to twelve digits.
- Docstring formulas for `A` and `B` in the cylindrical and both toroidal charts, which disagreed
  with the code they document.

### Added

- `orientation(::AnalyticField)` trait — `+1` for a right-handed chart, `-1` for a left-handed one —
  declared `-1` for the four left-handed charts, with a table of every chart's handedness in its
  docstring. `J` keeps its meaning and its value as the unsigned volume element.
- `test_curl` checks `B = ∇×A` in cartesian coordinates for all twenty equilibria, pulling each
  chart's covariant `A` onto the cartesian frame and central-differencing. This is the only
  statement about `B` that does not depend on the chart, so it is the one that catches an
  orientation error.
- `det(DF) ≈ orientation(equ) * J` alongside the existing `J ≈ sqrt(det(DF'DF))`, so a new chart
  that declares the wrong handedness fails rather than silently reversing its `B`.

### Changed

- The cartesian/cylindrical and cartesian/toroidal consistency checks lose their explicit minus
  sign: they had asserted `DF * B_cyl ≈ -B_car`, encoding the reversed field rather than flagging
  it. Per-equilibrium assertions go from 291 to 294.


## [0.6.3] - 2026-07-29

### Added

- **Common subexpression elimination in the generated field code.** `convert(Expr, ::Basic)` writes
  out SymEngine's expanded tree verbatim and SymEngine shares nothing, so a subexpression common to
  fifty branches of `∂ⱼ(Bᵢ / sqrt(gᵏˡBₖBₗ))` was emitted fifty times. A hash-consing pass now runs
  between `convert(Expr, f_expr)` and the `quote` wrapping each body, evaluating each distinct
  subexpression once into a local. Measured on `Solovev.ITER(xpoint = true)`:

  | function | `log` calls | speedup |
  |---|---|---|
  | `B` | 15 → 1 | 2.4× |
  | `b₁` | 21 → 1 | 3.3× |
  | `dBdx₁` | 60 → 1 | 6.3× |
  | `db₁dx₁` | 108 → 1 | 6.8× |
  | `d²b₁dx₁dx₁` | 367 → 1 | 12.8× |

  First-call latency falls with it, because there is less code to compile: `d²b₁dx₁dx₁` compiles in
  12 ms rather than 118 ms, and `load_equilibrium(Solovev.ITER(xpoint = true))` goes from 1.07 s to
  0.80 s end to end. The emitted code also gets smaller, not larger — sharing more than pays for the
  assignments.

  The pass only ever *names* subexpressions. It never reassociates, never reorders an operand, never
  folds a constant, so generated functions return bit-for-bit what they returned before. The test
  suite asserts this structurally for every generated function of every equilibrium it covers —
  8824 assertions.

- `cse = false` disables the pass, which is useful when reading generated code against a paper. It
  is accepted by `code`, `load_equilibrium` and the `@code` macros.

- **The `@code` macros forward options to `code`.** Macros cannot take keyword arguments, so options
  are written as `key = value` arguments the way `@testset` does. `output` comes along for free:

  ```julia
  ThetaPinch.@code cse = false            # defaults for the equilibrium's parameters
  ThetaPinch.@code(B₀, cse = false)       # or alongside them
  Solovev.@code_iter true                 # positional parameters still work
  code(equ; cse = false)                  # and directly, as before
  load_equilibrium(equ; cse = false)
  ```

- `SolovevFRC` is covered by the test suite for the first time, bringing the analytic testsets to
  twenty equilibria. It passes the full `@test_equilibrium` block at 291/291 like every other
  equilibrium.

### Changed

- **The default `PenningTrapBottle` and `PenningTrapAsymmetric` equilibria have changed.** Both
  modules' `init` functions take `(B₀, Bₚ, E₀)` but constructed the equilibrium as `(B₀, E₀, Bₚ)`,
  so the perturbation and electric field strengths were exchanged on the way in — for explicit calls
  and for the defaults alike:

  | | before | after |
  |---|---|---|
  | `PenningTrapBottle` | `Bₚ = 10`, `E₀ = 200` | `Bₚ = 200`, `E₀ = 10` |
  | `PenningTrapAsymmetric` | `Bₚ = 10`, `E₀ = 50` | `Bₚ = 50`, `E₀ = 10` |

  The constants were right and `init` was wrong: `DEFAULT_E₀ = 10.0` in all three Penning trap
  modules, and the swap was the only thing making the bottle and asymmetric traps disagree with
  `PenningTrapUniform` on the field strength. If you were relying on the 0.6.2 defaults, pass the
  parameters explicitly.

- Positional parameter defaults moved from the `@code` macro signatures to `init`, which already
  carried the same ones.

### Fixed

- `AxisymmetricTokamakToroidal.ITER()` and `@code_iter` raised `UndefVarError`. The module
  referenced `ITER_R₀`, `ITER_B₀` and `ITER_q₀`, which only the cylindrical and cartesian modules
  ever defined. Added, matching the cylindrical module's values.
- `Solovev.FRC()` and `@code_frc` never returned. The equilibrium had `R₀ = 0`, and `R₀` enters the
  metric as `g₁₁ = g₂₂ = R₀²`, so `inv(gmat)` was singular. `R₀ = 1` is the normalising major radius
  of Cerfon & Freidberg (*Phys. Plasmas* **17**, 032502, 2010), the paper the module already cites.
  `B₀ = 0` is left alone — a zero toroidal field is what makes it a field-reversed configuration,
  and `A₃ = ψ` does not involve `B₀`, so the poloidal field is untouched.
- `Solovev.@code_xpoint` no longer calls a nine-argument `SolovevXpointEquilibrium` that has no such
  method, and it now honours `doublex`.


## [0.6.2] - 2025-12-10

### Changed

- Equilibria that use logarithms now go through `NaNMath`: the axisymmetric tokamak cylindrical,
  toroidal and regularized charts, and the Solov'ev equilibria. Adds a dependency on `NaNMath`.
- Julia compat entry raised to 1.10, and the CI workflow updated to match.


## [0.6.1] - 2025-12-04

### Added

- Derivatives of the electric field, generated as `dEⁱdxⱼ`.
- Second derivatives of the metric coefficients, generated as `d²gᵢⱼdxₖdxₗ` and `d²gⁱʲdxₖdxₗ`.


## [0.6.0] - 2025-10-21

**This release is breaking.** `periodicity` is gone from the generated modules.

### Changed

- `periodicity` is replaced by `rangemin` and `rangemax`. Charts declare their bounds through
  `minx¹`/`minx²`/`minx³` and `maxx¹`/`maxx²`/`maxx³`, which default to `-Inf` and `+Inf`, and the
  generated module exports `rangemin` and `rangemax` in place of `periodicity`.


[0.7.1]: https://github.com/JuliaPlasma/ElectromagneticFields.jl/compare/v0.7.0...main
[0.7.0]: https://github.com/JuliaPlasma/ElectromagneticFields.jl/compare/v0.6.3...v0.7.0
[0.6.3]: https://github.com/JuliaPlasma/ElectromagneticFields.jl/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/JuliaPlasma/ElectromagneticFields.jl/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/JuliaPlasma/ElectromagneticFields.jl/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/JuliaPlasma/ElectromagneticFields.jl/compare/v0.5.3...v0.6.0

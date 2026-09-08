# Changelog

All notable changes to ElectromagneticFields.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Releases before v0.6.0 are
not covered here; see the git history for those.


## [Unreleased]

### Changed

- `src/analytic/analytic_equilibrium.jl` and `test/test_analytic.jl` are now Unicode
  NFC-normalised. They stored `ḡ` (21 times), `â` (3) and `ĉ` (3) as a base letter plus a combining
  mark, inherited from macOS rather than chosen, which makes no difference to the compiled code —
  Julia's parser normalises identifiers to NFC — but defeats every byte-matching tool: a `grep`
  pattern or an editor search typed in NFC matches nothing in such a file, silently.

  Worth naming, because it is the hazard this removes: the file already spelled the same glyph
  both ways. `functions["ḡ"]` was precomposed while the `ḡ` identifiers around it were decomposed,
  and that worked only because the parser normalises the identifier side and not the string. No
  non-NFC string literal remains under `src/`. `DF̄` is safe either way — a macron over `F` has no
  precomposed codepoint.

  A separate note for a reader of this file: `g̅` elsewhere in the same source is a *different*
  binding, `g` plus a combining overline rather than a macron. It is unchanged, and
  indistinguishable from `ḡ` on screen.

  The only bytes that change at runtime are the `"Dḡ"` and `"DDḡ"` labels passed to `symprint`,
  which prints them; string literals are not parser-normalised. Nothing reads that output back.

## [0.8.0] - 2026-08-10

### Changed

- **Plotting moved from Plots.jl to Makie**, and out of the package proper into a package
  extension. The twelve `RecipesBase.@recipe` definitions scattered through `src/analytic/` are
  replaced by `ext/ElectromagneticFieldsMakieExt.jl`, which is loaded as soon as `Makie` (or one
  of its backends) is. `RecipesBase` and `LaTeXStrings` are no longer dependencies of
  ElectromagneticFields.

  This is a breaking change: `plot(equ)` with Plots.jl no longer works. The replacement is

  ```julia
  using CairoMakie
  using ElectromagneticFields

  plot_equilibrium(Solovev.ITER())
  ```

  `plot_equilibrium(equ; size, kwargs...)` creates a `Figure` and returns it, while
  `plot_equilibrium!(position, equ; kwargs...)` draws into an existing one at any Makie grid
  position, e.g. `fig[1,2]`, and returns the `Axis` it created, or the vector of axes for the
  fields that draw more than one panel. The latter is what replaces Plots' `layout` for composing
  several equilibria into one figure, and it works uniformly for the single-panel fields and for
  those that draw a panel per vector potential component. The keyword arguments carry over
  unchanged (`nx`, `ny`, `levels`, `xlims`, `ylims`, plus `nτ` for Solov'ev and `ni` for ABC),
  except that ABC's `nl` is now spelled `levels` like everywhere else, and `aspect` replaces
  Plots' `aspect_ratio`. Contour panels take an opt-in `colorbar` keyword, and the Solov'ev
  equilibrium a `boundary` keyword to switch off the plasma boundary drawn on top of the flux
  surfaces.

  Equilibria without a plotting method — the three Penning traps — now report that in an
  `ArgumentError` instead of a `MethodError` on an internal helper.

- **`Documenter` is no longer a dependency.** `src/ElectromagneticFields.jl` carried a `using
  Documenter` that nothing used — `@doc raw` is Base, and `@ref` in a docstring is plain text until
  Documenter parses it at build time — so every downstream install pulled in Documenter and its
  tree. The inert `[targets] docs` entry goes with it; the doc environment declares Documenter
  itself.

- **Every generated function now returns the type of the coordinates it was given.** A body that
  does not mention the coordinates used to evaluate to whatever literal type SymEngine emitted, so
  the structurally constant components came out as `Int64` — `g₁₁`, the off-diagonal entries of `g`
  and `DF`, and `φ` and `E` of a purely magnetic equilibrium — while the rest were `Float64`. That
  made the matrix wrappers promote at runtime: `g` allocated 976 bytes and `DF` 1072 per call, where
  the 3×3 matrix they return is 144. Both are now 144, the scalar functions still allocate nothing,
  and `test_analytic.jl` asserts as much for every equilibrium. Values are unchanged; the one
  behavioural difference is that a coordinate-independent constant such as `B₃ = B₀ R₀` now takes the
  coordinate's type, so `Float32` coordinates give a `Float32` result where they previously gave
  `Float64` — the convention the `one(T)` / `zero(T)` chart traits already followed.

- The identity blocks in the documentation are wrapped in a hidden `@assert`, so an identity that
  stops holding fails the doc build instead of quietly rendering `false`.

- The `Documentation` workflow drops four redundant steps — a `Pkg.develop` that duplicates the
  `[sources]` entry in `docs/Project.toml`, a `Pkg.build`/`Pkg.precompile` pair covered by
  `Pkg.instantiate`, a `julia-buildpkg` for the main project the doc build never uses, and a
  standalone doctest run that repeats what `makedocs` already does — pins Julia to 1.12, which
  `[sources]` requires, and caches its depot.

- The stale Travis and Coveralls badges are gone, and the CI and Codecov badges in the documentation
  now match the working ones in the README, alongside a new Documentation badge in both.

### Fixed

- **Six of the contour plots were transposed.** Plots.jl expects the value matrix indexed as
  `z[j,i]` for `(x[i], y[j])`, and only the tokamak and Solov'ev recipes transposed accordingly;
  the recipes for ABC, Dipole, QuadraticPotentials, Singular, SymmetricQuadratic and ThetaPinch
  passed the matrix through as built and so plotted the mirror image about the diagonal. Makie
  uses the `z[i,j]` convention that the comprehensions already produce, and the ported code
  passes them straight through, so all twelve plots now show the field in the correct orientation.

  The six affected recipes are exactly the six with square default grids, `nx == ny == 100`,
  which is why Plots.jl never raised a dimension error; the two families that sample
  `nx = 100, ny = 120` are precisely the two that transposed correctly. The test suite now plots
  every field on a non-square grid as well, so a reintroduced transpose fails loudly.

- The quantity plotted as `|B|` for the ABC field was `|B|²`, and the one plotted as `B_z` for the
  symmetric quadratic field was `B₀ / (1 + x² + y²)` where the field is `B₀ (1 + x² + y²)`. Both
  helpers are used only for plotting — the generated evaluation routines were never affected.

- The logarithmic contour levels of the singular field are now derived from the finite magnitudes
  of the data. Previously they came from `maximum` and `minimum` directly, which failed on any
  grid that includes the singular line (`Inf` and `NaN` values, e.g. for odd `nx` and `ny`) and on
  any plot range in which a component of the vector potential does not change sign, such as a
  window that does not straddle the axis.

- The ABC field no longer evaluates `|B|` on the full three-dimensional grid to draw three
  mid-planes, which cost `O(nx³)` time and memory for `O(nx²)` values.

- **The tokamak and Solov'ev plots showed the wrong quantity.** They contoured `A₃ / R`, the
  *physical* toroidal component of the vector potential, and described it as the poloidal flux
  function. The flux function is the covariant component `A₃ = ψ` itself: contracting the magnetic
  field with `∇A₃` gives zero, while `B · ∇(A₃/R) ≠ 0`, so the contours drawn were not flux
  surfaces. For the cartesian tokamak, whose `A_y` at `y = 0` is that same physical component, the
  plot now shows `R A_y`. Only the `ψ = 0` contour was unaffected, which is why the red plasma
  boundary of the Solov'ev equilibria always looked right. Five figures in the documentation change.

- `plot_equilibrium` gives `size` and `figure` a defined precedence — the size chosen for the
  equilibrium, then `figure`, then an explicit `size` — instead of letting a size inside `figure`
  override the `size` argument as a side effect of the splat order. Documented and tested.

- The contour levels of the Solov'ev equilibria are anchored to the flux on the magnetic axis rather
  than spread evenly over the sampled range. `ψ` vanishes on the plasma boundary and grows without
  bound away from it, so an even spread spent nearly all its levels on the far field and left the
  flux surfaces of the plasma to four or five of them. They are now uniform in `ψ`, with one level
  exactly on the boundary and a quarter of them inside it — ten at the default `levels = 40`, down
  from `50`. Passing the levels themselves still bypasses this, as everywhere else.

- The default plot range of `SolovevSymmetric` is centred on its magnetic axis. The flux function
  depends on `x` through `(R₀ + x)⁴`, so the axis sits at `x = -R₀`, while the window was centred on
  `+R₀` — for any `R₀ ≠ 0` it showed a monotone ramp rather than closed flux surfaces. Its docstring
  described `R₀` as the position of the magnetic axis, which is off by a sign.

- The `φ` in the docstrings of the three Penning traps is missing its minus sign: the potential is
  `-E₀ (x²/2 + y²/2 - z²)`, which is what makes the stated `E = E₀ (x, y, -2z)` follow from
  `E = -∇φ` and what the code has always computed.

- The ABC page documented the wrong integrability criterion. The field lines are integrable when one
  of `a`, `b`, `c` vanishes; equal coefficients are not a special case, and `a = b = c = 1` — the
  package default — is the classic chaotic one.

- `singular.md` claimed the vector potential diverges like `r⁻³`. It diverges like `r⁻²`; only `|B|`
  goes as `r⁻³`.

### Removed

- The `examples/` directory. Its Jupyter notebooks predated the current API, carried no narrative,
  and one of them had been sitting in the repository with unresolved merge-conflict markers since
  2020. Everything worth keeping now lives in the documentation, one page per field, as executed
  `@example` blocks.


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

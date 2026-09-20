# Changelog

All notable changes to ElectromagneticFields.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Releases before v0.6.0 are
not covered here; see the git history for those.


## [Unreleased] — targeting 0.9.0

### Added

- **Aqua.jl runs as part of the test suite**, in `test/aqua_tests.jl`, which calls
  `Aqua.test_all`. This package carried no package-level quality check of any kind before, so all
  eight — method ambiguities, unbound type parameters, undefined exports, the agreement between
  `Project.toml` and `test/Project.toml`, stale dependencies, `[compat]` bounds, type piracy and
  persistent tasks — start from a clean pass rather than from a list of known failures.

  Two of them cover something this package does structurally. `test_piracies` covers the eight
  methods added to GeometricBase's `periodic`, `functions` and `parameters`, each of which
  dispatches on a type defined here; a new chart family brings a `periodic` method with it, so the
  property outlasts any list of the sites. `test_undefined_exports` covers an export list dominated
  by the musical-isomorphism accessors, which are defined across `src/analytic/` rather than beside
  the export.

- **The test suite has a `test/Project.toml` of its own**, carrying Aqua, CairoMakie,
  LinearAlgebra, SafeTestsets, StaticArrays and Test, each with a `[compat]` bound. The `[extras]`
  and `[targets]` sections of the package's `Project.toml` are gone with it, together with the
  `CairoMakie` bound that belonged to them. Test dependencies that only `[extras]` declared carried
  no bounds at all, which is what `Aqua.test_deps_compat` reports.

### Changed

- **The symbolic engine is replaced from SymEngine.jl to Symbolics.jl.** Dependencies change:
  SymEngine is removed, Symbolics, StaticArrays, GeometricBase and ConstructionBase are added.
  Julia floor stays 1.10.
  GeometricBase is there for the generic names `functions`, `parameters` and `periodic`, which
  this package now extends rather than defining its own — it previously exported a `periodicity`
  of its own, which collided with GeometricEquations' in any package using both.

  This change is **breaking** for any code that relied on the old API. Fields no longer inject
  generated code into a module. Instead, `FieldFunctions(equ, pert = ZeroPerturbation(); cse = true)`
  builds a struct holding the generated functions:

  ```julia
  using ElectromagneticFields

  field = FieldFunctions(SolovevEquilibriumITER())

  # every quantity is reached through a generic taking the field first
  B♭(field, 1.0, [1.05, 0.25, 0.5])

  # or all of them at once, as a NamedTuple
  functions(field)
  ```

  The generated functions are now RuntimeGeneratedFunctions built at runtime using
  `Symbolics.build_function(...; expression = Val{false})`, with optional common-subexpression
  elimination. A field is now a value: it can be constructed inside a function, stored, and passed
  to a vector field as part of the `parameters` NamedTuple of a GeometricEquations problem,
  without the world-age trap the old module-injection approach had.

  Building one is not free, and almost none of the cost is what one would guess. On the toroidal
  tokamak the symbolic trace takes 7.1 s and code generation 4.1 s, while compiling the generated
  code takes 0.6 s — 5% of the total. The rest is Julia compiling Symbolics' own machinery for the
  expression shapes a trace produces, which is paid per *shape* rather than per field.

  Being per shape, it is recoverable, and two things now recover it. The generated functions of a
  type already built are **cached** — keyed on the equilibrium and perturbation types rather than
  on the parameter values, which the code no longer contains — and a `PrecompileTools` workload
  **traces one field of every shipped type** during precompilation, so the specializations and
  the cache both land in the package image.

  The result is that constructing any equilibrium this package ships costs nothing: the first
  field in a session drops from 8.5 s to 0.00 s, and all twenty together from about 19 s to
  0.01 s, for any parameter values. The price is this package's own precompilation, 1.5 s →
  24.5 s, paid once per version. An equilibrium of your own is traced once per session and cached
  after that.

  The cache is keyed on types, not on the content of the methods behind them, so redefining an
  `A₁` or a metric coefficient in a running session leaves it stale. `clear_field_cache!()`
  discards it, and `FieldFunctions(equ; cache = false)` bypasses it for one call.

- **The generated code no longer has the equilibrium's parameters baked into it.** The symbolic
  trace runs against a copy of the equilibrium whose parameters are symbolic, so the code is
  written in terms of `R₀`, `B₀`, `q₀` … and takes their values as an argument. The equilibrium
  struct stays where the values live and the field passes them in on each call.

  What this changes for a caller is mostly that things are faster, but the property is worth
  knowing: the generated code depends on the equilibrium's *type* and not on its parameters, so
  two fields of the same type share their compiled functions exactly. That is what makes the
  precompiled workload above cover parameter values nobody anticipated, rather than only the ones
  it happened to name.

  It also fixes the meaning of `get_parameters`, which used to select which parameters to export
  as constants and now selects which are arguments of the generated code. Anything the `A₁`, `φ`
  or metric methods read and it does not list is frozen into the code as a literal. The default —
  every field of the struct but `name` — is what all the equilibria here want, so the Solov'ev
  override that excluded the derived coefficient vector `c` is gone: `A₃` reads `c`, so `c` has to
  travel with the rest. `parameters(field)` for a Solov'ev equilibrium therefore now includes it.

  A field type of your own must be constructible as `Type{T}(parameters...)` — the convention all
  of these follow, a parametric struct whose first member is `name` — so that the trace can build
  the symbolic copy. `symbolic_copy` takes a method for a type that departs from it.

- **A package can build its fields during its own precompilation** and pay nothing at load time.
  `@precompilable_fields` at the top level of a module, and `cache_module = @__MODULE__` passed to
  `FieldFunctions`, put the generated bodies in that module's cache so they survive into its
  package image; a `PrecompileTools.@compile_workload` over the accessors caches their compiled
  code as well. A field prepared this way evaluates in microseconds on the first call of a fresh
  session, with no trace, no code generation and no compilation. Without it the field still works
  and is simply rebuilt on each load.

- **The function API replaces ~450 scalar functions with ~41 generics returning tensors.** Each
  generic now takes the field as its first argument, followed by `t` (time) and `ξ` (coordinates).
  All old scalar values remain accessible as entries in the returned tensors.

  The new names use the musical isomorphisms: `♭` (`\flat`) lowers an index and denotes the
  covariant components, `♯` (`\sharp`) raises one and denotes the contravariant components, and
  `♮` (`\natural`) denotes the physical components. The bare letter is the magnitude, which has no
  index to raise or lower, so `B` is still `|B|`; two lowered indices, `B♭♭`, are the magnetic
  two-form. A leading `D` is a derivative with respect to the chart coordinates:

  | New name | Old names | Type |
  |---|---|---|
  | `B♭(field, t, ξ)` | `B₁, B₂, B₃` | SVector covariant components |
  | `B♯(field, t, ξ)` | `B¹, B², B³` | SVector contravariant components |
  | `B(field, t, ξ)` | `B` | magnitude |
  | `B♭♭(field, t, ξ)` | `B₁₁, B₁₂, ..., B₃₃` | SMatrix two-form |
  | `g♭(field, t, ξ)` | `g` | SMatrix covariant metric |
  | `g♯(field, t, ξ)` | `ḡ, g⁻¹` | SMatrix contravariant metric |
  | `DA♭(field, t, ξ)` | `dAᵢdxⱼ` | SMatrix first derivatives |
  | `DB(field, t, ξ)` | `dBdxᵢ` | SVector first derivatives of magnitude |
  | `DDB(field, t, ξ)` | `d²Bdxᵢdxⱼ` | SMatrix Hessian of the magnitude |
  | `DDA♭(field, t, ξ)` | `d²Aᵢdxⱼdxₖ` | SArray second derivatives |
  | `Dg♭`, `Dg♯`, `DDg♭`, `DDg♯` | `dgᵢⱼdxₖ`, `d²gᵢⱼdxₖdxₗ`, … | SArray metric derivatives |
  | `b♭`, `b♯`, `b♮` | `bᵢ`, `bⁱ`, `b₍ᵢ₎` | SVector unit magnetic field |
  | `a♭`, `a♯`, `a♮`, `c♭`, `c♯`, `c♮` | `a`, `a⃗`, `aₚ`, … | SVector perpendicular frame |
  | `E♭(field, t, ξ)`, `E♯` | `Eᵢ`, `Eⁱ` | SVector electric field |
  | `J(field, t, ξ)` | `J` | volume element `√\|g\| = \|det DF\|` |
  | `φ(field, t, ξ)` | `φ` | electrostatic potential |
  | `DF(field, t, ξ)`, `DF̄` | `DF`, `DF̄` | SMatrix tangent map and its inverse |
  | `to_cartesian`, `from_cartesian`, `rangemin`, `rangemax` | same | SVector, field first |

  Four things become data rather than functions: `parameters(field)` — which replaces the
  constants the old code spliced into the module — `coordinates(field)`, `periodic(field)`, and
  `orientation(field)`, which was a generated zero-argument function and is now a stored `Int`.
  `equilibrium(field)` and `perturbation(field)` return the objects the field was built from, and
  `functions(field)` returns all generated functions as a NamedTuple. Every one of the ~444 old
  scalar values is an entry of one of these tensors, verified against the SymEngine-0.8 output
  across all 20 equilibria and all 8844 scalar points — see `scripts/verify_against_symengine.jl`.

  A tensor computation shares subexpressions across all its entries and is therefore usually faster
  than the old per-component calls, and every accessor is allocation-free and type-stable. The test
  suite asserts both properties for all 20 equilibria.

- **Field constructors replace the old submodule `init()` pattern.** `Solovev.ITER()` is now
  `SolovevEquilibriumITER()`, exported from the package. The six Solov'ev alias modules
  (SolovevITER, SolovevNSTX, SolovevFRC, SolovevITERwXpoint, SolovevNSTXwXpoint,
  SolovevNSTXwDoubleXpoint) are superseded by `SolovevEquilibriumITER`, `SolovevEquilibriumNSTX`,
  `SolovevEquilibriumFRC`, `SolovevXpointEquilibriumITER`, `SolovevXpointEquilibriumNSTX` and
  `SolovevDoubleXpointEquilibriumNSTX` — names that already existed inside the Solovev module and
  are now exported from the package. `Solovev.ITER(xpoint = true)` and its siblings are gone. All other field modules follow
  the same pattern: `ThetaPinch.init()` → `ThetaPinchEquilibrium()`, `ABC.init()` → `ABCEquilibrium()`,
  `Dipole.init()` → `DipoleField()`, and likewise for Singular, QuadraticPotentials,
  SymmetricQuadratic, SolovevSymmetric, PenningTrap (with its three variants Uniform, Bottle and
  Asymmetric) and AxisymmetricTokamak (Cartesian, Cylindrical, Toroidal and the Toroidal
  Regularization gauge, each with an `AxisymmetricTokamak*ITER()` preset).
  The coordinate helpers (X, Y, Z, R, r, θ, ϕ, r²) are now package-level generics with one method
  per equilibrium type, not exported (to avoid collisions), and are accessible as
  `ElectromagneticFields.R` or via the `coordinates(field)` NamedTuple.

- **The type constructor is reached through `ConstructionBase.constructorof`** rather than
  `Base.typename(T).wrapper`, at the two places that strip a struct's type parameter: the symbolic
  copy the trace builds, and the generated-function cache key. Neither `Base.typename` nor the
  `.wrapper` field of `Core.TypeName` is public. ConstructionBase exists for this and defines the
  same thing, so the behaviour is unchanged for this package's types.

- Documentation restructured with new "Interface" and "Code Generation" pages; Usage, Coordinates,
  Fields, Plotting and the twelve analytic field pages rewritten for the new API.

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

- **The Solov'ev coefficient solve reaches Symbolics through a public API.** Imposing the boundary
  conditions means evaluating a symbolic expression at a point, and `substitute` folds arithmetic
  but not `log`, so a fully substituted `ψ₃` still carries an unevaluated `NaNMath.log(1.32)` and
  cannot be converted to a number. That step used `Symbolics.symbolic_to_float`, which is neither
  exported nor `Base.ispublic`, so any 7.x release may remove it. The expression is now built into
  a function with `Symbolics.build_function(…; expression = Val{false})` and called — the same
  machinery this package emits its field code with, applied one step earlier.

  The coefficients move, by less than the solve determines them. Four of the six shipped Solov'ev
  equilibria come out bit for bit identical; the largest change is `SolovevXpointEquilibriumITER`,
  at `6.1E-14` relative. The control for that figure is the same solve with every matrix entry
  perturbed by half an ULP, which is what two correctly rounded evaluators may disagree by: over
  20 draws it moves the same vector by a median `7.3E-14` and by up to `2.5E-13`. The SymEngine
  reference values are unchanged, 8844 of 8844 with no frame differences, and the test suite now
  asserts the boundary conditions the solve imposes, which nothing did before — every identity it
  checked holds for any coefficients whatever.

  It costs cold session time. Each distinct expression becomes a `RuntimeGeneratedFunction`
  that Julia compiles on its first call, and the six shipped Solov'ev equilibria make 702
  evaluations over 59 distinct expressions, so building all six in a fresh process roughly
  doubles, from about 0.37 s to about 0.7 s — median of five cold processes each, ratio
  between 1.8× and 2.0× across three runs, Julia 1.13 on aarch64 macOS. Compilation is
  the whole of it: building one function per distinct expression rather than one per call
  recovers none of the difference, and the second pass in one process is 0.018 s either way.

- **`periodicity` is replaced by `periodic`, which says which coordinates are periodic, per
  chart.** `periodicity` returned `(-Inf * ones(4), +Inf * ones(4))` for every equilibrium this
  package ships: four components for a three-coordinate chart, and no periodicity in any of them —
  including the two coordinates of the toroidal chart that are angles. That was the only
  definition of it anywhere in the package, and no equilibrium had ever overridden it, so the
  answer said nothing about any field. This package no longer answers `periodicity` at all, and no
  longer exports it.

  `periodic(equ)` and `periodic(field)` return an `SVector{3, Bool}`, one entry per coordinate of
  the chart: all `false` for the cartesian charts, `false, false, true` for the cylindrical ones,
  where the toroidal angle wraps, and `false, true, true` for the toroidal ones, where both angles
  do. It takes the equilibrium alone rather than a point and the equilibrium, because the answer
  does not depend on where in the chart it is asked.

  The name is `GeometricBase.periodic`, where `periodic(s::StateVariable)` already means one
  `Bool` per component and `isperiodic` means `any` of them. `GeometricBase.periodicity` is a
  different thing in this ecosystem: `GeometricEquations` answers it with an `(xmin, xmax)` tuple
  naming the periodic domain, and derives the per-component answer from those bounds under the
  name `getperiodicity`. Keeping the `Bool` vector on `periodicity` would have given one shared
  generic two incompatible shapes, and the mismatch is silent rather than loud —
  `per_lo, per_hi = periodicity(equ)` destructures a three-element `Bool` vector to
  `(false, false)` without error. A `GeometricEquations` problem still takes its own `periodicity`
  in its own form; `periodic(field)` is what says which components belong in it.

  It is a property of the chart, not of the equilibrium, so it is defined per chart family: five
  methods cover all twenty shipped equilibria. **There is no fallback.** A chart of your own that
  has not defined `periodic` raises a `MethodError` when a field is built from it. An
  all-`false` default was rejected because it is indistinguishable, at every call site, from a
  chart that genuinely has no periodic coordinate — which is exactly the failure being fixed.

  It cannot be read off `rangemin`/`rangemax`, which answer a different question: a bounded range
  does not imply that a coordinate wraps. The two agree on every chart here, since each bounds
  exactly its own angles, but a wall at `r = a` or a slab bounded in `z` would have a finite range
  in a coordinate that does not wrap.

  Two oddities go with the old shape. Four components for three coordinates fitted a `(q, t)` or
  phase-space convention that nothing here uses. And `FieldFunctions` obtained the value by
  calling the default with a hard-coded `zeros(3)`, so it came out `Float64` even for a `Float32`
  equilibrium; the new answer carries no element type from the equilibrium at all.

  This is breaking for any code that calls `periodicity` on anything from this package, whether an
  equilibrium, a perturbation or a field: the name is no longer answered here at all.

### Fixed

- **`LinearAlgebra` carries a `[compat]` bound.** It was the one entry in `[deps]` without one, so
  a resolve was free to pick a version this package had never been built against.
  `Aqua.test_deps_compat` reports exactly this, and fails without the entry.

- **The generated-code cache is keyed on the shape of the parameters, not only on their number.**
  `parameter_values` flattens a vector parameter entry by entry, and the generated code reads its
  parameter argument positionally, so which slot carries which meaning follows from the split. A
  key carrying only the total let two splits of one type share an entry: lengths `(2, 3)` and
  `(3, 2)` both flatten to five slots, and the second field was served the first one's code. The
  failure was silent — the `SVector` fits, so there was no bounds error, and `parameters(field)`
  read the right struct while every accessor computed from the wrong slots.

  No field this package ships can reach it. Exactly two structs carry a vector parameter,
  `SolovevEquilibrium` and `SolovevXpointEquilibrium`, and each carries exactly one;
  `ZeroPerturbation` has no parameters and `EzCosZPerturbation` one scalar. With a single vector
  per struct the total fixes that vector's length, and with it the split, which is why the old key
  held for every shipped field. An equilibrium of your own with two vector parameters reaches the
  collision with no perturbation involved.

  The shape carries `-1` for a scalar and the length for a vector, so a scalar cannot read as a
  vector of length zero. The number of slots follows from the shape, which makes the shape strictly
  finer than the number it replaces: no pair of fields the old key told apart is merged by the new
  one.

  A value that the generated code bakes in as a literal — anything the `A₁`, `φ` or metric methods
  read that `get_parameters` omits — stays invisible to the key, and no key over parameters can
  see it. The `get_parameters` docstring now says so: list the value as a parameter, or build the
  field with `cache = false`.

- **Contour plots of a field that diverges inside the plot window work again on Makie 0.24.15.**
  That release orders every traced contour line through `canonical_line_order`, which takes the
  smallest vertex of a closed line and then keeps the candidate rotations that equal it. A vertex
  holding a `NaN` equals nothing, so no candidate survives and the reduction over them throws
  `reducing over an empty collection is not allowed`.

  Such a vertex appears wherever the sampled grid meets a singular line. The tracer places each
  vertex at `(level - z₁) / (z₂ - z₁)` along a cell edge, so one non-finite sample makes every
  vertex of the lines through the cells around it `NaN`. `SingularEquilibrium` has `A₁ = A₂ = 0/0`
  and `|B| = 1/0` at the origin, which every grid with odd `nx` and `ny` over a window straddling
  the axis samples — `nx = 37, ny = 53` and `nx = 101, ny = 101` among them. All eight test jobs
  of the CI matrix were red on this, on every operating system and every Julia version.

  The panel values are now made finite before they reach `contour!`: a `NaN` becomes the lowest
  finite sample of the panel, which is the side the tracer already reads it on, and `±Inf` is
  clamped into the finite range. Finite samples pass through untouched, so no other plot changes,
  and the package no longer depends on Makie tolerating a `NaN` vertex. `[compat] Makie = "0.24"`
  is left as it is, because 0.24.14 and earlier were never affected.

- **A field rebuilt after `clear_field_cache!()` agrees with the one it replaces to a few ULP, not
  bit for bit.** Two symbolic traces of one equilibrium need not produce the same expression: the
  simplifier may choose any equivalent form, and SymbolicUtils 4.46.8 chooses a different one for
  about half of the 41 generated functions. The generated bodies then differ, and whether that
  reaches the result depends on what the platform contracts — on aarch64 macOS under Julia 1 it
  moves the last bit of 11 of them, on Windows under the Julia floor it moves none.

  Nothing about a field's accuracy changes; both forms evaluate the same quantity. What changes is
  the guarantee the test suite states. It compared the values of a rebuilt field bit for bit, which
  made it a lottery over the platform: red on aarch64 macOS, green elsewhere for no better reason
  than rounding. It now compares the parameters exactly, and the values at `rtol = 1e-12`, four
  orders tighter than `≈` alone and far tighter than any real error in a formula.

- **Both Penning trap docstrings disagreed with their code, and had done since the fields were
  added.** `PenningTrapBottleEquilibrium` had its two `Bₚ` terms swapped between the vector
  potential and the magnetic field, so both formulas were wrong; at `(0.3, 0.4, 0.5)` with the
  defaults the docstring's `A` gave `[-50, -25, 25]` where the code gives
  `[-37.87, 14.10, -12.00]`. `PenningTrapAsymmetricEquilibrium` gave the first component of `B`
  as `B₀/3` where the curl of its own `A` is `B₀/6`. Both docstrings also named the perturbation
  parameter `B₁` where the struct field is `Bₚ`. The code was right in every case; only the
  documentation changes.

- `contravariant_to_physical` carried `covariant_to_physical`'s docstring verbatim, describing a
  one-form where it takes a vector.

### Removed

- The `@code` macros — `@code` and the per-preset variants `@code_iter`, `@code_nstx`,
  `@code_frc`, `@code_xpoint`, `@code_iter_xpoint`, `@code_nstx_xpoint` and
  `@code_nstx_double_xpoint`. Code generation is now done by calling `FieldFunctions`, which
  returns a struct holding the functions.

- Internal helpers dropped in the SymEngine → Symbolics transition: the hand-rolled common-
  subexpression elimination pass (made redundant by `Symbolics.build_function(...; cse = true)`),
  `code_arguments`, `fnesc`, `replace_expr!`, `symprint`, and the dead `Γ` (Christoffel symbol)
  and `connection` routines.

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

# Usage

Working with ElectromagneticFields.jl has two steps. First an equilibrium is constructed, which
is a small immutable object holding nothing but the parameters of the field. Then Julia code for
evaluating that field is generated from it, and the generated functions are what the actual
computation uses.


## Constructing an Equilibrium

Every field lives in its own submodule and provides an `init` function with sensible defaults, so
the shortest way to get an equilibrium is

```@example usage
using ElectromagneticFields

equ = AxisymmetricTokamakCylindrical.init()
```

All parameters can be passed explicitly, here the major radius ``R_0``, the magnetic field
strength ``B_0`` at the magnetic axis, and the safety factor ``q_0``:

```@example usage
equ = AxisymmetricTokamakCylindrical.init(6.2, 5.3, 2.0)
```

Some fields ship named configurations in addition to `init`, e.g. the Solov'ev equilibrium comes
with parameter sets for ITER, NSTX and a field reversed configuration:

```@example usage
Solovev.ITER()
```


## Generating the Evaluation Code

An equilibrium object on its own cannot be evaluated. The evaluation routines are generated from
it symbolically, differentiated where necessary, and then spliced into a module as plain Julia
functions. There are three ways of doing that.

The most convenient one is the `@code` macro that every field module provides. It takes the same
arguments as `init`, and splices the generated functions into the current module:

```@example usage
AxisymmetricTokamakCylindrical.@code(6.2, 5.3, 2.0)
nothing # hide
```

The functions are now available directly. They all take the time as their first argument,
followed by the three coordinates of the evaluation point, which for this equilibrium are
``(R, Z, \phi)``:

```@example usage
t = 0.0
x = [6.5, 0.5, 0.0]

B(t, x)
```

Both a splatted and a vector call are defined, so `B(t, x)` and `B(t, x...)` are the same thing:

```@example usage
B(t, x...) == B(t, x)
```

### `load_equilibrium`

Where the generated functions should end up in a module of their own rather than in the current
one, use [`load_equilibrium`](@ref):

```julia
module Tokamak end

load_equilibrium(equ; target_module = Tokamak)

Tokamak.B(t, x)
```

There is a catch: the methods `load_equilibrium` defines live in a world age newer than the frame
that called it, so that frame cannot call them. At the top level this is invisible, because the
world age advances between statements, but inside a function the freshly defined `B` is "too new"
and Julia reports a `MethodError` on a name that plainly exists. For those cases pass a function,
which is run through `Base.invokelatest` and receives the target module:

```@example usage
module Tokamak end

load_equilibrium(equ; target_module = Tokamak) do mod
    mod.B(0.0, [6.5, 0.5, 0.0])
end
```

The do-block form returns whatever the function returns.

### `code`

Finally, [`code`](@ref) returns the generated code as an expression instead of evaluating it,
which is what to reach for when the code should be inspected or written to a file:

```@example usage
expr = code(equ)
typeof(expr)
```


## What Gets Generated

Rather more than just the magnetic field. The naming follows the usual conventions of
differential geometry: subscripts denote covariant components, superscripts contravariant ones,
and parenthesised subscripts the components in the physical (orthonormal) frame. See
[Coordinates](coordinates.md) for what those three representations are and
[Fields](fields.md) for how the field quantities below are derived from the vector potential.

| | |
|---|---|
| `A₁, A₂, A₃` / `A¹, A², A³` | components of the vector potential |
| `B` | absolute value of the magnetic field |
| `B₁, B₂, B₃` / `B¹, B², B³` / `B₍₁₎, B₍₂₎, B₍₃₎` | components of the magnetic field |
| `b₁, b₂, b₃`, … | components of the unit vector along the magnetic field |
| `E₁, E₂, E₃` / `E¹, E², E³` | components of the electric field |
| `φ` | electrostatic potential |
| `g`, `ḡ` and `gᵢⱼ`, `gⁱʲ` | metric and its inverse |
| `J` | Jacobian determinant of the coordinate transformation |
| `DF`, `DF̄` | Jacobian matrix of the chart and its inverse |
| `from_cartesian`, `to_cartesian` | coordinate transformations |
| `rangemin`, `rangemax` | bounds of the coordinate domain |
| `orientation` | handedness of the chart, `+1` or `-1` |

In addition, derivatives of most of these quantities are generated, e.g. `dBdx₁`, `dA₂dx₃` or
`d²Bdx₁dx₂`, and the parameters of the equilibrium are spliced in as constants, here `R₀`, `B₀`
and `q₀`:

```@example usage
R₀, B₀, q₀
```

The chart itself is available as well. `J` is the Jacobian determinant, `DF` the Jacobian matrix,
and the two are related through the orientation of the chart — this equilibrium uses a
left-handed ``(R, Z, \phi)`` chart, so the determinant of `DF` is `-J`:

```@example usage
using LinearAlgebra

orientation(), det(DF(t, x)) ≈ orientation() * J(t, x)
```


## Evaluating on a Grid

Since the generated functions are ordinary Julia functions, sampling a field is just a
comprehension:

```@example usage
Rgrid = LinRange(0.5 * R₀, 1.5 * R₀, 100)
Zgrid = LinRange(-0.5 * R₀, 0.5 * R₀, 120)

Bfield = [B(t, Rgrid[i], Zgrid[j], 0.0) for i in eachindex(Rgrid), j in eachindex(Zgrid)]

size(Bfield)
```

Note the index order: `Bfield[i,j]` holds the value at `(Rgrid[i], Zgrid[j])`, which is also the
convention Makie's `contour` expects, so such an array can be passed straight to a plotting call.
See [Plotting](plotting.md) for how to turn this into a figure.

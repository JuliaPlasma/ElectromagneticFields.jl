# Electromagnetic Fields

Every equilibrium in this package is defined by two inputs: the three covariant components of the
magnetic vector potential ``A_i`` and, where there is one, the electrostatic potential
``\varphi``. Everything else — the magnetic field, its magnitude, the unit vector along
it, the perpendicular frame, the electric field, and the derivatives of all of these — is derived
from them symbolically at code generation time.

This page describes what is derived and how, and what the generated functions are called. The three
coordinate representations that keep appearing below — covariant, contravariant and physical — are
the subject of [Coordinates](coordinates.md).

```@example fields
using ElectromagneticFields
using LinearAlgebra

field = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0))

t = 0.0
ξ = [6.5, 0.5, 0.25]
nothing # hide
```


## In the Language of Differential Geometry

The construction is the standard one. The vector potential is a one-form, the magnetic field is its
exterior derivative and hence a two-form, and the vector one usually calls the magnetic field is
what the Hodge star and the musical isomorphism make of that two-form:

```math
\begin{aligned}
A^1 &= A_i \, dx^i , \\
B^2 &= d A^1 = B_{ij} \, dx^i \wedge dx^j , \qquad B_{ij} = \frac{1}{2} (\partial_i A_j - \partial_j A_i ) , \\
B^1 &= \star B^2 = \det DF \, \epsilon_{klm} g^{ik} g^{jl} B_{ij} \, dx^m \equiv B_i \, dx^i , \qquad B_m = \det DF \, \epsilon_{klm} g^{ik} g^{jl} B_{ij} , \\
\vec{B} &= B^{1 \sharp} = g^{ij} B_i \partial_j \equiv B^i \partial_i , \qquad B^i = g^{ij} B_j , \\
| B | &= \sqrt{ i_{\vec{B}} B^1 } = \sqrt{ B^i B_i }
\end{aligned}
```

The determinant appearing in the Hodge star is the *signed* one, ``\det DF = \mathrm{orientation}
\cdot J``, and not the volume element ``J = \sqrt{|g|}`` that the generated function of that name
returns. The star is orientation-dependent, and using ``J`` in place of ``\det DF`` on one of the
left-handed charts reverses the magnetic field. See
[Volume Element and Orientation](@ref) for why several of the charts here are left-handed.

The electric field is the simpler half. There is no time dependence anywhere in the generated code,
so the induction term drops and the electric field is minus the gradient of the electrostatic
potential, a one-form like ``A^1``:

```math
E^1 = - d \varphi \equiv E_i \, dx^i , \qquad E_i = - \partial_i \varphi , \qquad E^i = g^{ij} E_j .
```


## Vector Potential

`A♭` gives the covariant components, and they are exactly what the equilibrium supplied.
For the cylindrical tokamak used here they are

```math
A (R, Z, \phi) = \frac{B_0}{2} \left( R_0 \frac{Z}{R} , \, - R_0 \ln \frac{R}{R_0} , \, \frac{r^2}{q_0} \right)^T .
```

```@example fields
A♭(field, t, ξ)
```

The contravariant components `A♯` are obtained by raising the index with the inverse metric. In
this chart ``g^{11} = g^{22} = 1``, so only the third component differs, by a factor ``1/R^2``:

```@example fields
A♯(field, t, ξ)
```

First and second derivatives are generated as `DA♭` and `DDA♭`, differentiated with respect to the
chart coordinates ``\xi`` rather than cartesian ones: `DA♭(field, t, ξ)[i,j]` is ``\partial_j
A_i``. There are no physical components of `A`.

The vector potential is only determined up to a gauge, and different gauges of the same field are
genuinely different equilibria here, since the generated code follows whatever was written down.
`AxisymmetricTokamakToroidalRegularization` exists for that reason: it carries the same magnetic
field as `AxisymmetricTokamakToroidal`, in a gauge whose poloidal vector potential is regular on the
magnetic axis.


## Magnetic Field

Two naming conventions are easy to trip over, so they are worth stating before the table.

The first is that `B` is a scalar. It is ``|B|``, the magnitude of the magnetic field, and the
components live under the decorated names. The second is that `B♭♭`, with two flats, is the
two-form ``B_{ij}`` and not the covariant vector `B♭` — two lowered indices rather than one.

| | |
|---|---|
| `B` | ``\|B\|``, the magnitude of the magnetic field |
| `B♭` | covariant components |
| `B♯` | contravariant components |
| `B♮` | physical components |
| `B♭♭` | components of the two-form ``B_{ij}`` |
| `DB♭` | derivatives of the covariant components |
| `DB`, `DDB` | first and second derivatives of ``\|B\|`` |

```@example fields
B(field, t, ξ)
```

The three representations, as columns — covariant, contravariant, physical:

```@example fields
[B♭(field, t, ξ) B♯(field, t, ξ) B♮(field, t, ξ)]
```

The covariant toroidal component is the outlier: it is ``B_0 R_0``, larger than the field strength
by very nearly a factor ``R``, because the toroidal coordinate basis vector has length ``R`` rather
than one. It falls short of exactly ``R`` only because ``|B|`` carries the poloidal field as well.
Contracting the covariant components with the contravariant ones nevertheless gives the right
answer, as does the euclidean norm of the physical ones:

```@example fields
Bcov = B♭(field, t, ξ)
Bcon = B♯(field, t, ξ)
Bphy = B♮(field, t, ξ)

magnitudes = (sqrt(Bcon' * Bcov) ≈ B(field, t, ξ), norm(Bphy) ≈ B(field, t, ξ))
@assert all(magnitudes) # hide
magnitudes
```

The two-form is antisymmetric by construction, and its entries are the curl of ``A`` before the
Hodge star has been applied:

```@example fields
B♭♭(field, t, ξ)
```


## The Frame along the Magnetic Field

Guiding-centre and gyrokinetic formulations need a frame adapted to the magnetic field, so one is
generated alongside it. `b` is the unit vector along ``B``,

```math
b = \frac{B}{|B|} ,
```

and `a` and `c` complete it to an orthonormal triad. `a` is the cross product of the first
coordinate basis vector not parallel to `b` with `b` itself, `c` is `b × a`, and both are then
normalised in the metric. Like every vector here they come in all three representations, as
`a♭`, `a♯` and `a♮`.

Orthonormality is only visible as such in the physical components, where the Gram matrix of the
triad is the identity:

```@example fields
F = [a♮(field, t, ξ) b♮(field, t, ξ) c♮(field, t, ξ)]

@assert F' * F ≈ I # hide
round.(F' * F; digits = 12)
```

Derivatives of the unit vector are generated as `Db♭` and `DDb♭` for the covariant components and
`Db♮` for the physical ones. No curvature or torsion is generated; those have to be assembled from
these derivatives.


## Electrostatic Potential and Electric Field

`φ` defaults to zero, so a purely magnetic equilibrium generates a vanishing potential and a
vanishing electric field rather than no function at all. Note that this is `φ` (U+03C6), while
several equilibria also export a coordinate function `ϕ` (U+03D5) for the toroidal angle; they are
different names.

The Penning traps and the quadratic potentials are the fields that define a potential. Taking the
uniform Penning trap, with

```math
\varphi (x,y,z) = - E_0 \left( \frac{x^2}{2} + \frac{y^2}{2} - z^2 \right) ,
```

a second field alongside the tokamak — which is no longer a difficulty, since each one is a value
of its own rather than a set of names spliced into the module:

```@example fields
trap = FieldFunctions(PenningTrapUniformEquilibrium(100.0, 10.0))

x = [0.5, 0.3, 0.2]

φ(trap, t, x)
```

```@example fields
E♭(trap, t, x)
```

The chart is cartesian, so here the covariant, contravariant and physical components all coincide.
Derivatives are generated as `DE♭`. Unlike the magnetic field, the electric field has neither
physical components nor a magnitude function.


## Independence of the Chart

The point of carrying all these representations is that the field itself does not depend on the
chart it is written in. The same tokamak equilibrium is available in cartesian and in cylindrical
coordinates, and at the same physical point the two agree — the contravariant components transform
with `DF`, the covariant ones with `DF̄'`, and the magnitude is invariant.

```@example fields
cylindrical = field
cartesian = FieldFunctions(AxisymmetricTokamakCartesianEquilibrium(6.2, 5.3, 2.0))
nothing # hide
```

```@example fields
x = to_cartesian(cylindrical, t, ξ)

Bcon_car = B♯(cartesian, t, x)
Bcov_car = B♭(cartesian, t, x)

transforms = (DF(cylindrical, t, ξ) * Bcon ≈ Bcon_car,
              DF̄(cylindrical, t, ξ)' * Bcov ≈ Bcov_car,
              B(cylindrical, t, ξ) ≈ B(cartesian, t, x))
@assert all(transforms) # hide
transforms
```

The physical components make the same point more directly. The cartesian chart is its own physical
frame, so the physical components of the cylindrical field are simply the components of the
cartesian one:

```@example fields
physical = B♮(cylindrical, t, ξ) ≈ Bcon_car
@assert physical # hide
physical
```

See [Usage](usage.md) for the code generation machinery used here, and
[Analytic Fields](analytic/index.md) for the fields the package ships.

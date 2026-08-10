# Coordinates

Most of the fields in this package are not most naturally written down in cartesian coordinates. A
tokamak equilibrium wants cylindrical or toroidal coordinates, a Solov'ev equilibrium wants
cylindrical coordinates normalised to the major radius, and only the simplest fields are genuinely
cartesian. ElectromagneticFields.jl therefore lets every equilibrium carry its own chart, and
generates all evaluation routines in the coordinates of that chart.

The price of that freedom is that a vector no longer has just three numbers. In a curvilinear chart
the coordinate basis is neither orthogonal nor normalised, so there are several inequivalent ways to
write the components of one and the same vector. This page describes the three that the package
uses, and how they show up in the generated code.


## Charts

A chart is a map from curvilinear coordinates ``\xi`` to the cartesian coordinates ``x`` of the
ambient space,

```math
F : \xi = (\xi^1, \xi^2, \xi^3) \mapsto x = (x^1, x^2, x^3) ,
```

together with its inverse. An equilibrium declares its chart by overloading the traits `x¹`, `x²`,
`x³` for ``F`` and `ξ¹`, `ξ²`, `ξ³` for ``F^{-1}``, plus the metric coefficients `g₁₁` … `g₃₃`, the
volume element `J`, the orientation, and optionally the bounds `minx¹` … `maxx³` of the coordinate
domain. Everything not overloaded falls back to a cartesian default: identity metric, infinite
range, positive orientation.

Five charts are in use:

| chart | ``(\xi^1, \xi^2, \xi^3)`` | metric ``g_{ij}`` | ``J`` | orientation |
|:--|:--|:--|:--|:--:|
| `CartesianEquilibrium` and its subtypes | ``(x, y, z)`` | ``\mathrm{diag}(1, 1, 1)`` | ``1`` | ``+1`` |
| `AxisymmetricTokamakCylindrical` | ``(R, Z, \phi)`` | ``\mathrm{diag}(1, 1, R^2)`` | ``R`` | ``-1`` |
| `AxisymmetricTokamakToroidal` | ``(r, \theta, \phi)`` | ``\mathrm{diag}(1, r^2, R^2)`` | ``r R`` | ``-1`` |
| `AxisymmetricTokamakToroidalRegularization` | ``(r, \theta, \phi)`` | ``\mathrm{diag}(1, r^2, R^2)`` | ``r R`` | ``-1`` |
| `AbstractSolovevEquilibrium` | ``(R/R_0, Z/R_0, \phi)`` | ``\mathrm{diag}(R_0^2, R_0^2, R^2)`` | ``R R_0^2`` | ``-1`` |

Note that `SolovevSymmetric` is a cartesian equilibrium despite its name, and that the Solov'ev
chart coordinates are dimensionless, which is where the factors of ``R_0^2`` in its metric come
from.

The examples on this page all use the cylindrical tokamak, which is a good representative: it is
curvilinear, its metric is non-trivial, and it is left-handed.

```@example coordinates
using ElectromagneticFields
using LinearAlgebra

AxisymmetricTokamakCylindrical.@code(6.2, 5.3, 2.0)

t = 0.0
ξ = [6.5, 0.5, 0.25]
nothing # hide
```

The chart itself is generated as `to_cartesian` and `from_cartesian`, and componentwise as `x¹`,
`x²`, `x³` and `ξ¹`, `ξ²`, `ξ³`:

```@example coordinates
to_cartesian(t, ξ)
```

```@example coordinates
roundtrip = from_cartesian(t, to_cartesian(t, ξ)) ≈ ξ
@assert roundtrip # hide
roundtrip
```

The bounds of the coordinate domain come as `rangemin` and `rangemax`. For this chart only the
toroidal angle is bounded:

```@example coordinates
rangemin(t, ξ), rangemax(t, ξ)
```


## Tangent Map and Metric

Differentiating the chart gives the tangent map ``DF``, and the metric is the pullback of the
cartesian one along it:

```math
{DF^i}_j = \frac{\partial x^i}{\partial \xi^j} ,
\qquad
g_{ij} = \sum_k {DF^k}_i \, {DF^k}_j ,
\qquad
g^{ij} = (g^{-1})_{ij} .
```

Both are generated in two forms: as individual components `DFᵢⱼ`, `DF̄ᵢⱼ`, `gᵢⱼ` and `gⁱʲ`, and as
matrix-valued wrappers `DF`, `DF̄`, `g` and `ḡ`. The naming is worth a moment: at component level a
raised index pair denotes the inverse metric, `g¹¹` … `g³³`, while at matrix level the inverse
carries an overbar, `ḡ`. The same overbar marks the inverse tangent map ``DF^{-1}``, whose
components are the derivatives of the inverse chart,

```math
{\bar{DF}^i}_j = \frac{\partial \xi^i}{\partial x^j} .
```

```@example coordinates
DF(t, ξ)
```

The three identities relating them hold pointwise:

```@example coordinates
checks = (DF̄(t, ξ) ≈ inv(DF(t, ξ)),
          DF(t, ξ)' * DF(t, ξ) ≈ g(t, ξ),
          ḡ(t, ξ) ≈ inv(g(t, ξ)))
@assert all(checks) # hide
checks
```

Derivatives of the metric are generated as well, `dgᵢⱼdxₖ` and `dgⁱʲdxₖ` for the first and
`d²gᵢⱼdxₖdxₗ` and `d²gⁱʲdxₖdxₗ` for the second, which is what a geometric integrator needs to
assemble Christoffel symbols.


## Covariant, Contravariant and Physical Components

A vector field and a one-form are different objects, but in the presence of a metric they carry the
same information, and the package generates both representations of every field it computes. The
metric raises and lowers the indices:

```math
v_i = g_{ij} \, v^j , \qquad v^i = g^{ij} \, v_j .
```

Components with a lower index are called *covariant*, components with an upper index
*contravariant*. In the code the two conversions are one-liners,
`covariant_to_contravariant` and `contravariant_to_covariant`, and the naming convention follows the
indices directly: `B₁, B₂, B₃` are covariant, `B¹, B², B³` contravariant.

Neither of them is what a measurement returns. Covariant and contravariant components are taken with
respect to the coordinate basis and its dual, and in a curvilinear chart those basis vectors are
neither unit length nor mutually orthogonal — the covariant toroidal component of ``B`` below is
larger than ``|B|`` by nearly a factor of ``R``, and carries different units than the other two
components.
What is usually quoted in the literature, and what this package calls the *physical* components, are
the components in an orthonormal frame:

```math
v_{(i)} = {DF^i}_j \, v^j = {\bar{DF}^j}_i \, v_j .
```

That is, physical components are obtained by pushing the vector forward to the ambient cartesian
frame. They all carry the same units, and their euclidean norm is the length of the vector. They are
written with parenthesised indices, `B₍₁₎, B₍₂₎, B₍₃₎`.

Some texts define physical components instead by normalising the coordinate basis vectors,
``v_{\langle i \rangle} = \sqrt{g_{ii}} \, v^i``. For an orthogonal chart that is also an
orthonormal frame, so it gives the same magnitude, but its components differ from the ones here by
the rotation relating the local frame to the cartesian one.

The unit vector along the magnetic field, `b`, is available in all three representations, as the
vector-valued wrappers `b` (covariant), `b⃗` (contravariant) and `bₚ` (physical):

```@example coordinates
[b(t, ξ) b⃗(t, ξ) bₚ(t, ξ)]
```

The three columns are the same vector. Lowering, raising and pushing forward take one into another:

```@example coordinates
conversions = (g(t, ξ) * b⃗(t, ξ) ≈ b(t, ξ),
               ḡ(t, ξ) * b(t, ξ) ≈ b⃗(t, ξ),
               bₚ(t, ξ) ≈ DF(t, ξ) * b⃗(t, ξ),
               bₚ(t, ξ) ≈ DF̄(t, ξ)' * b(t, ξ))
@assert all(conversions) # hide
conversions
```

Only two of the three give the length of the vector directly. Contracting a covariant with a
contravariant vector is a metric-free operation and needs no ``g``; the same contraction of the
physical components is an ordinary euclidean dot product. Contracting covariant with covariant, on
the other hand, is meaningless.

```@example coordinates
@assert b⃗(t, ξ)' * b(t, ξ) ≈ 1 && norm(bₚ(t, ξ)) ≈ 1 # hide
b⃗(t, ξ)' * b(t, ξ), norm(bₚ(t, ξ))
```

Summarising the notation:

| | |
|---|---|
| `Bᵢ`, e.g. `B₁` | covariant components |
| `Bⁱ`, e.g. `B¹` | contravariant components |
| `B₍ᵢ₎`, e.g. `B₍₁₎` | physical components |
| `b`, `a`, `c` | vector-valued wrapper, covariant |
| `b⃗`, `a⃗`, `c⃗` | vector-valued wrapper, contravariant |
| `bₚ`, `aₚ`, `cₚ` | vector-valued wrapper, physical |
| `g`, `DF` | matrix-valued metric and tangent map |
| `ḡ`, `DF̄` | their inverses |


## Volume Element and Orientation

`J` is the volume element ``\sqrt{|g|} = |\det DF|``. It is not computed from the metric but
declared by the chart, so that the generated expressions stay in the closed form the chart author
intended rather than passing through a symbolic square root.

Being an absolute value, `J` says nothing about handedness, and handedness matters: the Hodge star
that turns the magnetic two-form into a one-form, and the cross product used to build the
perpendicular frame, both need the *signed* determinant

```math
\det DF = \mathrm{orientation} \cdot J .
```

Four of the five charts above are left-handed, because the toroidal angle sits in the third slot
where the right-handed ordering would put the second poloidal coordinate — ``(R, Z, \phi)`` rather
than ``(R, \phi, Z)``. This is not a corner case, and getting it wrong has no visible symptom other
than a magnetic field pointing the wrong way. Each equilibrium therefore declares its handedness
explicitly, and the generated module exposes it as `orientation()`, the one generated function that
takes no arguments at all, since the sign depends on neither time nor position. See
[`orientation`](@ref ElectromagneticFields.orientation) for the full discussion.

```@example coordinates
orientation(), J(t, ξ), det(DF(t, ξ))
```

```@example coordinates
signed = det(DF(t, ξ)) ≈ orientation() * J(t, ξ)
@assert signed # hide
signed
```


## Adding a Chart

A new chart is defined by overloading, for the new equilibrium type,

* `x¹`, `x²`, `x³` — the chart map to cartesian coordinates,
* `ξ¹`, `ξ²`, `ξ³` — its inverse,
* `g₁₁` … `g₃₃` — the metric coefficients, in the chart's own coordinates,
* `J` — the volume element,
* `orientation` — `+1` or `-1`, and
* optionally `minx¹` … `maxx³` — the bounds of the coordinate domain.

Everything else — the inverse metric, the tangent map and its inverse, all derivatives, and the
whole tower of field representations described in [Fields](fields.md) — is derived symbolically from
these by the code generator. A perturbation added to an equilibrium must live on the same chart; the
generator asserts that the two agree on `J` and `g`.

The consistency of a new chart is checked by the test suite, which asserts
``J = \sqrt{\det(DF^T DF)}`` and ``\det DF = \mathrm{orientation} \cdot J`` for every equilibrium,
so a chart whose metric, volume element and handedness do not fit together fails immediately.

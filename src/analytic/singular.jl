
const DEFAULT_SINGULAR_B₀ = 1.0

@doc raw"""
Singular magnetic field in (x,y,z) coordinates with covariant components of
the vector potential given by
```math
A (x,y,z) = \frac{B_0}{\sqrt{(x^2 + y^2)}^3} \big( y , \, - x , \, 0 \big)^T
```
resulting in the magnetic field with covariant components
```math
B(x,y,z) = B_0 \, \begin{pmatrix}
0 \\
0 \\
(x^2 + y^2)^{-3/2} \\
\end{pmatrix}
```

Parameters: `B₀`
"""
struct SingularEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    B₀::T
    SingularEquilibrium{T}(B₀::T) where {T <: Number} = new("Singular Magnetic Field", B₀)
end

function SingularEquilibrium(B₀::T = DEFAULT_SINGULAR_B₀) where {T <: Number}
    SingularEquilibrium{T}(B₀)
end

function Base.show(io::IO, equ::SingularEquilibrium)
    print(io, equ.name)
end

r²(x::AbstractVector, equ::SingularEquilibrium) = X(x, equ)^2 + Y(x, equ)^2
r(x::AbstractVector, equ::SingularEquilibrium) = sqrt(r²(x, equ))
R(x::AbstractVector, equ::SingularEquilibrium) = r(x, equ)
θ(x::AbstractVector, equ::SingularEquilibrium) = atan(Y(x, equ), X(x, equ))
ϕ(x::AbstractVector, equ::SingularEquilibrium) = θ(x, equ)

A₁(x::AbstractVector, equ::SingularEquilibrium) = +equ.B₀ * x[2] / sqrt(x[1]^2 + x[2]^2)^3
A₂(x::AbstractVector, equ::SingularEquilibrium) = -equ.B₀ * x[1] / sqrt(x[1]^2 + x[2]^2)^3
A₃(x::AbstractVector, equ::SingularEquilibrium) = zero(eltype(x))

B(x::AbstractVector, equ::SingularEquilibrium) = equ.B₀ / sqrt(x[1]^2 + x[2]^2)^3

function get_functions(::SingularEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ, r² = r²)
end

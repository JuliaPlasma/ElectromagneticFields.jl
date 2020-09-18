@doc raw"""
Symmetric quadratic equilibrium in (x,y,z) coordinates with covariant components of
the vector potential given by
```math
A (x,y,z) = \frac{B_0}{4} \big( - y \, (2 + x^2 + y^2) , \, x \, (2 + x^2 + y^2) , \, 0 \big)^T
```
resulting in the magnetic field with covariant components
```math
B(x,y,z) = B_0 \, \begin{pmatrix}
0 \\
0 \\
1 + x^2 + y^2 \\
\end{pmatrix}
```

Parameters: `B₀`
"""
struct SymmetricQuadratic{T <: Number} <: CartesianEquilibrium
    name::String
    B₀::T
    SymmetricQuadratic{T}(B₀::T) where T <: Number = new("SymmetricQuadraticEquilibrium", B₀)
end

SymmetricQuadratic(B₀::T) where T <: Number = SymmetricQuadratic{T}(B₀)


function Base.show(io::IO, equ::SymmetricQuadratic)
    print(io, "Symmetric Quadratic Equilibrium")
end


r²(x::AbstractVector, equ::SymmetricQuadratic) = X(x,equ)^2 + Y(x,equ)^2
r(x::AbstractVector, equ::SymmetricQuadratic) = sqrt(r²(x,equ))
R(x::AbstractVector, equ::SymmetricQuadratic) = r(x,equ)
θ(x::AbstractVector, equ::SymmetricQuadratic) = atan(Y(x,equ), X(x,equ))
ϕ(x::AbstractVector, equ::SymmetricQuadratic) = θ(x,equ)

A₁(x::AbstractVector, equ::SymmetricQuadratic) = - equ.B₀ * x[2] * (2 + x[1]^2 + x[2]^2) / 4
A₂(x::AbstractVector, equ::SymmetricQuadratic) = + equ.B₀ * x[1] * (2 + x[1]^2 + x[2]^2) / 4
A₃(x::AbstractVector, equ::SymmetricQuadratic) = zero(eltype(x))

get_functions(::SymmetricQuadratic) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)

macro symmetric_quadratic_equilibrium()
    generate_equilibrium_code(SymmetricQuadratic(1.); output=false)
end

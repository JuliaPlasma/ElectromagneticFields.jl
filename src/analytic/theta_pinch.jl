@doc raw"""
θ-pinch equilibrium in (x,y,z) coordinates with covariant components of
the vector potential given by
```math
A (x,y,z) = \frac{B_0}{2} \big( - y , \, x , \, 0 \big)^T
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = \big( 0 , \, 0 , \, B_0 \big)^T .
```

Parameters:
    B₀: B-field at magnetic axis
"""
struct ThetaPinch{T <: Number} <: CartesianEquilibrium
    name::String
    B₀::T

    function ThetaPinch{T}(B₀::T) where T <: Number
        new("ThetaPinchEquilibrium", B₀)
    end
end

ThetaPinch(B₀::T=1.0) where T <: Number = ThetaPinch{T}(B₀)


function Base.show(io::IO, equ::ThetaPinch)
    print(io, "θ-Pinch Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  B₀ = ", equ.B₀)
end


R(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = r(x,equ)
r(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = sqrt(r²(x,equ))
r²(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = X(x,equ)^2 + Y(x,equ)^2
θ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = atan(Y(x,equ), X(x,equ))
ϕ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = θ(x,equ)

A₁(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = - equ.B₀ * Y(x,equ) / 2
A₂(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = + equ.B₀ * X(x,equ) / 2
A₃(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number} = zero(eltype(x))


get_functions(::ThetaPinch) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)


macro zpinch_equilibrium(R₀, B₀)
    generate_equilibrium_code(ThetaPinch(R₀, B₀); output=false)
end

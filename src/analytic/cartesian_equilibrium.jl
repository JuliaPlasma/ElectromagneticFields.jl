module AnalyticCartesianField

    import ..ElectromagneticFields
    import ..ElectromagneticFields: AnalyticEquilibrium, AnalyticPerturbation

    export CartesianField, CartesianEquilibrium, CartesianPerturbation

    abstract type CartesianEquilibrium <: AnalyticEquilibrium end
    abstract type CartesianPerturbation <: AnalyticPerturbation end

    CartesianField = Union{CartesianEquilibrium, CartesianPerturbation}

    X(x::AbstractVector, ::CartesianField) = x[1]
    Y(x::AbstractVector, ::CartesianField) = x[2]
    Z(x::AbstractVector, ::CartesianField) = x[3]

    ElectromagneticFields.J(x::AbstractVector, ::CartesianField) = one(eltype(x))

    ElectromagneticFields.x¹(ξ::AbstractVector, ::CartesianField) = ξ[1]
    ElectromagneticFields.x²(ξ::AbstractVector, ::CartesianField) = ξ[2]
    ElectromagneticFields.x³(ξ::AbstractVector, ::CartesianField) = ξ[3]

    ElectromagneticFields.ξ¹(x::AbstractVector, ::CartesianField) = x[1]
    ElectromagneticFields.ξ²(x::AbstractVector, ::CartesianField) = x[2]
    ElectromagneticFields.ξ³(x::AbstractVector, ::CartesianField) = x[3]

    ElectromagneticFields.g₁₁(x::AbstractVector, ::CartesianField) = one(eltype(x))
    ElectromagneticFields.g₂₂(x::AbstractVector, ::CartesianField) = one(eltype(x))
    ElectromagneticFields.g₃₃(x::AbstractVector, ::CartesianField) = one(eltype(x))

end

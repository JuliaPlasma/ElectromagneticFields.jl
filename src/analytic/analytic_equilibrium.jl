
using SymPy

simplify(x::Real) = x
simplify(x::SymPy.Sym) = x

abstract type AnalyticEquilibrium <: MagneticEquilibrium end
abstract type AnalyticPerturbation <: MagneticEquilibrium end

analyticA₁(x::AbstractArray{T,1}, equ::ET) where {T,ET} = error("analyticA₁() not implemented for ", ET)
analyticA₂(x::AbstractArray{T,1}, equ::ET) where {T,ET} = error("analyticA₂() not implemented for ", ET)
analyticA₃(x::AbstractArray{T,1}, equ::ET) where {T,ET} = error("analyticA₃() not implemented for ", ET)
analyticMetric(x::AbstractArray{T,1}, equ::ET) where {T,ET} = error("analyticMetric() not implemented for ", ET)
analyticHcoeffs(x::AbstractArray{T,1}, equ::ET) where {T,ET} = error("analyticHcoeffs() not implemented for ", ET)


struct ZeroPerturbation <: AnalyticPerturbation
    name::String
    ZeroPerturbation() = new("ZeroPerturbation")
end

analyticA₁(x::AbstractArray{T,1}, pert::ZeroPerturbation) where {T} = zero(T)
analyticA₂(x::AbstractArray{T,1}, pert::ZeroPerturbation) where {T} = zero(T)
analyticA₃(x::AbstractArray{T,1}, pert::ZeroPerturbation) where {T} = zero(T)


"""
Generate functions for evaluating analytic equilibria.
"""
function generate_equilibrium_functions(equ::AnalyticEquilibrium, pert::AnalyticPerturbation; output=0)
    # define symbols for coordinates x = (x₁, x₂, x₃),
    # positive=true is set so that sqrt(x^2) does not become |x^2|
    # x₁, x₂, x₃ = symbols("x₁, x₂, x₃", real=true, positive=true)
    # x = [x₁, x₂, x₃]
    x1, x2, x3 = symbols("x1, x2, x3", real=true, positive=true)
    x = [x1, x2, x3]

    # obtain metric
    g = analyticMetric(x, equ)
    symprint("g", g, output, 2)

    # compute metric coefficients for derivatives
    h = [simplify(sqrt(g[i,i])) for i in 1:3]
    symprint("h", h, output, 2)

    # compute Jacobian
    J = simplify(sqrt(det(g)))
    symprint("J", J, output, 2)

    # obtain vector potential
    A  = [analyticA₁(x, equ) + analyticA₁(x, pert),
          analyticA₂(x, equ) + analyticA₂(x, pert),
          analyticA₃(x, equ) + analyticA₃(x, pert)]
    symprint("A", A, output, 2)

    # compute Jacobian of vector potential A
    DA = [simplify(diff(A[i], x[j])) for i in 1:3, j in 1:3]
    symprint("DA", DA, output, 2)

    # compute second derivative of vector potential A
    DDA1 = [simplify(diff(diff(A[1], x[i]), x[j])) for i in 1:3, j in 1:3]
    symprint("DDA1", DDA1, output, 2)

    DDA2 = [simplify(diff(diff(A[2], x[i]), x[j])) for i in 1:3, j in 1:3]
    symprint("DDA2", DDA2, output, 2)

    DDA3 = [simplify(diff(diff(A[3], x[i]), x[j])) for i in 1:3, j in 1:3]
    symprint("DDA3", DDA3, output, 2)

    # compute magnetic field B = curl A
    B  = [( simplify(diff(h[3] * A[3], x[2]) - diff(h[2] * A[2], x[3])) * h[1] / J ),
          ( simplify(diff(h[1] * A[1], x[3]) - diff(h[3] * A[3], x[1])) * h[2] / J ),
          ( simplify(diff(h[2] * A[2], x[1]) - diff(h[1] * A[1], x[2])) * h[3] / J )]
    symprint("B", B, output, 2)

    # compute absolute value |B| of B
    Babs = simplify( sqrt(B[1]^2 + B[2]^2 + B[3]^2) )
    symprint("|B|", Babs, output, 2)

    # compute magnetic unit vector
    b  = [simplify( B[1] / Babs ),
          simplify( B[2] / Babs ),
          simplify( B[3] / Babs )]
    symprint("b", b, output, 2)

    # compute Jacobian of magnetic field B
    DB = [diff(B[i], x[j]) for i in 1:3, j in 1:3]
    symprint("DB", DB, output, 2)

    # compute Jacobian of magnetic unit vector b
    Db = [diff(b[i], x[j]) for i in 1:3, j in 1:3]
    symprint("Db", Db, output, 2)

    # compute second derivative of magnetic unit vector b
    DDb1 = [diff(diff(b[1], x[i]), x[j]) for i in 1:3, j in 1:3]
    symprint("DDb1", DDb1, output, 2)

    DDb2 = [diff(diff(b[2], x[i]), x[j]) for i in 1:3, j in 1:3]
    symprint("DDb2", DDb2, output, 2)

    DDb3 = [diff(diff(b[3], x[i]), x[j]) for i in 1:3, j in 1:3]
    symprint("DDb3", DDb3, output, 2)

    # compute first derivatives of absolute value of magnetic field
    DBabs = [diff(Babs, x[j]) for j in 1:3]
    symprint("D|B|", DBabs, output, 2)

    # compute second derivatives of absolute value of magnetic field
    DDBabs = [diff(diff(Babs, x[i]), x[j]) for i in 1:3, j in 1:3]
    symprint("DD|B|", DDBabs, output, 2)


    # collect all functions to generate code for
    functions = Dict{String,Any}()
    indices   = ["₁", "₂", "₃"]

    functions["J"] = J
    functions["B"] = Babs

    for i in 1:3
        functions["h" * indices[i]] = h[i]
        functions["A" * indices[i]] = A[i]
        functions["B" * indices[i]] = B[i]
        functions["b" * indices[i]] = b[i]
        functions["dBdx" * indices[i]] = DBabs[i]
    end

    for i in 1:3
        for j in 1:3
            functions["g"  * indices[i] * indices[j]] = g[i,j]
            functions["dA" * indices[i] * "dx" * indices[j]] = DA[i,j]
            functions["dB" * indices[i] * "dx" * indices[j]] = DB[i,j]
            functions["db" * indices[i] * "dx" * indices[j]] = Db[i,j]

            functions["d²A₁" * "dx" * indices[i] * "dx" * indices[j]] = DDA1[i,j]
            functions["d²A₂" * "dx" * indices[i] * "dx" * indices[j]] = DDA2[i,j]
            functions["d²A₃" * "dx" * indices[i] * "dx" * indices[j]] = DDA3[i,j]
            functions["d²b₁" * "dx" * indices[i] * "dx" * indices[j]] = DDb1[i,j]
            functions["d²b₂" * "dx" * indices[i] * "dx" * indices[j]] = DDb2[i,j]
            functions["d²b₃" * "dx" * indices[i] * "dx" * indices[j]] = DDb3[i,j]
            functions["d²B"  * "dx" * indices[i] * "dx" * indices[j]] = DDBabs[i,j]
        end
    end

    functions
end


"""
Generate code for evaluating analytic equilibria.
"""
function generate_equilibrium_code(equ, pert=ZeroPerturbation(); output=0)

    if output ≥ 1
        println()
        println(equ)
        println()
    end

    functions = generate_equilibrium_functions(equ; output=output)

    equ_code = :(  )

    # generate Julia code and export functions
    for (key, value) in functions
        f_symb = Symbol(key)
        f_expr = value

        output ≥ 1 ? println("Generating function ", key) : nothing

        # patch for removing Warnings with Julia v0.6 and SymPy v1.0
        f_str  = sympy_meth(:julia_code, f_expr)
        f_str  = replace(f_str, ".+" => " .+ ")
        f_str  = replace(f_str, ".-" => " .- ")
        f_str  = replace(f_str, ".*" => " .* ")
        f_str  = replace(f_str, "./" => " ./ ")
        f_str  = replace(f_str, ".^" => " .^ ")
        f_body = Meta.parse(f_str)
        # f_body = Meta.parse(sympy_meth(:julia_code, f_expr))
        output ? println("   ", f_body) : nothing

        f_code = quote
            export $(esc(f_symb))
            # function $(esc(f_symb))(x₁, x₂, x₃)
            function $(esc(f_symb))(x1, x2, x3)
                $f_body
            end
            function $(esc(f_symb))(t::Number, x::Vector)
                $(esc(f_symb))(x[1],x[2],x[3])
            end
        end

        # append f_code to equ_code
        push!(equ_code.args, f_code)
    end

    return equ_code
end


"""
Evaluate functions for evaluating analytic equilibria.
"""
function load_equilibrium(equ, pert=ZeroPerturbation(); target_module=Main, output=0)

    if output ≥ 1
        println()
        println(equ)
        println()
    end

    functions = generate_equilibrium_functions(equ, pert; output=output)

    # generate Julia code and export functions
    for (key, value) in functions
        f_symb = Symbol(key)
        f_expr = value

        output ≥ 1 ? println("Generating function ", key) : nothing

        f_str  = sympy_meth(:julia_code, f_expr)
        f_str  = replace(f_str, ".+" => " .+ ")
        f_str  = replace(f_str, ".-" => " .- ")
        f_str  = replace(f_str, ".*" => " .* ")
        f_str  = replace(f_str, "./" => " ./ ")
        f_str  = replace(f_str, ".^" => " .^ ")
        f_body = Meta.parse(f_str)
        # f_body = Meta.parse(sympy_meth(:julia_code, f_expr))

        output ≥ 2 ? println("   ", f_body) : nothing

        f_code = quote
            export $f_symb
            # function $f_symb(x₁, x₂, x₃)
            function $f_symb(x1, x2, x3)
                $f_body
            end
            function $f_symb(t::Number, x::AbstractArray{T,1}) where {T <: Number}
                $f_symb(x[1],x[2],x[3])
            end
        end

        Core.eval(target_module, f_code)
    end
end


function symprint(name, symexpr, output=1, detail_level=0)
    if output ≥ detail_level
        println(name, " = ", symexpr, "\n")
    end
end

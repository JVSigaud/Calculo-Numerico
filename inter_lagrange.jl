using BenchmarkTools
using LinearAlgebra
# using plots

struct init{T}
    x₀::T
    y₀::T
end

Base.@kwdef mutable struct interpol
    dados::init
    size_::Int32=size(dados.x₀,1)
    W::Vector = Vector{Float64}(undef,size_)
    Diff::Dict{Int16,Vector{Float64}} = Dict(0=> dados.y₀)
end

function diff_div!(pol::interpol)
    size::Int16 = pol.size_ - 1
    for i in 1:size        
        pol.Diff[i] = [(pol.Diff[i-1][j+1] - pol.Diff[i-1][j]) / (pol.dados.x₀[j+i] - pol.dados.x₀[j]) for j in 1:pol.size_-i]
    end
end

function lagrange!(P::interpol)
    for i in 1:P.size_
        P.W[i] = P.dados.y₀[i]/(prod(P.dados.x₀[i] .- deleteat!(P.dados.x₀[1:end],i)))
    end
end

function fit_lagrange(P::interpol,x::Float64)
    aprox = Vector{Float64}(undef,P.size_)
    for i in 1:P.size_
        aprox[i] = prod(x .- deleteat!(P.dados.x₀[1:end],i))
    end 
    return sum(P.W .* aprox)
end

# function fit_lagrange(P::interpol,x::LinRange{Float64,Int64})
#     aprox = Vector{Float64}(undef,P.size_)
#     y = similar(aprox)
#     for j in P.size_
#         for i in 1:P.size_
#             aprox[i] = prod(x .- deleteat!(P.dados.x₀[1:end],i))
#         end 
#     y[j] = sum(P.W .* aprox)
#     end
# end


b = interpol(dados = init{Vector}([1,4,5],[1,2,0]))
@btime lagrange!(b)
@btime fit_lagrange(b,2.0)


# println(fit_lagrange(b,3.0))

# X = LinRange(-1,1,1000)
# println(X)
# Y = 1/(1 .+ X .^ 2)
# println(typeof(Y))
# @btime lagrange!(interpol(dados = init(X::LinRange{Float64, Int64}, Y::Transpose{Float64, Vector{Float64}})))

# println(size(Z.W))


# inter = fit_lagrange(Z,X)
# @btime fit_lagrange(Z,X)
# println(b.dados.x₀)
# println(b.W)
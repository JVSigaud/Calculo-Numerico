using BenchmarkTools
using LinearAlgebra


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


#teste
b = interpol(dados = init{Vector}([1,4,5],[1,2,0]))
@btime lagrange!(b) 
@btime fit_lagrange(b,2.0) #valor de x na aproximação da interpolação



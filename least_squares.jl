using LinearAlgebra
using BenchmarkTools
using LinearSolve
using Statistics

struct init{T,t}
    x₀::T
    y₀::t
end

Base.@kwdef mutable struct Minimos_quadrados
    dados::init
    size_g::Int64
    size_::Int64=size(dados.x₀,1)
    equation_g::Vector{Function} = Vector{Function}(undef,size_g)
    A::Matrix{Float64} = zeros(Float64,size_g,size_g)
    B::Vector{Float64}= Vector{Float64}(undef,size_g)
    u::Vector{Float64}=Vector{Float64}(undef,size_g)

end

function g_vector!(P::Minimos_quadrados)
    @inbounds for i in 0:P.size_g-1
        f(x::Vector{Float64}) = x.^i
        P.equation_g[i+1] = f
    end
end
function A_matrix!(P::Minimos_quadrados)
    for j in 1:P.size_g, i in 1:P.size_g
        P.A[i,j] = dot(transpose(P.equation_g[i](P.dados.x₀)),P.equation_g[j](P.dados.x₀))
    end
end

function b_vector!(P::Minimos_quadrados)
    for i in 1:P.size_g
        P.B[i] = dot(transpose(P.equation_g[i](P.dados.x₀)),P.dados.y₀)
    end
end
function solve!(P::Minimos_quadrados)
    g_vector!(P)
    A_matrix!(P)
    b_vector!(P)
    prob = LinearProblem(P.A,P.B)  
    sol = solve(prob)
    P.u = sol.u
    corr = Statistics.cor(P.dados.x₀,P.dados.y₀)
    println("coeficientes:\na = $(b.u[2])\nb = $(b.u[1])\nCoeficiente de correlação:\nr = $(corr)\nCoeficiente de determinação:\nr² = $(corr^2)")
end
  
function fit(P::Minimos_quadrados,x::Vector)
return P.u[1] .* x .+ P.u[2]
end
  
a = [5.2,6.8,11.2,16.8,17.8,19.6,23.4,25.4,32.0,44.8]
c = [1.1,1.222,1.373,1.571,1.833,2.200,2.75,3.666,5.5,11.0]
c= log.(c)
println(c)
b = Minimos_quadrados(dados = init{typeof(a),typeof(c)}(a,c),size_g=2)
solve!(b) #obtenção dos coeficientes

x_fit = collect(range(0,10,length = 10))
y_fit = fit(b,x_fit) #Não rwsolvido para o grau N.


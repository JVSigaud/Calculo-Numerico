using Zygote
using CairoMakie

Base.@kwdef mutable struct init
    ϵ::Float64 = 10^(-16)
    x₀::Real = 1.5
    pontos = Vector{Real}()
    iter = 0.0
end

function Metodo_Newton!(l::init,f::Function)
    ψ = f(l.x₀)
    while abs(ψ) > l.ϵ
        l.x₀ = l.x₀ - ψ/f'(l.x₀)
        l.iter += 1
        ψ = f(l.x₀)
        push!(l.pontos,ψ)
    end
    return nothing
end

#Exemplo 

f(ϕ) = exp(-ϕ^2)-cos(ϕ) 
β = init()
Metodo_Newton!(β,f)

aprox = β.x₀
iter = β.iter

using Zygote
using CairoMakie
CairoMakie.activate!()

Base.@kwdef mutable struct init
    ϵ::Float64 = 10^(-16)
    x₀::Real = 1.5
    pontos = Vector{Real}()
    iter = 0.0
end

f(ϕ) = exp(-ϕ^2)-cos(ϕ) 
# f(ϕ) = sin(10*cos(ϕ)+ϕ^2)

function Metodo_Newton!(l::init)
    ψ = f(l.x₀)
    while abs(ψ) > l.ϵ
        l.x₀ = l.x₀ - ψ/f'(l.x₀)
        l.iter += 1
        ψ = f(l.x₀)
        push!(l.pontos,ψ)
    end
    return nothing
end

β = init()
Metodo_Newton!(β)

aprox = β.x₀
iter = β.iter
print("Aproximação = $aprox  iterações = $iter f(x) = $(f(aprox))")

fig = lines(1:β.iter,abs.(β.pontos))
for i in 0:10
    print(i)
end

x = [i for i in 0:10]
println(x)
println(x[1:11])
println(prod(x[2:4]))



# for i in 1:pol.size_-1
#     pol.Diff[i] = [(pol.Diff[i-1][j] - pol.Diff[i-1][j+1]) / (pol.dados.x₀[i+j] - pol.dados.x₀[i]) for j in 1:pol.size_-i]
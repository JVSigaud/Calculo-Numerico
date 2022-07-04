using BenchmarkTools
using CairoMakie
CairoMakie.activate!()

struct Bound{t}
    limite::t
end

#Apenas para Runge Kutta de segunda ordem
Base.@kwdef mutable struct RungeKutta
    bounds::Bound
    y₀::Vector{Float64}
    h::Float64
    n::Int64 = floor(Int,abs((bounds.limite[end]-bounds.limite[1])/h))
    xᵢ::Vector{Float64} = collect(range(bounds.limite[1],bounds.limite[2],length = n+1))
    yᵢ::Matrix{Float64} = zeros(Float64,2,n+1)
end

function solve!(p::RungeKutta,F)
    G(x,y,z)= [z;F(x,y,z)]
    p.yᵢ[:,1] .= p.y₀
    @inbounds for i in range(1,p.n)
        p.yᵢ[:,i+1] .= p.yᵢ[:,i]  .+ p.h .* G(p.xᵢ[i],p.yᵢ[:1,i],p.yᵢ[:2,i])
        p.yᵢ[:,i+1] .= p.yᵢ[:,i] .+ (p.h/2) .* (G(p.xᵢ[i],p.yᵢ[:1,i],p.yᵢ[:2,i]) .+ G(p.xᵢ[i+1],p.yᵢ[:1,i+1],p.yᵢ[:2,i+1]))
        # println("Iteração $i\nθ = $(p.yᵢ[:1,i])\nω = $(p.yᵢ[:2,i])")
        # println("")
    end
    return nothing
end


#exemplo
intervalo = (0,10.0)
H = 0.05
f(x,y,z) = -(9.81/1)*sin(y)
y = [π/3,0]

β = RungeKutta(bounds = Bound{typeof(intervalo)}(intervalo),y₀ = y, h = H)
@btime solve!(β,f)


lines(β.xᵢ,β.yᵢ[1,:]; label="Original", linewidth=2, linestyle=nothing,yticks=range(0,2),
    figure=(; figure_padding=5, resolution=(1080, 800), font="sans",
        backgroundcolor=:grey90, fontsize=16),
    axis=(; xlabel="t", title="RungeKuttaSegundaOrdem", xgridstyle=:dash,
        ygridstyle=:dash))

axislegend("legend")
current_figure()





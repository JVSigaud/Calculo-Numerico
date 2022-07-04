using BenchmarkTools
struct init{T,t}
    x₀::T
    y₀::t
end

Base.@kwdef mutable struct interpol
    dados::init
    size_::Int32=size(dados.x₀,1)
    Diff::Dict{Int32,Vector{Float64}} = Dict(0=> dados.y₀)
    W::Vector{Float64} = Vector{Float64}(undef,size_)

end

function diff_div!(pol::interpol)
    size::Int16 = pol.size_ - 1
    for i in 1:size        
        pol.Diff[i] = [(pol.Diff[i-1][j+1] - pol.Diff[i-1][j]) / (pol.dados.x₀[j+i] - pol.dados.x₀[j]) for j in 1:pol.size_-i]
    end
end

function mult!(Pol::interpol,x::Float64,vec::Vector{Float64})
    size::Int16 = Pol.size_ - 1
    prod::Float64 =1.0
    for i in 1:size
        prod = 1.0
        for j in 1:i
            prod *= x - Pol.dados.x₀[j]
        end
        vec[i] = Pol.Diff[i][1]*prod
    end
    return vec
end



function fit_newton(Pol::interpol,x::Float64)
    vec::Vector{Float64} = Vector{Float64}(undef,Pol.size_-1)
    vec = mult!(Pol,x,vec) 
    return  Pol.dados.y₀[1] + sum(vec)
end

# b = interpol(dados = init{Vector}([-1,0,1],[2,0,-1]))
# @btime diff_div!(b)
# println(b.Diff)
# @btime fit_newton(b,0.5)

using CairoMakie
CairoMakie.activate!()


x₀ = range(-5,5,length = 10000)
f(x) = 1 / (1 + x^2)
y₀ = map(f,x₀)

lines(x₀,y₀; label="Original", linewidth=2, linestyle=nothing,yticks=range(0,2),
    figure=(; figure_padding=5, resolution=(1080, 800), font="sans",
        backgroundcolor=:grey90, fontsize=16),
    axis=(; xlabel="x", title="Interpolação Newton", xgridstyle=:dash,
        ygridstyle=:dash))

for i in range(4,11,step = 2)
  x₂ = range(-5,5,i)
  y₂ = map(f,x₂)
  newton = interpol(dados = init{typeof(x₂),typeof(y₂)}(x₂,y₂))

  diff_div!(newton)

  test = Vector{Float64}(undef,size(y₀))
  function teste!(test::Vector{Float64})
      for (i,j) in enumerate(x₀)
          test[i] = fit_newton(newton,j)
      end
      
  end
  teste!(test)

  # scatterlines!(x₁,y₁; label="spline",color="red")
  lines!(x₀,test;label = "newton$i")
end
axislegend("legend")
current_figure()
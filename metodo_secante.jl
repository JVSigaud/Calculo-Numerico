Base.@kwdef mutable struct init
    ϵ::Float64 = 10^(-6)
    x₀::Real = 40
    x₁::Real = 45
    k::Float64 = 0.08
    g::Int32 = 32
    W::Real = 527
    D::Int64 = -300
    B::Int64 = 470
    F::Function = x ->  ((1/((k^2)*g))*(W*(W-B)*log(1+(k*x/(W-B)))-W*k*x) - D)
    # F::Function = x ->  sin(10*cos(x)+x^2)
    pontos = Vector{Real}()
    iter= 0.0

end

function metodo_secante!(C::init)
    ϕ = C.F(C.x₁)
    while abs(C.x₀-C.x₁) >= C.ϵ
        x₍ₖ₊₁₎= C.x₁ - (ϕ/((ϕ-C.F(C.x₀))/(C.x₁-C.x₀))) 
        C.x₀ = C.x₁
        C.x₁ = x₍ₖ₊₁₎
        C.iter += 1 
        ϕ = C.F(C.x₁)
    end
    
end

function testa_vel(resultado::init)
    if resultado.x₁ > 40
        println("O barril foi danificado!") 
    else
         println("O barril não foi danificado!")
    end
end


#1)
println("Questão 1 e 2")
β = init()

metodo_secante!(β)

println("Aproximação: $(β.x₁)")
println("Número de iterações: $(β.iter)")
println("f(aprox): $(β.F(β.x₁))")
testa_vel(β)
print("\n")

#2)
println("Questão 3")
β = init(D = -200) #Basta trocar as condições iniciais para novos testes.
#Método que implementa e altera os parametros iniciais necessarios do algoritmo até a convergencia.
metodo_secante!(β)

println("Aproximação: $(β.x₁)")
println("Número de iterações: $(β.iter)")
println("f(aprox): $(β.F(β.x₁))")
testa_vel(β)

# # ComponentArrays
# https://github.com/SciML/ComponentArrays.jl are array blocks that can be accessed through a named index.
# You can compose differential equations without a DSL like `ModelingToolkit`.
using ComponentArrays
using OrdinaryDiffEq
using SimpleUnPack: @unpack

# ## Composing Lorenz and Lotka-Volterra systems
# ### Lorenz system
function lorenz!(D, u, p, t; f=0.0)
    @unpack σ, ρ, β = p
    @unpack x, y, z = u

    D.x = σ*(y - x)
    D.y = x*(ρ - z) - y - f
    D.z = x*y - β*z
    return nothing
end

#---
tspan = (0.0, 20.0)
lorenz_p = (σ=10.0, ρ=28.0, β=8/3)
lorenz_ic = ComponentArray(x=1.0, y=0.0, z=0.0)
lorenz_prob = ODEProblem(lorenz!, lorenz_ic, tspan, lorenz_p)

# ### Lotka-Volterra system
function lotka!(D, u, p, t; f=0.0)
    @unpack α, β, γ, δ = p
    @unpack x, y = u

    D.x =  α*x - β*x*y + f
    D.y = -γ*y + δ*x*y
    return nothing
end

#---
lotka_p = (α=2/3, β=4/3, γ=1.0, δ=1.0)
lotka_ic = ComponentArray(x=1.0, y=1.0)
lotka_prob = ODEProblem(lotka!, lotka_ic, tspan, lotka_p)

# ### Composed system
function composed!(D, u, p, t)
    c = p.c
    @unpack lorenz, lotka = u

    lorenz!(D.lorenz, lorenz, p.lorenz, t, f=c*lotka.x)
    lotka!(D.lotka, lotka, p.lotka, t, f=c*lorenz.x)
    return nothing
end

#---
comp_p = (lorenz=lorenz_p, lotka=lotka_p, c=0.01)
comp_ic = ComponentArray(lorenz=lorenz_ic, lotka=lotka_ic)
comp_prob = ODEProblem(composed!, comp_ic, tspan, comp_p)

# ### Solve problem
# We can solve the composed system
@time comp_sol = solve(comp_prob, Tsit5())

# Extract a state variable
ts = 0:0.1:20
map(u -> u.lotka.x, comp_sol(ts).u)

# We can also unit test one of the component systems
lotka_sol = solve(lotka_prob, Tsit5())
lotka_sol(1.0).x

# ## Symbolic mapping in `ODEFunction`
# - `syms=keys(carray)`.
# - Only valid in non-nested `ComponentArrays`.

lotka_func = ODEFunction(lotka!; syms=keys(lotka_ic))
@time lab_sol = solve(ODEProblem(lotka_func, lotka_ic, tspan, lotka_p), Tsit5());

lab_sol(0.0:0.1:20, idxs=:x)

# ## Update a subset of the elements in ComponentArrays
using ForwardDiff
newval = ForwardDiff.Dual(3.0, 1.0)
v = ComponentVector(a = 1.0, b = 2.0, c = 3.0)

# This will throw an error because the new value has a different type than the existing values
try
    ComponentVector(v; a = newval)
catch e
    println("Error: ", e)
end

# Find the common type of existing values and new value
T = promote_type(eltype(newval), Float64)
# Convert the existing values to the common type
vT = ComponentVector{T}(v)
# Create a new ComponentVector with the updated value
ComponentVector(vT; a = newval)

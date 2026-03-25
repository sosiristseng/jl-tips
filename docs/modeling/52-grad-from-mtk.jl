# # Get vector fields from ModelingToolkit systems
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots

"Numpy's meshgrid clone"
function meshgrid(xrange, yrange)
    xx = [x for y in yrange, x in xrange]
    yy = [y for y in yrange, x in xrange]
    return xx, yy
end

"Get the gradients of the vector field of an ODE system at a grid of points in the state space."
function get_gradient(prob, xsym, ysym, xx, yy; t = nothing)
    ## The order of state variables (unknowns) in the ODE system is not guaranteed. So we may need to swap the order of x and y when calling ∂F.
    swap_or_not(x, y; xidx=1) = xidx == 1 ? [x, y] : [y, x]
    ∂F = prob.f
    ps = prob.p
    sys = prob.f.sys
    xidx = ModelingToolkit.variable_index(sys, xsym)
    yidx = ModelingToolkit.variable_index(sys, ysym)
    dx = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[xidx], xx, yy)
    dy = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[yidx], xx, yy)
    return (dx, dy)
end

"Normalize the gradients for better visualization. The `transform` function can be used to apply a non-linear transformation to the gradients."
function normalize_gradient(dx, dy, xscale, yscale; transform=cbrt)
    maxdx = maximum(abs, dx)
    maxdy = maximum(abs, dy)
    dx_norm = @. transform(dx / maxdx) * xscale
    dy_norm = @. transform(dy / maxdy) * yscale
    return (dx_norm, dy_norm)
end

# Example non-linear ODE system
function model(; tend = 1.5)
    @parameters k1=20 k2=5 k3=5 k4=5 k5=2 n=4
    @variables s1(t)=0.0 s2(t)=0.0
    eqs = [
        D(s1) ~ k1 / (1 + s2^n) - (k3 + k5) * s1
        D(s2) ~ k2 + k5 * s1 - k4 * s2
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time "Build problem" prob = model()
@unpack s1, s2 = prob.f.sys
xrange = 0:0.1:2
yrange = 0:0.1:2
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(prob, s1, s2, xx, yy)
dx_norm, dy_norm = normalize_gradient(dx, dy, step(xrange), step(yrange))
quiver(xx, yy, quiver=(dx_norm, dy_norm); aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), color=:gray)

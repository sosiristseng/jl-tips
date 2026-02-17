# # Get vector fields from ModelingToolkit systems
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots

# Get the 2D gradient from the ODE system
function get_gradient(prob, xsym, ysym; t = nothing, xrange=range(0, 2, 21), yrange=range(0, 2, 21))
    ## The order of unkowns in the ODE system is not guaranteed. So we may need to swap the order of x and y when calling ∂F.
    swap_or_not(x, y; xidx=1) = xidx == 1 ? [x, y] : [y, x]
    ∂F = prob.f
    ps = prob.p
    sys = prob.f.sys
    xidx = ModelingToolkit.variable_index(sys, xsym)
    yidx = ModelingToolkit.variable_index(sys, ysym)
    xx = [x for y in yrange, x in xrange]
    yy = [y for y in yrange, x in xrange]
    dx = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[xidx], xx, yy)
    dy = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[yidx], xx, yy)
    return (; xx, yy, dx, dy)
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
@time prob = model()
@unpack s1, s2 = prob.f.sys
@unpack xx, yy, dx, dy = get_gradient(prob, s1, s2)

# Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.1
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength

quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), color=:gray)

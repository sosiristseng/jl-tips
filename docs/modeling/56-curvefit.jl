#===
# Curve fitting

- [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl)
- [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl)
- [CurveFit.jl](https://github.com/SciML/CurveFit.jl)

## LsqFit

`LsqFit.jl` package is a small library that provides basic least-squares fitting in pure Julia.
===#

using LsqFit
@. model(x, p) = p[1] * exp(-x * p[2])

## Generate data
xdata = range(0, stop=10, length=20)
ydata = model(xdata, [1.0 2.0]) + 0.01 * randn(length(xdata))
## Initial guess
p0 = [0.5, 0.5]

## Fit the model
@time fit = curve_fit(model, xdata, ydata, p0; autodiff=:forwarddiff)

## The result should be close to `[1.0 2.0]`
coef(fit)

#===
## CurveFit

Linear, special function, and nonlinear curve fitting in Julia.

The algorithms used in `CurveFit.jl` are better suited for ill-conditioned nonlinear systems as stated in [JuliaCon 2025](https://youtu.be/mdcCjaYSNNc)

See a list of algorithms in the [documentation](https://docs.sciml.ai/CurveFit/stable/tutorials/getting_started/).
===#
using CurveFit

## Define a nonlinear function: y = a[1] + a[2] * x^a[3]
fn(a, x) = @. a[1] + a[2] * x^a[3]

## True parameters
true_params = [3.0, 2.0, 0.7]

## Generate sample data
x = collect(1.0:0.5:10.0)
y = fn(true_params, x)

## Create problem with initial guess for parameters
u0 = [0.5, 0.5, 0.5]
prob = NonlinearCurveFitProblem(fn, u0, x, y)
@time sol = solve(prob)

println("Fitted parameters: ", sol.u)
println("Prediction at x=5: ", sol(5.0))

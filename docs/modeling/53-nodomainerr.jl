#===
# Avoid DomainErrors

Some functions such as `sqrt(x)`, `log(x)`, and `pow(x)`, throw `DomainError` exceptions with negative `x`, interrupting differential equation solvers. The respective functions in https://github.com/JuliaMath/NaNMath.jl return `NaN` instead of throwing a `DomainError`. Then, the differential equation solvers will reject the solution and retry with a smaller time step.
===#
import NaNMath as nm
nm.sqrt(-1.0) ## returns NaN

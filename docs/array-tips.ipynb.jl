# # Julia array tips
# ## Convert a vector of vectors to a matrix
# Use the [stack](https://docs.julialang.org/en/v1/base/arrays/#Base.stack) function.

stack([[1,0,1],[0,0,1]], dims=1) == [ 1 0 1 ; 0 0 1 ]

# If you only need broadcast operations, [RecursiveArrayTools.jl](https://github.com/SciML/RecursiveArrayTools.jl) wraps a vector of vectors.
using RecursiveArrayTools

a = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
b = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
vA = VectorOfArray(a)
vB = VectorOfArray(b)

vA .* vB

# Convert a VectorOfArray into a matrix.
Matrix(vA .* vB)

# ## Convert a linear index to 2D indices
# Use `CartesianIndices((nrow, ncol))`: See https://discourse.julialang.org/t/julia-usage-how-to-get-2d-indexes-from-1d-index-when-accessing-a-2d-array/61440
x = rand((7, 10))
CI = CartesianIndices((7, 10))
for i in 1:length(x)
    r = CI[i][1]
    c = CI[i][2]
    @assert x[i] == x[r, c]
end

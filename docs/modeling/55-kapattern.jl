#===

# King-Altman method for enzyme kinetics

To calculate steady-state enzyme state fractions in a catalytic cycle using the [King-Altman method](https://link.springer.com/article/10.1186/1471-2105-10-238)

1. Calculate the kinetic matrix, where Ae = de (`A` is the matrix, `e` is the vector of enzyme states, and `de` is the time derivative of `e`)
2. To find the weight of enzyme state `i`, make a (n-1)x(n-1) matrix by skipping the `i`-th row and `i`-th column of the input matrix
3. Calculate the determinant of the matrix from step 2. That is the weight of enzyme state `i`.
4. The fraction of enzyme state `i` is the weight of state `i` over the sum of all states.

===#

using Symbolics
using LinearAlgebra: det

# Make a (n-1)x(n-1) matrix by skipping the i-th row and i-th column of the input matrix
function skip_colrow(mat, i::Int)
    return mat[[1:i-1; i+1:size(mat, 1)], [1:i-1; i+1:size(mat, 2)]]
end

# Accumulate rates into the transition rate (kinectic) matrix
function accumulate_rate!(mat, rate, src::Int, dst::Int)
	mat[dst, src] += rate
	mat[src, src] -= rate
	return mat
end

# Take a complicated symbolic matrix as an example (Complex I model in Gauthier, 2013)
mat_c1g = let
	@variables a12 a21 a65 a56 a61 a16 a23 a32 a34 a43 a47 a74 a57 a75 a42 a24
	mat = fill(Num(0), 7, 7)
	accumulate_rate!(mat, a12, 1, 2)
	accumulate_rate!(mat, a21, 2, 1)
	accumulate_rate!(mat, a65, 6, 5)
	accumulate_rate!(mat, a56, 5, 6)
	accumulate_rate!(mat, a61, 6, 1)
	accumulate_rate!(mat, a16, 1, 6)
	accumulate_rate!(mat, a23, 2, 3)
	accumulate_rate!(mat, a32, 3, 2)
	accumulate_rate!(mat, a34, 3, 4)
	accumulate_rate!(mat, a43, 4, 3)
	accumulate_rate!(mat, a47, 4, 7)
	accumulate_rate!(mat, a74, 7, 4)
	accumulate_rate!(mat, a57, 5, 7)
	accumulate_rate!(mat, a75, 7, 5)
	accumulate_rate!(mat, a42, 4, 2)
	accumulate_rate!(mat, a24, 2, 4)
end

# The weight of each state is
@time weights_c1g = [(-1)^(7-1) * det(skip_colrow(mat_c1g, i)) |> expand for i in 1:7]

# Total weight is the sum of all weights
total_weight_c1g = sum(weights_c1g)

# And the fraction of each state: fi = wi / wTotal
# For example, the fraction of state 1 is
fraction_state1_c1g = weights_c1g[1] / total_weight_c1g

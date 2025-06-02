using KohnShamResolution
using BenchmarkTools

original_stdout = stdout
output_file = open("tests/Performance/polynomial.txt", "w")
redirect_stdout(output_file)

P = RandPolynomial(30,-30)
Q = RandPolynomial(30,-30)
x = 10.0
n = 10

println("P = $P")
println("Q = $Q")
println("x = $x")
println("n = $n")

# Performance operation

println("Evaluation")
@btime P(x)

println("P + x")
@btime P + x
println("P * x")
@btime P * x
println("P ^ n")
@btime P ^ n

println("P + Q")
@btime P + Q
println("P * Q")
@btime P * Q


redirect_stdout(original_stdout)
close(output_file)

println("Performance finished")
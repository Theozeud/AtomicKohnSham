# Theoretical Eigenvalues
function theoretical_eigenvalue(n, z, l, T = Float64)
    -T(z)^2/(T(2)*T(n+l)^2)
end

function theoretical_eigenvalue(problem)
    theoretical_eigenvalue.(problem.opts.nλ, problem.z, problem.l, problem.T)
end

# Theoretical Eigenvectors
function theoretical_eigenvector(n, z, l, T = Float64)
    exp_part(x) = exp(-z*abs(x))*x
    C = zeros(T, n-l)
    C[1] = z^(3/2)/sqrt(π)
    for i in 2:(n - l)
        C[i] = -2/n * (n-l-i)/(i*(2*l+i+1)) * C[i - 1]
    end
    return x -> begin
        val = C[n - l]
        for j in (n - l - 1):-1:1
            val = val * x + C[j]
        end
        val * exp_part(x) * x^l
    end
end

function theoretical_eigenvector(problem)
    theoretical_eigenvector.(problem.opts.nU, problem.z, problem.l, problem.T)
end

using Pkg
Pkg.activate("./RK")
using RK

using Plots

function OutputSolution()
    # Output numerical solution in the end
    nz = length(history(tspan[1]))
    if nz == 1
        println("y(", tspan[2], ") = ", sol[2][end])
    else
        for i = 1 : nz
            println("y", i, "(", tspan[2], ") = ", (sol[2]')[end - nz + i])
        end
    end
end


# Speed test
stepsize = 1e-2;

# Example 1
idefun(t,y,z,i) = ((t - 1) * exp(t^2) * i) / (exp(-1) * y[1] - 1)
K(t,s,y)        = y * exp(-s * t)
delays_int(t)   = t - 1
history(t)      = exp(t)

tspan = [1.1 5]

print("Elapsed time:")
@time sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize)

OutputSolution()

# Example 2
idefun(t,y,z,i) = -y.^2 - t * exp(t^2) * z^4 * i
K(t,s,y)        = y * exp(s - s * t)
delays(t,y)     = t / 2
delays_int(t)   = t - 1
history(t)      = exp(-t)

print("Elapsed time:")
@time sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)

OutputSolution()

# Example 3
idefun(t,y,z,i) = exp(1) - exp(t^2) / (z[1]^2) * (i[1] - exp(-2 * t) * i[2]) * (t - 1)
K(t,s,y)        = [ y * exp(-s * t);
                    y * exp(t * (2 - s))]
delays(t,y)     = t - 1
delays_int(t)   = [ t - 1;
                    t - 2 ]
history(t)      = exp(t)

tspan = [0 5]

print("Elapsed time:")
@time sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)

OutputSolution()

# Example 4
idefun(t,y,z,i) = -z[1]^((t + 1) / 2) * z[2] * y.^2 * (1 + exp(t^2) * t * i[1]) / exp(0.5)
K(t,s,y)        = y * exp(s - s * t)
delays(t,y)     = [ (log(y[1]))^2 / (t + 1) - 0.5;
                    (t - 1) / 4 ]
delays_int(t)   = t / 2 - 1
history(t)      = exp(-t)

tspan = [0 5]

print("Elapsed time:")
@time sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)

OutputSolution()

# Example 5
idefun(t,y,z,i) = [ -2.5 * i[1];
                    -15 * i[2];
                    0.001 * exp(t)]
K(t,s,y)        = [ sin(y[1]);
                    sin(y[2])]
delays_int(t)   = [ t - 1;
                    t - 1 ]
history(t)      = [ 1.5;
                    1.5;
                    0.001 * exp(t)]

tspan = [0 10]

print("Elapsed time:")
@time sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize)

OutputSolution()
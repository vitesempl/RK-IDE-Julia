using Pkg
Pkg.activate("./")
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


# Example 6 (system of equations)
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
stepsize = 1e-2

sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize)

OutputSolution()

display(plot(sol, 
            linewidth=3, 
            xlabel="TIME", 
            ylabel="SOLUTION", 
            xticks=trunc(Int,tspan[1]):(tspan[2]-trunc(Int,tspan[1]))/10:tspan[2], 
            title="Example 6 (system of equations)", 
            legend=:outertopright))

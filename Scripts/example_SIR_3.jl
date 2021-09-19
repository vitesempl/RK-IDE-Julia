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


# Example (Yukihiko)
idefun(t,y,z,i) = -2.5 * i
K(t,s,y)        = sin(y)
delays_int(t)   = t - 1
history(t)      = 1.5

tspan = [0 25]
stepsize = 1e-2

sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize)

OutputSolution()

display(plot(sol, 
            linewidth=3, 
            xlabel="TIME", 
            ylabel="SOLUTION", 
            xticks=trunc(Int,tspan[1]):(tspan[2]-trunc(Int,tspan[1]))/10:tspan[2], 
            title="Example (Yukihiko)", 
            legend=:outertopright))

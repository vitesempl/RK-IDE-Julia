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


function OutputConvOrder(; dde=false)
    # Check convergence order of examples
    # "dde=true" - if equation is DDE (with dicrete delays)
    n = 8
    err = zeros(n)
    nsteps = zeros(n)

    for steppow = 1 : n
        stepsize = (2.0)^(-steppow)
        
        if dde
            sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)
        else
            sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize)
        end

        err[steppow] = abs(analytic_sol - sol[2][end])
        nsteps[steppow] = stepsize
    end

    println()
    println("Convergence order: ", (log10(err[end-1]) - log10(err[end])) / (log10((2.0)^(-n+1)) - log10((2.0)^(-n))))

    X_ticks = [(10.0)^(-x) for x in 1:15]
    Y_ticks = [(10.0)^(-x) for x in 0:0.5:3]
    display(plot(err, 
                nsteps, 
                linewidth=3, 
                xaxis=:log10, 
                yaxis=:log10, 
                xticks=X_ticks, 
                yticks=Y_ticks, 
                xlabel="ERROR", 
                ylabel="STEPSIZE",
                legend=false, 
                title="Convergence order"))
end


# Example 2 (integral+discrete delays)
idefun(t,y,z,i) = (1 + exp(-pi / 2)) * y - exp(-pi / 2) * z - 2 * exp(-2 * t) * i
K(t,s,y)        = y * exp(t + s)
delays(t,y)     = t - pi / 2
delays_int(t)   = t - pi / 2
history(t)      = cos(t)

tspan = [0 10]
stepsize = 1e-2

sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)

OutputSolution()

display(plot(sol, 
            linewidth=3, 
            xlabel="TIME", 
            ylabel="SOLUTION", 
            xticks=trunc(Int,tspan[1]):(tspan[2]-trunc(Int,tspan[1]))/10:tspan[2], 
            title="Example 2 (integral+discrete delays)", 
            legend=:outertopright))

# Check convergence order
fun(t) = cos(t)
analytic_sol = fun(tspan[end]);

# dde=true - if equation is DDE (with dicrete delays 'z')
OutputConvOrder(dde=true) 

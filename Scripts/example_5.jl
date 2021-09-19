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


# Example 5 (2 discrete delays + State-dependent)
idefun(t,y,z,i) = -z[1]^((t + 1) / 2) * z[2] * y.^2 * (1 + exp(t^2) * t * i[1]) / exp(0.5)
K(t,s,y)        = y * exp(s - s * t)
delays(t,y)     = [ (log(y[1]))^2 / (t + 1) - 0.5;
                    (t - 1) / 4 ]
delays_int(t)   = t / 2 - 1
history(t)      = exp(-t)

tspan = [0 5]
stepsize = 1e-2

sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)

OutputSolution()

display(plot(sol, 
            linewidth=3, 
            xlabel="TIME", 
            ylabel="SOLUTION", 
            xticks=trunc(Int,tspan[1]):(tspan[2]-trunc(Int,tspan[1]))/10:tspan[2], 
            title="Example 5 (2 discrete delays + State-dependent)", 
            legend=:outertopright))

# Check convergence order
fun(t) = exp(-t)
analytic_sol = fun(tspan[end]);

# dde=true - if equation is DDE (with dicrete delays 'z')
OutputConvOrder(dde=true) 

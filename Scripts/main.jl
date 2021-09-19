#using Pkg
#if !(haskey(Pkg.dependencies(), Base.UUID("7b8691fb-d5bd-47c8-98ff-aba88156f9b3")))
#    Pkg.add(path="./")
#end

#include("./RK.jl")
#using .RK
# using RK

using Pkg
Pkg.activate("./RK")
using RK

using Plots

# Output numerical solution in the end
function OutputSolution()
    nz = length(history(tspan[1]))
    if nz == 1
        println("y(",tspan[2],") = ",sol[2][end])
    else
        for i = 1 : nz
            println("y",i,"(",tspan[2],") = ",(sol[2]')[end-nz+i])
        end
    end
end

# Check 4-order of examples
# 1 - if equation has discrete delays z, 0 - if otherwise
function OutputConvOrder(; dde=false)
    n      = 8
    err    = zeros(n)
    nsteps = zeros(n)

    for steppow = 1:n
        stepsize = (2.0)^(-steppow)
        
        if dde
            sol = ide_solve(idefun,K,delays_int,history,tspan,stepsize,delays)
        else
            sol = ide_solve(idefun,K,delays_int,history,tspan,stepsize)
        end

        err[steppow] = abs(analytic_sol - sol[2][end])
        nsteps[steppow] = stepsize
    end

    println("Convergence order: ",(log10(err[end-1])-log10(err[end]))/(log10((2.0)^(-n+1))-log10((2.0)^(-n))))

    X_ticks = [(10.0)^(-x) for x in 1:10]
    display(plot(err, 
                nsteps, 
                xaxis=:log10, 
                yaxis=:log10, 
                xticks=X_ticks, 
                xlabel="ERROR", 
                ylabel="STEPSIZE",
                legend=false, 
                title="Convergence order"))
end

# Example 1 (only integral)
tspan = [1.1 5]
stepsize = 1e-2
idefun(t,y,z,i) = ((t - 1) * exp(t^2) * i) / (exp(-1) * y[1] - 1)
K(t,s,y)        = y * exp(-s * t)
delays_int(t)   = t - 1 # lower integration limit
history(t)      = exp(t)

sol = ide_solve(idefun,K,delays_int,history,tspan,stepsize)

OutputSolution()

display(plot(sol, 
            linewidth=2, 
            xlabel="TIME", 
            ylabel="SOLUTION", 
            xticks=trunc(Int,tspan[1]):(tspan[2]-trunc(Int,tspan[1]))/10:tspan[2], 
            title="Example 1 (Only integral)", 
            legend=:outertopright))

# Check convergence order of Example 1
fun(t) = exp(t)
analytic_sol = fun(tspan[end])

# dde=true - if equation is DDE (with dicrete delays)
OutputConvOrder() 


# Example 2 (integral+discrete delays)
tspan = [0 10];
stepsize = 1e-2
idefun(t,y,z,i) = (1+exp(-pi/2))*y-exp(-pi/2)*z-2*exp(-2*t)*i;
K(t,s,y)        = y*exp(t+s); 
delays(t,y)     = t-pi/2;
delays_int(t)   = t-pi/2;
history(t)      = cos(t);

sol = ide_solve(idefun,K,delays_int,history,tspan,stepsize,delays);

OutputSolution()

display(plot(sol, 
            linewidth=2, 
            xlabel="TIME", 
            ylabel="SOLUTION", 
            xticks = 0:tspan[2]/10:tspan[2], 
            title="Example 2 (integral+discrete delays)", 
            legend=:outertopright))

# Check convergence order of Example 2
fun(t) = cos(t)
analytic_sol = fun(tspan[end])

# dde=true - if equation is DDE (with dicrete delays)
OutputConvOrder(dde=true)

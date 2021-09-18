using Pkg
if !("RK" in keys(Pkg.installed()))
    Pkg.add(path="./")
end

using RK
using Plots

# Output numerical solution in the end
function Output()
    nz = length(history(tspan[1]));
    if nz == 1
        println("y(",tspan[2],") = ",sol[2][end]);
    else
        for i = 1 : nz
            println("y",i,"(",tspan[2],") = ",(sol[2]')[end-nz+i]);
        end
    end
end

# Check 4-order of examples
# 1 - if equation has discrete delays z, 0 - if otherwise
function err_calc(hasDelays)
    n      = 8;
    err    = zeros(n);
    nsteps = zeros(n);

    for steppow = 1:n
        stepsize = (2.0)^(-steppow);
        
        if hasDelays == 1
            sol = ide_delay_solve(idefun,delays,K,delays_int,history,tspan,stepsize);
        else
            sol = ide_solve(idefun,K,delays_int,history,tspan,stepsize);
        end

        err[steppow] = abs(true_sol - sol[2][end]);
        nsteps[steppow] = stepsize;
    end

    println("Convergence order: ",(log10(err[end-1])-log10(err[end]))/(log10((2.0)^(-n+1))-log10((2.0)^(-n))));

    display(plot(err, nsteps, xaxis=:log10, yaxis=:log10, xlabel="ERROR", ylabel="STEPSIZE",legend=false, title = "Convergence order"))
end

# Example 1 (only integral)
tspan = [1.1 5];
idefun(t,y,i) = ((t-1)*exp(t*t)*i)/(exp(-1)*y[1]-1);
K(t,s,y)      = y*exp(-s*t);
delays_int(t) = t-1; # delays of integrals
history(t)    = exp(t);

sol = ide_solve(idefun,K,delays_int,history,tspan,1e-2);

# Output vector of times t     #print(sol[1])
# Output vector of solutions y #print(sol[2])
Output();

display(plot(sol, linewidth=2, xlabel="TIME", ylabel="SOLUTION", title="Example 1 (Only integral)", legend=:outertopright))

# Check 4-order of Example 1
function fun(t)
    return exp(t);
end
true_sol = fun(tspan[end]);

# 1 - if equation has discrete delays z, 0 - if otherwise
err_calc(0) 

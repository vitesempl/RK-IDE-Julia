# RK-IDE-Julia

Package for solving Differential equations with Discrete and Distributed delays by the explicit Runge-Kutta method.

# Tools and packages

* Julia 1.6 (VERSION ≥ v"1.2")
* Plots
* Jupyter Notebook / Google Colaboratory

# Installation

Enter the `Pkg` REPL from the `Julia` REPL. To add a package, use `add`:
```
pkg> add https://github.com/vitesempl/RK-IDE-Julia
```
or add package from the `Julia` REPL:
```
using Pkg
Pkg.add(url="https://github.com/vitesempl/RK-IDE-Julia")
```

# Usage

```
using RK
```

```
sol = ide_solve(idefun, K, delays_int, history, tspan, stepsize, delays)
```

|  #  | Argument    | Description |
| --- | :---        |    :---     |
|  1  | idefun      | right-hand side function (*t* - time, *y* - solution, *z* - discrete delays, *i* - integrals) |
|  2  | K           | Kernel (integrated function) |
|  3  | delays_int  | distributed delays function (lower integration limit) |
|  4  | history     | history function |
|  5  | tspan       | solution interval |
|  6  | stepsize    | step of numerical Runge-Kutta method |
|  7  | delays      | (optional) discrete delays function (if idefun has 'z') |
|  8  | overlapping | (optional) if equation has overlapping in discrete delays. This option uses the 7-step method |


Examples of use can be found in the Notebook and in scripts from the folder `/Scripts`.

# Authors

* Aleksandr Lobaskin (Saint Petersburg State University)
* Alexey Eremin (Saint Petersburg State University)

# Special thanks

* Yukihiko Nakata (Shimane University) (examples from SIR-models)

# License

"RK-IDE-Julia" is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).
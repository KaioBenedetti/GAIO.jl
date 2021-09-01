using GAIO
using ForwardDiff
using StaticArrays
using BenchmarkTools

@inline function rk4(f, x, τ)
    τ½ = τ/2

    k = f(x)
    dx = @. k/6

    k = f(@. x+τ½*k)
    dx = @. dx + k/3

    k = f(@. x+τ½*k)
    dx = @. dx + k/3

    k = f(@. x+τ*k)
    dx = @. dx + k/6

    return @. x+τ*dx
end

function lorenz_v(x)
    s = 10.0
    rh = 28.0
    b = 0.4
    return s * (x[2] - x[1]), rh * x[1] - x[2] - x[1] * x[3], x[1] * x[2] - b * x[3]
end

@inline function lorenz_f(x)
    h = 0.01
    n = 20

    for i = 1:n
        x = rk4(lorenz_v, x, h)
    end

    return x
end

domain = Box((0.0, 0.0, 30.0), (30.0, 30.0, 40.0))
partition = RegularPartition(domain, 7)
B = partition[:]

no_of_points = 100
F = BoxMap(lorenz_f, domain, no_of_points=no_of_points)
trials = @benchmark F(B)
t = minimum(trials).time

length(B)*no_of_points*20*(4*9+45)/(t*1e-9)

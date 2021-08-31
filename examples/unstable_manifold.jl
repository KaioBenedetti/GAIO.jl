using GAIO

# the Lorenz system
const σ, ρ, β = 10.0, 28.0, 0.4
v((x,y,z)) = (σ*(y-x), ρ*x-y-x*z, x*y-β*z)
f(x) = rk4_flow_map(v, x)

center, radius = (0.0, 0.0, 25.0), (30.0, 30.0, 30.0)
Q = Box(center, radius)
depth = 18
P = RegularPartition(Q, depth)

F = BoxMap(f, Q)

x = (sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1)         # equilibrium
W = unstable_set!(F, P[x])

plot(W)

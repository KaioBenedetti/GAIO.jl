var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In the following, we will showcase some of the algorithms GAIO is capable of based on one example each.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"tip: Tip\nUsing SVector instead of Vector for higher-dimensional objects will in general reduce the computation time. However, in the GAIO environment, you are free to use any of the two.","category":"page"},{"location":"examples/#Table-of-Contents","page":"Examples","title":"Table of Contents","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Pages = [\"examples.md\"]","category":"page"},{"location":"examples/#Relative-Global-Attractor-of-the-Hénon-System","page":"Examples","title":"Relative Global Attractor of the Hénon System","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The Hénon map is the two-dimensional quadratic map ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign x_n+1 = 1-alpha x_n^2+y_n \ny_n+1 = beta x_nendalign ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"characterizing the Hénon system. The so called classical Hénon map has the parameters alpha = 14 and beta = 03.  In the following, we will demonstrate how to compute the relative global attractor for the classical Hénon map. Let us start by generating n points on each face of the square -11^2 subset mathbbR^2 as well as initializing the Hénon map f and the boxmap g corresponding to the dynamics of f","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"function henon()\n    generate_points = n -> [\n        [(x, -1.0) for x in LinRange(-1, 1, n)];\n        [(x,  1.0) for x in LinRange(-1, 1, n)];\n        [(-1.0, x) for x in LinRange(-1, 1, n)];\n        [( 1.0, x) for x in LinRange(-1, 1, n)];\n    ]\n\n    f = x -> SVector(1/2 - 2*1.4*x[1]^2 + x[2], 0.3*x[1])\n       \n    g = PointDiscretizedMap(f, generate_points(20))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"note: Note\nIn order to fit the relative global attractor into -11^2 we had to scale the first equation.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In order to create the boxset which is one of the input parameters of the function relative_attractor, we first need to initialize the domain, that is the box corrsponding to -11^2, and how/if it is already partitioned. It is natural to choose a regular partition which is initialized with the whole domain on depth 0, that is the partition is not yet subdivided.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    partition = RegularPartition(Box(SVector(0.0, 0.0), SVector(1.0, 1.0)))\n    boxset = boxset_full(partition)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Finally we are able to call relative_attractor by","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    return relative_attractor(boxset, g, 20)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"which will return the boxset corresponding to the relative global attractor we receive after subdividing the domain twenty times. We can review the result by plotting the attractor with the plot-command.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"[insert beautiful plot here]","category":"page"},{"location":"examples/#Chain-Recurrent-Set-for-the-Knotted-Flow-Map","page":"Examples","title":"Chain Recurrent Set for the Knotted Flow Map","text":"","category":"section"},{"location":"examples/#Root-Covering","page":"Examples","title":"Root Covering","text":"","category":"section"},{"location":"examples/#Unstable-Manifold-for-the-Lorenz-System","page":"Examples","title":"Unstable Manifold for the Lorenz System","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Let us consider the Lorenz System, which is the following three-dimensional continuous system ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign fracmathcaldxmathcaldt = s(y-x)  \nfracmathcaldymathcaldt = rx - y-xz \n fracmathcaldzmathcaldt = xy - bz  endalign","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example, we will choose the parameter values as s = 10 r = 28 b = 04 and we are looking for the unstable manifold through the equilibrium point x_0 = (sqrt(b*(r-1)) sqrt(b*(r-1)) r-1 ), which is a subset of the Lorenz attractor.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In order to compute this, we first need to initialize the function lorenz_f, which we will achieve by solving the Lorenz System with a Runge-Kutta 4th order ODE solver.  By choosing an equally spaced covering of the cube -11^3 we can initialize the boxmap g:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    grid = LinRange(-1, 1, 7)\n    points = collect(Iterators.product(grid, grid, grid))\n    g = PointDiscretizedMap(lorenz_f, points)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now we can choose the domain, which has to contain the fixed point and for this algorithm needs to be subdivided already. We decide to subdivide the domain 24 times, which will give us a regular partition with 2^24 boxes.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    domain = Box(SVector(0.0, 0.0, 27.0), SVector(30.0, 30.0, 40.0))\n    partition = RegularPartition(domain, 24)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We then define the starting point, which is the equilibrium point we mentioned earlier and initialize the target set, which at the moment contains only the small box containing the fixed point but will store all the boxes we collect in the course of the algorithm and will eventually contain the unstable manifold.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    rh = 28.0\n    b = 0.4\n    x0 = (sqrt(b*(rh-1)), sqrt(b*(rh-1)), rh-1)\n\n    boxset = partition[x0]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The only thing left to do now is call the function unstable_set","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    unstable_set!(boxset, g)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"[insert breathtaking plot here]","category":"page"},{"location":"examples/#Transition-Matrix-for-the-Lorenz-System","page":"Examples","title":"Transition Matrix for the Lorenz System","text":"","category":"section"},{"location":"general/#General-usage","page":"General usage","title":"General usage","text":"","category":"section"},{"location":"general/","page":"General usage","title":"General usage","text":"The base of the numerical set oriented methods of this framework are BoxSet and BoxMap, thus, in the following, we will have a closer look at the two and other useful things to know when working with GAIO.","category":"page"},{"location":"general/#Table-of-Contents","page":"General usage","title":"Table of Contents","text":"","category":"section"},{"location":"general/","page":"General usage","title":"General usage","text":"Pages = [\"general.md\"]","category":"page"},{"location":"general/#BoxSet","page":"General usage","title":"BoxSet","text":"","category":"section"},{"location":"general/#BoxMap","page":"General usage","title":"BoxMap","text":"","category":"section"},{"location":"general/#The-Subdivision-Algorithm","page":"General usage","title":"The Subdivision Algorithm","text":"","category":"section"},{"location":"general/#Regular-vs.-tree-partition","page":"General usage","title":"Regular vs. tree partition","text":"","category":"section"},{"location":"general/#Plots","page":"General usage","title":"Plots","text":"","category":"section"},{"location":"#GAIO-Global-Analysis-of-Invariant-Objects","page":"Home","title":"GAIO - Global Analysis of Invariant Objects","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"GAIO is a framework establishing the global numerical analysis of Dynamical Systems based on set oriented techniques.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It allows the set oriented computation and/or visualization of ","category":"page"},{"location":"","page":"Home","title":"Home","text":"invariant sets (e.g. periodic points, global attractor, chain recurrent set) of arbitrary dimension or topology\ninvariant manifolds (stable/unstable manifolds of an arbitrary invariant set)\ninvariant measures, almost invariant sets, coherent sets","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Just ","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"GAIO\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"from the Julia REPL","category":"page"},{"location":"algorithms/#Algorithms-and-Mathematical-Background","page":"Algorithms","title":"Algorithms and Mathematical Background","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"note: Note\nIn the following, f will always refer to the map describing the dynamics of a system, while g will be the corresponding BoxMap.","category":"page"},{"location":"algorithms/#The-Relative-Global-Attractor","page":"Algorithms","title":"The Relative Global Attractor","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Consider a time discrete dynamical system induced by the map f mathbbR^d to mathbbR^d. Let Q subset mathbbR^d compact. The set ","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"A_Q = bigcap_k geq 0 f^k(Q)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"is called the global attractor relative to Q. The relative global attractor can be seen as the set which is eventually approached by every orbit originating in Q. In particular, A_Q contains each invariant set in Q and therefore all the potentially interesting dynamics. Thus it is of great interest to be able to compute a relative global attractor numerically.  The idea of the algorithm is to cover the relative global attractor with boxes and recursively tighten the covering by refining appropriately selected boxes.","category":"page"},{"location":"algorithms/#Mathematical-Background-of-the-Algorithm","page":"Algorithms","title":"Mathematical Background of the Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Mathematically, the algorithm to compute the global attractor relative to Q takes two input arguments: a compact set Q as well as a map f, which describes the dynamics. Now in each iteration, two steps happen:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"subdivision-step: The domain B_k-1 is subdivided once, i.e. every box is bisected along one axis, which gives rise to a new partition of the domain, hatB_k, with double the amount of boxes.\nselection_step: For each box B of the new partition we check, if there is another Box B that is mapped into B under g, i.e. if f(B) cap B neq emptyset. If not, we remove B from the domain.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"After removing every non-hit box, we arrive at the new domain B_k, and as k to infty, the collection of boxes B_k converges to the relative global attractor A_Q.","category":"page"},{"location":"algorithms/#Implementation-of-the-Algorithm","page":"Algorithms","title":"Implementation of the Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"So, how was taken care of the discretization, which is necessary for the implementation? In other words, how did we translate this Algorithm into GAIO?","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"function relative_attractor(boxset::BoxSet, g::BoxMap, depth::Int)\n    for k = 1:depth\n        boxset = subdivide(boxset)\n        boxset = g(boxset; target=boxset)\n    end\n\n    return boxset\nend","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The first thing one notices is that the implementation has a third input parameter depth, which describes the level of approximation. Since in each step of the algorithm the initial domain is divided in half, the final partition after depth many steps will contain n = 2^textdepth boxes, i.e. every box in the final covering is frac1n times the size of the initial box. Besides depth we have the input parameters boxset and g, which is the BoxMap describing the dynamics, that is a function which maps boxes to boxes. g is going to be implemented as a PointDiscretizedMap, which is a struct containing the underlying pointwise map corresponding to the dynamics as well as the set of reference points, implemented as an Array of Tuples, that will be the discretization of a box.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"struct PointDiscretizedMap{F,P} <: BoxMap\n    f::F\n    points::P\nend","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"boxset is struct, which carries both a set and the partition of the set, that is the information about the size of each single box in the set.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"struct BoxSet{P <: BoxPartition,S <: AbstractSet}\n    partition::P\n    set::S\nend","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"For this algorithm the partition needs to be the full regular partition of the initial set (i.e. the starting domain). Now, for each step in the algorithm, the current set is subdivided via the subdivision algorithm:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"boxset = subdivide(boxset)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"And it is checked, if a box is hit by another box under the dynamics:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"function map_boxes_with_target(g, source::BoxSet, target::BoxSet)\n    result = boxset_empty(target.partition)\n\n    for (_, hit) in ParallelBoxIterator(g, source, target.partition)\n        if hit !== nothing # check that point was inside domain\n            if hit in target.set\n                push!(result.set, hit)\n            end\n        end\n    end\n\n    return result\nend","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"where implementationally we start with an empty boxset and store each box that is hit in it.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"","category":"page"},{"location":"algorithms/#unstable_set","page":"Algorithms","title":"unstable_set","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"In the following we are presenting the algorithm to cover invariant manifolds within some domain Q (which has to contain a fixed point).","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"note: Note\nFor simplicity, we will explain the algorithm for the case of the unstable manifold. However one can compute the stable manifold as well by considering the boxmap describing the inverse map f^-1 as input argument for the algorithm.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Usually, the computation of the unstable manifold is relatively simple in 1D, but the higher the dimension, the more complicated it becomes. GAIO is able to compute the unstable manifold for arbitrary dimension.","category":"page"},{"location":"algorithms/#Mathematical-Background-of-the-Algorithm-2","page":"Algorithms","title":"Mathematical Background of the Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The unstable manifold is defined as","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"W^U(x_0) = x lim_k to - infty f^k(x) = x_0 ","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"where x_0 is a fixed point of f.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The idea behind the algorithm to compute the unstable manifold can be explained in two steps. Before starting we need to identify a hyperbolic fixed point and the region Q, which we are going to compute the manifold in. The region Q needs to be already partitioned into small boxes.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Since a fixed point is always part of the unstable manifold, we need to identify a small region/box containing this fixed point.\nThe small box containing the fixed point is then mapped forward under the dynamics defined by f and the boxes that are hit under the image are added to the box collection. Then those newly included boxes are mapped forward and the procedure is repeated. ","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"With these two steps we obtain a covering of part of the global unstable manifold.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"warning: Warning\nOne might not be able to compute the parts of the unstable manifold whose preimage lies outside the domain Q. Thus, it is important to choose Q large enough.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"How was this algorithm translated into the language of GAIO?","category":"page"},{"location":"algorithms/#Implementation-of-the-Algorithm-2","page":"Algorithms","title":"Implementation of the Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"function unstable_set!(boxset::BoxSet, g::BoxMap)\n    boxset_new = boxset\n\n    while !isempty(boxset_new)\n        boxset_new = g(boxset_new)\n\n        setdiff!(boxset_new, boxset)\n        union!(boxset, boxset_new)\n    end\n\n    return boxset\nend","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Let us start with the input arguments for the algorithm: Again, like in the relative attractor Algorithm from above, g is the BoxMap describing the underlying dynamics f. It thus stores f and a set of reference points necessary for the discretization of the boxes. The other input argument boxset includes two things:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"The domain Q we are going to compute the unstable manifold in (Q can be implemented as a large Box) and the underlying partition of the domain. Unlike in the previous algorithm, the domain will not be subdivided along the algorithms course, but we need to pass a partition which is already subdivided to the depth d (and therefore the level of accuracy) we want our final boxcovering to have.\nSince boxset is going to store all the new boxes we aquire in every iteration of the algorithm, it has to be initialized containing no other box than the single box of size frac12^d around the fixed point that is part of the unstable manifold we intend to compute.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Note: This algorithm works with two mutable sets of boxes: boxset, which collects the boxes we aquire in each iteration and will eventually cover part of the unstable manifold, and boxset_new, which will be overwritten in each iteration and contains only the boxes which will be newly added to our collection. To initialize, we set","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"boxset_new = boxset","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Now we repeat the following steps:  First, we map the newly aquired boxes one step forward in time","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"boxset_new = g(boxset_new)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Note: Mapping only the newly acquired boxes from the previous step saves memory and computation time since we already computed the images of the old boxes in previous steps and thus those boximages are already part of the collector boxset. Now we need to update boxset_new and boxset. As mentioned prior, we only want to consider boxes in each iteration step, that have not been 'hit' under g by any boxes we acquired in a previous iteration step, because that would mean that this box image already is part of our box collection. To differ between truly new boxes and boxes we already added, we take the setdifference between the images of boxes boxset_new and the whole boxcollection boxset:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"setdiff!(boxset_new, boxset)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"Now boxset_new contains nothing but the truly new boxes. boxset is then updated by adding those new image boxes to our collection of boxes by forming the union with the already existing collection:","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"union!(boxset, boxset_new)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"We repeat these steps as long as","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"while !isempty(boxset_new)","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"is true. Thus, the iteration will end when no new boxes can be added to the boxcollection, because we e.g. got so close to the border of the domain Q, that every further image of boxes lies beyond the border, or the unstable manifold oscillates so strongly, that our chosen level of accuracy can no longer distinguish between the oscillations.","category":"page"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"","category":"page"},{"location":"algorithms/#transition_matrix","page":"Algorithms","title":"transition_matrix","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"","category":"page"},{"location":"algorithms/#chain*recurrent*set","page":"Algorithms","title":"chainrecurrentset","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"","category":"page"},{"location":"algorithms/#root_covering","page":"Algorithms","title":"root_covering","text":"","category":"section"}]
}

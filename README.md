# balaban-segment-intersections

Automatically exported from code.google.com/p/balaban-segments-intersections

Implementation for segment's intersection algorithms by I.J.Balaban.
Ivan J. Balaban, An optimal algorithm for finding segments intersections,
Proceedings of the eleventh annual symposium on Computational geometry,
p.211-219, June 05-07, 1995, Vancouver, British Columbia,Canada.

Both suboptimal and optimal algorithms are implemented.
Also parallel version of suboptimal algorithm is given.
Suboptimal algorithm called "fast" in the implementation because it is faster than optimal one in most practical cases.
For comparison trivial algorithm, simple sweep line and Bentley \& Ottmann algorithms are presented too.

The aim of the project is to show practical behavior of the algorithms.

Project was compiles under MS VisualStudio 2019. Hopefully it can be compiled in any other C++ enviroment.
More information on [wiki](https://github.com/ivvaan/balaban-segments-intersections/wiki) pages.

**No degenerate case processing is implemented in the master branch. If you need to process them, see the degenerate\_cases branch: CIntegerSegmentCollection class**


The project contains implementations of several segment‑intersection algorithms. The main focus is on the suboptimal and optimal algorithms from Ivan J. Balaban: "An optimal algorithm for finding segments intersections", Proceedings of the 11th Annual Symposium on Computational Geometry, p.211–219, June 05–07, 1995, Vancouver, British Columbia, Canada.

There are two implementations: [Old](https://github.com/ivvaan/balaban-segments-intersections/wiki/Old-implementation) and [New](https://github.com/ivvaan/balaban-segments-intersections/wiki/New-implementation). They differ in the way segments are handled — i.e., they use different segment‑handling functions (SHF).

The project implements the following algorithms:
1. [Trivial](https://github.com/ivvaan/balaban-segments-intersections/wiki/Trivial-algorithm) — Old and New implementations
2. [Simple plane sweep](https://github.com/ivvaan/balaban-segments-intersections/wiki/Simple-plane-sweep) — Old and New implementations
3. [Balaban suboptimal](https://github.com/ivvaan/balaban-segments-intersections/wiki/Suboptimal-algorithm) — Old and New implementations (CHECK OUT THIS FIRST!)
4. [Balaban optimal](https://github.com/ivvaan/balaban-segments-intersections/wiki/Optimal-algorithm) — Old and New implementations  
5. [Balaban suboptimal parallel](https://github.com/ivvaan/balaban-segments-intersections/wiki/Suboptimal-parallel) — Old and New implementations
6. [Bentley & Ottmann sweep line](https://github.com/ivvaan/balaban-segments-intersections/wiki/Bentley-&-Ottmann) — Old implementation only (no corresponding SHF in the New implementation)
7. [Special faster version of Balaban suboptimal for line segments (intersecting pairs only)](https://github.com/ivvaan/balaban-segments-intersections/wiki/Intersecting-pairs-only-algorithm)
8. [Balaban suboptimal with smaller memory usage](https://github.com/ivvaan/balaban-segments-intersections/wiki/Suboptimal-with-smaller-memory-usage) — Old and New implementations

Notes on segments and linearity
- The algorithms operate more efficiently with straight line segments than with non‑linear segments (arcs). A mixed set containing any non‑linear segment must be treated as a non‑linear set (i.e., `SegmentsColl::is_line_segments == false`).
- Supported segment types:
  1. Line segments represented by endpoints
  2. Line segments represented by coefficients `y = a*x + b` with X interval `[x1, x2]`
  3. Arc segments (circular arcs)
  4. Graph (cycle) — edges of a cycle (graph) distribution

Distributions for testing
- [Random](https://github.com/ivvaan/balaban-segments-intersections/wiki/Random-distribution)
- [Long](https://github.com/ivvaan/balaban-segments-intersections/wiki/Long-distribution)
- [Short](https://github.com/ivvaan/balaban-segments-intersections/wiki/Short-distribution)
- [Mixed](https://github.com/ivvaan/balaban-segments-intersections/wiki/Mixed-distribution)
- [Param](https://github.com/ivvaan/balaban-segments-intersections/wiki/Param-distribution)
- [Circle](https://github.com/ivvaan/balaban-segments-intersections/wiki/Circle-distribution)

Distribution notes
- All distributions apply to line segments. For arc segments only Random, Short and Param distributions are applicable.
- Graph (cycle) distribution is supported only in the New implementation.

[Registration (result collection) options](https://github.com/ivvaan/balaban-segments-intersections/wiki/Intersection-registrars)
1. Simple total intersection counting
2. Per‑segment intersection counting
3. Storing intersecting pairs and intersection points (vector storage) — be careful with memory, K can be O(N²)

When using per‑segment registration you can obtain either:
- the total number of intersections, or
- the maximum number of intersections found for any single segment  
(Internally the program also keeps the full per‑segment counts when this registrator is selected.)

Implementation and API pointers
- Test harness and collection creation: `test_coll.cpp` (`create_test_collection`, `delete_test_collection`, `find_intersections` glue)
- Old implementation (function pointers SHF): `intersection_finder.h` / `intersection_finder.cpp` and `test_coll.cpp` (class `SegmentFunctions<T>` adapters)
- New implementation (template SHF / SegmentCollection interface): `new_int_finder.h`, collection classes `lines1_collection.h`, `lines2_collection.h`, `arcs_collection.h`, `graph_collection.h`
- Algorithms (new implementation): `new_int_finder.h`, `fast_finder.h`, `optimal_finder.h`
- Low-level helpers and contracts:
  - [SetCurSeg-functions](https://github.com/ivvaan/balaban-segments-intersections/wiki/SetCurSeg-functions) — how `SetCurSeg*` variants must clip and prepare segments for intersection routines
  - [Sentinels and guards](https://github.com/ivvaan/balaban-segments-intersections/wiki/Sentinels-and-guards) — sentinel/guard optimization contract and examples
  - [Cut dispatching functions](https://github.com/ivvaan/balaban-segments-intersections/wiki/Cut-dispatching-functions) — strip splitting control flow and dispatch logic
  - [Intersection registrars](https://github.com/ivvaan/balaban-segments-intersections/wiki/Intersection-registrars) — registration system architecture, `reg_type` dispatch, and available registrar implementations
- Visualization / SVG output: `test_coll.cpp` (`write_SVG`) and `seg_int.cpp` (HTML template & SVG insertion)

Reference
- Ivan J. Balaban, "An optimal algorithm for finding segments intersections", Proc. 11th ACM Symposium on Computational Geometry, 1995. Available from ACM Digital Library.




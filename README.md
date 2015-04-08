# balaban-segments-intersections
Automatically exported from code.google.com/p/balaban-segments-intersections

Implementation for segment's intersection algorithms by I.J.Balaban.
Ivan J. Balaban, An optimal algorithm for finding segments intersections, 
Proceedings of the eleventh annual symposium on Computational geometry, 
p.211-219, June 05-07, 1995, Vancouver, British Columbia,Canada.

Both suboptimal and optimal algorithms are implemented.
Also parallel version of suboptimal algorithm is given. 
Suboptimal algorithm called "fast" in the implementation because it is faster than optimal one in most practical cases. 
For comparison trivial algorithm, simple sweep line and Bentley & Ottmann algorithms are presented too.

The aim of the project is to show practical behavior of the algorithms. 
No any degenerate cases processing is implemented. 
So only segments in common situation are processed correctly. 
For example all segment ends must have different x coordinates.

Project was compiles under MS VisualStudio 2008-2010. Hopefully it can be compiled in any other C++ enviroment.

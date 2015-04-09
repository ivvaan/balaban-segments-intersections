# balaban-segments-intersections
Contains robust implementation of the suboptimal algorithm. 
The implementation uses only integers and performs "logical" calculations.

All float point operations moved out from the implementation to CRobustSegmentCollection. 
So if you implement the CRobustSegmentCollection functions robust way the whole procedure would be robust. 

Function balaban_fast works with logical N segments represented by integers 0,1..N-1. 
The ends for the segment s represented by integers 2s (begin point) and 2s+1 (end point). 
The function assumes that all ends logically different and begin point lies "on the left" from end point.
Ends also represents logical x-coordinates from 0 to 2N-1.

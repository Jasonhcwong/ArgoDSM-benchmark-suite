Single Source Shortest Path
===========================

You can run ```make``` to create executables for each program or use the following commands below

**Input Graph from a File**

To run with P number of threads and from source node S:
  ```./sssp P <input_file> <output_file> S```
  
  It will then ask for the input file, enter:
  sample.txt
  OR any other file such as road networks from the SNAP datasets (e.g. roadNet-CA)
  https://snap.stanford.edu/data/#road

**Notes**

This version of sssp parallelizes the O(V^2) version of Dijkstra's Algorithm, given by Yen's Optimization.
Paper: J.Y.Yen, "An algorithm for finding shortest routes from all source nodes to a given destination in general networks", Quarterly of Applied Mathematics 01/1970.

SSSP has a parameter P_max that is specified by ```cuberoot(N)*3```, which represents the number of iterations for the outer loop.

The executable then outputs the time in seconds that the program took to run.

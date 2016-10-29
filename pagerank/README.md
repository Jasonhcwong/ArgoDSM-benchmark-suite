ArgoDSM PageRank Benchmark
==========================

Run ```make``` to generate executables, then use the syntax explained below

The first argument to the executable specifies whether you want to read the graph from a file (1), or generate a synthetic one internally (0).

**Input Graph from File**

To run with P number of threads over all nodes, and an input file,
    ```./pagerank P <input_file> <output_file>```

**Notes**

The executable then outputs the time in seconds that the program took to run.
It also outputs a file that contains the pageranks (normalized to 1).

Some very small differences in pageranks might occur due to floating point round offs within the program.

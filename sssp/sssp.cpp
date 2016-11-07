/*
 *
    Distributed Under the MIT license
    Uses the Range based of Bellman-Ford/Dijkstra Algorithm to find shortest path distances
    For more info, see Yen's Optimization for the Bellman-Ford Algorithm, or the Pannotia Benchmark Suite.
    Programs by Masab Ahmad (UConn)
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <pthread.h>
#include <string>

#include "argo.hpp"


// Single Source Shortest Path Constants
int    MAX_VERTICES     = 2000000;
int    MAX_DEGREE       = 32;
int    MAX_THREADS      = 256;
int    MAX_DISTANCE     = 100000000;


struct Thread_args {
   int       global_thread_id;
   int*      distances;
   int       vertices;
};


// node local variables
int global_num_threads;
int local_num_threads;
argo::globallock::global_tas_lock *lock; // single lock
std::vector<argo::globallock::global_tas_lock*>* locks; //[2097152]; //change the number of locks to approx or greater N


// ArgoDSM global variables
Thread_args* arguments;
bool*        lock_flag;
bool*        lock_flags;
int*         global_size;
int*         distances;
bool*        exists;
int*         weights;
int*         graph;
bool*        terminate;
int*         global_range;       // starting range set to 1



void* do_work(void* argptr) {

    Thread_args* args = static_cast<Thread_args*>(argptr);

    int  global_thread_id = args->global_thread_id;
    int  vertices         = args->vertices;

    int count = 0;
    int start = 0;
    int stop = 1;
    int neighbor = 0;
    int v = 0;

    argo::barrier(local_num_threads);

    while (!(*terminate)) {
        while (!(*terminate)) {
            for(v = start; v < stop; v++) {
                if (!exists[v])
                    continue;
                for (int i = 0; i < MAX_DEGREE; i++) {
                    if(v < vertices) neighbor = graph[v*MAX_DEGREE+i];
                    if(neighbor >= vertices) break;

                    locks->at(neighbor)->lock();

                    if (distances[graph[v*MAX_DEGREE+i]] > (distances[v] + weights[v*MAX_DEGREE+i]))    //relax, update distance
                        distances[graph[v*MAX_DEGREE+i]] = distances[v] + weights[v*MAX_DEGREE+i];

                    locks->at(neighbor)->unlock();
                }
            }

            argo::barrier(local_num_threads);

            if(global_thread_id == 0) {
                int range = (*global_range) * MAX_DEGREE; //change this for range heuristic e.g. range = range+DEG;
                if(range >= vertices)
                    range = vertices;
                *global_range = range;
            }

            argo::barrier(local_num_threads);

            int range = *global_range;
            start = (int) ((range / global_num_threads) * global_thread_id);
            stop  = (int) ((range / global_num_threads) * (global_thread_id+1));

            if(stop > range) stop = range;

            if(start == vertices || v > vertices-1) *terminate = true;

            argo::barrier(local_num_threads);
        }

        argo::barrier(local_num_threads);

        if(global_thread_id == 0) {
            count++;
            if(count < MAX_THREADS) {
                *terminate = false;
                *global_range = 1;
            }
        }

        start = 0;
        stop = 1;

        argo::barrier(local_num_threads);
    }
    return NULL;
}



void read_graph_from_file(std::string filename) {

    std::ifstream input(filename);

    int int_max = std::numeric_limits<int>::max();

    for(int j = 0; j < MAX_VERTICES; j++) {
        for(int i = 0; i < MAX_DEGREE; i++) {
            weights[j*MAX_DEGREE+i] = int_max;
            graph[j*MAX_DEGREE+i] = int_max;
        }
        exists[j] = false;
    }

    int number0, number1, number2;
    int previous_node = -1;
    int inter = -1;
    int vertex_cnt = 0;

    std::string line;
    while (input >> number0 >> number1 >> number2) {
        if (number0 >= MAX_VERTICES) {
            std::cout << "Node " << number0 << " exceeds maximum graph size of " << MAX_VERTICES << "\n";
            exit (EXIT_FAILURE);
        }
        exists[number0] = true;
        exists[number1] = true;
        if(number0 == previous_node)
           inter++;
        else inter = 0;

        if (inter >= MAX_DEGREE) {
            std::cout << "Node " << number0 << " exceeds maximum maximum degree of " << MAX_DEGREE << "\n";
            exit (EXIT_FAILURE);
        }

        bool is_existing = false;
        for (int i = 0; i < inter; ++i) {
           if (graph[number0*MAX_DEGREE+i] == number1) {
              is_existing = true;
              break;
           }
        }

        if (!is_existing) {
           graph[number0*MAX_DEGREE+inter] = number1;
           weights[number0*MAX_DEGREE+inter] = number2;
           previous_node = number0;
        }
        if(number0 > vertex_cnt) vertex_cnt = number0;
        if(number1 > vertex_cnt) vertex_cnt = number1;
    }

    *global_size = ++vertex_cnt;

    input.close();
}



void write_sssp_to_file(std::string filename, int* distances, int size) {
    std::ofstream output(filename);
    output << std::fixed;
    output << std::setprecision(6);
    output << "distances:\n";
    for(int i = 0; i < size; i++)
        if(distances[i] < MAX_DISTANCE)
            output << "distance(" << i << ") = " << distances[i] << "\n";
    output.close();
}



int main(int argc, char** argv) {

    argo::init(0.5*1024*1024*1024UL);

    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " <thread-count> <input-file> <output-file>" << std::endl;
        return 1;
    }

    global_num_threads = atoi(argv[1]);
    std::string input_filename = argv[2];
    std::string output_filename = argv[3];
    int source_vertex = atoi(argv[4]);

    if (!global_num_threads) {
        std::cout << "Thread count must be a valid integer greater than 0." << std::endl;
        printf ("Error:  ");
        return 1;
    }

    Thread_args*  arguments    = argo::conew_array<Thread_args>(global_num_threads);
    terminate    = argo::conew_<bool>(false);
    global_size  = argo::conew_<int>(0);
    global_range = argo::conew_<int>(1);
    lock_flag    = argo::conew_<bool>(false);
    lock_flags   = argo::conew_array<bool>(MAX_VERTICES);
    distances    = argo::conew_array<int>(MAX_VERTICES);
    exists       = argo::conew_array<bool>(MAX_VERTICES);
    weights      = argo::conew_array<int>(MAX_VERTICES*MAX_DEGREE);
    graph        = argo::conew_array<int>(MAX_VERTICES*MAX_DEGREE);

    argo::barrier();

    lock  = new argo::globallock::global_tas_lock(lock_flag);

    if (argo::node_id() == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    int vertices = *global_size;

    locks = new std::vector<argo::globallock::global_tas_lock*>(vertices);

    for (int i = 0; i < vertices; i++)
        locks->at(i) = new argo::globallock::global_tas_lock(&lock_flags[i]);

    // Initialize ArgoDSM global variables
    if (argo::node_id() == 0) {
        for(int i = 0; i < MAX_VERTICES; i++) {
            lock_flags[i] = false;
            distances[i] = MAX_DISTANCE;
        }
        distances[source_vertex] = 0;

        for(int j = 0; j < global_num_threads; j++)
            arguments[j].global_thread_id = j;
    }

    argo::barrier();

    // Divide the work as equal as possible among global_num_threads
    local_num_threads = global_num_threads / argo::number_of_nodes();

    int node_threads_begin = argo::node_id() * local_num_threads;
    std::vector<pthread_t> threads(local_num_threads);

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].vertices = vertices;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    argo::barrier();

    if (argo::node_id() == 0)
        write_sssp_to_file(output_filename, distances, vertices);

    delete lock;
    delete locks;

    argo::codelete_(terminate);
    argo::codelete_(global_range);
    argo::codelete_(global_size);
    argo::codelete_array(lock_flags);
    argo::codelete_array(distances);
    argo::codelete_array(exists);
    argo::codelete_array(weights);
    argo::codelete_array(graph);

    argo::finalize();

    return 0;
}

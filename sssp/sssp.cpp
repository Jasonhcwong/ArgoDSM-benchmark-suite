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
#include "synchronization/cohort_lock.cpp"



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
int node_id;
int number_of_nodes;
int global_num_threads;
int local_num_threads;
argo::globallock::cohort_lock *lock; // single lock
std::vector<argo::globallock::cohort_lock*>* locks; //[2097152]; //change the number of locks to approx or greater N


// ArgoDSM global variables
Thread_args* arguments;
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

                    if (distances[neighbor] > (distances[v] + weights[v*MAX_DEGREE+i]))    //relax, update distance
                        distances[neighbor] = distances[v] + weights[v*MAX_DEGREE+i];

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
            double tid_d = global_thread_id;
            double P_d = global_num_threads;
            double start_d = (range/P_d) * tid_d;
            double stop_d = (range/P_d) * (tid_d+1.0);
            start = start_d;
            stop = stop_d;

            if(stop > range) stop = range;

            //std::cout << "Range: " << range << " > Thread " << global_thread_id << ": " <<  start << " - " << stop << "\n";

            if(start == vertices || v >= vertices) {
                lock->lock();
                *terminate = true;
                lock->unlock();
            }

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

    std::vector<bool>* exists_tmp = new std::vector<bool>(MAX_VERTICES);
    std::vector<int>* weights_tmp = new std::vector<int>(MAX_VERTICES * MAX_DEGREE);
    std::vector<int>* graph_tmp = new std::vector<int>(MAX_VERTICES * MAX_DEGREE);

    std::vector<int>* degree = new std::vector<int>(MAX_VERTICES);

    std::ifstream input(filename);

    int int_max = std::numeric_limits<int>::max();

    for(int j = 0; j < MAX_VERTICES; j++) {
        for(int i = 0; i < MAX_DEGREE; i++) {
            weights_tmp->at(j*MAX_DEGREE+i) = int_max;
            graph_tmp->at(j*MAX_DEGREE+i) = int_max;
        }
        degree->at(j) = 0;
        exists_tmp->at(j) = false;
    }

    int number0, number1, number2;
    int vertex_cnt = 0;

    std::string line;
    while (input >> number0 >> number1 >> number2) {

        if (number0 >= MAX_VERTICES) {
            std::cout << "Node " << number0 << " exceeds maximum graph size of " << MAX_VERTICES << "\n";
            exit (EXIT_FAILURE);
        }

        if (degree->at(number0) >= MAX_DEGREE) {
            std::cout << "Node " << number0 << " exceeds maximum maximum degree of " << MAX_DEGREE << "\n";
            exit (EXIT_FAILURE);
        }

        bool is_existing = false;
        for (int i = 0; i < degree->at(number0); ++i) {
           if (graph_tmp->at(number0*MAX_DEGREE+i) == number1) {
              is_existing = true;
              break;
           }
        }

        exists_tmp->at(number0) = true;
        exists_tmp->at(number1) = true;

        if (!is_existing) {
           graph_tmp->at(number0*MAX_DEGREE+degree->at(number0)) = number1;
           weights_tmp->at(number0*MAX_DEGREE+degree->at(number0)) = number2;
           degree->at(number0)++;
        }

        if(number0 > vertex_cnt) vertex_cnt = number0;
        if(number1 > vertex_cnt) vertex_cnt = number1;
    }

    *global_size = ++vertex_cnt;

    for (int j = 0; j < vertex_cnt; j++) {
	exists[j] = exists_tmp->at(j);
	for (int i = 0; i < degree->at(j); i++) {
	    weights[j] = weights_tmp->at(j*MAX_DEGREE+i);
	    graph[j] = graph_tmp->at(j*MAX_DEGREE+i);
	}
    }

    delete degree;
    delete exists_tmp;
    delete weights_tmp;
    delete graph_tmp;

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

    auto start = std::chrono::system_clock::now();

    argo::init(0.5*1024*1024*1024UL);

    if (node_id == 0)
    	std::cout << "argo::init: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() / 1000 << " s" << std::endl;

    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " <thread-count> <input-file> <output-file> <source-vertex>" << std::endl;
        return 1;
    }

    node_id = argo::node_id();
    number_of_nodes = argo::number_of_nodes();

    local_num_threads = atoi(argv[1]);
    global_num_threads = local_num_threads * number_of_nodes;
    std::string input_filename = argv[2];
    std::string output_filename = argv[3];
    int source_vertex = atoi(argv[4]);

    if (!global_num_threads) {
        std::cout << "Thread count must be a valid integer greater than 0." << std::endl;
        printf ("Error:  ");
        return 1;
    }

    start = std::chrono::system_clock::now();

    arguments    = argo::conew_array<Thread_args>(global_num_threads);
    terminate    = argo::conew_<bool>(false);
    global_size  = argo::conew_<int>(0);
    global_range = argo::conew_<int>(1);
    distances    = argo::conew_array<int>(MAX_VERTICES);
    exists       = argo::conew_array<bool>(MAX_VERTICES);
    weights      = argo::conew_array<int>(MAX_VERTICES*MAX_DEGREE);
    graph        = argo::conew_array<int>(MAX_VERTICES*MAX_DEGREE);

    argo::barrier();

    if (node_id == 0)
    	std::cout << "Global declarations: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() / 1000 << " s" << std::endl;

    start = std::chrono::system_clock::now();

    lock  = new argo::globallock::cohort_lock();

    if (node_id == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    if (node_id == 0)
    	std::cout << "Reading graph: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() / 1000 << " s" << std::endl;

    start = std::chrono::system_clock::now();

    int vertices = *global_size;

    locks = new std::vector<argo::globallock::cohort_lock*>(vertices);

    for (int i = 0; i < vertices; i++)
        locks->at(i) = new argo::globallock::cohort_lock();

    // Initialize ArgoDSM global variables
    if (node_id == 0) {
        for(int i = 0; i < MAX_VERTICES; i++)
            distances[i] = MAX_DISTANCE;
        distances[source_vertex] = 0;

        for(int j = 0; j < global_num_threads; j++)
            arguments[j].global_thread_id = j;
    }

    argo::barrier();

    if (node_id == 0)
    	std::cout << "Global initializations: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() / 1000 << " s" << std::endl;

    // Divide the work as equal as possible among global_num_threads

    int node_threads_begin = node_id * local_num_threads;
    std::vector<pthread_t> threads(local_num_threads);

    start = std::chrono::system_clock::now();

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].vertices = vertices;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    argo::barrier();

    if (node_id == 0)
    	std::cout << "Parallel runtime: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() / 1000 << " s" << std::endl;

    if (node_id == 0) {
        write_sssp_to_file(output_filename, distances, vertices);
        std::cout << "\nSingle Source Shortest Path";
        std::cout << "\nArgo nodes: " << number_of_nodes;
        std::cout << "\nGlobal threads: " << global_num_threads;
        std::cout << "\nLocal threads: " << local_num_threads;
	std::cout << "\nGraph: " << input_filename;
	std::cout << "\nVertices: " << vertices << "\n";
    }

    for (int i = 0; i < vertices; i++)
        delete locks->at(i);

    delete lock;
    delete locks;

    argo::codelete_(terminate);
    argo::codelete_(global_range);
    argo::codelete_(global_size);
    argo::codelete_array(distances);
    argo::codelete_array(exists);
    argo::codelete_array(weights);
    argo::codelete_array(graph);

    argo::finalize();

    return 0;
}

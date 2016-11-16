#include <cassert>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <pthread.h>
#include <string>

#include "argo.hpp"

// Pagerank Constants
int    MAX_VERTICES     = 2000000;
int    MAX_DEGREE       = 32;
double DAMPING_FACTOR   = 0.85;
int    ITERATIONS       = 40;


struct Thread_args {
    int global_thread_id;
    int local_thread_id;
    int data_begin;
    int data_end;
    int vertices;
    pthread_barrier_t* barrier;
};

struct Chunk {
    int begin;
    int size;
    std::vector<int> inlinks;
    std::vector<bool> exists;
    std::vector<int> outlinks;
    std::vector<double> pagerank;
};


// ArgoDSM node local variables
int node_id;
int number_of_nodes;
std::vector<double>* dp_local;
double dp_node;
int global_num_threads;
int local_num_threads;
Chunk* chunk;


// ArgoDMS global variables
int*         global_size;
Thread_args* arguments;
double*      dp;
double*	     dp_global;
double*      pagerank;
int*         inlinks;
bool*        exists;
int*         outlinks;
int*         graph;


void* do_work(void* argptr) {
    Thread_args* args = static_cast<Thread_args*>(argptr);

    int global_thread_id        = args->global_thread_id;
    int local_thread_id         = args->local_thread_id;
    int data_begin              = args->data_begin;
    int data_end                = args->data_end;
    int vertices                = args->vertices;
    pthread_barrier_t* barrier  = args->barrier;

    int iterations              = ITERATIONS;

    argo::barrier(local_num_threads);

    while(iterations > 0) {

        dp_local->at(local_thread_id) = 0;
        for(int i = data_begin; i < data_end; i++)
            if(chunk->outlinks[i] == 0)
                dp_local->at(local_thread_id) += DAMPING_FACTOR*(pagerank[chunk->begin+i]/vertices);

        pthread_barrier_wait(barrier);

        if (local_thread_id == 0) {
	    double sum = 0;
	    for (int i = 0; i < local_num_threads; i++)
            	sum += dp_local->at(i); 
	    dp_global[node_id] = sum;
        }

        argo::barrier(local_num_threads);

        if (global_thread_id == 0) {
	    double sum = 0;
	    for (int i = 0; i < number_of_nodes; i++)
            	sum += dp_global[i];
	    *dp = sum;
        }

        argo::barrier(local_num_threads);

        for(int i = data_begin; i < data_end; i++) {
            if (chunk->exists[i]) {
                chunk->pagerank[i] = ((1-DAMPING_FACTOR)/vertices) + (*dp);
                for(int j = 0; j < chunk->inlinks[i]; j++)
                    chunk->pagerank[i] += DAMPING_FACTOR * (pagerank[graph[(chunk->begin+i)*MAX_DEGREE+j]] / outlinks[graph[(chunk->begin+i)*MAX_DEGREE+j]]);
            }
        }

        argo::barrier(local_num_threads);

        for(int i = data_begin; i < data_end; i++)
            if(chunk->exists[i])
                pagerank[chunk->begin+i] = chunk->pagerank[i];

        argo::barrier(local_num_threads);

        iterations--;
	std::cout << "Iteration: " << ITERATIONS - iterations << "\n";
    }
    return NULL;
}



void read_graph_from_file(std::string filename) {

    std::ifstream input(filename);

    for(int i = 0; i < MAX_VERTICES; i++) {
        inlinks[i]  = 0;
        exists[i]   = false;
        outlinks[i] = 0;
        for(int j = 0; j < MAX_DEGREE; j++)
            graph[i*MAX_DEGREE+j] = std::numeric_limits<int>::max();
    }

    int number0, number1, inter;
    int vertex_cnt = 0;

    std::string line;
    while (input >> number0 >> number1) {
        inter = inlinks[number1];
        graph[number1*MAX_DEGREE+inter] = number0;
        inlinks[number1]++;
        outlinks[number0]++;
        exists[number0] = true;
        exists[number1] = true;
        if(number0 > vertex_cnt) vertex_cnt = number0;
        if(number1 > vertex_cnt) vertex_cnt = number1;
    }
    *global_size = ++vertex_cnt;

    input.close();
}



void write_pagerank_to_file(std::string filename, double* pagerank, int size) {
    std::ofstream output(filename);
    //output << std::fixed;
    //output << std::setprecision(6);
    for(int i = 0; i < size; i++) {
        if(exists[i])
            output << "pr(" << i << ") = " << pagerank[i] << "\n";
    }
    output.close();
}



int main(int argc, char* argv[]) {

    // Startup ArgoDSM
    argo::init(0.5*1024*1024*1024UL);

    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <thread-count> <input-file> <output-file>" << std::endl;
        return 1;
    }

    node_id = argo::node_id();
    number_of_nodes = argo::number_of_nodes();
    global_num_threads = atoi(argv[1]);
    local_num_threads = global_num_threads / number_of_nodes;

    std::string input_filename = argv[2];
    std::string output_filename = argv[3];

    if (!global_num_threads) {
        std::cout << "Thread count must be a valid integer greater than 0." << std::endl;
        return 1;
    }

    // Declare ArgoDSM global memory variables
    arguments   = argo::conew_array<Thread_args>(global_num_threads);
    dp          = argo::conew_<double>(0);
    dp_local    = new std::vector<double>(local_num_threads);
    dp_global   = argo::conew_array<double>(number_of_nodes);

    global_size = argo::conew_<int>(MAX_VERTICES);
    pagerank    = argo::conew_array<double>(MAX_VERTICES);
    inlinks     = argo::conew_array<int>(MAX_VERTICES);
    exists      = argo::conew_array<bool>(MAX_VERTICES);
    outlinks    = argo::conew_array<int>(MAX_VERTICES);
    graph       = argo::conew_array<int>(MAX_VERTICES*MAX_DEGREE);

    argo::barrier();

    // Read in graph from file by one node
    if (node_id == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    int vertices = *global_size;

    // Initialize ArgoDSM global memory variables
    if (node_id == 0)
        for(int i = 0; i < vertices; i++)
            pagerank[i] = 1.0/vertices; //  1.0-DAMPING_FACTOR;

    argo::barrier();

    // Divide the work as equal as possible among global_num_threads
    if (node_id == 0) {
        int node_chunk = 0;
        int thread_chunk = vertices / global_num_threads;
        int rem = vertices % global_num_threads;
        for (int i = 0; i < global_num_threads; i++) {
            if (i % local_num_threads == 0) node_chunk = 0;
            arguments[i].global_thread_id = i;
            if (rem != 0 && i < rem) {
                arguments[i].data_begin = node_chunk;
                arguments[i].data_end = arguments[i].data_begin + (thread_chunk+1);
                node_chunk += (thread_chunk+1);
            } else {
                arguments[i].data_begin = node_chunk;
                arguments[i].data_end = arguments[i].data_begin + thread_chunk;
                node_chunk += thread_chunk;
            }
        }
    }

    argo::barrier();

    // Copy relevant graph chunk to local memory of the argo node
    chunk = new Chunk();
    int node_threads_begin = node_id * local_num_threads;
    chunk->size = arguments[node_threads_begin + local_num_threads - 1].data_end;
    chunk->begin = 0;
    for (int i = 0; i < node_threads_begin; i++)
        chunk->begin += arguments[i].data_end - arguments[i].data_begin;

    chunk->pagerank = std::vector<double>(chunk->size);
    chunk->inlinks     = std::vector<int>(chunk->size);
    chunk->exists   = std::vector<bool>(chunk->size);
    chunk->outlinks = std::vector<int>(chunk->size);

    for (int j = 0; j < chunk->size; j++) {
        chunk->inlinks[j]     = inlinks[chunk->begin+j];
        chunk->exists[j]   = exists[chunk->begin+j];
        chunk->outlinks[j] = outlinks[chunk->begin+j];
    }

    std::vector<pthread_t> threads(local_num_threads);
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, local_num_threads);

    auto start = std::chrono::system_clock::now();

    // Startup threads for computation
    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].local_thread_id = i;
        arguments[j].vertices = vertices;
        arguments[j].barrier = &barrier;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    // Join threads
    for (auto &t : threads)
        pthread_join(t, nullptr);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    argo::barrier();

    if (node_id == 0) {
        write_pagerank_to_file(output_filename, pagerank, vertices);
        std::cout << "\nPagerank\n";
        std::cout << "Argo nodes: " << number_of_nodes;
        std::cout << "\nGlobal threads: " << global_num_threads;
        std::cout << "\nLocal threads: " << local_num_threads;
        std::cout << "\nGraph: " << input_filename;
        std::cout << "\nVertices: " << vertices;
        std::cout << "\nIterations: " << ITERATIONS << "\n";
        std::cout << "Runtime: " << elapsed.count() << " ms\n" << std::endl;
    }

    delete dp_local;
    delete chunk;

    argo::codelete_array(dp_global);
    argo::codelete_array(inlinks);
    argo::codelete_array(inlinks);
    argo::codelete_array(exists);
    argo::codelete_array(outlinks);
    argo::codelete_array(graph);
    argo::codelete_array(pagerank);
    argo::codelete_(dp);
    argo::codelete_array(arguments);

    argo::finalize();
}

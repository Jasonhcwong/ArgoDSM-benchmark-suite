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
int    MAX_DEGREE       = 16;
double DAMPING_FACTOR   = 0.85;
int    ITERATIONS       = 40;


struct Thread_args {
    int global_thread_id;
    int data_begin;
    int data_end;
    int vertices;
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
int global_num_threads;
int local_num_threads;
argo::globallock::global_tas_lock *lock; // single lock
Chunk* chunk;


// ArgoDMS global variables
int*         global_size;
Thread_args* arguments;
bool*        lock_flag;
double*      dp;
double*      dp_threads;
double*      pagerank;
int*         inlinks;
bool*        exists;
int*         outlinks;
int*         graph;


void* do_work(void* argptr) {
    Thread_args* args = static_cast<Thread_args*>(argptr);

    int global_thread_id = args->global_thread_id;
    int data_begin       = args->data_begin;
    int data_end         = args->data_end;
    int vertices         = args->vertices;
    int iterations       = ITERATIONS;

    argo::barrier(local_num_threads);

    while(iterations > 0) {
        // think about this
        if(global_thread_id == 0)
            *dp = 0;

        argo::barrier(local_num_threads);

        dp_threads[global_thread_id] = 0;
        for(int i = data_begin; i < data_end; i++)
            if(chunk->outlinks[i] == 0)
                dp_threads[global_thread_id] += DAMPING_FACTOR*(pagerank[chunk->begin+i]/vertices);

        lock->lock();
        *dp = *dp + dp_threads[global_thread_id];
        lock->unlock();

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
        std::cout << "Wrong number of arguments" << std::endl;
        return 0;
    }

    global_num_threads = atoi(argv[1]);
    std::string input_filename = argv[2];
    std::string output_filename = argv[3];

    // Declare ArgoDSM global memory variables
    arguments   = argo::conew_array<Thread_args>(global_num_threads);
    dp          = argo::conew_<double>(0);
    dp_threads  = argo::conew_array<double>(global_num_threads);
    lock_flag   = argo::conew_<bool>(false);
    lock        = new argo::globallock::global_tas_lock(lock_flag);

    global_size = argo::conew_<int>(MAX_VERTICES);
    pagerank    = argo::conew_array<double>(MAX_VERTICES);
    inlinks     = argo::conew_array<int>(MAX_VERTICES);
    exists      = argo::conew_array<bool>(MAX_VERTICES);
    outlinks    = argo::conew_array<int>(MAX_VERTICES);
    graph       = argo::conew_array<int>(MAX_VERTICES*MAX_DEGREE);

    argo::barrier();

    // Read in graph from file by hopefully only one node
    if (argo::node_id() == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    int vertices = *global_size;

    // Initialize ArgoDSM global memory variables
    if (argo::node_id() == 0)
        for(int i = 0; i < vertices; i++)
            pagerank[i] = 1.0/vertices; //  1.0-DAMPING_FACTOR;

    argo::barrier();

    // Divide the work as equal as possible among global_num_threads
    local_num_threads = global_num_threads / argo::number_of_nodes();

    if (argo::node_id() == 0) {
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
    int node_threads_begin = argo::node_id() * local_num_threads;
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

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].vertices = vertices;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    argo::barrier();

    if (argo::node_id() == 0)
        write_pagerank_to_file(output_filename, pagerank, vertices);

    delete chunk;
    delete lock;

    argo::codelete_array(inlinks);
    argo::codelete_array(exists);
    argo::codelete_array(outlinks);
    argo::codelete_array(graph);
    argo::codelete_array(pagerank);
    argo::codelete_(lock_flag);
    argo::codelete_array(dp_threads);
    argo::codelete_(dp);
    argo::codelete_array(arguments);

    argo::finalize();
}

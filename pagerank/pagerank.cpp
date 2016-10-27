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

struct Thread_args {
    int global_thread_id;
    int data_begin;
    int data_end;
    int vertices;
};

struct Chunk {
    int begin;
    int size;
    std::vector<int> test;
    std::vector<bool> exists;
    std::vector<bool> dangling;
    std::vector<int> outlinks;
    std::vector<std::vector<int>> W_index;
    std::vector<double> pagerank;
};


// ArgoDSM node local variables
int global_num_threads;
int local_num_threads;
argo::globallock::global_tas_lock *lock; // single lock
Chunk* chunk;

// Pagerank Constants
int    MAX_VERTICES     = 2000000;
int    MAX_DEGREE       = 16;
double INITIAL_PAGERANK = 0.15;
double DAMPING_FACTOR   = 0.85;
int    ITERATIONS       = 1;


// ArgoDMS global variables
int*         io_node_id;
int*         global_size;
Thread_args* arguments;
bool*        lock_flag;
double*      dp;
double*      dp_threads;
double*      pagerank;
int*         test;
bool*        exists;
bool*        dangling;
int*         outlinks;
int**        W_index;



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

        for(int i = data_begin; i < data_end; i++)
            if(chunk->dangling[i])
                dp_threads[global_thread_id] += DAMPING_FACTOR*(pagerank[chunk->begin+i]/vertices);

        lock->lock();
        *dp = *dp + dp_threads[global_thread_id];
        lock->unlock();

        argo::barrier(local_num_threads);

        for(int i = data_begin; i < data_end; i++) {
            if (chunk->exists[i]) {
                chunk->pagerank[i] = 0.15;
                for(int j = 0; j < chunk->test[i]; j++)
                    chunk->pagerank[i] += DAMPING_FACTOR * (pagerank[chunk->W_index[chunk->begin+i][j]] / outlinks[chunk->W_index[chunk->begin+i][j]]);
            }
            if (chunk->pagerank[i] >= 1.0)
                chunk->pagerank[i] = 1.0;
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
    if (!input.good()) {
        input.close();
        return;
    }
    lock->lock();
    *io_node_id = argo::node_id();

    for(int i = 0; i < MAX_VERTICES; i++) {
        test[i]     = 0;
        exists[i]   = false;
        dangling[i] = false;
        outlinks[i] = 0;
        for(int j = 0; j < MAX_DEGREE; j++)
            W_index[i][j] = std::numeric_limits<int>::max();
    }

    int number0, number1, inter;
    int vertex_cnt = 0;

    std::string line;
    while (input >> number0 >> number1) {
        inter = test[number1];
        W_index[number1][inter] = number0;
        test[number1]++;
        outlinks[number0]++;
        exists[number0] = true;
        exists[number1] = true;
        dangling[number1] = true;
        if(number0 > vertex_cnt) vertex_cnt = number0;
        if(number1 > vertex_cnt) vertex_cnt = number1;
    }
    for(int i = 0; i < MAX_VERTICES; i++)
        if(exists[i] && dangling[i])
            dangling[i] = false;
    *global_size = ++vertex_cnt;

    lock->unlock();
    input.close();
}



void write_pagerank_to_file(std::string filename, double* pagerank, int size) {
    std::ofstream output(filename);
    output << std::fixed;
    output << std::setprecision(6);
    for(int i = 0; i < size; i++) {
        if(exists[i])
            output << "pr(" << i << ") = " << pagerank[i] << "\n";
    }
    output.close();
}



int main(int argc, char* argv[]) {

    // Startup ArgoDSM
    argo::init(0.5*1024*1024*1024UL);

    global_num_threads = atoi(argv[1]);
    std::string filename = argv[2];

    // Declare ArgoDSM global memory variables
    arguments   = argo::conew_array<Thread_args>(global_num_threads);
    dp          = argo::conew_<double>(0);
    dp_threads  = argo::conew_array<double>(global_num_threads);
    lock_flag   = argo::conew_<bool>(false);
    lock        = new argo::globallock::global_tas_lock(lock_flag);

    global_size = argo::conew_<int>(MAX_VERTICES);
    io_node_id  = argo::conew_<int>(0);
    pagerank    = argo::conew_array<double>(MAX_VERTICES);
    test        = argo::conew_array<int>(MAX_VERTICES);
    exists      = argo::conew_array<bool>(MAX_VERTICES);
    dangling    = argo::conew_array<bool>(MAX_VERTICES);
    outlinks    = argo::conew_array<int>(MAX_VERTICES);
    W_index     = argo::conew_array<int*>(MAX_VERTICES);

    for (int i = 0; i < MAX_VERTICES; i++)
        W_index[i] = argo::conew_array<int>(MAX_DEGREE);

    // Read in graph from file by hopefully only one node
    read_graph_from_file(filename);

    argo::barrier();

    int vertices = *global_size;

    // Initialize ArgoDSM global memory variables
    if (argo::node_id() == 0)
        for(int i = 0; i < vertices; i++)
            pagerank[i] = INITIAL_PAGERANK;

    argo::barrier();

    // Divide the work as equal as possible among global_num_threads
    local_num_threads = global_num_threads / argo::number_of_nodes();

    if (argo::node_id() == 0) {
        int node_chunk = 0;
        int thread_chunk = vertices / global_num_threads;
        int rem = vertices % global_num_threads;
        for (int i = 0; i < global_num_threads; i++) {
            arguments[i].global_thread_id = i;
            if (node_chunk % local_num_threads == 0)
                node_chunk = 0;
            if (rem != 0 && arguments[i].global_thread_id < rem) {
                arguments[i].data_begin = node_chunk;
                arguments[i].data_end = arguments[i].data_begin + (thread_chunk+1);
                node_chunk += (thread_chunk+1);
                rem--;
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
    chunk->test     = std::vector<int>(chunk->size);
    chunk->exists   = std::vector<bool>(chunk->size);
    chunk->dangling = std::vector<bool>(chunk->size);
    chunk->outlinks = std::vector<int>(chunk->size);
    chunk->W_index  = std::vector<std::vector<int>>(chunk->size);

    for (int j = 0; j < chunk->size; j++) {
        chunk->test[j]     = test[chunk->begin+j];
        chunk->exists[j]   = exists[chunk->begin+j];
        chunk->dangling[j] = dangling[chunk->begin+j];
        chunk->outlinks[j] = outlinks[chunk->begin+j];
        chunk->W_index[j]  = std::vector<int>(MAX_DEGREE);
        for (int i = 0; i < MAX_DEGREE; i++)
            chunk->W_index[j][i] = W_index[chunk->begin+j][i];
    }

    std::vector<pthread_t> threads(local_num_threads);

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[i].vertices = vertices;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    argo::barrier();

    write_pagerank_to_file("file.txt", pagerank, vertices);

    delete chunk;
    delete lock;

    argo::codelete_array(test);
    argo::codelete_array(exists);
    argo::codelete_array(dangling);
    argo::codelete_array(outlinks);
    argo::codelete_array(W_index);
    argo::codelete_array(pagerank);
    argo::codelete_(lock_flag);
    argo::codelete_array(dp_threads);
    argo::codelete_(dp);
    argo::codelete_array(arguments);

    argo::finalize();
}

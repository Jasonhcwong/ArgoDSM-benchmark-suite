#include <pthread.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <mutex>

#include "argo.hpp"

// Connected Components Constants
int    MAX_VERTICES     = 2000000;
int    MAX_DEGREE       = 32;
int    MAX_COMPONENTS   = 100000000;


struct Thread_args {
    int local_thread_id;
    int data_begin;
    int data_end;
    int vertices;
    pthread_barrier_t* barrier;
};

struct Chunk {
    int begin;
    int size;
    std::vector<bool> exists;
    std::vector<int> graph;
    std::vector<int> edges;
};


// ArgoDSM node local variables
int node_id;
int number_of_nodes;
int global_num_threads;
int local_num_threads;
std::vector<int>* triangles_local;
Chunk* chunk;

argo::globallock::global_tas_lock *lock;
std::vector<argo::globallock::global_tas_lock*>* locks; //[2097152]; //change the number of locks to approx or greater N

// ArgoDMS global variables
Thread_args* arguments;
int*         global_size;
bool*        lock_flag;
bool*        lock_flags;
bool*        exists;
int*         edges;
int*         graph;
int*         triangle_edges;
int*         triangles_global;



void* do_work(void* argptr) {
    Thread_args* args = static_cast<Thread_args*>(argptr);

    int local_thread_id         = args->local_thread_id;
    int data_begin              = args->data_begin;
    int data_end                = args->data_end;
    pthread_barrier_t* barrier  = args->barrier;

    argo::barrier(local_num_threads);

    // First phase
    for (int j = data_begin; j < data_end; j++) {
        if (chunk->exists[j]) {
            for (int i = 0; i < chunk->edges[j]; i++) {
                int neighbor = chunk->graph[j*MAX_DEGREE+i];
                locks->at(neighbor)->lock();
                triangle_edges[neighbor]++;
                locks->at(neighbor)->unlock();
            }
        }
    }

    argo::barrier(local_num_threads);

    //Find triangles for each thread
    for (int i = data_begin; i < data_end; i++) {
        if (chunk->exists[i]) {
            int triangles = (triangle_edges[chunk->begin+i] / 3);
            if (triangles > 0)
                triangles_local->at(local_thread_id) += triangles;
        }
    }

    pthread_barrier_wait(barrier);

    if (local_thread_id == 0) {
        int local_sum = 0;
        for (int i = 0; i < local_num_threads; i++)
            local_sum += triangles_local->at(i);
        triangles_global[node_id] = local_sum;
    }

    argo::barrier(local_num_threads);

    return NULL;
}



void read_graph_from_file(std::string filename) {

    std::ifstream input(filename);

    for(int j = 0; j < MAX_VERTICES; j++) {
        for(int i = 0; i < MAX_DEGREE; i++) {
            graph[j*MAX_DEGREE+i] = MAX_COMPONENTS;
        }
        exists[j] = false;
        edges[j]  = 0;
    }

    int number0, number1;
    int vertex_cnt = 0;

    std::string line;
    while (input >> number0 >> number1) {

        if (number0 >= MAX_VERTICES) {
            std::cout << "Node " << number0 << " exceeds maximum graph size of " << MAX_VERTICES << "\n";
            exit (EXIT_FAILURE);
        }

        if (edges[number0] >= MAX_DEGREE) {
            std::cout << "Node " << number0 << " exceeds maximum maximum degree of " << MAX_DEGREE << "\n";
            exit (EXIT_FAILURE);
        }

        graph[number0*MAX_DEGREE+edges[number0]] = number1;
        edges[number0]++;

        exists[number0] = true;
        exists[number1] = true;

        if(number0 > vertex_cnt) vertex_cnt = number0;
        if(number1 > vertex_cnt) vertex_cnt = number1;
    }

    *global_size = ++vertex_cnt;

    /*
    std::cout << "Graph:\n";
    for (int j = 0; j < vertex_cnt; j++) {
        std::cout << "["<<j<<"]> ";
        for (int i = 0; i < MAX_DEGREE; i++) {
            if (graph[j*MAX_DEGREE+i] != MAX_COMPONENTS)
                std::cout << "[" << graph[j*MAX_DEGREE+i] << "] ";
        }
        std::cout << "\n";
    }
    */
}



int main(int argc, char** argv) {

    // Startup ArgoDSM
    argo::init(0.5*1024*1024*1024UL);

    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <thread-count> <input-file>" << std::endl;
        return 1;
    }

    node_id = argo::node_id();
    number_of_nodes = argo::number_of_nodes();
    global_num_threads = atoi(argv[1]);
    local_num_threads = global_num_threads / number_of_nodes;

    std::string input_filename = argv[2];

    if (!global_num_threads) {
        std::cout << "Thread count must be a valid integer greater than 0." << std::endl;
        return 1;
    }

    // Declare ArgoDSM global memory variables
    arguments           = argo::conew_array<Thread_args>(global_num_threads);
    global_size         = argo::conew_<int>(0);
    lock_flag           = argo::conew_<bool>(false);
    lock_flags          = argo::conew_array<bool>(MAX_VERTICES);
    triangle_edges      = argo::conew_array<int>(MAX_VERTICES);
    triangles_global    = argo::conew_array<int>(number_of_nodes);
    exists              = argo::conew_array<bool>(MAX_VERTICES);
    edges               = argo::conew_array<int>(MAX_VERTICES);
    graph               = argo::conew_array<int>((MAX_VERTICES)*(MAX_DEGREE));

    argo::barrier();

    // Read in graph from file by one node
    if (node_id == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    int vertices = *global_size;

    // Initialize ArgoDSM global variables
    if (node_id == 0)
        for(int i = 0; i < vertices; i++)
            triangle_edges[i] = 0;

    lock             = new argo::globallock::global_tas_lock(lock_flag);
    locks            = new std::vector<argo::globallock::global_tas_lock*>(vertices);
    triangles_local  = new std::vector<int>(local_num_threads);

    for (int i = 0; i < vertices; i++)
        locks->at(i) = new argo::globallock::global_tas_lock(&lock_flags[i]);

    argo::barrier();

    // Divide the work as equal as possible among global_num_threads
    if (node_id == 0) {
        int node_chunk = 0;
        int thread_chunk = vertices / global_num_threads;
        int rem = vertices % global_num_threads;
        for (int i = 0; i < global_num_threads; i++) {
            if (i % local_num_threads == 0) node_chunk = 0;
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

    chunk->exists = std::vector<bool>(chunk->size);
    chunk->edges  = std::vector<int>(chunk->size);
    chunk->graph  = std::vector<int>(chunk->size*MAX_DEGREE);

    for (int j = 0; j < chunk->size; j++) {
        chunk->exists[j] = exists[chunk->begin+j];
        chunk->edges[j]  = edges[chunk->begin+j];
        for (int i = 0; i < MAX_DEGREE; i++)
            chunk->graph[j*MAX_DEGREE+i] = graph[(chunk->begin+j)*MAX_DEGREE+i];
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

    int total_triangles = 0;
    if (node_id == 0) {
        for (int i = 0; i < number_of_nodes; i++)
            total_triangles += triangles_global[i];
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    argo::barrier();

    if (node_id == 0) {
        std::cout << "\nTriangle counting";
        std::cout << "\nArgo nodes: " << number_of_nodes;
        std::cout << "\nGlobal threads: " << global_num_threads;
        std::cout << "\nLocal threads: " << local_num_threads;
        std::cout << "\nGraph: " << input_filename;
        std::cout << "\nTriangles: " << total_triangles;
        std::cout << "\nRuntime: " << elapsed.count() << " ms\n" << std::endl;
    }

    for (int i = 0; i < vertices; i++)
        delete locks->at(i);

    delete triangles_local;
    delete locks;
    delete lock;
    delete chunk;

    argo::codelete_array(arguments);
    argo::codelete_(global_size);
    argo::codelete_(lock_flag);
    argo::codelete_(lock_flags);
    argo::codelete_array(exists);
    argo::codelete_array(edges);
    argo::codelete_array(graph);
    argo::codelete_array(triangle_edges);
    argo::codelete_array(triangles_global);

    argo::finalize();

    return 0;
}

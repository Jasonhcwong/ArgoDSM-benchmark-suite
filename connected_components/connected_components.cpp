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
    std::vector<bool> exists;
    std::vector<int> graph;
    std::vector<int> edges;
};


// ArgoDSM node local variables
int number_of_nodes;
int node_id;
int global_num_threads;
int local_num_threads;
argo::globallock::global_tas_lock *lock; // single lock
std::mutex mutex;
int not_changed_local;
Chunk* chunk;


// ArgoDMS global variables
Thread_args* arguments;
int*         not_changed_global;
int*         global_size;
bool*        lock_flag;
int*         components;
bool*        exists;
int*         edges;
int*         graph;


void* do_work(void* argptr) {
    Thread_args* args = static_cast<Thread_args*>(argptr);

    int global_thread_id        = args->global_thread_id;
    int local_thread_id         = args->local_thread_id;
    int data_begin              = args->data_begin;
    int data_end                = args->data_end;
    pthread_barrier_t* barrier  = args->barrier;

    bool mod                    = true;
    int iterations              = 0;

    argo::barrier(local_num_threads);

    //Each component is its own, first phase
    for(int i = data_begin; i < data_end; i++)
        components[chunk->begin+i] = chunk->begin+i;

    argo::barrier(local_num_threads);

    //start connecting, second phase
    while(*not_changed_global < number_of_nodes) {
        mod = false;
        iterations++;
        for (int j = data_begin; j < data_end; j++) {
            if (chunk->exists[j]) {
                for (int i = 0; i < chunk->edges[j]; i++) {
                    int neighbor = chunk->graph[j*MAX_DEGREE+i];
                    if((components[chunk->begin+j] < components[neighbor]) && (components[neighbor] == components[components[neighbor]])) {
                        components[components[neighbor]] = components[chunk->begin+j];
                        mod = true;
                    }
                }
            }
        }

        if (local_thread_id == 0)
            not_changed_local = 0;

        if (global_thread_id == 0)
            *not_changed_global = 0;

        argo::barrier(local_num_threads);

        if (!mod) {
            mutex.lock();
            not_changed_local++;
            mutex.unlock();
        }

        pthread_barrier_wait(barrier);

        //For termination Condition
        if (local_thread_id == 0) {
            if (not_changed_local == local_num_threads) {
                lock->lock();
                (*not_changed_global)++;
                lock->unlock();
            }
        }

        //Third phase, assign components
        for(int i = data_begin; i < data_end; i++) {
            while(components[chunk->begin+i] != components[components[chunk->begin+i]]) {
                components[chunk->begin+i] = components[components[chunk->begin+i]];
            }
        }

        argo::barrier(local_num_threads);
    }
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



void write_connected_components_to_file(std::string filename, int* components, int size) {
    std::ofstream output(filename);
    for(int i = 0; i < size; i++)
        if(exists[i]) {
            output << "Vertex: " << i << " Component: " << components[i] << "\n";
        }
    output.close();
}



int main(int argc, char** argv) {

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
        printf ("Error:  ");
        return 1;
    }

    Thread_args*  arguments = argo::conew_array<Thread_args>(global_num_threads);
    global_size         = argo::conew_<int>(0);
    not_changed_global  = argo::conew_<int>(0);
    lock_flag           = argo::conew_<bool>(false);
    components          = argo::conew_array<int>(MAX_VERTICES);
    exists              = argo::conew_array<bool>(MAX_VERTICES);
    edges               = argo::conew_array<int>(MAX_VERTICES);
    graph               = argo::conew_array<int>((MAX_VERTICES)*(MAX_DEGREE));

    argo::barrier();

    lock  = new argo::globallock::global_tas_lock(lock_flag);

    if (node_id == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    int vertices = *global_size;


    // Initialize ArgoDSM global variables
    if (node_id == 0)
        for(int j = 0; j < global_num_threads; j++)
            arguments[j].global_thread_id = j;

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

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].local_thread_id = i;
        arguments[j].vertices = vertices;
        arguments[j].barrier = &barrier;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    argo::barrier();

    if (node_id == 0) {
        write_connected_components_to_file(output_filename, components, vertices);
        std::cout << "\nConnected components\n";
        std::cout << "Argo nodes: " << number_of_nodes << "\nGlobal threads: " << global_num_threads << "\nLocal threads: " << local_num_threads << "\nGraph: " << input_filename << "\n";
        std::cout << "Runtime: " << elapsed.count() << " ms\n" << std::endl;
    }

    delete lock;

    argo::codelete_array(arguments);
    argo::codelete_(global_size);
    argo::codelete_(not_changed_global);
    argo::codelete_(lock_flag);
    argo::codelete_array(components);
    argo::codelete_array(exists);
    argo::codelete_array(edges);
    argo::codelete_array(graph);

    argo::finalize();

    return 0;
}

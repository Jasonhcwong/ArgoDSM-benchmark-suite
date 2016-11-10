#include <pthread.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "argo.hpp"

#define BILLION 1E9

// Connected Components Constants
int    MAX_VERTICES     = 2000000;
int    MAX_DEGREE       = 16;
int    MAX_COMPONENTS   = 100000000;


struct Thread_args {
    int global_thread_id;
    int vertices;
};


// ArgoDSM node local variables
int global_num_threads;
int local_num_threads;
argo::globallock::global_tas_lock *lock; // single lock


// ArgoDMS global variables
Thread_args* arguments;
int*         change;
int*         global_size;
bool*        lock_flag;
int*         components;
bool*        exists;
int*         edges;
int*         graph;



void* do_work(void* argptr) {
    Thread_args* args = static_cast<Thread_args*>(argptr);

    int global_thread_id     = args->global_thread_id;
    double vertices_d        = args->vertices;

    bool mod                 = true;
    int v                    = 0;
    int iterations           = 0;
    int start =  0;
    int stop  = 0;
    double tid_d = global_thread_id;
    double P_d = global_num_threads;

    //Chunk work for threads via double precision
    double start_d = (tid_d) * (vertices_d/P_d);
    double stop_d  = (tid_d+1.0) * (vertices_d/P_d);
    start =  start_d;
    stop  =  stop_d;

    argo::barrier(local_num_threads);

    //Each component is its own, first phase
    for(v = start; v < stop; v++)
        components[v] = v;

    argo::barrier(local_num_threads);

    //start connecting, second phase
    while(*change < global_num_threads) {
        mod = false;
        iterations++;
        for (v = start; v < stop; v++) {
            if (exists[v]) {
                for (int i = 0; i < edges[v]; i++) {
                    int neighbor = graph[v*MAX_DEGREE+i];
                    if((components[v] < components[neighbor]) && (components[neighbor] == components[components[neighbor]])) {
                        mod = true;
                        components[components[neighbor]] = components[v];
                    }
                }
            }
        }

        if (global_thread_id == 0)
            *change = 0;

        argo::barrier(local_num_threads);

        //Third phase, assign components
        for(v = start; v < stop; v++) {
            while(components[v] != components[components[v]]) {
                components[v] = components[components[v]];
            }
        }

        //For termination Condition
        if (!mod) {
            lock->lock();
            (*change)++;
            lock->unlock();
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

    global_num_threads = atoi(argv[1]);
    std::string input_filename = argv[2];
    std::string output_filename = argv[3];

    if (!global_num_threads) {
        std::cout << "Thread count must be a valid integer greater than 0." << std::endl;
        printf ("Error:  ");
        return 1;
    }

    Thread_args*  arguments = argo::conew_array<Thread_args>(global_num_threads);
    global_size  = argo::conew_<int>(0);
    change       = argo::conew_<int>(0);
    lock_flag    = argo::conew_<bool>(false);
    components   = argo::conew_array<int>(MAX_VERTICES);
    exists       = argo::conew_array<bool>(MAX_VERTICES);
    edges        = argo::conew_array<int>(MAX_VERTICES);
    graph        = argo::conew_array<int>((MAX_VERTICES)*(MAX_DEGREE));

    argo::barrier();

    lock  = new argo::globallock::global_tas_lock(lock_flag);

    if (argo::node_id() == 0)
        read_graph_from_file(input_filename);

    argo::barrier();

    int vertices = *global_size;


    // Initialize ArgoDSM global variables
    if (argo::node_id() == 0)
        for(int j = 0; j < global_num_threads; j++)
            arguments[j].global_thread_id = j;

    argo::barrier();

    // Divide the work as equal as possible among global_num_threads
    local_num_threads = global_num_threads / argo::number_of_nodes();

    int node_threads_begin = argo::node_id() * local_num_threads;
    std::vector<pthread_t> threads(local_num_threads);

    auto start = std::chrono::system_clock::now();

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].vertices = vertices;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    argo::barrier();

    if (argo::node_id() == 0) {
        write_connected_components_to_file(output_filename, components, vertices);
        std::cout << "\nConnected components\n";
        std::cout << "Argo nodes: " << argo::number_of_nodes() << "\nGlobal threads: " << global_num_threads << "\nLocal threads: " << local_num_threads << "\nGraph: " << input_filename << "\n";
        std::cout << "Runtime: " << elapsed.count() << " ms\n" << std::endl;
    }

    delete lock;

    argo::codelete_array(arguments);
    argo::codelete_(global_size);
    argo::codelete_(change);
    argo::codelete_(lock_flag);
    argo::codelete_array(components);
    argo::codelete_array(exists);
    argo::codelete_array(edges);
    argo::codelete_array(graph);

    argo::finalize();

    return 0;
}

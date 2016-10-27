#include <cassert>
#include <limits>
#include <iostream>
#include <vector>

#include <pthread.h>
#include <string.h>

#include "argo.hpp"

struct thread_args {
    int                global_thread_id;
    int                local_thread_id;
    int                data_begin;
    int                data_end;
    int                vertices;
    pthread_barrier_t* barrier;
};


// ArgoDMS global variables
thread_args* arguments;

bool *lock_flag;
int*  test;          // test variable arrays for graph storage
bool* exists;
bool* dangling;
int*  outlinks;      // array for outlinks

double* dp;
double* dp_threads;

double* pagerank;

double** W;
int** W_index;




// ArgoDSM node local variables
argo::globallock::global_tas_lock *lock; // single lock
int local_data_begin;
std::vector<int>*  test_local;          // test variable arrays for graph storage
std::vector<bool>* exists_local;
std::vector<bool>* dangling_local;
std::vector<int>*  outlinks_local;      // array for outlinks
std::vector<double>* pagerank_local;      // array for outlinks


void initialize_single_source(int vertices, double initial_rank);
void init_weights(int vertices, int degree);


void* do_work(void* argptr) {
    thread_args* args = static_cast<thread_args*>(argptr);

    // Let's declutter the code a bit, while also avoiding unnecessary global memory accesses
    int       global_thread_id  = args->global_thread_id;
    //int       local_thread_id   = args->local_thread_id;
    const int vertices          = args->vertices;
    int       data_begin        = args->data_begin;
    int       data_end          = args->data_end;
    pthread_barrier_t* barrier  = args->barrier;

    double    damping_factor    = 0.85;
    int       iterations        = 1;

    pthread_barrier_wait(barrier);
    //argo::barrier();

    while(iterations > 0) {
        // think about this
        if(global_thread_id == 0)
            *dp = 0;

        pthread_barrier_wait(barrier);
        argo::barrier();

        for(int i = data_begin; i < data_end; i++)
            if(dangling_local->at(i))
                dp_threads[global_thread_id] += damping_factor*(pagerank[local_data_begin+i]/vertices);

        lock->lock();
        *dp = *dp + dp_threads[global_thread_id];
        lock->unlock();

        pthread_barrier_wait(barrier);
        argo::barrier();

        for(int i = data_begin; i < data_end; i++) {
            if (exists_local->at(i)) {
                pagerank_local->at(i) = 0.15;
                for(int j = 0; j < test_local->at(i); j++)
                    pagerank_local->at(i) += damping_factor * pagerank[W_index[local_data_begin+i][j]] / outlinks_local->at(W_index[local_data_begin+i][j]);
            }
            if (pagerank_local->at(i) >= 1.0)
                pagerank_local->at(i) = 1.0;
        }

        pthread_barrier_wait(barrier);
        //argo::barrier();

        for(int i = data_begin; i < data_end; i++)
            if(exists_local->at(i)) {
                pagerank[local_data_begin+i] = pagerank_local->at(i);
                printf("%d\n", local_data_begin+i);
            }

        pthread_barrier_wait(barrier);
        argo::barrier();

        iterations--;

    }

    return NULL;
}



int main(int argc, char* argv[]) {

    FILE *file0 = NULL;
    FILE *f = NULL;
    int vertices = 0;                  // Total vertices
    int degree = 0;                       // Edges per vertex
    const int select = atoi(argv[1]);  // 0 for synthetic, 1 for file read
    char filename[100];

    // generate graph input
    if (select == 0) {
        vertices = atoi(argv[3]);
        degree = atoi(argv[4]);
        printf("\nGraph with Parameters: Vertices:%d degree:%d\n", vertices, degree);
    }
    // graph file input
    if (select == 1) {
        vertices = 2097152; // 4194304; // can be read from file if needed, this is a default upper limit
        degree = 16;           // also can be read from file if needed, upper limit here again
        strcpy(filename,argv[3]);
        f = fopen(filename,"r");
    }

    const int num_threads = atoi(argv[2]);  // number of threads

    if (degree > vertices) {
        fprintf(stderr, "Degree of graph cannot be grater than number of Vertices\n");
        exit(EXIT_FAILURE);
    }


    // Initialize ArgoDSM with 1 GB of Storage
    argo::init(0.5*1024*1024*1024UL);

    // ArgoDSM memory allocations
    arguments = argo::conew_array<thread_args>(num_threads);

    dp         = argo::conew_<double>(0);
    dp_threads = argo::conew_array<double>(num_threads);

    pagerank = argo::conew_array<double>(vertices);
    test     = argo::conew_array<int>(vertices);
    exists   = argo::conew_array<bool>(vertices);
    dangling = argo::conew_array<bool>(vertices);
    outlinks = argo::conew_array<int>(vertices);

    lock_flag = argo::conew_<bool>(false);       // ArgoDSM pthread_mutex_init(&lock, NULL);
    lock = new argo::globallock::global_tas_lock(lock_flag);

    W_index = argo::conew_array<int*>(vertices);

    for(int i = 0; i < vertices; i++) {
        W_index[i] = argo::conew_array<int>(degree);
    }

    // ArgoDSM memory initialization done by a single node
    if (argo::node_id() == 0) {

        for(int i = 0; i < vertices; i++) {
            for(int j = 0; j < degree; j++) {
                W_index[i][j] = std::numeric_limits<int>::max();
            }
            test[i]     = 0;
            exists[i]   = false;
            dangling[i] = false;
            outlinks[i] = 0;
        }

        // if graph read from file
        if (select == 1) {
            int lines = 0;
            int nodecount = 0;
            int lines_to_check = 0;      // file processing variables
            char c;
            int number0;
            int number1;
            int inter = -1;

            for(c = getc(f); c != EOF; c = getc(f)) {
                if(c=='\n')
                    lines++;
            }
            fclose(f);

            file0 = fopen(filename,"r");

            for(c = getc(file0); c != EOF; c = getc(file0)) {
                if(c=='\n')
                    lines_to_check++;

                if(lines_to_check>3 && lines_to_check<lines) {
                    int f0 = fscanf(file0, "%d %d", &number0,&number1);

                    if(f0 != 2 && f0 != EOF) {
                        printf ("Error: Read %d values, expected 2. Parsing failed.\n",f0);
                        exit (EXIT_FAILURE);
                    }

                    inter = test[number1];
                    W_index[number1][inter] = number0;
                    test[number1]++;
                    outlinks[number0]++;
                    exists[number0] = true;
                    exists[number1] = true;
                    dangling[number1] = true;
                    if(number0 > nodecount) nodecount = number0;
                    if(number1 > nodecount) nodecount = number1;
                }
            }
            nodecount++;
            for(int i = 0; i < vertices; i++) {
                if(exists[i] && dangling[i])
                    dangling[i] = false;
            }
            printf("\nLargest Vertex: %d", nodecount);
            vertices = nodecount;
        }

        // graph to be generated synthetically
        if(select == 0) {
            int div = vertices;
            if(degree >= vertices)
                div = degree;
            init_weights(vertices, degree);
            for(int i = 0; i < vertices; i++) {
                outlinks[i] = rand()%(div); //random outlinks
                if(outlinks[i] == 0)
                    outlinks[i] = vertices;
                exists[i] = true;
                if (i % 100 == 0)
                    dangling[i] = true;
                test[i] = degree;
            }
        }

        // initialize pageranks
        initialize_single_source(vertices, 0.15);
        printf("\nInitialization Done");
    }

    argo::barrier();

    int local_num_threads = num_threads / argo::number_of_nodes();

    if (argo::node_id() == 0) {
        int node_chunk = 0;
        int thread_chunk = vertices / num_threads;
        int rem = vertices % num_threads;
        for (int i = 0; i < num_threads; i++) {
            arguments[i].global_thread_id = i;
            arguments[i].vertices = vertices;
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

    int node_threads_begin = argo::node_id() * local_num_threads;
    int node_chunk = arguments[node_threads_begin + local_num_threads - 1].data_end;

    test_local     = new std::vector<int>(node_chunk);
    exists_local   = new std::vector<bool>(node_chunk);
    dangling_local = new std::vector<bool>(node_chunk);
    outlinks_local = new std::vector<int>(node_chunk);
    pagerank_local = new std::vector<double>(node_chunk);

    local_data_begin = 0;
    for (int i = 0; i < node_threads_begin; i++)
        local_data_begin += arguments[i].data_end - arguments[i].data_begin;

    for (int i = 0; i < node_chunk; i++) {
        test_local->at(i)     = test[local_data_begin+i];
        exists_local->at(i)   = exists[local_data_begin+i];
        dangling_local->at(i) = dangling[local_data_begin+i];
        outlinks_local->at(i) = outlinks[local_data_begin+i];
    }

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, local_num_threads);
    std::vector<pthread_t> threads(local_num_threads);

    for (int i = 0; i < local_num_threads; i++) {
        int j = node_threads_begin + i;
        arguments[j].local_thread_id = i;
        arguments[j].barrier = &barrier;
        pthread_create(&threads[i], nullptr, do_work, &arguments[j]);
    }

    for (auto &t : threads)
        pthread_join(t, nullptr);

    argo::barrier();

    if (argo::node_id() == 0) {
        //Print pageranks to file
        FILE *f1 = fopen("file.txt", "w");

        for(int i = 0; i < vertices; i++) {
            if(exists[i])
                fprintf(f1,"pr(%d) = %f\n", i,pagerank[i]);
        }
        printf("\n");
        fclose(f1);
    }

    delete test_local;
    delete exists_local;
    delete dangling_local;
    delete outlinks_local;
    delete pagerank_local;
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



void initialize_single_source(int vertices, double initial_rank) {
    for(int i = 0; i < vertices; i++)
        pagerank[i] = initial_rank;
}



void init_weights(int vertices, int degree) {
    for(int i = 0; i < vertices; i++)
        for(int j = 0; j < degree; j++)
            W_index[i][j]= -1;

    for(int i = 0; i < vertices; i++) {
        for(int j = 0; j < degree; j++) {
            if(W_index[i][j] == -1) {
                int neighbor = rand() % (i + degree * 2);
                if(neighbor < j)
                    W_index[i][j] = neighbor;
                else
                    W_index[i][j] = vertices-1;
            }
            if(W_index[i][j] >= vertices)
                W_index[i][j] = vertices-1;
        }
    }
}

#include <cassert>
#include <limits>
#include <iostream>
#include <vector>

#include <pthread.h>
#include <string.h>

#include "argo.hpp"

struct thread_args {
    int         tid;
    int         data_begin;
    int         data_end;
    int         vertices;
};


// global variables
bool *lock_flag;
argo::globallock::global_tas_lock *lock; // single lock

int *test;          // test variable arrays for graph storage
int *exist;
int *test2;
int *dangling;
int *outlinks;      // array for outlinks

double *dp;      // dangling pointer variable
double *dp_tid;

double *pr_tmp;      // temporary pageranks
double *pagerank;

double** W;
int** W_index;

void initialize_single_source(int vertices, double initial_rank);
void init_weights(int vertices, int degree);


void* do_work(void* argptr) {
    // get the arguments
    thread_args* args = static_cast<thread_args*>(argptr);

    // Let's declutter the code a bit, while also avoiding unnecessary global
    // memory accesses
    int tid                  = args->tid;
    const int vertices       = args->vertices;
    int i_start              = args->data_begin;
    int i_stop               = args->data_end;

    double d                = 0.85;          // damping coefficient

    // pagerank iteration count
    int iterations = 1;

    argo::barrier(); // ArgoDSM pthread_barrier_wait(arg->barrier);

    while(iterations > 0) {
        if(tid == 0)
            *dp = 0;

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        // for no outlinks, dangling pages calculation
        for(int v = i_start; v < i_stop; v++)
            if(dangling[v] == 1)
                dp_tid[tid] = dp_tid[tid] + d*(pagerank[v]/vertices);

        lock->lock();           // ArgoDSM pthread_mutex_lock(&lock);
        *dp = *dp + dp_tid[tid];
        lock->unlock();         // ArgoDSM pthread_mutex_unlock(&lock);

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        // calculate pageranks
        for(int v = i_start; v < i_stop; v++) {
            if (exist[v] == 1) { // if vertex exists
                pr_tmp[v] = 0.15;
                for(int j = 0; j < test[v]; j++) // if inlink
                    pr_tmp[v] = pr_tmp[v] + (d * pagerank[W_index[v][j]] / outlinks[W_index[v][j]]);
                    /* replace d with dp if dangling pageranks required */
            }
            if (pr_tmp[v] >= 1.0) pr_tmp[v] = 1.0;
        }

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        // put temporary pageranks into final pageranks
        for(int v = i_start; v < i_stop; v++)
            if(exist[v] == 1)
                pagerank[v] = pr_tmp[v];

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

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

    printf("debug before argoDSM\n");

    // ArgoDSM memory allocations
    //
    dp       = argo::conew_<double>(0);
    dp_tid   = argo::conew_array<double>(num_threads);

    pagerank = argo::conew_array<double>(vertices);
    pr_tmp   = argo::conew_array<double>(vertices);
    test     = argo::conew_array<int>(vertices);
    exist    = argo::conew_array<int>(vertices);
    test2    = argo::conew_array<int>(vertices);
    dangling = argo::conew_array<int>(vertices);
    outlinks = argo::conew_array<int>(vertices);

    lock_flag = argo::conew_<bool>(false);       // ArgoDSM pthread_mutex_init(&lock, NULL);
    lock = new argo::globallock::global_tas_lock(lock_flag);

    W = argo::conew_array<double*>(vertices);
    W_index = argo::conew_array<int*>(vertices);

    for(int i = 0; i < vertices; i++) {
        W[i] = argo::conew_array<double>(degree);
        W_index[i] = argo::conew_array<int>(degree);
    }

    printf("debug before single node argoDSM\n");

    // ArgoDSM memory initialization done by a single node
    if (argo::node_id() == 0) {

        for(int i = 0; i < vertices; i++) {
            for(int j = 0; j < degree; j++) {
                W[i][j] = std::numeric_limits<double>::max();
                W_index[i][j] = std::numeric_limits<int>::max();
            }
            test[i]     = 0;
            exist[i]    = 0;
            test2[i]    = 0;
            dangling[i] = 0;
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
                    W[number0][inter] = 0; //drand48();
                    W_index[number1][inter] = number0;
                    test[number1]++;
                    outlinks[number0]++;
                    exist[number0] = 1;
                    exist[number1] = 1;
                    test2[number0] = 1;
                    dangling[number1] = 1;
                    if(number0 > nodecount) nodecount = number0;
                    if(number1 > nodecount) nodecount = number1;
                }
            }
            nodecount++;
            for(int i = 0; i < vertices; i++) {
                if(test2[i] == 1 && dangling[i] == 1)
                    dangling[i] = 0;
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
                exist[i] = 1;
                if (i % 100 == 0)
                    dangling[i] = 1;
                test[i] = degree;
            }
        }

        // initialize pageranks
        initialize_single_source(vertices, 0.15);
        printf("\nInitialization Done");
    }

    argo::barrier();

    // Start the threads
    int local_num_threads = num_threads / argo::number_of_nodes();
    int chunk = vertices / num_threads;
    std::vector<pthread_t> threads(local_num_threads);
    std::vector<thread_args> args(local_num_threads);
    for (int i = 0; i < local_num_threads; ++i) {
        int global_tid = (argo::node_id() * local_num_threads) + i;
        args[i].tid = global_tid;
        args[i].data_begin = global_tid * chunk;
        args[i].data_end = args[i].data_begin + chunk;
        args[i].vertices = vertices;
        pthread_create(&threads[i], nullptr, do_work, &args[i]);
    }
    // Join the threads
    for (auto &t : threads)
        pthread_join(t, nullptr);

    // Make sure everyone is done and get the changes
    argo::barrier();

    if (argo::node_id() == 0) {
        //Print pageranks to file
        FILE *f1 = fopen("file.txt", "w");

        for(int i = 0; i < vertices; i++) {
            if(exist[i] == 1)
                fprintf(f1,"pr(%d) = %f\n", i,pagerank[i]);
        }
        printf("\n");
        fclose(f1);
    }

    argo::finalize();
}



void initialize_single_source(int vertices, double initial_rank) {
   for(int i = 0; i < vertices; i++) {
      pagerank[i] = initial_rank;
      pr_tmp[i] = initial_rank;
   }
}


void init_weights(int vertices, int degree)
{
   // Initialize to -1
   for(int i = 0; i < vertices; i++)
      for(int j = 0; j < degree; j++)
         W_index[i][j]= -1;

   // Populate Index Array
   for(int i = 0; i < vertices; i++)
   {
      int max = degree;
      for(int j = 0; j < degree; j++)
      {
         if(W_index[i][j] == -1)
         {
            int neighbor = rand()%(i+max*2);
            if(neighbor<j)
               W_index[i][j] = neighbor;
            else
               W_index[i][j] = vertices-1;
         }
         else
         {
         }
         if(W_index[i][j]>=vertices)
         {
            W_index[i][j] = vertices-1;
         }
      }
   }

   // Populate Cost Array
   for(int i = 0; i < vertices; i++)
   {
      for(int j = 0; j < degree; j++)
      {
        if(W_index[i][j] == i)
         W[i][j] = 0;

         else
            W[i][j] = 0;//(double) (v) + 1;
      }
   }
}

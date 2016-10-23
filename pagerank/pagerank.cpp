#include <cassert>
#include <limits>
#include <iostream>
#include <vector>

#include <pthread.h>

#include "argo.hpp"

struct thread_args {
    int*      local_min;
    int*      global_min;
/*
    int*      Q;
    double*   PR;
    double**  W;
    int**     W_index;
    int       tid;
*/
    int       P;
    int       N;
    int       DEG;
};


// global variables
bool *lock_flag;
argo::globallock::global_tas_lock *lock; // single lock

argo::globallock::global_tas_lock *locks;

int* local_min_buffer;
double* dp_tid;
int* global_min_buffer;

int terminate = 0;  // terminate variable

int *test;          // test variable arrays for graph storage
int *exist;
int *test2;
int *dangling;
int *outlinks;      // array for outlinks

double* PR;
int* Q;

double dp = 0;      // dangling pointer variable
double *pgtmp;      // temporary pageranks
int nodecount = 0;

int initialize_single_source(double* PR, int* Q, int source, int N, double initial_rank);
void init_weights(int N, int DEG, double** W, int** W_index);


void* do_work(void* args) {
    // get the arguments
    thread_args* args = static_cast<thread_args*>(argptr);

    // Let's declutter the code a bit, while also avoiding unnecessary global
    // memory accesses
    int tid                    = args->tid;
    double* PR                 = args->PR;
    int** W_index              = args->W_index;
    const int N                = args->N;
    int v                      = 0;             // variable for current vertex
    double r                   = 0.15;          // damping coefficient
    double d                   = 0.85;          // damping coefficient
    double N_real              = N;
    double tid_d = tid;
    double P_d = args->P;

    // allocate work among threads
    double start_d = (tid_d) * (N_real/P_d);
    double stop_d = (tid_d+1.0) * (N_real/P_d);
    int i_start = start_d;
    int i_stop  = stop_d;

    // pagerank iteration count
    int iterations = 1;

    // barrier before starting work
    pthread_barrier_wait(arg->barrier);

    while(iterations > 0) {
        if(tid==0)
            dp=0;

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        // for no outlinks, dangling pages calculation
        for(v=i_start;v<i_stop;v++)
            if(dangling[v]==1)
                dp_tid[tid] = dp_tid[tid] + d*(PR[v]/N_real);

        lock->lock();           // ArgoDSM pthread_mutex_lock(&lock);
        dp = dp + dp_tid[tid];
        lock->unlock();         // ArgoDSM pthread_mutex_unlock(&lock);

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        v = 0;

        // calculate pageranks
        for(v = i_start; v < i_stop; v++) {
            if (exist[v]==1) { // if vertex exists
                pgtmp[v] = r;
                for(int j = 0; j < test[v]; j++) // if inlink
                    pgtmp[v] = pgtmp[v] + (d * PR[W_index[v][j]] / outlinks[W_index[v][j]]);
                    /* replace d with dp if dangling PRs required */
            }
            if (pgtmp[v] >= 1.0) pgtmp[v] = 1.0;
        }

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        // put temporary pageranks into final pageranks
        for(v = i_start; v < i_stop; v++)
            if(exist[v] == 1)
                PR[v] = pgtmp[v];

        argo::barrier();        // ArgoDSM pthread_barrier_wait(arg->barrier);

        iterations--;
    }

    return NULL;
}



int main(int argc, char* argv[]) {

    FILE *file0 = NULL;
    FILE *f = NULL;
    int N = 0;                         //Total vertices
    int DEG = 0;                       //Edges per vertex
    const int select = atoi(argv[1]);  //0 for synthetic, 1 for file read
    char filename[100];

    int lines_to_check = 0;      //file processing variables
    char c;
    int number0;
    int number1;
    int inter = -1;

    // generate graph input
    if (select == 0) {
        N = atoi(argv[3]);
        DEG = atoi(argv[4]);
        printf("\nGraph with Parameters: N:%d DEG:%d\n",N,DEG);
    }
    // graph file input
    if (select == 1) {
        N = 2097152;  // 4194304; can be read from file if needed, this is a default upper limit
        DEG = 16;     // also can be read from file if needed, upper limit here again
        strcpy(filename,argv[3]);
        f = fopen(filename,"r");
    }
    const int P = atoi(argv[2]);  // number of threads

    if (DEG > N) {
        fprintf(stderr, "Degree of graph cannot be grater than number of Vertices\n");
        exit(EXIT_FAILURE);
    }

    // Initialize ArgoDSM with 1 GB of Storage
    argo::init(1*1024*1024*1024UL);

    //Memory allocations
    PR       = argo::conew_array<double>(N);
    Q        = argo::conew_array<int>(N);
    test     = argo::conew_array<int>(N);
    exist    = argo::conew_array<int>(N);
    test2    = argo::conew_array<int>(N);
    dangling = argo::conew_array<int>(N);
    pgtmp    = argo::conew_array<double>(N);
    outlinks = argo::conew_array<int>(N);


    double** W = (double**) malloc(N*sizeof(double*));
    int** W_index = (int**) malloc(N*sizeof(int*));

    for(int i = 0; i < N; i++)
    {
      //W[i] = (int *)malloc(sizeof(int)*N);
      int ret = posix_memalign((void**) &W[i], 64, DEG*sizeof(double));
      int re1 = posix_memalign((void**) &W_index[i], 64, DEG*sizeof(int));
      if (ret != 0 || re1!=0)
      {
         fprintf(stderr, "Could not allocate memory\n");
         exit(EXIT_FAILURE);
      }
    }

    //Memory initialization
    for(int i=0;i<N;i++)
    {
      for(int j=0;j<DEG;j++)
      {
         W[i][j] = 1000000000;
         W_index[i][j] = INT_MAX;
      }
      test[i]=0;
      exist[i]=0;
      test2[i]=0;
      dangling[i]=0;
      outlinks[i]=0;
    }

    //If graph read from file
    if(select==1)
    {
      int lines=0;
      for(c=getc(f); c!=EOF; c=getc(f))
      {
         if(c=='\n')
            lines++;
      }
      fclose(f);

      file0 = fopen(filename,"r");

      for(c=getc(file0); c!=EOF; c=getc(file0))
      {
         if(c=='\n')
            lines_to_check++;

         if(lines_to_check>3 && lines_to_check<lines)
         {
            int f0 = fscanf(file0, "%d %d", &number0,&number1);
            if(f0 != 2 && f0 != EOF)
            {
               printf ("Error: Read %d values, expected 2. Parsing failed.\n",f0);
               exit (EXIT_FAILURE);
            }

            inter = test[number1];

            W[number0][inter] = 0; //drand48();
            W_index[number1][inter] = number0;
            test[number1]++;
            outlinks[number0]++;
            exist[number0]=1; exist[number1]=1;
            test2[number0]=1;
            dangling[number1]=1;
            if(number0 > nodecount)
               nodecount = number0;
            if(number1 > nodecount)
               nodecount = number1;
         }
      }
      nodecount++;
      for(int i=0;i<N;i++)
      {
         if(test2[i]==1 && dangling[i]==1)
            dangling[i]=0;
      }
      //printf("\n\nFile Read %d",lines_to_check);

      // Calculate total nodes, in order to calculate an initial weight.
      /*for(int i=0;i<N;i++)
        {
        if (test1[i]==1)
        nodecount++;
        }*/
      printf("\nLargest Vertex: %d",nodecount);
      N = nodecount;
    }

    //If graph to be generated synthetically
    if(select==0)
    {
         int div = N;
         if(DEG>=N)
             div = DEG;
      init_weights(N, DEG, W, W_index);
      for(int i=0;i<N;i++)
      {
         outlinks[i] = rand()%(div); //random outlinks
         if(outlinks[i]==0)
            outlinks[i] = N;
      }
    }

    //Synchronization parameters
    pthread_barrier_init(&barrier, NULL, P);
    pthread_mutex_init(&lock, NULL);

    for(int i=0; i<N; i++)
    {
      if(select==0)
      {
         exist[i]=1;
         if(i%100==0)
         {
            dangling[i]=1;
         }
         test[i] = DEG;
      }
      pthread_mutex_init(&locks[i], NULL);
    }

    //Initialize PageRanks
    initialize_single_source(PR, Q, 0, N, 0.15);
    printf("\nInitialization Done");

    //Thread arguments
    for(int j = 0; j < P; j++) {
      thread_arg[j].local_min  = local_min_buffer;
      thread_arg[j].global_min = &global_min_buffer;
      thread_arg[j].Q          = Q;
      thread_arg[j].PR         = PR;
      thread_arg[j].W          = W;
      thread_arg[j].W_index    = W_index;
      thread_arg[j].tid        = j;
      thread_arg[j].P          = P;
      thread_arg[j].N          = N;
      thread_arg[j].DEG        = DEG;
      thread_arg[j].barrier    = &barrier;
    }

    //Start CPU clock
    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);

    // Enable Graphite performance and energy models
    //CarbonEnableModels();

    //Spawn Threads
    for(int j = 1; j < P; j++) {
      pthread_create(thread_handle+j,
            NULL,
            do_work,
            (void*)&thread_arg[j]);
    }
    do_work((void*) &thread_arg[0]);  //Spawn main

    //Join threads
    for(int j = 1; j < P; j++) { //mul = mul*2;
      pthread_join(thread_handle[j],NULL);
    }

    // Disable Graphite performance and energy models
    //CarbonDisableModels();

    //Read clock and print time
    clock_gettime(CLOCK_REALTIME, &requestEnd);
    double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;
    printf( "\nTime:%lf seconds\n", accum );

    //printf("\ndistance:%d \n",D[N-1]);

    //Print pageranks to file
    FILE *f1 = fopen("file.txt", "w");

    for(int i = 0; i < N; i++) {
      if(exist[i]==1)
         fprintf(f1,"pr(%d) = %f\n", i,PR[i]);
    }
    printf("\n");
    fclose(f1);

    return 0;
}


void init_weights(int N, int DEG, double** W, int** W_index)
{
   // Initialize to -1
   for(int i = 0; i < N; i++)
      for(int j = 0; j < DEG; j++)
         W_index[i][j]= -1;

   // Populate Index Array
   for(int i = 0; i < N; i++)
   {
      int max = DEG;
      for(int j = 0; j < DEG; j++)
      {
         if(W_index[i][j] == -1)
         {
            int neighbor = rand()%(i+max*2);
            if(neighbor<j)
               W_index[i][j] = neighbor;
            else
               W_index[i][j] = N-1;
         }
         else
         {
         }
         if(W_index[i][j]>=N)
         {
            W_index[i][j] = N-1;
         }
      }
   }

   // Populate Cost Array
   for(int i = 0; i < N; i++)
   {
      for(int j = 0; j < DEG; j++)
      {
         /*if(v > 0.8 || W_index[i][j] == -1)
           {       W[i][j] = MAX;
           W_index[i][j] = -1;
           }

           else*/ if(W_index[i][j] == i)
         W[i][j] = 0;

         else
            W[i][j] = 0;//(double) (v) + 1;
      }
   }
}


/*
int main(int argc, char* argv[]) {
    int data_length = 160000;
    int num_threads = 16;
    int local_num_threads;

    // We totally need 10GB for this application
    argo::init(1*1024*1024*1024UL);

    local_num_threads = num_threads / argo::number_of_nodes();

    // Initialize the lock
    lock_flag = argo::conew_<bool>(false);
    lock = new argo::globallock::global_tas_lock(lock_flag);
    // Allocate the array
    data = argo::conew_array<int>(data_length);
    max = argo::conew_<int>(std::numeric_limits<int>::min());

    // Initialize the input data
    if (argo::node_id() == 0) {
        for (int i = 0; i < data_length; ++i)
            data[i] = i * 11 + 3;
    }

    // Make sure initialization is done and distribute the changes
    argo::barrier();

    // Start the threads
    int chunk = data_length / num_threads;
    std::vector<pthread_t> threads(local_num_threads);
    std::vector<thread_args> args(local_num_threads);
    for (int i = 0; i < local_num_threads; ++i) {
        int global_tid = (argo::node_id() * local_num_threads) + i;
        args[i].data_begin = global_tid * chunk;
        args[i].data_end = args[i].data_begin + chunk;
        pthread_create(&threads[i], nullptr, parmax, &args[i]);
    }
    // Join the threads
    for (auto &t : threads)
        pthread_join(t, nullptr);

    // Make sure everyone is done and get the changes
    argo::barrier();

    argo::codelete_array(data);
        delete lock;
    argo::codelete_(lock_flag);
    // Print the result
    if (argo::node_id() == 0)
        std::cout << "Max found to be " << *max << std::endl;
    assert(*max == ((data_length - 1) * 11 + 3));

    argo::finalize();
}
*/



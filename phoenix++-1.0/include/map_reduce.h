/* Copyright (c) 2007-2011, Stanford University
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Stanford University nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY STANFORD UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL STANFORD UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef MAP_REDUCE_H_
#define MAP_REDUCE_H_

#include <assert.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <cstring>

#include "stddefines.h"
#include "processor.h"
#include "scheduler.h"
#include "task_queue.h"
#include "combiner.h"
#include "container.h"
#include "locality.h"
#include "thread_pool.h"

#include "argo/argo.hpp"

template<typename Impl, typename D, typename K, typename V,
         class Container = hash_container<K, V, buffer_combiner> >
class MapReduce
{
public:

    /* Standard data types for the function arguments and results */
    typedef D data_type;
    typedef V value_type;
    typedef K key_type;
    typedef Container container_type;

    typedef typename container_type::input_type map_container;
    typedef typename container_type::output_type reduce_iterator;

    struct keyval
    {
        key_type key;
        value_type val;
    };

protected:

    // Parameters.
    uint64_t num_threads;               // # of threads to run.
    uint64_t thread_offset;             // cores to skip when assigning threads.

    thread_pool* threadPool;            // Thread pool.
    task_queue* taskQueue;              // Queues of tasks.

    container_type container;
    std::vector<keyval>* final_vals;    // Array to send to merge task.

    uint64_t num_map_tasks;
    uint64_t num_reduce_tasks;


    virtual void run_map(data_type* data, uint64_t len);
    virtual void run_reduce();
    virtual void run_merge();

    virtual void map_worker(
        thread_loc const& loc, double& time, double& user_time, int& tasks);
    virtual void reduce_worker(
        thread_loc const& loc, double& time, double& user_time, int& tasks);
    virtual void merge_worker(
        thread_loc const& loc, double& time, double& user_time, int& tasks);

    // Data passed to the callback functions.
    struct thread_arg_t
    {
        // in
        MapReduce* mr;
        // out
        double user_time;
        double time;
        int tasks;
    };

    static void map_callback(void* arg, thread_loc const& loc) {
        thread_arg_t* t = (thread_arg_t*)arg;
        t->mr->map_worker(loc, t->time, t->user_time, t->tasks);
    }
    static void reduce_callback(void* arg, thread_loc const& loc) {
        thread_arg_t* t = (thread_arg_t*)arg;
        t->mr->reduce_worker(loc, t->time, t->user_time, t->tasks);
    }
    static void merge_callback(void* arg, thread_loc const& loc) {
        thread_arg_t* t = (thread_arg_t*)arg;
        t->mr->merge_worker(loc, t->time, t->user_time, t->tasks);
    }
    void start_workers (void (*callback)(void*, thread_loc const&),
                        int num_threads, char const* stage);

    // the default split function...
    int split(data_type &a) { return 0; }

    // the default map function...
    void map(data_type const& a, map_container& m) const {}

    // the default reduce function...
    void reduce(key_type const& key, reduce_iterator const& values,
                std::vector<keyval>& out) const {
        value_type val;
        while (values.next(val))
        {
            keyval kv = {key, val};
            out.push_back(kv);
        }
    }

    // the default locator function...
    void* locate(data_type* data, uint64_t) const {
        return (void*)data;
    }

public:

    MapReduce() : threadPool(NULL), taskQueue(NULL) {
        // Determine the number of threads to use.
        // First check for an environment variable, then use the
        // number of processors
        int threads = atoi(GETENV("MR_NUMTHREADS"));
        //setThreads(threads > 0 ? threads : proc_get_num_cpus(), 0);
        setThreads(8, 0);
    }

    virtual ~MapReduce() {
        if (this->threadPool != NULL) delete this->threadPool;
        if (this->taskQueue != NULL) delete this->taskQueue;
    }

    // override the default thread offset and thread count.
    MapReduce& setThreads(int num_threads, sched_policy const* policy = NULL) {
        this->num_threads = (num_threads > 0) ? num_threads : this->num_threads;

        if (this->threadPool != NULL) delete this->threadPool;
        if (this->taskQueue != NULL) delete this->taskQueue;

        // Create thread pool and task queue
        sched_policy_strand_fill default_policy(0);
        this->threadPool = new thread_pool(
            num_threads, policy == NULL ? &default_policy : policy);
        this->taskQueue = new task_queue(num_threads, num_threads);

        return *this;
    }

    /* The main MapReduce engine. This is the function called by the
     * application. It is responsible for creating and scheduling all map
     * and reduce tasks, and also organizes and maintains the data which is
     * passed from application to map tasks, map tasks to reduce tasks, and
     * reduce tasks back to the application. Results are stored in result.
     * A return value less than zero represents an error. This function is
     * not thread safe.
     */
    int run(data_type *data, uint64_t count, std::vector<keyval>& result);

    // This version assumes that the split function is provided.
    int run(std::vector<keyval>& result);

    void emit_intermediate(typename container_type::input_type& i,
                           key_type const& k, value_type const& v) const {
        i[k].add(v);
    }
};

template<typename Impl, typename D, typename K, typename V, class Container>
int MapReduce<Impl, D, K, V, Container>::
run (std::vector<keyval>& result)
{
    timespec begin;
    std::vector<D> data;

    uint64_t *count;
    D chunk;
    D *argo_data;

    uint64_t number_of_nodes = argo::number_of_nodes();
    uint64_t node_id = argo::node_id();

    // Run splitter to generate chunks
    int i = 0;
    get_time (begin);
    count = argo::conew_<uint64_t>(0);
    if (argo::node_id() == 0) {
        while (static_cast<Impl *>(this)->split(chunk)) {
            data.push_back(chunk);
        }
        *count = data.size();
    }
    argo::barrier();
    argo_data = argo::conew_array<D>(*count);
    if (argo::node_id() == 0) {
        std::memcpy(argo_data, &data[0], sizeof(D) * (*count));
    }
    argo::barrier();
    print_time_elapsed("split phase", begin);
    int chunk_size = (int)ceil((double)(*count) / number_of_nodes);
    int index = chunk_size * node_id;
    if (index + chunk_size > *count) chunk_size = *count - index;
    printf("node%d: count: %d, index: %d, chunk_size: %d\n", node_id, *count, index, chunk_size);
    return run(&argo_data[index], chunk_size, result);
}

template<typename Impl, typename D, typename K, typename V, class Container>
int MapReduce<Impl, D, K, V, Container>::
run (D *data, uint64_t count, std::vector<keyval>& result)
{
    timespec begin;
    timespec run_begin = get_time();
    // Initialize library
    get_time (begin);


    // Compute task counts (should make this more adjustable) and then
    // allocate storage
    this->num_map_tasks = std::min(count, this->num_threads) * 16;
    this->num_reduce_tasks = this->num_threads;
    printf ("num_map_tasks = %d, node_id = %d\n", num_map_tasks, argo::node_id());
    printf ("num_reduce_tasks = %d, node_id = %d\n", num_reduce_tasks, argo::node_id());

    container.init(this->num_threads, this->num_reduce_tasks);
    this->final_vals = new std::vector<keyval>[this->num_threads];
    for (uint64_t i = 0; i < this->num_threads; i++) {
        // Try to avoid a reallocation. Very costly on Solaris.
        this->final_vals[i].reserve(100);
    }
    printf("1 argo node: %d, finalvals size: %d\n", argo::node_id(), this->final_vals->size());

    print_time_elapsed("library init", begin);
    // Run map tasks and get intermediate values
    get_time (begin);
    run_map(&data[0], count);
    argo::barrier();
    print_time_elapsed("map phase", begin);

    dprintf("In scheduler, all map tasks are done, now scheduling reduce tasks\n");

    // Run reduce tasks and get final values
    get_time (begin);
    run_reduce();
    argo::barrier();
    print_time_elapsed("reduce phase", begin);

    dprintf("In scheduler, all reduce tasks are done, now scheduling merge tasks\n");

    get_time (begin);
    run_merge();
    print_time_elapsed("merge phase", begin);

    result.swap(*this->final_vals);

    // Delete structures
    delete [] this->final_vals;

    print_time_elapsed("run time", run_begin);

    return 0;
}

/**
 * Run map tasks and get intermediate values
 */
template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::
run_map (data_type * data, uint64_t count)
{

    printf("node%i: run_map count: %d\n", argo::node_id(), count);
    // Compute map task chunk size
    uint64_t chunk_size =
        std::max(1, (int)ceil((double)count / this->num_map_tasks));

    // Generate tasks by splitting input data and add to queue.
    //for(uint64_t i = node_number; i < this->num_map_tasks; i +=number_of_nodes)
    for (uint64_t i = 0; i < this->num_map_tasks; i++)
    {
        uint64_t start = chunk_size * i;

        if (start < count)
        {
            uint64_t len = std::min(chunk_size, count - start);
            int lgrp = loc_mem_to_lgrp (
                           static_cast<Impl *>(this)->locate(data + start, len));
            task_queue::task_t task =
                // For debugging, last element is normally padding
            {    i, len, (uint64_t)(data + start), lgrp };

            this->taskQueue->enqueue_seq (task, this->num_map_tasks, lgrp);
            //printf("node%i: i: %d, start: %d, chunk_size: %d, len: %d, data+start: %p, lgrp: %d\n", argo::node_id(), i, start, chunk_size, len, (data + start), lgrp);
        }

    }

    this->start_workers (&map_callback, std::min(num_map_tasks, num_threads), "map");
}

/**
 * Dequeue the latest task and run it
 */
template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::
map_worker(thread_loc const & loc, double & time, double & user_time, int& tasks)
{
    timespec begin = get_time();
    typename container_type::input_type t = container.get(loc.thread);
    task_queue::task_t task;
    while (taskQueue->dequeue (task, loc)) {
        tasks++;
        timespec user_begin = get_time();
        for (data_type* data = (data_type*)task.data;
                data < (data_type*)task.data + task.len; ++data) {
            static_cast<Impl *>(this)->map(*data, t);
        }
        user_time += time_elapsed(user_begin);
    }

    container.add(loc.thread, t);
    time += time_elapsed(begin);
}

/**
 * Run reduce tasks and get final values.
 */
template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::run_reduce ()
{
    // Create tasks and enqueue...
    for (uint64_t i = 0; i < this->num_reduce_tasks; ++i) {
        task_queue::task_t task = {    i, 0, i, 0 };
        this->taskQueue->enqueue_seq(task, this->num_reduce_tasks);
    }

    this->start_workers (&reduce_callback,
                         std::min(this->num_reduce_tasks, num_threads), "reduce");
}

/**
 * Dequeue next reduce task and do it
 */
template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::reduce_worker (
    thread_loc const & loc, double & time, double & user_time, int& tasks)
{
    timespec begin = get_time();

    task_queue::task_t task;
    while (taskQueue->dequeue (task, loc)) {
        tasks++;

        typename container_type::iterator i = container.begin(task.data);

        timespec user_begin = get_time();
        K key;
        reduce_iterator values;

        while (i.next(key, values))
        {
            if (values.size() > 0)
                static_cast<Impl *>(this)->reduce(
                    key, values, this->final_vals[loc.thread]);
        }
        user_time += time_elapsed(user_begin);
    }

    time += time_elapsed(begin);
}

/**
 * Merge all reduced data
 */
template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::run_merge ()
{
    size_t total = 0;
    for (size_t i = 0; i < num_threads; i++) {
        total += this->final_vals[i].size();
    }

    std::vector<keyval>* final = new std::vector<keyval>[1];
    final[0].reserve(total);

    for (size_t i = 0; i < num_threads; i++) {
        final[0].insert(final[0].end(), this->final_vals[i].begin(),
                        this->final_vals[i].end());
    }
    delete [] this->final_vals;
    this->final_vals = final;

    // argo's merge
    // number of keyval pairs
    int *result_count;
    result_count = argo::conew_array<int>(argo::number_of_nodes());
    result_count[argo::node_id()] = this->final_vals->size();
    argo::barrier();
    // allocate argo data
    int sum = 0;
    int start = 0;
    keyval *argo_result_tmp;
    for (int i = 0; i < argo::number_of_nodes(); i++) {
        sum += result_count[i];
        if (i < argo::node_id()) start += result_count[i];
    }
    argo_result_tmp = argo::conew_array<keyval>(sum);
    // echo node move its own result to argo memory
    std::memcpy(&argo_result_tmp[start], this->final_vals->data(), sizeof(keyval) * this->final_vals->size());
    argo::barrier();
    // merger results from all nodes using a map container
    if (0 == argo::node_id()) {
        std::map<K, keyval> mymap;
        std::pair<typename std::map<K, keyval>::iterator, bool> ret;
        for (int i = 0; i < sum; i++) {
            ret = mymap.insert(std::pair<K, keyval>(argo_result_tmp[i].key, argo_result_tmp[i]));
            if (ret.second == false) {
                argo_result_tmp[i].val += ret.first->second.val;
                mymap.erase(ret.first);
                mymap.insert(std::pair<K, keyval>(argo_result_tmp[i].key, argo_result_tmp[i]));
            }
        }
        this->final_vals->clear();

        // copy merged result to a vector container
        for (auto& iter : mymap) {
            this->final_vals->push_back(iter.second);
        }
    }
    argo::barrier();
    argo::codelete_array(result_count);
    argo::codelete_array(argo_result_tmp);

}

template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::
merge_worker (thread_loc const & loc, double & time, double & user_time, int& tasks)
{
    // do nothing at all unless it turns out to be a bottleneck to merge in serial.
}

template<typename Impl, typename D, typename K, typename V, class Container>
void MapReduce<Impl, D, K, V, Container>::start_workers (void (*func)(void*, thread_loc const&), int num_threads, char const * stage)
{
    thread_arg_t* th_arg_array = new thread_arg_t[num_threads];
    thread_arg_t** th_arg_ptrarray = new thread_arg_t*[num_threads];

    thread_arg_t args = { this, 0, 0, 0 };
    for (int thread = 0; thread < num_threads; ++thread)
    {
        th_arg_array[thread] = args;
        th_arg_ptrarray[thread] = &(th_arg_array[thread]);
    }

    CHECK_ERROR (threadPool->set(func, (void **)th_arg_ptrarray, num_threads));
    // Start worker threads
    CHECK_ERROR (threadPool->begin());
    dprintf("Status: All %d threads have been created\n", num_threads);
    // Barrier, wait for all threads to finish.
    CHECK_ERROR (threadPool->wait());

#ifdef TIMING
    double user_time = 0, work_time = 0, max_user_time = 0,
           min_user_time = std::numeric_limits<double>::max(), max_work_time = 0,
           min_work_time = std::numeric_limits<double>::max();
    for (int thread = 0; thread < num_threads; ++thread)
    {
        dprintf("Thread %d: ran %d in %.3f\n", thread,
                th_arg_array[thread].tasks, th_arg_array[thread].time);
        user_time += th_arg_array[thread].user_time;
        min_user_time = std::min(min_user_time, th_arg_array[thread].user_time);
        max_user_time = std::max(max_user_time, th_arg_array[thread].user_time);
        work_time += th_arg_array[thread].time;
        min_work_time = std::min(min_work_time, th_arg_array[thread].time);
        max_work_time = std::max(max_work_time, th_arg_array[thread].time);
    }
    if (max_user_time > 0)
        fprintf (stderr, "%s avg user time: %.3f    (%.3f, %.3f)\n",
                 stage, user_time / num_threads, min_user_time, max_user_time);
    if (max_work_time > 0)
        fprintf (stderr, "%s avg thread time: %.3f    (%.3f, %.3f)\n",
                 stage, work_time / num_threads, min_work_time, max_work_time);
#endif

    delete [] th_arg_ptrarray;
    delete [] th_arg_array;

    dprintf("Status: All tasks have completed\n");
}

template<typename Impl, typename D, typename K, typename V,
         class Container = hash_container<K, V, buffer_combiner> >
class MapReduceSort : public MapReduce<Impl, D, K, V, Container>
{
public:
    typedef typename MapReduce<Impl, D, K, V, Container>::keyval keyval;

protected:

    // default sorting order is by key. User can override.
    bool sort(keyval const& a, keyval const& b) const { return a.key < b.key; }

    struct sort_functor {
        MapReduceSort* mrs;
        sort_functor(MapReduceSort* mrs) : mrs(mrs) {}
        bool operator()(keyval const& a, keyval const& b) const {
            return static_cast<Impl const*>(mrs)->sort(a, b);
        }
    };

    virtual void run_merge ()
    {
        // how many lists to merge in a single task.
        static const int merge_factor = 2;
        int merge_queues = this->num_threads;

        // First sort each queue in place
        for (int i = 0; i < merge_queues; i++)
        {
            task_queue::task_t task =
            { i, 0, (uint64_t)&this->final_vals[i], 0 };
            this->taskQueue->enqueue_seq(task, merge_queues);
        }
        this->start_workers(&this->merge_callback, this->num_threads, "merge");

        // Then merge
        std::vector<keyval>* merge_vals;
        while (merge_queues > 1) {
            uint64_t resulting_queues =
                (uint64_t)std::ceil(merge_queues / (double)merge_factor);

            // swap queues
            merge_vals = this->final_vals;
            this->final_vals = new std::vector<keyval>[resulting_queues];

            // distribute tasks into task queues using locality information
            // if provided.
            int queue_index = 0;
            for (uint64_t i = 0; i < resulting_queues; i++)
            {
                int actual = std::min(merge_factor, merge_queues - queue_index);
                task_queue::task_t task
                = { i,
                    (uint64_t)actual,
                    (uint64_t)&merge_vals[queue_index],
                    0
                  };
                int lgrp = loc_mem_to_lgrp (&merge_vals[queue_index][0]);
                // For debugging, normally this is padding.
                task.pad = lgrp;
                this->taskQueue->enqueue_seq (task, resulting_queues, lgrp);
                queue_index += actual;
            }

            // Run merge tasks and get merge values.
            this->start_workers (&this->merge_callback,
                                 std::min(resulting_queues, this->num_threads), "merge");

            delete [] merge_vals;
            merge_queues = resulting_queues;
        }

        assert(merge_queues == 1);

        // argo's merge
        // number of keyval pairs
        int *result_count;
        result_count = argo::conew_array<int>(argo::number_of_nodes());
        result_count[argo::node_id()] = this->final_vals->size();
        argo::barrier();
        // allocate argo data
        int sum = 0;
        int start = 0;
        keyval *argo_result_tmp;
        for (int i = 0; i < argo::number_of_nodes(); i++) {
            sum += result_count[i];
            if (i < argo::node_id()) start += result_count[i];
        }
        argo_result_tmp = argo::conew_array<keyval>(sum);
        // echo node move its own result to argo memory
        std::memcpy(&argo_result_tmp[start], this->final_vals->data(), sizeof(keyval) * this->final_vals->size());
        argo::barrier();
        // merger results from all nodes using a map container
        if (0 == argo::node_id()) {
	    std::map<K, keyval> mymap;
            std::pair<typename std::map<K, keyval>::iterator, bool> ret;
            for (int i = 0; i < sum; i++) {
                ret = mymap.insert(std::pair<K, keyval>(argo_result_tmp[i].key, argo_result_tmp[i]));
                if (ret.second == false) {
                    argo_result_tmp[i].val += ret.first->second.val;
                    mymap.erase(ret.first);
                    mymap.insert(std::pair<K, keyval>(argo_result_tmp[i].key, argo_result_tmp[i]));
                }
            }
            this->final_vals->clear();

            // copy merged result to a vector container
            for (auto& iter : mymap) {
                this->final_vals->push_back(iter.second);
            }
            // sort merged result
            std::stable_sort(this->final_vals->begin(), this->final_vals->end(), sort_functor(this));
        }
        argo::barrier();
        argo::codelete_array(result_count);
        argo::codelete_array(argo_result_tmp);

    }

    virtual void merge_worker (thread_loc const& loc, double& time,
                               double& user_time, int& tasks)
    {
        timespec begin = get_time();
        task_queue::task_t task;
        while (this->taskQueue->dequeue (task, loc)) {
            tasks++;
            std::vector<keyval>* vals = (std::vector<keyval>*)task.data;
            uint64_t length = task.len;
            uint64_t out_index = task.id;

            if (length == 0)
            {
                // this case really just means sort my list in place.
                // stable_sort ensures that the order of same keyvals with
                // the same key emitted in reduce remains the same in sort
                std::stable_sort(vals->begin(), vals->end(), sort_functor(this));
            }
            else if (length == 1)
            {
                // if only one list, don't merge, just move over.
                (*vals).swap(this->final_vals[out_index]);
            }
            else if (length == 2)
            {
                // stl merge is nice and fast for 2.
                this->final_vals[out_index].resize(vals[0].size() + vals[1].size());
                std::merge(vals[0].begin(), vals[0].end(),
                           vals[1].begin(), vals[1].end(),
                           this->final_vals[out_index].begin(), sort_functor(this));
            }
            else
            {
                // for more, do a multiway merge.
                assert(0);
            }
        }
        time += time_elapsed(begin);
    }
};

#endif // MAP_REDUCE_H_

// vim: ts=8 sw=4 sts=4 smarttab smartindent

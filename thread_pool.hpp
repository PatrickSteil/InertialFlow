#pragma once

#include <condition_variable>
#include <cstddef>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

// A small fixed-size work-queue thread pool that supports fork-join style
// use: a task may itself submit further tasks (e.g. a recursive
// algorithm's "now spawn my children" step), and wait_idle() blocks until
// the entire resulting task graph -- not just the tasks submitted before
// wait_idle() was called -- has finished running.
//
// Each task is a std::function<void(size_t)>; the size_t argument is the
// id (0..size()-1) of the worker thread executing it, stable for that
// worker's whole lifetime. It's meant to index into per-worker scratch
// state (e.g. Partitioner::FlowWorkspace) so concurrently running tasks
// never contend over the same scratch buffers, even though which task
// lands on which worker is otherwise unpredictable (any idle worker can
// pick up any queued task).
class ThreadPool {
 public:
  // Throws std::invalid_argument if num_threads < 1.
  explicit ThreadPool(std::size_t num_threads);

  // Stops and joins every worker. Any tasks still queued at this point are
  // simply dropped; callers that care about every submitted task having
  // run should call wait_idle() first.
  ~ThreadPool();

  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

  // Enqueues a task. Safe to call both from outside the pool (e.g. the
  // owner kicking off the first task) and from inside a running task (to
  // spawn more work) -- in the latter case the newly submitted task counts
  // towards the same outstanding-work total that wait_idle() waits on, so
  // a task spawning children before it returns can never cause wait_idle()
  // to wake up prematurely.
  void submit(std::function<void(std::size_t)> task);

  // Blocks the calling thread until every task submitted so far --
  // including tasks submitted by other tasks while this call is blocked --
  // has finished running. Must be called from outside the pool (i.e. not
  // from within a task running on one of this pool's own workers, which
  // would deadlock).
  void wait_idle();

  std::size_t size() const { return workers.size(); }

 private:
  void worker_loop(std::size_t id);

  std::vector<std::thread> workers;
  std::deque<std::function<void(std::size_t)>> queue;
  std::mutex mtx;
  std::condition_variable cv_task;  // a task was pushed, or stop was set
  std::condition_variable cv_idle;  // pending dropped to zero
  // Number of tasks submitted but not yet finished running (queued or
  // currently executing). wait_idle() blocks until this reaches zero.
  std::size_t pending = 0;
  bool stop = false;
};

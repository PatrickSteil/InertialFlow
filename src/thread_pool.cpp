#include "thread_pool.hpp"

#include <stdexcept>

ThreadPool::ThreadPool(std::size_t num_threads) {
  if (num_threads < 1) {
    throw std::invalid_argument("ThreadPool needs at least 1 thread");
  }
  workers.reserve(num_threads);
  for (std::size_t id = 0; id < num_threads; ++id) {
    workers.emplace_back([this, id]() { worker_loop(id); });
  }
}

ThreadPool::~ThreadPool() {
  {
    std::lock_guard<std::mutex> lock(mtx);
    stop = true;
  }
  cv_task.notify_all();
  for (auto &t : workers) {
    if (t.joinable())
      t.join();
  }
}

void ThreadPool::submit(std::function<void(std::size_t)> task) {
  {
    std::lock_guard<std::mutex> lock(mtx);
    queue.push_back(std::move(task));
    ++pending;
  }
  cv_task.notify_one();
}

void ThreadPool::wait_idle() {
  std::unique_lock<std::mutex> lock(mtx);
  cv_idle.wait(lock, [&] { return pending == 0; });
}

void ThreadPool::worker_loop(std::size_t id) {
  while (true) {
    std::function<void(std::size_t)> task;
    {
      std::unique_lock<std::mutex> lock(mtx);
      cv_task.wait(lock, [&] { return stop || !queue.empty(); });
      if (stop)
        return;
      task = std::move(queue.front());
      queue.pop_front();
    }

    task(id);

    {
      std::lock_guard<std::mutex> lock(mtx);
      --pending;
      if (pending == 0)
        cv_idle.notify_all();
    }
  }
}

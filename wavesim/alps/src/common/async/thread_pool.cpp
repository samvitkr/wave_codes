//
// Created by xuanx004 on 3/4/24.
//

#include "thread_pool.h"

#include <common/base/logging.h>

namespace alps::async {

thread_pool::thread_pool(const concurrency_t   thread_count_,
                         std::function<void()> init_func)
  : thread_count(thread_count_.value() != 0u
                   ? thread_count_
                   : concurrency_t{std::thread::hardware_concurrency()})
  , threads(std::make_unique<std::thread[]>(thread_count.value()))
  , init_func_(std::move(init_func))
  , logger(alps::get_logger("thread_pool"))
{
  logger->debug("Creating thread pool with {} threads", thread_count.value());
  create_threads();
}

thread_pool::~thread_pool()
{
  logger->trace("Destroying thread pool");
  wait_for_tasks();
  destroy_threads();
}

void thread_pool::reset(std::function<void()> init_func)
{
  reset(thread_count, std::move(init_func));
}

void thread_pool::reset(const concurrency_t   thread_count_,
                        std::function<void()> init_func)
{
  const bool was_paused = paused;
  paused                = true;
  wait_for_tasks();
  destroy_threads();
  thread_count = thread_count_.value() != 0 ? thread_count_ : thread_count;
  threads      = std::make_unique<std::thread[]>(thread_count.value());
  paused       = was_paused;
  init_func_   = std::move(init_func);

  logger->debug("Resetting thread pool with {} threads", thread_count.value());
  create_threads();
}

void thread_pool::wait_for_tasks()
{
  waiting = true;
  std::unique_lock<std::mutex> tasks_lock(tasks_mutex);
  task_done_cv.wait(tasks_lock, [this] {
    return (tasks_total == (paused ? tasks.size() : 0));
  });
  waiting = false;
}

void thread_pool::create_threads()
{
  running = true;
  for (std::size_t i = 0; i < thread_count.value(); ++i) {
    threads[i] = std::thread(&thread_pool::worker, this);
  }
}

void thread_pool::destroy_threads()
{
  running = false;
  task_available_cv.notify_all();
  for (std::size_t i = 0; i < thread_count.value(); ++i) {
    threads[i].join();
  }
}

void thread_pool::worker()
{
  init_func_();
  while (running) {
    std::unique_lock<std::mutex> tasks_lock(tasks_mutex);
    task_available_cv.wait(tasks_lock,
                           [&] { return !tasks.empty() || !running; });
    if (running && !paused) {
      auto task = std::move(tasks.front());
      tasks.pop();
      tasks_lock.unlock();
      task();
      --tasks_total;
      if (waiting) task_done_cv.notify_one();
    }
  }
}

} // namespace alps::async

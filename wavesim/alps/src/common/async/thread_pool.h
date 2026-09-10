#pragma once

/**
 * @file thread_pool.h
 * @author Barak Shoshany (baraksh@gmail.com) (http://baraksh.com)
 * @version 3.0.0
 * @date 2022-05-30
 * @copyright Copyright (c) 2022 Barak Shoshany. Licensed under the MIT license.
 * If you use this library in software of any kind, please provide a link to the
 * GitHub repository https://github.com/bshoshany/thread-pool in the source code
 * and documentation. If you use this library in published research, please cite
 * it as follows: Barak Shoshany, "A C++17 Thread Pool for High-Performance
 * Scientific Computing", doi:10.5281/zenodo.4742687, arXiv:2105.00613 (May
 * 2021)
 *
 * @brief alps::thread_pool: a fast, lightweight, and easy-to-use C++17 thread
 * pool library.
 */

#include <atomic>             // std::atomic
#include <condition_variable> // std::condition_variable
#include <exception>          // std::current_exception
#include <future>             // std::future, std::promise
#include <memory> // std::make_shared, std::make_unique, std::shared_ptr, std::unique_ptr
#include <mutex>  // std::mutex, std::scoped_lock, std::unique_lock
#include <queue>  // std::queue
#include <thread> // std::thread
#include <type_traits> // std::common_type_t, std::decay_t, std::is_void_v, std::invoke_result_t
#include <utility> // std::move, std::swap
#include <vector>  // std::vector

#include <common/base/logging_fwd.h>

namespace alps::async {
class concurrency_t
{
 private:
  using value_t =
    std::invoke_result_t<decltype(std::thread::hardware_concurrency)>;

 public:
  template<typename IntType>
  concurrency_t(IntType n,
                std::enable_if_t<std::is_unsigned_v<IntType>>* = nullptr)
    : n_{static_cast<value_t>(n)}
  {}

  auto value() const { return n_; }

 private:
  value_t n_{0};
};

/**
 * @brief A helper class to facilitate waiting for and/or getting the results of
 * multiple futures at once.
 */
template<typename T>
class [[nodiscard]] multi_future
{
 public:
  /**
   * @brief Construct a multi_future object with the given number of futures.
   *
   * @param num_futures_ The desired number of futures to store.
   */
  explicit multi_future(const size_t num_futures_ = 0)
    : futures(num_futures_)
  {}

  /**
   * @brief Get the results from all the futures stored in this multi_future
   * object, rethrowing any stored exceptions.
   *
   * @return If the futures return void, this function returns void as well.
   * Otherwise, it returns a vector containing the results.
   */
  [[nodiscard]] std::conditional_t<std::is_void_v<T>, void, std::vector<T>>
  get()
  {
    if constexpr (std::is_void_v<T>) {
      for (size_t i = 0; i < futures.size(); ++i) {
        futures[i].get();
      }
      return;
    } else {
      std::vector<T> results(futures.size());
      for (size_t i = 0; i < futures.size(); ++i) {
        results[i] = futures[i].get();
      }
      return results;
    }
  }

  /**
   * @brief Get a reference to one of the futures stored in this multi_future
   * object.
   *
   * @param i The index of the desired future.
   * @return The future.
   */
  [[nodiscard]] std::future<T>& operator[](const size_t i)
  {
    return futures[i];
  }

  /**
   * @brief Append a future to this multi_future object.
   *
   * @param future The future to append.
   */
  void push_back(std::future<T> future)
  {
    futures.push_back(std::move(future));
  }

  /**
   * @brief Get the number of futures stored in this multi_future object.
   *
   * @return The number of futures.
   */
  [[nodiscard]] size_t size() const { return futures.size(); }

  /**
   * @brief Wait for all the futures stored in this multi_future object.
   */
  void wait() const
  {
    for (size_t i = 0; i < futures.size(); ++i) {
      futures[i].wait();
    }
  }

 private:
  /**
   * @brief A vector to store the futures.
   */
  std::vector<std::future<T>> futures;
};

/**
 * @brief A fast, lightweight, and easy-to-use C++17 thread pool class.
 */
class [[nodiscard]] thread_pool
{
 public:
  /**
   * @brief Construct a new thread pool.
   *
   * @param thread_count_ Number of threads. The default value is the total
   * number of hardware threads available, as reported by the implementation.
   * @param init_func An optional initialization function to be executed at the
   * creation of the threads; e.g. to set the Cuda device.
   */
  explicit thread_pool(
    concurrency_t         thread_count_ = std::thread::hardware_concurrency(),
    std::function<void()> init_func     = [] {});

  /**
   * @brief Destruct the thread pool.
   */
  ~thread_pool();

  /**
   * @brief Number of tasks currently waiting in the queue
   */
  auto num_tasks_queued() const
  {
    const std::scoped_lock tasks_lock(tasks_mutex);
    return tasks.size();
  }

  /**
   * @brief Number of unfinished tasks: queued or running in a thread.
   */
  auto num_tasks_total() const { return tasks_total.load(); }

  /**
   * @brief Get the number of threads in the pool.
   */
  concurrency_t num_threads() const { return thread_count; }

  /**
   * @brief Push a function with no return value into the task queue.
   */
  template<typename F>
  void push_task(F&& task)
  {
    static_assert(std::is_void_v<std::invoke_result_t<std::decay_t<F>>>,
                  "push_task() can only be used with functions that return "
                  "void.");
    {
      const std::scoped_lock tasks_lock(tasks_mutex);
      tasks.emplace(std::forward<F>(task));
    }
    ++tasks_total;
    task_available_cv.notify_one();
  }

  /**
   * @brief Reset the threads in the pool.
   *
   * Waits for all currently running tasks to be completed, then destroys all
   * threads in the pool and creates a new thread pool. Any tasks that were
   * waiting in the queue before the pool was reset will then be executed by the
   * new threads. If the pool was paused before resetting it, the new pool will
   * be paused as well.
   * An optional initialization function can be provided to
   * be execute at the creation of the threads; e.g. to set Cuda device. The
   * thread count can also be changed by supplying a second argument.
   */
  void reset(std::function<void()> init_func = [] {});

  void reset(
    concurrency_t         thread_count_ = 0u,
    std::function<void()> init_func     = [] {});

  /**
   * @brief Submit a functor with signature f() into the task queue and get a
   * future object.
   */
  template<typename F, typename R = std::invoke_result_t<std::decay_t<F>>>
  [[nodiscard]] std::future<R> submit(F&& func)
  {
    auto task_promise = std::make_shared<std::promise<R>>();
    push_task([task = std::forward<F>(func), promise = task_promise] {
      try {
        if constexpr (std::is_void_v<R>) {
          task();
          promise->set_value();
        } else {
          promise->set_value(task());
        }
      } catch (...) {
        try {
          promise->set_exception(std::current_exception());
        } catch (...) {}
      }
    });
    return task_promise->get_future();
  }

  /**
   * @brief Wait for tasks to be completed.
   */
  void wait_for_tasks();

  /**
   * @brief Pause the pool.
   */
  void pause() { paused = true; }

  /**
   * @brief Unpause the pool.
   */
  void unpause() { paused = false; }

 private:
  void create_threads();

  void destroy_threads();

  void worker();

  concurrency_t thread_count{0u};

  // A smart pointer to manage the memory allocated for the threads.
  std::unique_ptr<std::thread[]> threads = nullptr;

  std::function<void()> init_func_ = [] {};

  // A queue of tasks to be executed by the threads.
  std::queue<std::function<void()>> tasks = {};

  // A condition variable used to notify worker() of a new task
  std::condition_variable task_available_cv = {};

  // A condition variable used to notify wait_for_tasks() that a tasks is done
  std::condition_variable task_done_cv = {};

  mutable std::mutex tasks_mutex = {};

  // An atomic variable to keep track of the total number of unfinished tasks
  std::atomic<size_t> tasks_total = 0;

  alps::Logger logger;

  // An atomic variable indicating that wait_for_tasks() is active and expects
  // to be notified whenever a task is done.
  std::atomic<bool> waiting = false;

  // An atomic variable indicating whether the workers should pause.
  std::atomic<bool> paused = false;

  // An atomic variable indicating to the workers to keep running.
  std::atomic<bool> running = false;
};

} // namespace alps::async

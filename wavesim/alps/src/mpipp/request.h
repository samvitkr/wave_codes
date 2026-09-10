#pragma once

#include <mpipp/config.h>
#include <mpipp/status.h>

#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace mpipp {

/// \brief Represents a non-blocking communication request.
class irequest
{
 protected:
  MPI_Request           req;
  std::function<void()> post_hook;

 public:
  irequest() = delete;

  explicit irequest(MPI_Request mpi_req) noexcept
    : req{mpi_req}
  {}

  ~irequest();

  irequest(irequest&& other) noexcept
    : req{other.req}
    , post_hook{std::move(other.post_hook)}
  {
    other.req       = MPI_REQUEST_NULL;
    other.post_hook = nullptr;
  }

  // disable copy semantics to avoid double free
  irequest(const irequest&)       = delete;
  void operator=(const irequest&) = delete;
  // disable move assignment to avoid implicit MPI_Request_free calls
  irequest& operator=(irequest&& other) = delete;

  void cancel();

  std::optional<status> test();

  status wait();

  std::optional<status> get_status() const;

  friend class irequest_pool;
};

//------------------------------------------------------------------

/// \brief Container for managing a list of non-blocking communication requests.
class irequest_pool
{
 protected:
  std::vector<MPI_Request>           reqs{};
  std::vector<status>                stats{};
  std::vector<std::function<void()>> post_hooks{};

 public:
  using size_type = std::vector<MPI_Request>::size_type;

  irequest_pool() = default;

  irequest_pool(const irequest_pool&) = delete;

  irequest_pool(irequest_pool&& other) noexcept
    : reqs(std::move(other.reqs))
    , stats(std::move(other.stats))
    , post_hooks(std::move(other.post_hooks))
  {}

  ~irequest_pool();

  void operator=(const irequest_pool&) = delete;

  irequest_pool& operator=(irequest_pool&& other) = delete;

  void reserve(size_type count);

  size_type size() const noexcept { return reqs.size(); }

  bool empty() const noexcept { return reqs.empty(); }

  const status& get_status(size_type i) const { return stats.at(i); }

  void cancel(size_type i);

  void cancelall();

  void push(irequest_pool&& other);

  void push(irequest&& other);

  void push(irequest&& other, std::function<void()> post_hook);

  void push(MPI_Request req, std::function<void()> post_hook = nullptr);

  /// Returns the index of the completed request, or std::nullopt if none
  /// completed or all requests are inactive.
  std::optional<size_type> waitany();

  /// Returns the index of the completed request, or std::nullopt if none
  /// completed or all requests are inactive/completed.
  std::optional<size_type> testany();

  void waitall();

  bool testall();

  std::vector<int> waitsome();

  std::vector<int> testsome();
};

} // namespace mpipp

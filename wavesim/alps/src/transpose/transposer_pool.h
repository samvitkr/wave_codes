#pragma once

#include "transposer_base.h"

#include <utility>
#include <vector>

namespace alps::transpose {

// Key type for transposers (dimension sizes)
struct TransposerKey
{
  int      n0;
  int      n1;
  int      np;
  MPI_Comm comm;

  bool operator==(const TransposerKey& rhs) const noexcept
  {
    if (this->n0 != rhs.n0 || this->n1 != rhs.n1 || this->np != rhs.np)
      return false;
    int  result{MPI_UNEQUAL};
    auto error = MPI_Comm_compare(comm, rhs.comm, &result);
    return (error == MPI_SUCCESS) && (result == MPI_IDENT);
  }
};

/// A singleton managing a pool of Transposers
/** Contains Transposers of different transpose sizes but the same number
 * type (T) and the same execution device (ExecSpace).
 **/
template<typename T, typename ExecSpace>
class TransposerPool
{
 private:
  using transposer_t = std::unique_ptr<TransposerBase<T, ExecSpace>>;
  using element_t    = std::pair<TransposerKey, transposer_t>;

  static constexpr TransposerOptions DefaultOptions = []() constexpr {
    TransposerOptions options;
    options.method          = TransposeMethod::Autotune;
    options.tune_tile_sizes = true;
    return options;
  }();

 public:
  /// Return a Singleton instance of TransposerPool
  static TransposerPool& get_instance()
  {
    static TransposerPool instance;
    return instance;
  }

  /// Return a Transposer with the specified dimensions
  const TransposerBase<T, ExecSpace>&
  get_transposer(const mpipp::communicator& comm,
                 int                        n0,
                 int                        n1,
                 int                        nz,
                 TransposerOptions          options = DefaultOptions)
  {
    auto const key = TransposerKey{n0, n1, comm.size(), comm.raw_handle()};
    for (auto const& item : plans) {
      if (item.first == key) return *(item.second);
    }
    return *(emplace(comm, n0, n1, nz, std::move(options)));
  }

  /// Create a Transposer with the specified dimensions
  const transposer_t& emplace(const mpipp::communicator& comm,
                              int                        n0,
                              int                        n1,
                              int                        max_nz,
                              TransposerOptions          options);

  ~TransposerPool();

  // Delete copy and move constructors and assignment operators
  TransposerPool(const TransposerPool&)            = delete;
  TransposerPool(TransposerPool&&)                 = delete;
  TransposerPool& operator=(const TransposerPool&) = delete;
  TransposerPool& operator=(TransposerPool&&)      = delete;

 private:
  TransposerPool();

  std::vector<element_t> plans{};
};

} // namespace alps::transpose

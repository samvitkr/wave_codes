#pragma once

#include <atomic>
#include <vector>

#include <mpipp/config.h>
#include <mpipp/constants.h>

namespace mpipp {

//--------------------------------------------------------------------
/**
 * @brief A reference-counted wrapper of a MPI communicator
 *
 * @note Reference counting is disabled for MPI_COMM_WORLD, MPI_COMM_SELF,
 * and MPI_COMM_NULL.
 * @note Reference counting is only done on the host.
 */
class communicator
{
 private:
  using ref_count_value_t = long;

 public:
  /**
   * @brief Create a communicator by assuming the ownership of an existing MPI
   * communicator.
   */
  [[nodiscard]] static communicator create_by_take_over(MPI_Comm comm);

  /**
   * @brief Create a communicator by duplicating an existing MPI communicator.
   */
  [[nodiscard]] static communicator create_by_duplicate(MPI_Comm comm);

  /// @brief Create a communicator from an existing MPI communicator but do not
  /// manage its lifetime
  explicit communicator(MPI_Comm mpi_comm) noexcept
    : comm_{mpi_comm}
  {}

  communicator() = default;

  communicator(const communicator& other) noexcept;

  communicator(communicator&& other) noexcept;

  communicator& operator=(const communicator& other) noexcept;

  communicator& operator=(communicator&& other) noexcept;

  [[nodiscard]] int size() const noexcept
  {
    int result{};
    MPI_Comm_size(comm_, &result);
    return result;
  }

  [[nodiscard]] int rank() const noexcept
  {
    int result{};
    MPI_Comm_rank(comm_, &result);
    return result;
  }

  void barrier() const noexcept { MPI_Barrier(comm_); }

  /// @brief Duplicate the communicator
  [[nodiscard]] communicator duplicate() const;

  /// @brief Split the communicator
  [[nodiscard]] communicator split(int color, int key = 0) const;

  /// @brief Split the communicator based on shared memory
  [[nodiscard]] communicator split_shared(int key = 0) const;

  bool operator==(const communicator& other) const noexcept;

  bool operator!=(const communicator& other) const noexcept;

  [[nodiscard]] equality_type compare(const communicator& other) const noexcept;

  [[nodiscard]] equality_type compare(MPI_Comm other) const noexcept;

  [[nodiscard]] bool is_valid() const noexcept
  {
    return comm_ != MPI_COMM_NULL;
  }

  [[nodiscard]] auto use_count() const noexcept
  {
    return (this->ref_count_ != nullptr)
           ? ref_count_->load(std::memory_order_relaxed)
           : static_cast<ref_count_value_t>(0);
  }

  void reset() noexcept;

  [[nodiscard]] MPI_Comm raw_handle() const noexcept { return comm_; }

  ~communicator() noexcept { reset(); }

 private:
  communicator(MPI_Comm comm, ref_count_value_t initial_ref_count);

  void increment_ref_count() const noexcept;

  ref_count_value_t decrement_ref_count() const noexcept;

  bool is_mpi_const_comm() const noexcept
  {
    return comm_ == MPI_COMM_WORLD || comm_ == MPI_COMM_SELF
        || comm_ == MPI_COMM_NULL;
  }

  void free_comm();

  // mutable because the underlying object may be changed by the MPI library
  MPI_Comm comm_{MPI_COMM_NULL};

  std::atomic<ref_count_value_t>* ref_count_{nullptr};
};

inline communicator COMM_WORLD() noexcept
{
  return communicator{MPI_COMM_WORLD};
}

inline communicator COMM_SELF() noexcept
{
  return communicator{MPI_COMM_SELF};
}

/**
 * @brief Translate ranks from one communicator to another.
 *
 * @param from The source communicator.
 * @param ranks The ranks in the source communicator to translate.
 * @param to The target communicator.
 * @return The ranks in the target communicator. MPI_UNDEFINED is used for
 * ranks that are not in the target communicator.
 */
std::vector<int> translate_ranks(const communicator&     from,
                                 const std::vector<int>& ranks,
                                 const communicator&     to);

} // namespace mpipp

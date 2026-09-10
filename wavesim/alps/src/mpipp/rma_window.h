#pragma once

#include <atomic>
#include <cstddef>
#include <type_traits>

#include <nonstd/span.hpp>

#include <mpipp/comm.h>
#include <mpipp/config.h>
#include <mpipp/datatype.h>
#include <mpipp/environment.h>
#include <mpipp/operator.h>
#include <mpipp/utility.h>

namespace mpipp {

enum class lock_type
{
  exclusive,
  shared
};

//--------------------------------------------------------------------
// Forward declarations
template<typename T>
class rma_window;

//--------------------------------------------------------------------
/**
 * @brief Base class: untyped RMA window wrapper (rma_window<void>)
 *
 * This is the base implementation containing all the core functionality
 * for MPI window management. The typed version rma_window<T> inherits
 * from this and adds type-safe interfaces.
 *
 * @note Reference counting is disabled for MPI_WIN_NULL.
 * @note Reference counting is only done on the host.
 */
template<>
class rma_window<void>
{
 protected:
  using ref_count_value_t = long;

 public:
  /**
   * @brief Create a window by assuming the ownership of an existing MPI window.
   */
  [[nodiscard]] static rma_window<void> create_by_take_over(MPI_Win win);

  /**
   * @brief Create a window using MPI_Win_create (low-level interface)
   *
   * @param base Base pointer to the memory region
   * @param size Size of the memory region in bytes
   * @param disp_unit Displacement unit (typically sizeof(element))
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  rma_window(void*               base,
             std::size_t         size,
             int                 disp_unit,
             const communicator& comm,
             MPI_Info            info = MPI_INFO_NULL);

  /**
   * @brief Create a window using MPI_Win_create with span (untyped)
   *
   * @tparam U Element type of the span
   * @tparam N Extent of the span
   * @param data Span of data to expose
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  template<typename U, std::size_t N>
  rma_window(nonstd::span<U, N>  data,
             const communicator& comm,
             MPI_Info            info = MPI_INFO_NULL)
    : rma_window(data.data(), data.size_bytes(), sizeof(U), comm, info)
  {}

  /**
   * @brief Allocate a window using MPI_Win_allocate (untyped)
   *
   * @param size Size of the memory region in bytes
   * @param disp_unit Displacement unit (typically sizeof(element))
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  rma_window(std::size_t         size,
             int                 disp_unit,
             const communicator& comm,
             MPI_Info            info = MPI_INFO_NULL);

  /**
   * @brief Create a window from an existing MPI window but do not manage its
   * lifetime
   */
  explicit rma_window(MPI_Win win) noexcept;

  rma_window() = default;
  rma_window(const rma_window& other) noexcept;
  rma_window(rma_window&& other) noexcept;
  rma_window& operator=(const rma_window& other) noexcept;
  rma_window& operator=(rma_window&& other) noexcept;

  [[nodiscard]] bool is_valid() const noexcept { return win_ != MPI_WIN_NULL; }

  [[nodiscard]] MPI_Win raw_handle() const noexcept { return win_; }

  [[nodiscard]] auto use_count() const noexcept
  {
    return (ref_count_ != nullptr) ? ref_count_->load(std::memory_order_relaxed)
                                   : static_cast<ref_count_value_t>(0);
  }

  [[nodiscard]] void* base_address() const noexcept { return base_ptr_; }

  [[nodiscard]] std::size_t size_bytes() const noexcept { return size_; }

  [[nodiscard]] int disp_unit() const noexcept { return disp_unit_; }

  /**
   * @brief Access local memory as a typed span (untyped window, specify type)
   *
   * @tparam U Element type to cast to
   * @return nonstd::span<U> Span of local memory
   */
  template<typename U>
  [[nodiscard]] nonstd::span<U> local_memory() const
  {
    if ((int)sizeof(U) != disp_unit_) {
      throw std::runtime_error("Element size does not match displacement unit");
    }
    return nonstd::span<U>{static_cast<U*>(base_ptr_), size_ / sizeof(U)};
  }

  // Epoch management
  void fence(int assert = 0) const noexcept { MPI_Win_fence(assert, win_); }
  void lock(lock_type type, int rank, int assert = 0) const noexcept
  {
    int mpi_lock_type =
      (type == lock_type::exclusive) ? MPI_LOCK_EXCLUSIVE : MPI_LOCK_SHARED;
    MPI_Win_lock(mpi_lock_type, rank, assert, win_);
  }
  void unlock(int rank) const noexcept { MPI_Win_unlock(rank, win_); }
  void lock_all(int assert = 0) const noexcept
  {
    MPI_Win_lock_all(assert, win_);
  }
  void unlock_all() const noexcept { MPI_Win_unlock_all(win_); }

  void reset() noexcept;
  ~rma_window() noexcept;

 protected:
  rma_window(MPI_Win win, std::atomic<ref_count_value_t>* ref_count) noexcept;
  void              increment_ref_count() const noexcept;
  ref_count_value_t decrement_ref_count() const noexcept;
  void              free_win() noexcept;

  template<typename U>
  void check_disp_unit() const
  {
    if (is_valid() && (int)sizeof(U) != disp_unit_) {
      throw std::runtime_error("Element size does not match displacement unit");
    }
  }

  MPI_Win win_{MPI_WIN_NULL};

  void*       base_ptr_{nullptr};
  std::size_t size_{0};
  int         disp_unit_{1};

  std::atomic<ref_count_value_t>* ref_count_{nullptr};
};

//--------------------------------------------------------------------
/**
 * @brief Typed RMA window wrapper
 *
 * @tparam T The type of elements in the window
 *
 * Inherits from rma_window<void> and adds type-safe interfaces for
 * accessing local memory and creating typed windows.
 */
template<typename T>
class rma_window : public rma_window<void>
{
 public:
  /**
   * @brief Create a window by assuming the ownership of an existing MPI window.
   */
  [[nodiscard]] static rma_window<T> create_by_take_over(MPI_Win win)
  {
    return rma_window<T>{rma_window<void>::create_by_take_over(win)};
  }

  /**
   * @brief Create a window using MPI_Win_create with a memory span
   *
   * @tparam N Extent of the span
   * @param data Span of data to expose
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  template<std::size_t N>
  rma_window(nonstd::span<T, N>  data,
             const communicator& comm,
             MPI_Info            info = MPI_INFO_NULL)
    : rma_window<void>(data.data(), data.size_bytes(), sizeof(T), comm, info)
  {}

  /**
   * @brief Allocate a window using MPI_Win_allocate (typed)
   *
   * @param count Number of elements to allocate
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  rma_window(std::size_t         count,
             const communicator& comm,
             MPI_Info            info = MPI_INFO_NULL)
    : rma_window<void>(count * sizeof(T), sizeof(T), comm, info)
  {}

  explicit rma_window(MPI_Win win)
    : rma_window<void>(win)
  {
    check_disp_unit<T>();
  }

  /**
   * @brief Access local memory as a typed span (typed window)
   *
   * @return nonstd::span<T> Span of local memory
   */
  [[nodiscard]] nonstd::span<T> local_memory() const noexcept
  {
    return nonstd::span<T>{static_cast<T*>(this->base_ptr_),
                           this->size_ / sizeof(T)};
  }

 protected:
  explicit rma_window(rma_window<void> other)
    : rma_window<void>(std::move(other))
  {
    check_disp_unit<T>();
  }
};

//--------------------------------------------------------------------
// RMA Operations (free functions)
//--------------------------------------------------------------------

/**
 * @brief Put a single element to a remote window
 *
 * @tparam T Element type (must be an MPI type)
 * @tparam W Window type (must be T or void)
 * @param origin Origin data to send
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param win RMA window
 */
template<
  typename T,
  typename W,
  std::enable_if_t<is_mpi_type_v<T>
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void put(const T&             origin,
         int                  target_rank,
         std::size_t          target_disp,
         const rma_window<W>& win)
{
  MPI_Put(&origin,
          1,
          get_type<T>(),
          target_rank,
          size_t_to_mpi_aint(target_disp),
          1,
          get_type<T>(),
          win.raw_handle());
}

/**
 * @brief Get a single element from a remote window
 *
 * @tparam T Element type (must be an MPI type)
 * @tparam W Window type (must be T or void)
 * @param result Result buffer
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param win RMA window
 */
template<
  typename T,
  typename W,
  std::enable_if_t<(!std::is_const_v<T> && is_mpi_type_v<T>)
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void get(T&                   result,
         int                  target_rank,
         std::size_t          target_disp,
         const rma_window<W>& win)
{
  MPI_Get(&result,
          1,
          get_type<T>(),
          target_rank,
          size_t_to_mpi_aint(target_disp),
          1,
          get_type<T>(),
          win.raw_handle());
}

/**
 * @brief Put a span of elements to a remote window
 *
 * @tparam T Element type (must be an MPI type)
 * @tparam N Extent of the span
 * @tparam W Window type (must be T or void)
 * @param origin Origin data to send
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param win RMA window
 */
template<typename T,
         std::size_t N,
         typename W,
         std::enable_if_t<is_mpi_type_v<std::remove_const_t<T>>
                            && (std::is_same_v<W, std::remove_const_t<T>>
                                || std::is_same_v<W, void>),
                          int> = 0>
void put(nonstd::span<T, N>   origin,
         int                  target_rank,
         std::size_t          target_disp,
         const rma_window<W>& win)
{
  using value_type = std::remove_const_t<T>;
  MPI_Put(origin.data(),
          static_cast<int>(origin.size()),
          get_type<value_type>(),
          target_rank,
          size_t_to_mpi_aint(target_disp),
          static_cast<int>(origin.size()),
          get_type<value_type>(),
          win.raw_handle());
}

/**
 * @brief Get a span of elements from a remote window
 *
 * @tparam T Element type (must be an MPI type, non-const)
 * @tparam N Extent of the span
 * @tparam W Window type (must be T or void)
 * @param result Result buffer
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param win RMA window
 */
template<
  typename T,
  std::size_t N,
  typename W,
  std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void get(nonstd::span<T, N>   result,
         int                  target_rank,
         std::size_t          target_disp,
         const rma_window<W>& win)
{
  MPI_Get(result.data(),
          static_cast<int>(result.size()),
          get_type<T>(),
          target_rank,
          size_t_to_mpi_aint(target_disp),
          static_cast<int>(result.size()),
          get_type<T>(),
          win.raw_handle());
}

/**
 * @brief Accumulate a single element to a remote window
 *
 * @tparam T Element type (must be an MPI type)
 * @tparam Op Operation type (predefined op from mpipp/operator.h)
 * @tparam W Window type (must be T or void)
 * @param origin Origin data to accumulate
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param op_obj MPI operation
 * @param win RMA window
 */
template<
  typename T,
  typename Op,
  typename W,
  std::enable_if_t<is_mpi_type_v<T>
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void accumulate(const T&             origin,
                int                  target_rank,
                std::size_t          target_disp,
                const Op&            op_obj,
                const rma_window<W>& win)
{
  static_assert(is_accumulate_compatible_v<T, Op>,
                "Data type or operation cannot be used with accumulate.");
  MPI_Accumulate(&origin,
                 1,
                 get_type<T>(),
                 target_rank,
                 size_t_to_mpi_aint(target_disp),
                 1,
                 get_type<T>(),
                 op_obj.raw_handle(),
                 win.raw_handle());
}

/**
 * @brief Get and accumulate a single element
 *
 * @tparam T Element type (must be an MPI type, non-const)
 * @tparam Op Operation type (predefined op from mpipp/operator.h)
 * @tparam W Window type (must be T or void)
 * @param origin Origin data to accumulate
 * @param result Result buffer for the original value
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param op_obj MPI operation
 * @param win RMA window
 */
template<
  typename T,
  typename Op,
  typename W,
  std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void get_accumulate(const T&             origin,
                    T&                   result,
                    int                  target_rank,
                    std::size_t          target_disp,
                    const Op&            op_obj,
                    const rma_window<W>& win)
{
  static_assert(is_get_accumulate_compatible_v<T, Op>,
                "Data type or operation cannot be used with get_accumulate.");
  MPI_Get_accumulate(&origin,
                     1,
                     get_type<T>(),
                     &result,
                     1,
                     get_type<T>(),
                     target_rank,
                     size_t_to_mpi_aint(target_disp),
                     1,
                     get_type<T>(),
                     op_obj.raw_handle(),
                     win.raw_handle());
}

/**
 * @brief Atomic fetch-and-op operation
 *
 * @tparam T Element type (must be an MPI type)
 * @tparam Op Operation type (predefined op from mpipp/operator.h)
 * @tparam W Window type (must be T or void)
 * @param origin Origin data
 * @param result Result buffer for the original value
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param op_obj MPI operation
 * @param win RMA window
 */
template<
  typename T,
  typename Op,
  typename W,
  std::enable_if_t<is_mpi_type_v<T>
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void fetch_and_op(const T&             origin,
                  T&                   result,
                  int                  target_rank,
                  std::size_t          target_disp,
                  const Op&            op_obj,
                  const rma_window<W>& win)
{
  static_assert(is_fetch_and_op_compatible_v<T, Op>,
                "Data type or operation cannot be used with fetch_and_op.");
  MPI_Fetch_and_op(&origin,
                   &result,
                   get_type<T>(),
                   target_rank,
                   size_t_to_mpi_aint(target_disp),
                   op_obj.raw_handle(),
                   win.raw_handle());
}

/**
 * @brief Atomic compare-and-swap operation
 *
 * @tparam T Element type (must be an MPI type)
 * @tparam W Window type
 * @param origin Origin data (new value)
 * @param compare Compare value
 * @param result Result buffer for the original value
 * @param target_rank Target rank
 * @param target_disp Target displacement (in units of disp_unit)
 * @param win RMA window
 */
template<
  typename T,
  typename W,
  std::enable_if_t<is_cas_compatible_v<T>
                     && (std::is_same_v<W, T> || std::is_same_v<W, void>),
                   int> = 0>
void compare_and_swap(const T&             origin,
                      const T&             compare,
                      T&                   result,
                      int                  target_rank,
                      std::size_t          target_disp,
                      const rma_window<W>& win)
{
  static_assert(is_cas_compatible_v<T>,
                "Data type cannot be used with compare_and_swap");
  MPI_Compare_and_swap(&origin,
                       &compare,
                       &result,
                       get_type<T>(),
                       target_rank,
                       size_t_to_mpi_aint(target_disp),
                       win.raw_handle());
}
} // namespace mpipp

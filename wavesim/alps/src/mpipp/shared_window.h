#pragma once

#include <cstddef>

#include <nonstd/span.hpp>

#include <mpipp/comm.h>
#include <mpipp/config.h>
#include <mpipp/rma_window.h>

namespace mpipp {

//--------------------------------------------------------------------
// Forward declarations
template<typename T>
class shared_window;

//--------------------------------------------------------------------
/**
 * @brief Base class: untyped shared-memory window wrapper (shared_window<void>)
 *
 * Inherits from rma_window<void> and adds shared_query functionality.
 */
template<>
class shared_window<void> : public rma_window<void>
{
 public:
  /**
   * @brief Create a shared window by assuming ownership of an existing MPI
   * window.
   */
  [[nodiscard]] static shared_window<void> create_by_take_over(MPI_Win win);

  /**
   * @brief Allocate a shared window using MPI_Win_allocate_shared (low-level)
   *
   * @param size Size of the memory region in bytes
   * @param disp_unit Displacement unit (typically sizeof(element))
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  shared_window(std::size_t         size,
                int                 disp_unit,
                const communicator& comm,
                MPI_Info            info = MPI_INFO_NULL);

  explicit shared_window(MPI_Win win);

  /**
   * @brief Query shared memory region of a specific rank (untyped window,
   * specify type)
   *
   * @tparam U Element type to cast to
   * @param rank Rank to query
   * @return nonstd::span<U> Span of shared memory for the specified rank
   */
  template<typename U>
  [[nodiscard]] nonstd::span<U> shared_query(int rank) const
  {
    MPI_Aint size{};
    int      disp_unit{};
    void*    base_ptr{nullptr};
    MPI_Win_shared_query(this->win_, rank, &size, &disp_unit, (void*)&base_ptr);
    if (disp_unit != (int)sizeof(U)) {
      throw std::runtime_error(
        "Displacement unit does not match the size of the specified type");
    }
    return nonstd::span<U>{static_cast<U*>(base_ptr), size_t(size) / sizeof(U)};
  }

  // Inherit local_memory from rma_window
  using rma_window<void>::local_memory;

 protected:
  shared_window(MPI_Win win, std::atomic<ref_count_value_t>* ref_count);
};

//--------------------------------------------------------------------
/**
 * @brief Typed shared-memory window wrapper
 *
 * @tparam T The type of elements in the window
 *
 * Inherits from shared_window<void> and adds type-safe interfaces.
 */
template<typename T>
class shared_window final : public shared_window<void>
{
 public:
  /**
   * @brief Create a shared window by assuming ownership of an existing MPI
   * window.
   */
  [[nodiscard]] static shared_window<T> create_by_take_over(MPI_Win win)
  {
    return shared_window<T>{shared_window<void>::create_by_take_over(win)};
  }

  /**
   * @brief Allocate a shared window using MPI_Win_allocate_shared (typed)
   *
   * @param count Number of elements to allocate
   * @param comm Communicator for the window
   * @param info MPI info object (default: MPI_INFO_NULL)
   */
  shared_window(std::size_t         count,
                const communicator& comm,
                MPI_Info            info = MPI_INFO_NULL)
    : shared_window<void>(count * sizeof(T), sizeof(T), comm, info)
  {}

  explicit shared_window(MPI_Win win)
    : shared_window<void>(win)
  {
    check_disp_unit<T>();
  }

  /**
   * @brief Query shared memory region of a specific rank (typed window)
   *
   * @param rank Rank to query
   * @return nonstd::span<T> Span of shared memory for the specified rank
   */
  [[nodiscard]] nonstd::span<T> shared_query(int rank) const noexcept
  {
    MPI_Aint size{};
    int      disp_unit{};
    void*    base_ptr{nullptr};
    MPI_Win_shared_query(
      this->win_, rank, &size, &disp_unit, (void*)(&base_ptr));
    return nonstd::span<T>{static_cast<T*>(base_ptr), size_t(size) / sizeof(T)};
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

 private:
  explicit shared_window(shared_window<void> base_window)
    : shared_window<void>{std::move(base_window)}
  {
    check_disp_unit<T>();
  }
};
} // namespace mpipp

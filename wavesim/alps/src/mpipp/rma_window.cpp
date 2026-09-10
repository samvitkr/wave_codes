#include "rma_window.h"

#include <stdexcept>

#include <fmt/format.h>

#include <mpipp/utility.h>

namespace mpipp {

namespace {
void* get_window_base(MPI_Win win)
{
  void* ptr{nullptr};
  int   flag{0};
  if (win != MPI_WIN_NULL) {
    MPI_Win_get_attr(win, MPI_WIN_BASE, (void*)(&ptr), &flag);
  }
  return (flag != 0) ? ptr : nullptr;
}

std::size_t get_window_size(MPI_Win win)
{
  MPI_Aint* psize{nullptr};
  int       flag{0};
  if (win != MPI_WIN_NULL) {
    MPI_Win_get_attr(win, MPI_WIN_SIZE, &psize, &flag);
  }
  return (flag != 0 && psize != nullptr) ? (size_t)(*psize) : 0;
}

int get_window_disp_unit(MPI_Win win)
{
  int* pdisp{nullptr};
  int  flag{0};
  if (win != MPI_WIN_NULL) {
    MPI_Win_get_attr(win, MPI_WIN_DISP_UNIT, &pdisp, &flag);
  }
  return (flag != 0 && pdisp != nullptr) ? *pdisp : 1;
}
} // anonymous namespace

rma_window<void> rma_window<void>::create_by_take_over(MPI_Win win)
{
  if (win == MPI_WIN_NULL) {
    return rma_window<void>{win, nullptr}; // no reference counting
  }
  return rma_window<void>{win, new std::atomic<ref_count_value_t>{1}};
}

rma_window<void>::rma_window(void*               base,
                             std::size_t         size,
                             int                 disp_unit,
                             const communicator& comm,
                             MPI_Info            info)
  : win_([=] {
    if (disp_unit <= 0) {
      throw std::invalid_argument(
        "mpipp::rma_window<void>::create: disp_unit must be positive");
    }
    if (size > 0 && base == nullptr) {
      throw std::invalid_argument("mpipp::rma_window<void>::create: non-zero "
                                  "size requires non-null base");
    }

    MPI_Win win{MPI_WIN_NULL};
    auto    result = MPI_Win_create(
      base, size_t_to_mpi_aint(size), disp_unit, info, comm.raw_handle(), &win);
    if (result != MPI_SUCCESS) {
      throw std::runtime_error(
        fmt::format("MPI_Win_create failed with error code {}", result));
    }
    return win;
  }())
  , base_ptr_{base}
  , size_{size}
  , disp_unit_{disp_unit}
  , ref_count_{is_valid() ? new std::atomic<ref_count_value_t>{1} : nullptr}
{}

rma_window<void>::rma_window(std::size_t         size,
                             int                 disp_unit,
                             const communicator& comm,
                             MPI_Info            info)
  : win_([=] {
    if (disp_unit <= 0) {
      throw std::invalid_argument(
        "mpipp::rma_window<void>::allocate: disp_unit must be positive");
    }

    MPI_Win win{MPI_WIN_NULL};
    void*   base_ptr{nullptr};
    auto    result = MPI_Win_allocate(size_t_to_mpi_aint(size),
                                   disp_unit,
                                   info,
                                   comm.raw_handle(),
                                   (void*)(&base_ptr),
                                   &win);
    if (result != MPI_SUCCESS) {
      throw std::runtime_error(
        fmt::format("MPI_Win_allocate failed with error code {}", result));
    }
    return win;
  }())
  , base_ptr_{get_window_base(win_)}
  , size_{size}
  , disp_unit_{disp_unit}
  , ref_count_{is_valid() ? new std::atomic<ref_count_value_t>{1} : nullptr}
{}

rma_window<void>::rma_window(MPI_Win                         win,
                             std::atomic<ref_count_value_t>* ref_count) noexcept
  : win_{win}
  , base_ptr_{get_window_base(win)}
  , size_{get_window_size(win)}
  , disp_unit_{get_window_disp_unit(win)}
  , ref_count_{ref_count}
{}

rma_window<void>::rma_window(MPI_Win win) noexcept
  : rma_window(win, nullptr) // no reference counting
{}

rma_window<void>::rma_window(const rma_window& other) noexcept
  : win_{other.win_}
  , base_ptr_{other.base_ptr_}
  , size_{other.size_}
  , disp_unit_{other.disp_unit_}
  , ref_count_{other.ref_count_}
{
  if (ref_count_ != nullptr) {
    increment_ref_count();
  }
}

rma_window<void>::rma_window(rma_window&& other) noexcept
  : win_{other.win_}
  , base_ptr_{other.base_ptr_}
  , size_{other.size_}
  , disp_unit_{other.disp_unit_}
  , ref_count_{other.ref_count_}
{
  other.win_       = MPI_WIN_NULL;
  other.base_ptr_  = nullptr;
  other.size_      = 0;
  other.disp_unit_ = 1;
  other.ref_count_ = nullptr;
}

rma_window<void>& rma_window<void>::operator=(const rma_window& other) noexcept
{
  if (this != &other) {
    reset();
    win_       = other.win_;
    base_ptr_  = other.base_ptr_;
    size_      = other.size_;
    disp_unit_ = other.disp_unit_;
    ref_count_ = other.ref_count_;
    if (ref_count_ != nullptr) {
      increment_ref_count();
    }
  }
  return *this;
}

rma_window<void>& rma_window<void>::operator=(rma_window&& other) noexcept
{
  if (this != &other) {
    reset();
    win_             = other.win_;
    other.win_       = MPI_WIN_NULL;
    base_ptr_        = other.base_ptr_;
    other.base_ptr_  = nullptr;
    size_            = other.size_;
    other.size_      = 0;
    disp_unit_       = other.disp_unit_;
    other.disp_unit_ = 1;
    ref_count_       = other.ref_count_;
    other.ref_count_ = nullptr;
  }
  return *this;
}

void rma_window<void>::reset() noexcept
{
  if (ref_count_ == nullptr) {
    win_ = MPI_WIN_NULL;
    return;
  }

  auto const count = decrement_ref_count();
  if (count == 0) {
    free_win();
  }
  win_       = MPI_WIN_NULL;
  base_ptr_  = nullptr;
  size_      = 0;
  disp_unit_ = 1;
  ref_count_ = nullptr;
}

void rma_window<void>::free_win() noexcept
{
  if (is_valid()) {
    // The window would already be freed if MPI is finalized
    if (!::mpipp::finalized()) {
      MPI_Win_free(&win_);
    }
  }
  win_ = MPI_WIN_NULL;
  delete ref_count_;
  ref_count_ = nullptr;
}

rma_window<void>::~rma_window() noexcept
{
  reset();
}

void rma_window<void>::increment_ref_count() const noexcept
{
  ref_count_->fetch_add(1, std::memory_order_relaxed);
}

typename rma_window<void>::ref_count_value_t
rma_window<void>::decrement_ref_count() const noexcept
{
  return ref_count_->fetch_sub(1, std::memory_order_acq_rel) - 1;
}
} // namespace mpipp

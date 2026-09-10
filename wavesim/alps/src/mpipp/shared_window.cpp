#include "shared_window.h"

#include <fmt/format.h>

#include <mpipp/utility.h>

namespace mpipp {

namespace {
bool check_shared_window_flavor(MPI_Win win)
{
  int* pflavor{nullptr};
  int  flag{};
  MPI_Win_get_attr(win, MPI_WIN_CREATE_FLAVOR, &pflavor, &flag);
  return flag != 0 && pflavor != nullptr && *pflavor == MPI_WIN_FLAVOR_SHARED;
}
} // namespace

shared_window<void> shared_window<void>::create_by_take_over(MPI_Win win)
{
  if (win == MPI_WIN_NULL) {
    return shared_window<void>{win, nullptr};
  }
  return shared_window<void>{win, new std::atomic<ref_count_value_t>(1)};
}

shared_window<void>::shared_window(std::size_t         size,
                                   int                 disp_unit,
                                   const communicator& comm,
                                   MPI_Info            info)
  : shared_window<void>(
      [=] {
        if (disp_unit <= 0) {
          throw std::invalid_argument(
            "mpipp::shared_window<void>::allocate: disp_unit must be positive");
        }

        MPI_Win win{MPI_WIN_NULL};
        void*   base_ptr{nullptr};
        auto    result = MPI_Win_allocate_shared(size_t_to_mpi_aint(size),
                                              disp_unit,
                                              info,
                                              comm.raw_handle(),
                                              (void*)(&base_ptr),
                                              &win);
        if (result != MPI_SUCCESS) {
          throw std::runtime_error(fmt::format(
            "MPI_Win_allocate_shared failed with error code {}", result));
        }
        return win;
      }(),
      is_valid() ? new std::atomic<ref_count_value_t>(1) : nullptr)
{}

shared_window<void>::shared_window(MPI_Win win)
  : shared_window{win, nullptr}
{}

shared_window<void>::shared_window(MPI_Win                         win,
                                   std::atomic<ref_count_value_t>* ref_count)
  : rma_window<void>{win, ref_count}
{
  if (this->win_ != MPI_WIN_NULL && !check_shared_window_flavor(this->win_)) {
    throw std::runtime_error(
      "MPI_Win must be created with MPI_Win_allocate_shared or have the "
      "MPI_WIN_CREATE_FLAVOR attribute");
  }
}

} // namespace mpipp

#pragma once

#include <mpipp/config.h>
#include <mpipp/datatype.h>

namespace mpipp {

struct status_ignore_t
{};

static constexpr status_ignore_t status_ignore{};

class status : public MPI_Status
{
 public:
  [[nodiscard]] int source() const { return MPI_SOURCE; }

  [[nodiscard]] int tag() const { return MPI_TAG; }

  [[nodiscard]] int error() const { return MPI_ERROR; }

  [[nodiscard]] bool is_cancelled() const
  {
    int result{};
    MPI_Test_cancelled(this, &result);
    return result != 0;
  }

  [[nodiscard]] bool is_canceled() const { return is_cancelled(); }

  template<typename T,
           typename = std::enable_if_t<is_mpi_type_v<std::remove_const_t<T>>>>
  [[nodiscard]] int count() const
  {
    using val_t = std::remove_const_t<T>;
    int result{};
    MPI_Get_count(this, get_type<val_t>(), &result);
    return result;
  }

  status()
    : MPI_Status{}
  {
    MPI_SOURCE = MPI_ANY_SOURCE;
    MPI_TAG    = MPI_ANY_TAG;
    MPI_ERROR  = MPI_SUCCESS;
  }
};

} // namespace mpipp

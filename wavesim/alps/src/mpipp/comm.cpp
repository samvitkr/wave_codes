#include "comm.h"

#include "environment.h"

#include <fmt/format.h>

namespace mpipp {

communicator::communicator(MPI_Comm comm, ref_count_value_t initial_ref_count)
  : comm_{comm}
  , ref_count_{(initial_ref_count < 1 || is_mpi_const_comm())
                 ? nullptr
                 : new std::atomic<ref_count_value_t>{initial_ref_count}}
{}

communicator communicator::create_by_take_over(MPI_Comm comm)
{
  if (comm == MPI_COMM_NULL || comm == MPI_COMM_WORLD
      || comm == MPI_COMM_SELF) {
    return {comm, 0}; // no reference counting
  }
  return {comm, 1};
}

communicator communicator::create_by_duplicate(MPI_Comm comm)
{
  if (comm == MPI_COMM_NULL) {
    return {MPI_COMM_NULL, 0}; // cannot duplicate MPI_COMM_NULL
  }
  MPI_Comm new_comm{MPI_COMM_NULL};
  MPI_Comm_dup(comm, &new_comm);
  return {new_comm, 1};
}

communicator communicator::duplicate() const
{
  if (is_mpi_const_comm()) {
    // the communicator is a constant, such as MPI_COMM_WORLD or MPI_COMM_SELF
    return {this->comm_, 0};
  }
  return create_by_duplicate(this->comm_);
}

communicator communicator::split(int color, int key) const
{
  MPI_Comm new_comm{MPI_COMM_NULL};
  MPI_Comm_split(comm_, color, key, &new_comm);
  return create_by_take_over(new_comm);
}

communicator communicator::split_shared(int key) const
{
  MPI_Comm new_comm{MPI_COMM_NULL};
  MPI_Comm_split_type(
    comm_, MPI_COMM_TYPE_SHARED, key, MPI_INFO_NULL, &new_comm);
  return create_by_take_over(new_comm);
}

communicator::communicator(const communicator& other) noexcept
  : comm_{other.comm_}
  , ref_count_{other.ref_count_}
{
  // Increment the reference count if on the host
  if (ref_count_ != nullptr) {
    increment_ref_count();
  }
}

communicator::communicator(communicator&& other) noexcept
  : comm_{other.comm_}
  , ref_count_{other.ref_count_}
{
  other.comm_      = MPI_COMM_NULL;
  other.ref_count_ = nullptr;
}

communicator& communicator::operator=(const communicator& other) noexcept
{
  if (this != &other) {
    reset();
    comm_      = other.comm_;
    ref_count_ = other.ref_count_;
    // Increment the reference count if on the host
    if (ref_count_ != nullptr) {
      increment_ref_count();
    }
  }
  return *this;
}

communicator& communicator::operator=(communicator&& other) noexcept
{
  if (this != &other) {
    reset();
    comm_            = other.comm_;
    other.comm_      = MPI_COMM_NULL;
    ref_count_       = other.ref_count_;
    other.ref_count_ = nullptr;
  }
  return *this;
}

equality_type communicator::compare(const communicator& other) const noexcept
{
  int result;
  MPI_Comm_compare(comm_, other.comm_, &result);
  return static_cast<equality_type>(result);
}

equality_type communicator::compare(MPI_Comm other) const noexcept
{
  int result;
  MPI_Comm_compare(comm_, other, &result);
  return static_cast<equality_type>(result);
}

void communicator::reset() noexcept
{
  if (ref_count_ == nullptr) {
    // The communicator is not reference counted (not managed by this class)
    comm_ = MPI_COMM_NULL;
    return;
  }

  auto const count = decrement_ref_count();
  if (count == 0) {
    free_comm();
  }
  comm_      = MPI_COMM_NULL;
  ref_count_ = nullptr;
}

void communicator::free_comm()
{
  if (is_valid() && !is_mpi_const_comm()) {
    // The communicator would already be freed if MPI is finalized
    if (!::mpipp::finalized()) {
      MPI_Comm_free(&comm_);
    }
  }
  comm_ = MPI_COMM_NULL;
  delete ref_count_;
  ref_count_ = nullptr;
}

void communicator::increment_ref_count() const noexcept
{
  ref_count_->fetch_add(1, std::memory_order_relaxed);
}

communicator::ref_count_value_t
communicator::decrement_ref_count() const noexcept
{
  return ref_count_->fetch_sub(1, std::memory_order_acq_rel) - 1;
}

bool communicator::operator==(const communicator& other) const noexcept
{
  int result{};
  MPI_Comm_compare(comm_, other.comm_, &result);
  return result == MPI_IDENT;
}

bool communicator::operator!=(const communicator& other) const noexcept
{
  int result{};
  MPI_Comm_compare(comm_, other.comm_, &result);
  return result != MPI_IDENT;
}

std::vector<int> translate_ranks(const communicator&     from,
                                 const std::vector<int>& ranks,
                                 const communicator&     to)
{
  if (ranks.empty()) {
    return {};
  }

#ifndef NDEBUG
  const int from_size = from.size();
  for (const int r : ranks) {
    throw std::out_of_range(fmt::format(
      "Rank {} is out of range for the source communicator of size {}",
      r,
      from_size));
  }
#endif

  MPI_Group group_from{MPI_GROUP_NULL}, group_to{MPI_GROUP_NULL};
  MPI_Comm_group(from.raw_handle(), &group_from);
  MPI_Comm_group(to.raw_handle(), &group_to);

  std::vector<int> result(ranks.size());
  MPI_Group_translate_ranks(group_from,
                            static_cast<int>(ranks.size()),
                            ranks.data(),
                            group_to,
                            result.data());

  MPI_Group_free(&group_from);
  MPI_Group_free(&group_to);

  return result;
}

} // namespace mpipp

#pragma once

#include <mpipp/comm.h>
#include <mpipp/config.h>
#include <mpipp/datatype.h>
#include <mpipp/operator.h>
#include <mpipp/request.h>
#include <mpipp/utility.h>

#include <nonstd/span.hpp>

#include <type_traits>

namespace mpipp {

// === collective ==================================================
// === barrier ===
// --- blocking barrier ---
inline void barrier(const communicator& comm)
{
  MPI_Barrier(comm.raw_handle());
}

// --- nonblocking barrier ---
inline irequest ibarrier(const communicator& comm)
{
  MPI_Request req{MPI_REQUEST_NULL};
  MPI_Ibarrier(comm.raw_handle(), &req);
  return irequest(req);
}

// === broadcast ===
// --- blocking broadcast ---
template<typename T,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
void bcast(T& data, int root, const communicator& comm)
{
  MPI_Bcast(&data, 1, get_type<T>(), root, comm.raw_handle());
}

template<typename T,
         std::size_t N,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
void bcast(nonstd::span<T, N> data, int root, const communicator& comm)
{
  auto count = std::size(data);
  MPI_Bcast(
    data.data(), to_int_size(count), get_type<T>(), root, comm.raw_handle());
}

// --- nonblocking broadcast ---
template<typename T,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
irequest ibcast(T& data, int root, const communicator& comm)
{
  MPI_Request req{MPI_REQUEST_NULL};
  MPI_Ibcast(&data, 1, get_type<T>(), root, comm.raw_handle(), &req);
  return irequest(req);
}

template<typename T,
         std::size_t N,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
irequest ibcast(nonstd::span<T, N> data, int root, const communicator& comm)
{
  MPI_Request req{MPI_REQUEST_NULL};
  auto        count = std::size(data);
  MPI_Ibcast(data.data(),
             to_int_size(count),
             get_type<T>(),
             root,
             comm.raw_handle(),
             &req);
  return irequest(req);
}

// === gather ===
// === root gets a single value from each rank and stores in contiguous memory
// --- blocking gather ---

// root or non-root
template<typename T, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void gather(const T& senddata, T* recvdata, int root, const communicator& comm)
{
  if (comm.rank() == root) {
    MPI_Gather(&senddata,
               1,
               get_type<T>(),
               recvdata,
               1,
               get_type<T>(),
               root,
               comm.raw_handle());
  } else {
    MPI_Gather(&senddata,
               1,
               get_type<T>(),
               nullptr,
               0,
               MPI_DATATYPE_NULL,
               root,
               comm.raw_handle());
  }
}

// non-root
template<typename T, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void gather(const T& senddata, int root, const communicator& comm)
{
#if !defined(NDEBUG)
  detail::check_non_root(root, comm);
#endif
  MPI_Gather(&senddata,
             1,
             get_type<T>(),
             nullptr,
             0,
             MPI_DATATYPE_NULL,
             root,
             comm.raw_handle());
}

// root or non-root
template<typename SendT,
         typename RecvT,
         std::size_t SendN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void gather(nonstd::span<SendT, SendN> senddata,
            RecvT*                     recvdata,
            int                        root,
            const communicator&        comm)
{
  auto count = to_int_size(senddata.size());
  if (comm.rank() == root) {
    MPI_Gather(senddata.data(),
               count,
               get_type<std::remove_const_t<SendT>>(),
               recvdata,
               count,
               get_type<RecvT>(),
               root,
               comm.raw_handle());
  } else {
    MPI_Gather(senddata.data(),
               count,
               get_type<std::remove_const_t<SendT>>(),
               nullptr,
               0,
               MPI_DATATYPE_NULL,
               root,
               comm.raw_handle());
  }
}

// non-root
template<typename SendT,
         std::size_t SendN,
         std::enable_if_t<is_mpi_type_v<std::remove_const_t<SendT>>, int> = 0>
void gather(nonstd::span<SendT, SendN> senddata,
            int                        root,
            const communicator&        comm)
{
#if !defined(NDEBUG)
  detail::check_non_root(root, comm);
#endif
  MPI_Gather(senddata.data(),
             to_int_size(senddata.size()),
             get_type<std::remove_const_t<SendT>>(),
             nullptr,
             0,
             MPI_DATATYPE_NULL,
             root,
             comm.raw_handle());
}

template<typename SendT,
         std::size_t SendN,
         typename RecvT,
         typename RecvCountT,
         std::size_t RecvCountN,
         typename RecvDisplT,
         std::size_t RecvDisplN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void gatherv(nonstd::span<SendT, SendN>           senddata,
             RecvT*                               recvdata,
             nonstd::span<RecvCountT, RecvCountN> recv_counts,
             nonstd::span<RecvDisplT, RecvDisplN> recv_displacements,
             int                                  root,
             const communicator&                  comm)
{
  static_assert(
    std::is_same_v<int, typename nonstd::span<RecvCountT>::value_type>,
    "gatherv: count types must be int");
  static_assert(
    std::is_same_v<int, typename nonstd::span<RecvDisplT>::value_type>,
    "gatherv: displacement types must be int");

  if (comm.rank() == root) {
    MPI_Gatherv(senddata.data(),
                to_int_size(senddata.size()),
                get_type<std::remove_const_t<SendT>>(),
                recvdata,
                recv_counts.data(),
                recv_displacements.data(),
                get_type<RecvT>(),
                root,
                comm.raw_handle());
  } else {
    MPI_Gatherv(senddata.data(),
                to_int_size(senddata.size()),
                get_type<std::remove_const_t<SendT>>(),
                nullptr,
                nullptr,
                nullptr,
                MPI_DATATYPE_NULL,
                root,
                comm.raw_handle());
  }
}

template<typename SendT,
         std::size_t SendN,
         typename RecvT,
         typename RecvCountT,
         std::size_t RecvCountN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void gatherv(nonstd::span<SendT, SendN>           senddata,
             RecvT*                               recvdata,
             nonstd::span<RecvCountT, RecvCountN> recv_counts,
             int                                  root,
             const communicator&                  comm)
{
  static_assert(
    std::is_same_v<int, typename nonstd::span<RecvCountT>::value_type>,
    "gatherv: count types must be int");

  if (comm.rank() == root) {
    auto recv_displacements = generate_displacements(recv_counts);
    gatherv(senddata,
            recvdata,
            recv_counts,
            nonstd::span(recv_displacements),
            root,
            comm);
  } else {
    gatherv(
      senddata, recvdata, nonstd::span<int>(), nonstd::span<int>(), root, comm);
  }
}

template<typename SendT,
         std::size_t SendN,
         typename RecvT,
         typename RecvCountT,
         std::size_t RecvCountN,
         typename RecvDisplT,
         std::size_t RecvDisplN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void allgatherv(nonstd::span<SendT, SendN>           senddata,
                RecvT*                               recvdata,
                nonstd::span<RecvCountT, RecvCountN> recv_counts,
                nonstd::span<RecvDisplT, RecvDisplN> recv_displacements,
                const communicator&                  comm)
{
  static_assert(
    std::is_same_v<int, typename nonstd::span<RecvCountT>::value_type>,
    "gatherv: count types must be int");
  static_assert(
    std::is_same_v<int, typename nonstd::span<RecvDisplT>::value_type>,
    "gatherv: displacement types must be int");

  MPI_Allgatherv(senddata.data(),
                 to_int_size(senddata.size()),
                 get_type<std::remove_const_t<SendT>>(),
                 recvdata,
                 recv_counts.data(),
                 recv_displacements.data(),
                 get_type<RecvT>(),
                 comm.raw_handle());
}

template<typename SendT,
         std::size_t SendN,
         typename RecvT,
         typename RecvCountT,
         std::size_t RecvCountN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void allgatherv(nonstd::span<SendT, SendN>           senddata,
                RecvT*                               recvdata,
                nonstd::span<RecvCountT, RecvCountN> recv_counts,
                const communicator&                  comm)
{
  static_assert(
    std::is_same_v<int, typename nonstd::span<RecvCountT>::value_type>,
    "gatherv: count types must be int");

  auto recv_displacements = generate_displacements(recv_counts);
  allgatherv(
    senddata, recvdata, recv_counts, nonstd::span(recv_displacements), comm);
}

// --- nonblocking gather ---
// -- not implemented --

// === allgather ===
// === get a single value from each rank and stores in contiguous memory
// --- blocking allgather ---
template<typename T, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void allgather(const T& senddata, T* recvdata, const communicator& comm)
{
  MPI_Allgather(
    &senddata, 1, get_type<T>(), recvdata, 1, get_type<T>(), comm.raw_handle());
}

template<typename SendT,
         typename RecvT,
         std::size_t SendN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void allgather(nonstd::span<SendT, SendN> senddata,
               RecvT*                     recvdata,
               const communicator&        comm)
{
  using val_t = std::remove_const_t<SendT>;
  auto count  = senddata.size();
  MPI_Allgather(senddata.data(),
                to_int_size(count),
                get_type<val_t>(),
                recvdata,
                to_int_size(count),
                get_type<val_t>(),
                comm.raw_handle());
}

// === all-to-all ===
// === each rank sends a single value to each rank
// --- blocking all-to-all ---
template<typename SendT,
         std::size_t SendN,
         typename CountT,
         typename RecvT,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int>                                            = 0,
         std::enable_if_t<std::is_integral_v<CountT>, int> = 0>
void alltoall(nonstd::span<SendT, SendN> senddata,
              CountT                     count,
              RecvT*                     recvdata,
              const communicator&        comm)
{
#if !defined(NDEBUG)
  detail::check_gather_size(count, std::size(senddata), comm);
#endif
  MPI_Alltoall(senddata.data(),
               to_int_size(count),
               get_type<std::remove_const_t<SendT>>(),
               recvdata,
               to_int_size(count),
               get_type<RecvT>(),
               comm.raw_handle());
}

// --- nonblocking all-to-all ---
template<typename SendT,
         std::size_t SendN,
         typename CountT,
         typename RecvT,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int>                                            = 0,
         std::enable_if_t<std::is_integral_v<CountT>, int> = 0>
irequest ialltoall(nonstd::span<SendT, SendN> senddata,
                   CountT                     count,
                   RecvT*                     recvdata,
                   const communicator&        comm)
{
#if !defined(NDEBUG)
  detail::check_gather_size(count, std::size(senddata), comm);
#endif

  MPI_Request req{MPI_REQUEST_NULL};
  MPI_Ialltoall(senddata.data(),
                to_int_size(count),
                get_type<std::remove_const_t<SendT>>(),
                recvdata,
                to_int_size(count),
                get_type<RecvT>(),
                comm.raw_handle(),
                &req);
  return irequest(req);
}

template<typename SendT,
         std::size_t SendN,
         typename SendCountT,
         std::size_t SendCountN,
         typename RecvT,
         typename RecvCountT,
         std::size_t RecvCountN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void alltoallv(nonstd::span<SendT, SendN>           senddata,
               nonstd::span<SendCountT, SendCountN> send_counts,
               RecvT*                               recvdata,
               nonstd::span<RecvCountT, RecvCountN> recv_counts,
               const communicator&                  comm)
{
  static_assert(
    std::is_same_v<int, typename nonstd::span<SendCountT>::value_type>
      && std::is_same_v<int, typename nonstd::span<RecvCountT>::value_type>,
    "alltoallv: count types must be int");

  auto send_displacements = generate_displacements(send_counts);
  auto recv_displacements = generate_displacements(recv_counts);
  MPI_Alltoallv(senddata.data(),
                send_counts.data(),
                send_displacements.data(),
                get_type<std::remove_const_t<SendT>>(),
                recvdata,
                recv_counts.data(),
                recv_displacements.data(),
                get_type<RecvT>(),
                comm.raw_handle());
}

template<typename SendT,
         std::size_t SendN,
         typename SendCountT,
         std::size_t SendCountN,
         typename SendDisplT,
         std::size_t SendDisplN,
         typename RecvT,
         typename RecvCountT,
         std::size_t RecvCountN,
         typename RecvDisplT,
         std::size_t RecvDisplN,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void alltoallv(nonstd::span<SendT, SendN>           senddata,
               nonstd::span<SendCountT, SendCountN> send_counts,
               nonstd::span<SendDisplT, SendDisplN> send_displacements,
               RecvT*                               recvdata,
               nonstd::span<RecvCountT, RecvCountN> recv_counts,
               nonstd::span<RecvDisplT, RecvDisplN> recv_displacements,
               const communicator&                  comm)
{
  static_assert(
    std::is_same_v<int, typename nonstd::span<SendCountT>::value_type>
      && std::is_same_v<int, typename nonstd::span<RecvCountT>::value_type>,
    "alltoallv: count types must be int");
  static_assert(
    std::is_same_v<int, typename nonstd::span<SendDisplT>::value_type>
      && std::is_same_v<int, typename nonstd::span<RecvDisplT>::value_type>,
    "alltoallv: displacement types must be int");

  MPI_Alltoallv(senddata.data(),
                send_counts.data(),
                send_displacements.data(),
                get_type<std::remove_const_t<SendT>>(),
                recvdata,
                recv_counts.data(),
                recv_displacements.data(),
                get_type<RecvT>(),
                comm.raw_handle());
}

// === reduce ===
// --- blocking reduce ---
template<typename T, typename Op, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void reduce(const T&            senddata,
            T&                  recvdata,
            const Op&           op,
            int                 root,
            const communicator& comm)
{
  static_assert(has_mpi_op_raw_handle_v<Op>, "An MPI operation is required");
  MPI_Reduce(&senddata,
             &recvdata,
             1,
             get_type<T>(),
             op.raw_handle(),
             root,
             comm.raw_handle());
}

template<typename SendT,
         std::size_t SendN,
         typename RecvT,
         std::size_t RecvN,
         typename Op,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void reduce(nonstd::span<SendT, SendN> senddata,
            nonstd::span<RecvT, RecvN> recvdata,
            const Op&                  op,
            int                        root,
            const communicator&        comm)
{
#if !defined(NDEBUG)
  if (senddata.size() != recvdata.size()) {
    throw std::runtime_error("reduce: send and recv data sizes do not match");
  }
#endif
  static_assert(has_mpi_op_raw_handle_v<Op>, "An MPI operation is required");
  MPI_Reduce(senddata.data(),
             recvdata.data(),
             to_int_size(senddata.size()),
             get_type<RecvT>(),
             op.raw_handle(),
             root,
             comm.raw_handle());
}

// === all-reduce ===
// --- blocking all-reduce ---
template<typename T, typename Op, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void allreduce(const T&            senddata,
               T&                  recvdata,
               const Op&           op,
               const communicator& comm)
{
  static_assert(has_mpi_op_raw_handle_v<Op>, "An MPI operation is required");
  MPI_Allreduce(
    &senddata, &recvdata, 1, get_type<T>(), op.raw_handle(), comm.raw_handle());
}

template<typename SendT,
         std::size_t SendN,
         typename RecvT,
         std::size_t RecvN,
         typename Op,
         std::enable_if_t<
           is_mpi_type_v<std::remove_const_t<SendT>> && !std::is_const_v<RecvT>
             && std::is_same_v<std::remove_const_t<SendT>, RecvT>,
           int> = 0>
void allreduce(nonstd::span<SendT, SendN> senddata,
               nonstd::span<RecvT, RecvN> recvdata,
               const Op&                  op,
               const communicator&        comm)
{
#if !defined(NDEBUG)
  if (senddata.size() != recvdata.size()) {
    throw std::runtime_error(
      "allreduce: send and recv data sizes do not match");
  }
#endif
  static_assert(has_mpi_op_raw_handle_v<Op>, "An MPI operation is required");
  MPI_Allreduce(senddata.data(),
                recvdata.data(),
                to_int_size(senddata.size()),
                get_type<RecvT>(),
                op.raw_handle(),
                comm.raw_handle());
}

} // namespace mpipp

#pragma once

#include <mpipp/comm.h>
#include <mpipp/config.h>
#include <mpipp/datatype.h>
#include <mpipp/request.h>
#include <mpipp/status.h>
#include <mpipp/utility.h>

#include <nonstd/span.hpp>

#include <type_traits>
#include <utility>

namespace mpipp {

// --- blocking send ---
template<typename T, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void send(const T& data, int dest, int tag, const communicator& comm)
{
  MPI_Send(&data, 1, get_type<T>(), dest, tag, comm.raw_handle());
}

template<typename T,
         std::size_t N,
         std::enable_if_t<is_mpi_type_v<std::remove_const_t<T>>, int> = 0>
void send(nonstd::span<T, N> data, int dest, int tag, const communicator& comm)
{
  using val_t = std::remove_const_t<T>;
  auto count  = std::size(data);
  MPI_Send(count > 0 ? data.data() : nullptr,
           to_int_size(count),
           get_type<val_t>(),
           dest,
           tag,
           comm.raw_handle());
}

// --- nonblocking send ---
template<typename T, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
[[nodiscard]] irequest
isend(const T& data, int dest, int tag, const communicator& comm)
{
  MPI_Request req{MPI_REQUEST_NULL};
  MPI_Isend(&data, 1, get_type<T>(), dest, tag, comm.raw_handle(), &req);
  return irequest(req);
}

template<typename T,
         std::size_t N,
         std::enable_if_t<is_mpi_type_v<std::remove_const_t<T>>, int> = 0>
[[nodiscard]] irequest
isend(nonstd::span<T, N> data, int dest, int tag, const communicator& comm)
{
  using val_t = std::remove_const_t<T>;
  MPI_Request req{MPI_REQUEST_NULL};
  const auto  count{std::size(data)};
  MPI_Isend(count > 0 ? data.data() : nullptr,
            to_int_size(count),
            get_type<val_t>(),
            dest,
            tag,
            comm.raw_handle(),
            &req);
  return irequest(req);
}

// === receive ===
// --- blocking receive ---
template<typename T,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
void recv(T& data, int source, int tag, const communicator& comm, status& s)
{
  MPI_Recv(&data, 1, get_type<T>(), source, tag, comm.raw_handle(), &s);
}

template<typename T,
         std::size_t N,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
void recv(nonstd::span<T, N>  data,
          int                 source,
          int                 tag,
          const communicator& comm,
          status&             s)
{
  auto count = std::size(data);
  MPI_Recv(count > 0 ? data.data() : nullptr,
           to_int_size(count),
           source,
           tag,
           comm.raw_handle(),
           &s);
}

template<typename T,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
void recv(T&                  data,
          int                 source,
          int                 tag,
          const communicator& comm,
          status_ignore_t /*s*/)
{
  MPI_Recv(
    &data, 1, get_type<T>(), source, tag, comm.raw_handle(), MPI_STATUS_IGNORE);
}

template<typename T,
         std::size_t N,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
void recv(nonstd::span<T, N>  data,
          int                 source,
          int                 tag,
          const communicator& comm,
          status_ignore_t /*s*/)
{
  auto count = std::size(data);
  MPI_Recv(count > 0 ? data.data() : nullptr,
           to_int_size(count),
           get_type<T>(),
           source,
           tag,
           comm.raw_handle(),
           MPI_STATUS_IGNORE);
}

// --- nonblocking receive ---
template<typename T,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
irequest irecv(T& data, int source, int tag, const communicator& comm)
{
  MPI_Request req{MPI_REQUEST_NULL};
  MPI_Irecv(&data, 1, get_type<T>(), source, tag, comm.raw_handle(), &req);
  return irequest(req);
}

template<typename T,
         std::size_t N,
         std::enable_if_t<!std::is_const_v<T> && is_mpi_type_v<T>, int> = 0>
irequest
irecv(nonstd::span<T, N> data, int source, int tag, const communicator& comm)
{
  MPI_Request req{MPI_REQUEST_NULL};
  auto        count = std::size(data);
  MPI_Irecv(count > 0 ? data.data() : nullptr,
            to_int_size(count),
            get_type<T>(),
            source,
            tag,
            comm.raw_handle(),
            &req);
  return irequest(req);
}

// === probe ===
// --- blocking probe ---
inline status probe(int source, int tag, const communicator& comm)
{
  status s;
  MPI_Probe(source, tag, comm.raw_handle(), &s);
  return s;
}

// --- nonblocking probe ---
inline std::pair<bool, status>
iprobe(int source, int tag, const communicator& comm)
{
  int    result;
  status s;
  MPI_Iprobe(source, tag, comm.raw_handle(), &result, &s);
  return std::make_pair(static_cast<bool>(result), s);
}

// === send and receive ===
// --- send and receive ---
template<typename T, std::enable_if_t<is_mpi_type_v<T>, int> = 0>
void sendrecv(const T&            senddata,
              int                 dest,
              int                 sendtag,
              T&                  recvdata,
              int                 source,
              int                 recvtag,
              const communicator& comm)
{
  MPI_Sendrecv(&senddata,
               1,
               get_type<T>(),
               dest,
               sendtag,
               &recvdata,
               1,
               get_type<T>(),
               source,
               recvtag,
               comm.raw_handle(),
               MPI_STATUS_IGNORE);
}

template<typename SendT,
         typename RecvT,
         std::size_t SendN,
         std::size_t RecvN,
         std::enable_if_t<is_mpi_type_v<std::remove_const_t<SendT>>
                            && !std::is_const_v<RecvT> && is_mpi_type_v<RecvT>,
                          int> = 0>
void sendrecv(nonstd::span<SendT, SendN> senddata,
              int                        dest,
              int                        sendtag,
              nonstd::span<RecvT, RecvN> recvdata,
              int                        source,
              int                        recvtag,
              const communicator&        comm)
{
  auto send_n = std::size(senddata);
  auto recv_n = std::size(recvdata);
  MPI_Sendrecv(send_n > 0 ? senddata.data() : nullptr,
               to_int_size(send_n),
               get_type<std::remove_const_t<SendT>>(),
               dest,
               sendtag,
               recv_n > 0 ? recvdata.data() : nullptr,
               to_int_size(recv_n),
               get_type<RecvT>(),
               source,
               recvtag,
               comm.raw_handle(),
               MPI_STATUS_IGNORE);
}

} // namespace mpipp

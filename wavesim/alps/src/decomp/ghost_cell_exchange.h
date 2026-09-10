#pragma once

#include "block_partition.h"
#include "mdcomm.h"
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/manager.h>

#include <mpipp/point2point.h>

#include <cstdint>
#include <utility>

namespace alps {

/**  Tag prefixes used for halo-cell exchange MPI messages.
 *
 * The full MPI tag for a message is intended to be constructed by combining a
 * one-byte prefix with an 8-bit message-local tag.
 */
struct HaloMsgTagPrefix
{
  static constexpr uint8_t Z{128};
};

/// Perform the operation of MPI_Sendrecv
/**
 * If MPI is not GPU aware, use host buffers to transfer data
 * @param sendData data to be sent
 * @param dstRank destination rank to send data to
 * @param recvData data to be received
 * @param srcRank source rank to receive data from
 * @param comm MPI communicator
 * @param tag message tag
 */
template<class... DT1, class... DT2>
void exchange(const Kokkos::View<DT1...>& sendData,
              const int                   dstRank,
              const Kokkos::View<DT2...>& recvData,
              const int                   srcRank,
              const mpipp::communicator&  comm,
              int                         tag)
{
  if (!sendData.span_is_contiguous()) {
    throw std::runtime_error("exchange: sendData is not contiguous");
  }
  if (!recvData.span_is_contiguous()) {
    throw std::runtime_error("exchange: recvData is not contiguous");
  }

  if constexpr (alps::mpi_can_access_v<
                  typename Kokkos::View<DT1...>::memory_space>
                && alps::mpi_can_access_v<
                  typename Kokkos::View<DT2...>::memory_space>) {
    mpipp::sendrecv(nonstd::span(sendData.data(), sendData.span()),
                    dstRank,
                    tag,
                    nonstd::span(recvData.data(), recvData.span()),
                    srcRank,
                    tag,
                    comm);
  } else {
    auto send_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                 PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                                 sendData);
    if (dstRank != MPI_PROC_NULL) {
      Kokkos::deep_copy(send_host, sendData);
    }

    auto recv_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                 PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                                 recvData);

    mpipp::sendrecv(nonstd::span(send_host.data(), send_host.span()),
                    dstRank,
                    tag,
                    nonstd::span(recv_host.data(), recv_host.span()),
                    srcRank,
                    tag,
                    comm);

    if (srcRank != MPI_PROC_NULL) {
      Kokkos::deep_copy(recvData, recv_host);
    }
  }
}

/// Perform the operation of MPI_Sendrecv using MPI_Irecv and MPI_Isend
/**
 * If MPI is not GPU aware, fall back to blocking exchange
 * @param requests MPI request pool to append the requests to
 * @param sendData data to be sent
 * @param dstRank destination rank to send data to
 * @param recvData data to be received
 * @param srcRank source rank to receive data from
 * @param comm MPI communicator
 * @param tag message tag
 */
template<class... DT1, class... DT2>
[[nodiscard]] mpipp::irequest_pool
async_exchange(const Kokkos::View<DT1...>& sendData,
               const int                   dstRank,
               const Kokkos::View<DT2...>& recvData,
               const int                   srcRank,
               const mpipp::communicator&  comm,
               int                         tag)
{
  if (!sendData.span_is_contiguous()) {
    throw std::runtime_error("exchange: sendData is not contiguous");
  }
  if (!recvData.span_is_contiguous()) {
    throw std::runtime_error("exchange: recvData is not contiguous");
  }

  mpipp::irequest_pool requests;
  if constexpr (alps::mpi_can_access_v<
                  typename Kokkos::View<DT1...>::memory_space>
                && alps::mpi_can_access_v<
                  typename Kokkos::View<DT2...>::memory_space>) {
    requests.push(mpipp::irecv(
      nonstd::span(recvData.data(), recvData.span()), srcRank, tag, comm));
    requests.push(mpipp::isend(
      nonstd::span(sendData.data(), sendData.span()), dstRank, tag, comm));
  } else {
    auto recv_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                 PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                                 recvData);
    requests.push(
      mpipp::irecv(
        nonstd::span(recv_host.data(), recv_host.span()), srcRank, tag, comm),
      [recv_host,
       recvData,
       srcRank,
       stream = RuntimeManager::instance().host_to_device_stream()] {
        if (srcRank != MPI_PROC_NULL) {
          Kokkos::deep_copy(stream, recvData, recv_host);
          stream.fence();
        }
      });

    auto send_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                 PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                                 sendData);
    if (dstRank != MPI_PROC_NULL) {
      auto stream = RuntimeManager::instance().device_to_host_stream();
      Kokkos::deep_copy(stream, send_host, sendData);
      stream.fence();
    }

    requests.push(
      mpipp::isend(
        nonstd::span(send_host.data(), send_host.span()), dstRank, tag, comm),
      [send_host = std::move(send_host)] { /* extend send_host lifetime */ });
  }

  return requests;
}

//-----------------------------------------------------------------------------

/** Update the upper halo cells in the z-direction using non-blocking operations
 *
 * @note The @p msg_tag is used to construct the MPI message tag for
 * differentiating messages. The full tag is a 16-bit integer with the format
 * [prefix (8 bits)][msg_tag (8 bits)].
 */
template<class D, class... RP>
[[nodiscard]] mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&     plan,
                          const HaloView<D, RP...>& data,
                          uint8_t                   msg_tag,
                          uint8_t prefix = HaloMsgTagPrefix::Z)
{
  if (plan.comm.dims[2] == 1) return {};

  auto upper_ghost =
    trailing_subview(data, index_range(local_end(data, 2), end(data, 2)))
      .view();

  auto lower_data =
    trailing_subview(data, index_range(0, -begin(data, 2))).view();
  auto upper_id = plan.comm.next_proc_on_axis[2];
  auto lower_id = plan.comm.prev_proc_on_axis[2];

  return async_exchange(lower_data,
                        lower_id,
                        upper_ghost,
                        upper_id,
                        plan.comm,
                        ((uint16_t)prefix << 8) | msg_tag);
}

/** Update the lower halo cells in the z-direction using non-blocking operations
 *
 * @note The @p msg_tag is used to construct the MPI message tag for
 * differentiating messages. The full tag is a 16-bit integer with the format
 * [prefix (8 bits)][msg_tag (8 bits)].
 */
template<class D, class... RP>
[[nodiscard]] mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&     plan,
                          const HaloView<D, RP...>& data,
                          uint8_t                   msg_tag,
                          uint8_t prefix = HaloMsgTagPrefix::Z)
{
  if (plan.comm.dims[2] == 1) return {};

  auto lower_ghost =
    trailing_subview(data, index_range(begin(data, 2), 0)).view();

  auto upper_data =
    trailing_subview(
      data,
      index_range(local_end(data, 2) + begin(data, 2), local_end(data, 2)))
      .view();

  auto upper_id = plan.comm.next_proc_on_axis[2];
  auto lower_id = plan.comm.prev_proc_on_axis[2];

  return async_exchange(upper_data,
                        upper_id,
                        lower_ghost,
                        lower_id,
                        plan.comm,
                        ((uint16_t)prefix << 8) | msg_tag);
}

/** Update the upper and lower halo cells in the z-direction using non-blocking
 * operations
 *
 * @note The @p msg_tag is used to construct the MPI message tag for
 * differentiating messages. The full tag is a 16-bit integer with the format
 * [prefix (8 bits)][msg_tag (8 bits)].
 */
template<class D, class... RP>
[[nodiscard]] mpipp::irequest_pool
async_update_halo_z(const BlockPartition&     plan,
                    const HaloView<D, RP...>& data,
                    uint8_t                   msg_tag,
                    uint8_t                   prefix = HaloMsgTagPrefix::Z)
{
  auto requests = async_update_halo_upper_z(plan, data, msg_tag, prefix);
  requests.push(async_update_halo_lower_z(plan, data, msg_tag, prefix));
  return requests;
}

//-----------------------------------------------------------------------------
// External template instantiations
extern template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&      plan,
                          const HaloView<double***>& data,
                          uint8_t                    msg_tag,
                          uint8_t                    prefix);

extern template mpipp::irequest_pool async_update_halo_upper_z(
  const BlockPartition&                                                   plan,
  const HaloView<double***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t msg_tag,
  uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&                           plan,
                          const HaloView<double***, default_memory_pool>& data,
                          uint8_t msg_tag,
                          uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&     plan,
                          const HaloView<float***>& data,
                          uint8_t                   msg_tag,
                          uint8_t                   prefix);

extern template mpipp::irequest_pool async_update_halo_upper_z(
  const BlockPartition&                                                  plan,
  const HaloView<float***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t msg_tag,
  uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&                          plan,
                          const HaloView<float***, default_memory_pool>& data,
                          uint8_t msg_tag,
                          uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&      plan,
                          const HaloView<double***>& data,
                          uint8_t                    msg_tag,
                          uint8_t                    prefix);

extern template mpipp::irequest_pool async_update_halo_lower_z(
  const BlockPartition&                                                   plan,
  const HaloView<double***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t msg_tag,
  uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&                           plan,
                          const HaloView<double***, default_memory_pool>& data,
                          uint8_t msg_tag,
                          uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&     plan,
                          const HaloView<float***>& data,
                          uint8_t                   msg_tag,
                          uint8_t                   prefix);

extern template mpipp::irequest_pool async_update_halo_lower_z(
  const BlockPartition&                                                  plan,
  const HaloView<float***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t msg_tag,
  uint8_t prefix);

extern template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&                          plan,
                          const HaloView<float***, default_memory_pool>& data,
                          uint8_t msg_tag,
                          uint8_t prefix);

//-----------------------------------------------------------------------------
/// Update the upper halo cells in the z-direction
template<class D, class... RP>
void update_halo_upper_z(const BlockPartition&     plan,
                         const HaloView<D, RP...>& data,
                         uint8_t                   msg_tag,
                         uint8_t                   prefix = HaloMsgTagPrefix::Z)
{
  async_update_halo_upper_z(plan, data, msg_tag, prefix).waitall();
}

/// Update the lower halo cells in the z-direction
template<class D, class... RP>
void update_halo_lower_z(const BlockPartition&     plan,
                         const HaloView<D, RP...>& data,
                         uint8_t                   msg_tag,
                         uint8_t                   prefix = HaloMsgTagPrefix::Z)
{
  async_update_halo_lower_z(plan, data, msg_tag, prefix).waitall();
}

/// Update the upper and lower halo cells in the z-direction
template<class D, class... RP>
void update_halo_z(const BlockPartition&     plan,
                   const HaloView<D, RP...>& data,
                   uint8_t                   msg_tag,
                   uint8_t                   prefix = HaloMsgTagPrefix::Z)
{
  update_halo_upper_z(plan, data, msg_tag, prefix);
  update_halo_lower_z(plan, data, msg_tag, prefix);
}

} // namespace alps

//
// Created by xuanx004 on 7/12/24.
//

#include "ghost_cell_exchange.h"

namespace alps {

// Explicit instantiation of the template functions
template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&      plan,
                          const HaloView<double***>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool async_update_halo_upper_z(
  const BlockPartition&                                                   plan,
  const HaloView<double***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t,
  uint8_t);

template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&                           plan,
                          const HaloView<double***, default_memory_pool>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&     plan,
                          const HaloView<float***>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool async_update_halo_upper_z(
  const BlockPartition&                                                  plan,
  const HaloView<float***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t,
  uint8_t);

template mpipp::irequest_pool
async_update_halo_upper_z(const BlockPartition&                          plan,
                          const HaloView<float***, default_memory_pool>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&      plan,
                          const HaloView<double***>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool async_update_halo_lower_z(
  const BlockPartition&                                                   plan,
  const HaloView<double***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t,
  uint8_t);

template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&                           plan,
                          const HaloView<double***, default_memory_pool>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&     plan,
                          const HaloView<float***>& data,
                          uint8_t,
                          uint8_t);

template mpipp::irequest_pool async_update_halo_lower_z(
  const BlockPartition&                                                  plan,
  const HaloView<float***, Kokkos::DefaultExecutionSpace::memory_space>& data,
  uint8_t,
  uint8_t);

template mpipp::irequest_pool
async_update_halo_lower_z(const BlockPartition&                          plan,
                          const HaloView<float***, default_memory_pool>& data,
                          uint8_t,
                          uint8_t);
} // namespace alps

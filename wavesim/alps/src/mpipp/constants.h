#pragma once

#include <mpipp/config.h>

namespace mpipp {

enum class equality_type
{
  ident     = MPI_IDENT,
  congruent = MPI_CONGRUENT,
  similar   = MPI_SIMILAR,
  unequal   = MPI_UNEQUAL
};
} // namespace mpipp

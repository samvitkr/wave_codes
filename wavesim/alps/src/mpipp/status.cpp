//
// Created by xuanx004 on 3/14/24.
//

#include "status.h"

namespace mpipp {

static_assert(sizeof(status) == sizeof(MPI_Status), "status size mismatch");
static_assert(alignof(status) == alignof(MPI_Status),
              "status alignment mismatch");

} // namespace mpipp

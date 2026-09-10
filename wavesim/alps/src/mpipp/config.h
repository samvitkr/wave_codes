#pragma once

#if !defined(c_plusplus) && !defined(__cplusplus)
/* Force MPICH not to define SEEK_SET, SEEK_CUR, and SEEK_END, which
   conflict with the versions in <stdio.h> and <cstdio>. */
#define MPICH_IGNORE_CXX_SEEK 1
#endif

/* We do not want to link in the OpenMPI CXX stuff */
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#define __MPIPP_OMPI_SKIP_MPICXX_DEFINED__
#endif

#include <mpi.h>

// undefine OMPI_SKIP_MPICXX to avoid conflicts with other libraries defining
// OMPI_SKIP_MPICXX
#if defined(__MPIPP_OMPI_SKIP_MPICXX_DEFINED__)
#undef OMPI_SKIP_MPICXX
#undef __MPIPP_OMPI_SKIP_MPICXX_DEFINED__
#endif

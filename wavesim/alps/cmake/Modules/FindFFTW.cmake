# - Find the FFTW library
#
# Original version of this file:
#   Copyright (c) 2015, Wenzel Jakob
#   https://github.com/wjakob/layerlab/blob/master/cmake/FindFFTW.cmake, commit 4d58bfdc28891b4f9373dfe46239dda5a0b561c6
# Modifications:
#   Copyright (c) 2017, Patrick Bos
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] [COMPONENTS component1 ... componentX] )
#
# It sets the following variables:
#   FFTW_FOUND                  ... true if fftw is found on the system
#   FFTW_<COMPONENT>_FOUND      ... true if the component is found on the system (<COMPONENT> name is upper-case)
#   FFTW_LIBRARIES              ... full paths to all found fftw libraries
#   FFTW_<COMPONENT>_LIBRARY    ... full path to one of the components (<COMPONENT> name is upper-case)
#   FFTW_INCLUDE_DIRS           ... fftw include directory paths
#
# The following variables will be checked by the function
#   FFTW_USE_STATIC_LIBS        ... if true, only static libraries are found, otherwise both static and shared.
#   FFTW_ROOT                   ... if set, the libraries are exclusively searched
#                                   under this path
#   FFTW_USE_MKL:               ... if set, MKL implementation is used for FFT
#
# This package supports the following components:
#   float
#   double
#   float_threads_lib
#   double_threads_lib
#   float_openmp_lib
#   double_openmp_lib
#

set(FFTW_VALID_LIBRARIES_TYPES double float
                               double_openmp float_openmp
                               double_threads float_threads)
set(FFTW_LIBRARIES_DOUBLE_NAMES "fftw3")
set(FFTW_LIBRARIES_DOUBLE_OPENMP_NAMES "fftw3" "fftw3_omp")
set(FFTW_LIBRARIES_DOUBLE_THREADS_NAMES "fftw3" "fftw3_threads")
set(FFTW_LIBRARIES_FLOAT_NAMES "fftw3f")
set(FFTW_LIBRARIES_FLOAT_OPENMP_NAMES "fftw3f" "fftw3f_omp")
set(FFTW_LIBRARIES_FLOAT_THREADS_NAMES "fftw3f" "fftw3f_threads")

macro(_fftw_set_cache_if_unset name value)
  if(NOT ${name})
    set(${name} "${value}" CACHE STRING "" FORCE)
  endif()
endmacro()

if(NOT FFTW_FIND_COMPONENTS)
  set(FFTW_SEARCH_TYPES DOUBLE)
  set(FFTW_REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDE_DIRS)
else()
  foreach(component ${FFTW_FIND_COMPONENTS})
    list(FIND FFTW_VALID_LIBRARIES_TYPES ${component} component_location)
    if(NOT component_location EQUAL -1)
      string(TOUPPER ${component} UPPERCOMPONENT)
      list(APPEND FFTW_SEARCH_TYPES ${UPPERCOMPONENT})
      list(APPEND FFTW_REQUIRED_VARS FFTW_${UPPERCOMPONENT}_LIBRARY)
      set(FFTW_${UPPERCOMPONENT}_LIBRARY)
    else()
      message(WARNING "Skip unknown FFTW component ${component}.")
    endif()
  endforeach()
  list(APPEND FFTW_REQUIRED_VARS FFTW_INCLUDE_DIRS)
endif()
set(FFTW_FOUND False)
set(FFTW_LIBRARIES)
set(FFTW_INCLUDE_DIRS)

if( FFTW_USE_MKL )
# ------------------------------------------------------------------------
#  Use FFTW interface from MKL
# ------------------------------------------------------------------------
  if ( NOT MKL_FOUND )
    find_package( MKL )
  endif()

  if ( MKL_FOUND )
    find_path(FFTW_INCLUDE_DIRS
      NAMES "fftw3.h"
      PATHS ${MKL_INCLUDE_DIRS}/fftw
    )
    foreach( _comp ${FFTW_SEARCH_TYPES} )
      set( FFTW_${_comp}_LIBRARY "${MKL_LIBRARIES}" )
    endforeach()
  else()
    message(SEND_ERROR "FindFFTW: FFTW_USE_MKL enabled, but MKL was not found.")
  endif()

else()

# ------------------------------------------------------------------------
#  Search regular FFTW library
# ------------------------------------------------------------------------

  if( NOT FFTW_ROOT )
    if ( DEFINED FFTW_DIR )
      set( FFTW_ROOT ${FFTW_DIR} )
    elseif ( DEFINED ENV{FFTW_DIR} )
      set( FFTW_ROOT $ENV{FFTW_DIR} )
    elseif ( DEFINED ENV{FFTW_ROOT} )
      set( FFTW_ROOT $ENV{FFTW_ROOT} )
    endif()
  endif()

  if( FFTW_ROOT ) # On cc[a|b|t] FFTW_DIR is set to the lib directory :(
      get_filename_component(_dirname ${FFTW_ROOT} NAME)
      if( _dirname MATCHES "lib" )
          set( FFTW_ROOT "${FFTW_ROOT}/.." )
      endif()
  endif()


  if( NOT FFTW_ROOT )
      # Check if we can use PkgConfig
      find_package(PkgConfig)

      #Determine from PKG
      if( PKG_CONFIG_FOUND AND NOT FFTW_ROOT )
          pkg_check_modules( PKG_FFTW QUIET "fftw3" )
      endif()
  endif()

# Check whether to search static or dynamic libs
  if( FFTW_USE_STATIC_LIBS )
    set( _FFTW_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
  endif()

  if( FFTW_ROOT )
      set( _default_paths NO_DEFAULT_PATH )
      set( _lib_paths ${FFTW_ROOT} )
      set( _include_paths ${FFTW_ROOT} )
  else()
      set( _lib_paths ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR} )
      set( _include_paths ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR} )
  endif()

# find includes
  if( NOT FFTW_INCLUDE_DIRS )
    find_path(FFTW_INCLUDE_DIRS
      NAMES "fftw3.h"
      PATHS ${_include_paths} ENV FFTW_INC
      PATH_SUFFIXES "include"
      ${_default_paths}
    )
  endif()

# find libs
  foreach(_comp ${FFTW_SEARCH_TYPES})
    if( NOT FFTW_${_comp}_LIBRARY )
      foreach(LIB ${FFTW_LIBRARIES_${_comp}_NAMES})
        find_library(FFTW_${LIB}_LIB
          NAMES ${LIB}
          PATHS ${_lib_paths} ENV LD_LIBRARY_PATH
          PATH_SUFFIXES "lib" "lib64"
          ${_default_paths}
        )
        if(FFTW_${LIB}_LIB)
          list(APPEND FFTW_${_comp}_LIBRARY ${FFTW_${LIB}_LIB})
        endif()
      endforeach()
      if(FFTW_${_comp}_LIBRARY)
        list( REMOVE_DUPLICATES FFTW_${_comp}_LIBRARY )
        _fftw_set_cache_if_unset(FFTW_${_comp}_LIBRARY ${FFTW_${_comp}_LIBRARY})
      endif()
    endif()
  endforeach()

  if( FFTW_USE_STATIC_LIBS )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${_FFTW_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES} )
  endif()

endif()

# ------------------------------------------------------------------------
#  End finding FFTW libraries
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
#  Call FPHSA helper, see https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
# ------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)

foreach( _comp ${FFTW_SEARCH_TYPES} )
  string( TOLOWER ${_comp} _lowercomp )
  if (FFTW_${_comp}_LIBRARY AND FFTW_INCLUDE_DIRS)
    # Define FFTW_<component>_FOUND as needed by the component handler
    set(FFTW_${_comp}_FOUND TRUE)
    set(FFTW_${_lowercomp}_FOUND TRUE)

    list(APPEND FFTW_LIBRARIES ${FFTW_${l}_LIBRARIES})
  else()
    set(FFTW_${_lowercomp}_FOUND FALSE)
  endif()
endforeach()
list( REMOVE_DUPLICATES FFTW_LIBRARIES )

find_package_handle_standard_args(FFTW
  REQUIRED_VARS ${FFTW_REQUIRED_VARS}
  HANDLE_COMPONENTS
)

# ------------------------------------------------------------------------
#  Add imported targets
# ------------------------------------------------------------------------
if(FFTW_FOUND)
  add_library(FFTW::fftw INTERFACE IMPORTED)
  foreach( _comp ${FFTW_SEARCH_TYPES} )
    string( TOLOWER ${_comp} _lowercomp )
    if(FFTW_${_lowercomp}_FOUND)
      set( _target FFTW::${_lowercomp} )
      add_library( ${_target} INTERFACE IMPORTED )
      target_link_libraries( ${_target} INTERFACE ${FFTW_${_comp}_LIBRARY} )
      target_include_directories( ${_target} INTERFACE ${FFTW_INCLUDE_DIRS} )

      target_link_libraries(FFTW::fftw INTERFACE ${_target})
    endif()
  endforeach()
endif()

mark_as_advanced(
  FFTW_INCLUDE_DIRS
  FFTW_LIBRARIES
)

if (FFTW_FIND_DEBUG)
  message(STATUS "FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")
  message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
endif()


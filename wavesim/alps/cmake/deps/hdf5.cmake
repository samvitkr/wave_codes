if(NOT ${PROJECT_NAME_UPPERCASE}_BUILD_HDF5)
  find_package(HDF5 COMPONENTS C)
else()
  CPMAddPackage(
    NAME hdf5
    VERSION 1.14.6
    URL https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
    URL_HASH SHA256=e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b
    OPTIONS "HDF5_BUILD_EXAMPLES OFF"
            "BUILD_TESTING OFF"
            "HDF5_ENABLE_SZIP_SUPPORT OFF"
            "HDF5_ENABLE_SZIP_ENCODING OFF"
            "HDF5_ENABLE_Z_LIB_SUPPORT OFF"
            "BUILD_SHARED_LIBS OFF"
            "HDF5_ENABLE_PARALLEL ON"
            "HDF5_ENABLE_THREADSAFE OFF"
            "HDF5_BUILD_HL_LIB OFF"
            "HDF5_BUILD_CPP_LIB OFF"
            "HDF5_BUILD_FORTRAN OFF"
            "HDF5_BUILD_TOOLS OFF"
            "HDF5_EXTERNALLY_CONFIGURED ON"
  )
  if(TARGET hdf5-static)
    add_library(HDF5::HDF5 ALIAS hdf5-static)
  else()
    message(FATAL_ERROR "HDF5 build is requested but the target hdf5-static is not found.")
  endif()
endif()

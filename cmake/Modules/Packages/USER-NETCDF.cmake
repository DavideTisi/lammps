if(PKG_USER-NETCDF)
  find_package(NetCDF REQUIRED)
  include_directories(${NETCDF_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${NETCDF_LIBRARIES})
  add_definitions(-DLMP_HAS_NETCDF -DNC_64BIT_DATA=0x0020)
endif()
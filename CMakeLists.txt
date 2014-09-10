generate_sourcelists("./")
set(RELATIVE_PATH "")
include(auto.cmake)

if(WITH_MPI)
  add_executable(gameoflife_adcirc ${SOURCES})
  target_link_libraries(gameoflife_adcirc ${LOCAL_LIBGEODECOMP_LINK_LIB})
endif(WITH_MPI)

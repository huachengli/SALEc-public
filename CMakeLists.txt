cmake_minimum_required(VERSION 3.13)
project(SALEc C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

LINK_LIBRARIES(m z)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_GZDIR")

#control whether to compile and link tracer module.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSALEC_USE_TRACER")
#control whether to compile and link checkpoint
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSALEC_USE_CHECKPOINT")



FIND_PACKAGE(MPI REQUIRED)
if(MPI_FOUND)
    set(CMAKE_C_COMPILER mpicc)
elseif(NOT MPI_FOUND)
    message("MPI UNFOUND")
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)
set(SRCS ${PROJECT_SOURCE_DIR}/src/Variable.c
         ${PROJECT_SOURCE_DIR}/src/Setup.c
         ${PROJECT_SOURCE_DIR}/src/Variable.c
         ${PROJECT_SOURCE_DIR}/src/Iteration.c
         ${PROJECT_SOURCE_DIR}/src/State.c
         ${PROJECT_SOURCE_DIR}/src/Output.c
         ${PROJECT_SOURCE_DIR}/src/Ghost.c
         ${PROJECT_SOURCE_DIR}/src/Error.c
         ${PROJECT_SOURCE_DIR}/src/Parallelism.c
         ${PROJECT_SOURCE_DIR}/src/Reconstruction.c
         ${PROJECT_SOURCE_DIR}/src/InputParser.c
         ${PROJECT_SOURCE_DIR}/src/Tracer.c)

add_executable(SALEc ${SRCS} ${PROJECT_SOURCE_DIR}/src/SALEc.c)
add_executable(WriteVtm ${SRCS} ${PROJECT_SOURCE_DIR}/src/WriteVtm.c)

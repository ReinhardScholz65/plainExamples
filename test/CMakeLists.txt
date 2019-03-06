cmake_minimum_required (VERSION 2.8)

if (NOT DEFINED LEMONADE_DIR)
message("LEMONADE_DIR is not provided. If build fails, use -DLEMONADE_DIR=/path/to/LeMonADE/")
endif()

SET (LEMONADE_INCLUDE_DIR ${LEMONADE_DIR}include/)
SET (LEMONADE_LIBRARY_DIR ${LEMONADE_DIR}lib/)

message( ${LEMONADE_DIR} )

include_directories("../")

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -O2 -std=c++11")
SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG} -O2 ")

## ###############  Tests linked to main  ############# ##
add_executable(testsLeMonADE test_main.cpp test_updaterExample.cpp)
target_link_libraries(testsLeMonADE LeMonADE )

cmake_minimum_required (VERSION 2.8)

# hard code include directories, better use 
#   cmake -DLEMONADE_DIR=/yourpath/to/Lemonade
# this is Martins lemonde library location, uncomment this line
#set (LEMONADE_DIR "/scratch/localuser/lemonade/lemonadeInstall/")
# here your lemonde library location should be inserted:

# set (LEMONADE_DIR "<your directory>")

if (NOT DEFINED LEMONADE_DIR)
message("LEMONADE_DIR is not provided. If build fails, use -DLEMONADE_DIR=/path/to/LeMonADE/")
endif()

SET (LEMONADE_INCLUDE_DIR ${LEMONADE_DIR}include/)
SET (LEMONADE_LIBRARY_DIR ${LEMONADE_DIR}lib/)

message( ${LEMONADE_INCLUDE_DIR} )

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY "./bin/")

include_directories("../../")

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

#define flags
SET (CMAKE_CXX_FLAGS " -O3 -msse2 -mssse3 -fexpensive-optimizations -std=c++11")
SET (CMAKE_C_FLAGS "-O3 -msse2 -mssse3 -fexpensive-optimizations -std=c++11")

## ###############  add your executables here ############# ##
add_executable(linearChainIdeal lemonade_example_main.cpp)
target_link_libraries(linearChainIdeal LeMonADE)

add_executable(linearChainReal lemonade_example_main_exclvol.cpp)
target_link_libraries(linearChainReal LeMonADE)

add_executable(linearChainSolvent lemonade_example_main_nninter.cpp)
target_link_libraries(linearChainSolvent LeMonADE)

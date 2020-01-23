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

# add_executable(starReal_17_64stars lemonade_star_main_exclvol.cpp)
# target_link_libraries(starReal_17_64stars LeMonADE)

# add_executable(starReal_65_16stars lemonade_star_main_exclvol.cpp)
# target_link_libraries(starReal_65_16stars LeMonADE)

# add_executable(starReal_65_128stars lemonade_star_main_exclvol.cpp)
# target_link_libraries(starReal_65_128stars LeMonADE)

# add_executable(starReal_65_128stars_AB lemonade_star_main_exclvol_AB.cpp)
# target_link_libraries(starReal_65_128stars_AB LeMonADE)

add_executable(starReal_65_1024stars_AB lemonade_star_main_exclvol_AB.cpp)
target_link_libraries(starReal_65_1024stars_AB LeMonADE)

# add_executable(star_readbfm lemonade_star_main_readbfm.cpp)
# target_link_libraries(star_readbfm LeMonADE)

# add_executable(star_make_reactive lemonade_star_main_read_make_reactive.cpp)
# target_link_libraries(star_make_reactive LeMonADE)

# add_executable(star_make_reactive lemonade_star_main_read_make_reactive.cpp)
# target_link_libraries(star_make_reactive LeMonADE)


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

# add_executable(starReal_65_1024stars_box128_AB lemonade_star_main_exclvol_AB.cpp)
# target_link_libraries(starReal_65_1024stars_box128_AB LeMonADE)

# check functionality of code for star networks:
# 'star_main_exclvol_AB.cpp' (August 2020) is a copy of
# 'lemonade_star_main_exclvol_AB.cpp' (January 2020)

# add_executable(starReal_65_1024stars_box128_AB star_main_exclvol_AB.cpp)
# target_link_libraries(starReal_65_1024stars_box128_AB LeMonADE)

# add_executable(starReal_65_1024stars_box128_AB_no_AA_no_BB star_main_exclvol_AB_no_AA_no_BB.cpp)
# target_link_libraries(starReal_65_1024stars_box128_AB_no_AA_no_BB LeMonADE)

# add_executable(starReal_65_1024stars_box128_AB_nn lemonade_star_main_exclvol_AB_nn_repulsion.cpp)
# target_link_libraries(starReal_65_1024stars_box128_AB_nn LeMonADE)

# add_executable(chainAB_128_128_box64_nn lemonade_AB_chain_exclvol_nn_repulsion.cpp)
# target_link_libraries(chainAB_128_128_box64_nn LeMonADE)

# add_executable(chainAB_64_256_box64_nn lemonade_AB_chain_exclvol_nn_repulsion.cpp)
# target_link_libraries(chainAB_64_256_box64_nn LeMonADE)

# add_executable(chainAB_48_341_box64_nn lemonade_AB_chain_exclvol_nn_repulsion.cpp)
# target_link_libraries(chainAB_48_341_box64_nn LeMonADE)

# add_executable(chainAB_32_512_box64_nn lemonade_AB_chain_exclvol_nn_repulsion.cpp)
# target_link_libraries(chainAB_32_512_box64_nn LeMonADE)

# add_executable(chainAB_16_1024_box64_nn_338 lemonade_AB_chain_exclvol_nn_repulsion.cpp)
# target_link_libraries(chainAB_16_1024_box64_nn_338 LeMonADE)

# add_executable(chainAB_24_683_box64_nn_020 lemonade_AB_chain_exclvol_nn_repulsion_mcs_1e5.cpp)
# target_link_libraries(chainAB_24_683_box64_nn_020 LeMonADE)

# add_executable(chainAB_40_410_box64_nn_011 lemonade_AB_chain_exclvol_nn_repulsion_mcs_1e5.cpp)
# target_link_libraries(chainAB_40_410_box64_nn_011 LeMonADE)

# add_executable(chainAB_48_341_box64_nn_009 lemonade_AB_chain_exclvol_nn_repulsion_mcs_1e5.cpp)
# target_link_libraries(chainAB_48_341_box64_nn_009 LeMonADE)

# add_executable(chainAB_56_293_box64_nn_0076 lemonade_AB_chain_exclvol_nn_repulsion_mcs_1e5.cpp)
# target_link_libraries(chainAB_56_293_box64_nn_0076 LeMonADE)

# add_executable(chainAB_64_256_box64_nn_0064 lemonade_AB_chain_exclvol_nn_repulsion_mcs_1e5.cpp)
# target_link_libraries(chainAB_64_256_box64_nn_0064 LeMonADE)

###########################################
#
# box 256, block copolymers
#

# add_executable(chainA12B45_16_65536_box256_nn_0330 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_16_65536_box256_nn_0330 LeMonADE)

# add_executable(chainA12B45_24_43691_box256_nn_0200 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_24_43691_box256_nn_0200 LeMonADE)

# add_executable(chainA12B45_32_32768_box256_nn_0146 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_32_32768_box256_nn_0146 LeMonADE)

# add_executable(chainA12B45_40_26214_box256_nn_0116 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_40_26214_box256_nn_0116 LeMonADE)

# add_executable(chainA12B45_48_21845_box256_nn_0096 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_48_21845_box256_nn_0096 LeMonADE)

# add_executable(chainA12B45_56_18725_box256_nn_0080 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_56_18725_box256_nn_0080 LeMonADE)

# add_executable(chainA12B45_64_16384_box256_nn_0070 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_64_16384_box256_nn_0070 LeMonADE)

###########################################
#
# box 128, block copolymers
#

# add_executable(chainA12B45_16_8192_box128_nn_0330 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_16_8192_box128_nn_0330 LeMonADE)

# add_executable(chainA12B45_24_5461_box128_nn_0208 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_24_5461_box128_nn_0208 LeMonADE)

# add_executable(chainA12B45_40_3277_box128_nn_0112 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_40_3277_box128_nn_0112 LeMonADE)

# add_executable(chainA12B45_48_2731_box128_nn_0087 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_48_2731_box128_nn_0087 LeMonADE)

# add_executable(chainA12B45_56_2341_box128_nn_0080 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_56_2341_box128_nn_0080 LeMonADE)

# add_executable(chainA12B45_64_2048_box128_nn_0066 lemonade_A12_B45_excvol_nn.cpp)
# target_link_libraries(chainA12B45_64_2048_box128_nn_0066 LeMonADE)

####################################################

# add_executable(chainAB_32_2048_box128_nn lemonade_AB_chain_exclvol_nn_repulsion.cpp)
# target_link_libraries(chainAB_32_2048_box128_nn LeMonADE)

####################################################

# add_executable(starReal_65_2048stars_box128_AB lemonade_star_main_exclvol_AB.cpp)
# target_link_libraries(starReal_65_2048stars_box128_AB LeMonADE)

# add_executable(star_readbfm lemonade_star_main_readbfm.cpp)
# target_link_libraries(star_readbfm LeMonADE)

# add_executable(chain_readbfm lemonade_star_main_readbfm.cpp)
# target_link_libraries(chain_readbfm LeMonADE)

# add_executable(chainAB_16_1024_box64_nn_200_restart lemonade_AB_chain_exclvol_nn_repuls_restart.cpp)
# target_link_libraries(chainAB_16_1024_box64_nn_200_restart LeMonADE)

# add_executable(star_make_reactive lemonade_star_main_read_make_reactive.cpp)
# target_link_libraries(star_make_reactive LeMonADE)

# add_executable(star_make_reactive lemonade_star_main_read_make_reactive.cpp)
# target_link_libraries(star_make_reactive LeMonADE)

# add_executable(chainAB_16_8192_box128_nn_restart lemonade_AB_chain_exclvol_nn_repuls_restart.cpp)
# target_link_libraries(chainAB_16_8192_box128_nn_restart LeMonADE)

# add_executable(chainAB_32_4096_box128_nn_restart lemonade_AB_chain_exclvol_nn_repuls_restart.cpp)
# target_link_libraries(chainAB_32_4096_box128_nn_restart LeMonADE)


# add_executable(star_A12_B45_80_box128 lemonade_star_A12_B45_nn_rep.cpp)
# target_link_libraries(star_A12_B45_80_box128 LeMonADE)

# add_executable(star_A12_B45_100stars_box128_AA_BB_AB lemonade_star_A12_B45_nn_rep_only_AB.cpp)
# target_link_libraries(star_A12_B45_100stars_box128_AA_BB_AB LeMonADE)

add_executable(star_A12_B45_100stars_box128_only_AB lemonade_star_A12_B45_nn_rep_only_AB.cpp)
target_link_libraries(star_A12_B45_100stars_box128_only_AB LeMonADE)



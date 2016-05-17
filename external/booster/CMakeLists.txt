# Sub-Libraries
add_subdirectory(bwtlcp)
add_subdirectory(huffrle)

# Booster libraries
add_library(booster
 "ac_defs.c" "arit2_defs.c" "dummy_defs.c" "gdict_defs.c" "huffrle_defs.c" "if_defs.c" "multihuff_defs.c" "range.c" 
 "rc_defs.c" "st_visit.c" "wv_defs.c"
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Binaries
add_executable(boost_compress booster.c)
target_link_libraries(boost_compress booster bwtlcp huffrle m)

add_executable(boost_decompress unboost.c)
target_link_libraries(boost_decompress booster bwtlcp huffrle m)
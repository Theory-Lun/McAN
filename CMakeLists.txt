cmake_minimum_required(VERSION 3.7)
project(McAN C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

option(MUTILTHREADING "enable mutil-threading" ON)
option(TEST_SPLIT_DELETION "enable splitting of deletions" OFF)
option(TEST_REVERSE_CHRONOLOGICAL_ORDER "enable sorting of haplotypes in reverse chronological order" OFF)

add_compile_options(-O3)
if(MUTILTHREADING)
    add_definitions(-DMUTILTHREADING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
if(TEST_SPLIT_DELETION)
    add_definitions(-DTEST_SPLIT_DELETION)
endif()
if(TEST_REVERSE_CHRONOLOGICAL_ORDER)
    add_definitions(-TEST_REVERSE_CHRONOLOGICAL_ORDER)
endif()

include_directories(src/struct/ src/tool hashmap/include cJSON)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_TARGET_EXPORT "Enable exporting of CMake targets. Disable when it causes problems!" OFF)
option(ENABLE_CJSON_VERSION_SO "Enables cJSON so version" OFF)
option(ENABLE_CJSON_TEST "Enable building cJSON test" OFF)
option(ENABLE_CJSON_UNINSTALL "Enable creating uninstall target" OFF)
option(ENABLE_LOCALES "Enable the use of locales" OFF)
add_subdirectory(hashmap)
add_subdirectory(cJSON)

#haplotype network
set(SOURCES
    src/calcDist.c
    src/countryDistribution.c
    src/findHap.c
    src/hapNetAll.c
    src/hapNetDirected.c
    src/main.c
    src/main_siteMask.c
    src/networkAnalysis.c
    src/opt.c
    src/output.c
    src/writeJson.c
    src/compare/compare.c
    src/IO/anaMetaMut.c
    src/IO/readMeta_acc2line_hash.c
    src/IO/readMetadata.c
    src/IO/readMut_acc2line_hash.c
    src/IO/readMutations.c
    src/IO/show_all.c
    src/struct/edges.c
    src/struct/haplotype.c
    src/struct/mutation.c
    src/struct/virus.c
    src/tool/calcDate.c
    src/tool/chomp.c
    src/tool/emptyStr.c
    src/tool/findSubstr.c
    src/tool/queue.c
    src/tool/unionFind.c
    src/tool/optionparser.cpp
    src/vcf2mut.cpp
    )
if(MUTILTHREADING)
    list(APPEND SOURCES src/calcDist_thread.cpp)
endif()

add_executable(McAN ${SOURCES})
target_link_libraries(McAN m hashmap cjson)
if(MUTILTHREADING)
    target_link_libraries(McAN pthread)
endif()

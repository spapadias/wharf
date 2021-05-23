# -------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING WITH METROPOLIS HASTINGS ----------------------------
# ------------------------------------------------ SET UP GOOGLE TESTS -------------------------------------------------

# link to threads library first
find_package(Threads REQUIRED)

# enable ExternalProject cmake module
include(ExternalProject)

# download and install GoogleTest
ExternalProject_Add(
    gtest
    URL https://github.com/google/googletest/archive/release-1.8.1.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # disable install step
    INSTALL_COMMAND ""
)

# get GoogleTest source and binary directories from cmake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest gtest_main)

# set libgtest properties
set_target_properties(
    libgtest PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# create a libgtest_main target to be used as a dependency by test programs
add_library(libgtest_main IMPORTED STATIC GLOBAL)
add_dependencies(libgtest_main gtest_main)

# set libgtest_main properties
set_target_properties(
    libgtest_main PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest_main.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# create a libgmock target to be used as a dependency by test programs
add_library(libgmock IMPORTED STATIC GLOBAL)
add_dependencies(libgmock gtest)

# set libgmock properties
set_target_properties(
    libgmock PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# set include directories
include_directories(
    "${source_dir}/googletest/include"
    "${source_dir}/googlemock/include"
)

# set GoogleTest library variable
set(GTEST_LIBRARIES libgtest libgtest_main libgmock)

# ---------------------------------------------------- TEST DATA -------------------------------------------------------

file(COPY ../experiments/data DESTINATION .)
file(COPY ../experiments/scripts/SNAPtoAdj DESTINATION .)

# -------------------------------------------------- COMPILE FLAGS -----------------------------------------------------

# build options
option(CLANG        "???"               OFF)
option(CILK         "???"               OFF)
option(OPENMP       "???"               OFF)
option(HOMEGROWN    "???"               OFF)
option(DOCK_DEBUG   "Dock debug mode"   OFF)
option(SERIAL       "???"               OFF)

#set(SERIAL ON)
set(DOCK_DEBUG ON)

# build variables
set(OMPFLAGS         "-DOPENMP -fopenmp")
set(CILKFLAGS        "-DCILK -fcilkplus")
set(HGFLAGS          "-DHOMEGROWN -pthread")
set(OPTIMIZATION     "-O3 -g")
set(CFLAGS           "-DEDGELONG -mcx16 ${OPTIMIZATION} -ldl -std=c++17 -march=native -Wall -Wno-subobject-linkage -DUSEMALLOC -DNDEBUG")
set(DOCK_DEBUG_FLAGS "-DEDGELONG -mcx16 ${OPTIMIZATION} -ldl -std=c++17 -march=native -Wall -Wno-subobject-linkage -DUSEMALLOC")

if(CLANG)
    set(CMAKE_CXX_COMPILER "clang++")
    set(PFLAGS ${CILKFLAGS})
elseif(CILK)
    set(CMAKE_CXX_COMPILER "g++")
    set(PFLAGS ${CILKFLAGS})
elseif(OPENMP)
    set(CMAKE_CXX_COMPILER "g++")
    set(PFLAGS ${OMPFLAGS})
elseif(HOMEGROWN)
    set(CMAKE_CXX_COMPILER "g++")
    set(PFLAGS ${HGFLAGS})
elseif(SERIAL)
    set(CMAKE_CXX_COMPILER "g++")
    set(PFLAGS "")
else()
    set(CMAKE_CXX_COMPILER "g++")
    set(PFLAGS ${HGFLAGS})
endif()

# set compile flags
if (DOCK_DEBUG)
    set(CMAKE_CXX_FLAGS "${DOCK_DEBUG_FLAGS} ${PFLAGS} -DDOCK_DEBUG -DDOCK_TIMER -DVERTEXLONG") # -DDOCK_DEBUG -DDOCK_TIMER
else()
    set(CMAKE_CXX_FLAGS "${CFLAGS} ${PFLAGS}")
endif()

# ------------------------------------------------------ TESTS ---------------------------------------------------------

# Dock tests
add_executable(dock dock.cpp)
target_link_libraries(dock ${GTEST_LIBRARIES})
add_test(NAME wharf COMMAND dock)

# Pairing function tests
add_executable(pairings pairings.cpp)
target_link_libraries(pairings ${GTEST_LIBRARIES})
add_test(NAME pairings COMMAND pairings)

# Sampler tests
add_executable(sampler sampler.cpp)
target_link_libraries(sampler ${GTEST_LIBRARIES})
add_test(NAME sampler COMMAND sampler)

# -------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING WITH METROPOLIS HASTINGS ----------------------------


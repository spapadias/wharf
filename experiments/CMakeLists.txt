# -------------------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING -----------------------------------------

option(CLANG     "???" OFF)
option(CILK      "???" OFF)
option(OPENMP    "???" OFF)
option(HOMEGROWN "???" OFF)
option(SERIAL    "???" OFF)

# compile flags options
set(OMPFLAGS     "-DOPENMP -fopenmp")
set(CILKFLAGS    "-DCILK -fcilkplus")
set(HGFLAGS      "-DHOMEGROWN -pthread")
set(OPTIMIZATION "-O3 -g")
set(CFLAGS       "-DEDGELONG  -mcx16 ${OPTIMIZATION} -ldl -std=c++17 -march=native -Wall -Wno-subobject-linkage -DUSEMALLOC ") #

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
set(CMAKE_CXX_FLAGS "${CFLAGS} ${PFLAGS} -DVERTEXLONG")

file(COPY ../experiments/data DESTINATION .)
file(COPY ../experiments/scripts/SNAPtoAdj DESTINATION .)
file(COPY ../experiments/scripts/SNAPtoAdj DESTINATION .)

# MEMORY FOOTPRINT EXPERIMENT EXECUTABLE
#add_executable(memory_footprint src/memory_footprint.cpp)
#add_executable(run_batch_updates src/run_batch_updates.cpp)
add_executable(learn_embeddings src/learn_embeddings.cpp ../training/word2vec.c)
add_executable(compression src/compression.cpp)
# -------------------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING -----------------------------------------
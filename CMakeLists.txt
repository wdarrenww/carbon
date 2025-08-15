cmake_minimum_required(VERSION 3.16)
project(carbon VERSION 1.0.0 LANGUAGES CXX)

# require modern c++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# compiler warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/W4)
endif()

# debug/release flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# main executable
add_executable(carbon
    src/main.cpp
)

target_include_directories(carbon
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)

# enable testing
enable_testing()
add_subdirectory(tests)

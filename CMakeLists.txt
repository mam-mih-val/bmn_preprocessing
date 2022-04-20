cmake_minimum_required(VERSION 3.11)
project(bmn_preprocessing)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(Boost REQUIRED COMPONENTS program_options)

include(cmake_modules/AnalysisTree.cmake)

add_subdirectory(src)

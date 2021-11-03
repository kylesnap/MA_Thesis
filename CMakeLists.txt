cmake_minimum_required(VERSION 3.20)
project(MA_Thesis)

set (CMAKE_CXX_STANDARD 20)

if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    conan_basic_setup()
else()
    message(WARNING "The file conanbuildinfo.cmake doesn't exist, you have to run conan install first")
endif()

add_executable(MA_Thesis main.cpp SimCell.cpp SimCell.h LmOLS.cpp LmOLS.h tests/LmOLS_Test.cpp DesignMat.cpp DesignMat.h tests/DesignMat_Tests.cpp tests/SimCell_Tests.cpp)
target_link_libraries(MA_Thesis ${CONAN_LIBS})
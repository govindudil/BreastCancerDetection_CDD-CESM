cmake_minimum_required(VERSION 3.25)
project( Dissertation )

set(CMAKE_CXX_STANDARD 17)

find_package( OpenCV REQUIRED )
# find_package( jsoncpp REQUIRED )

# Provide the path to the jsoncpp library
set(JSONCPP_INCLUDE_DIRS /usr/include/jsoncpp)
set(JSONCPP_LIBRARIES /usr/lib/x86_64-linux-gnu/libjsoncpp.so)

include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories( ${JSONCPP_INCLUDE_DIRS} )

set( CMAKE_CXX_FLAGS "-g -fopenmp" )

add_executable( Dissertation main.cpp glcmMain.cpp glcmMain.h)
target_link_libraries( Dissertation ${OpenCV_LIBS} )
target_link_libraries( Dissertation ${JSONCPP_LIBRARIES} )
target_link_libraries( Dissertation -fopenmp )
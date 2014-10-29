cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

SET(AppName fc)

project(EvenSimProject)

# FIND_PACKAGE(boost)
INCLUDE_DIRECTORIES("~/Library/boost_1_55_0") 

add_executable(${AppName} cluster.h cluster.cc findcluster.h findcluster.cc input.h input.cc main.cc)

target_link_libraries(${AppName} m pthread)
      

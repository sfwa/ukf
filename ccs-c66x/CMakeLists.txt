CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(c66ukf)

INCLUDE_DIRECTORIES(../include ../c)

set(CMAKE_C_FLAGS "-O3 -Weverything -Werror -Wno-documentation -Wno-padded -Wno-unknown-pragmas -fPIC")

ADD_LIBRARY(c66ukf SHARED cukf.c)
ADD_LIBRARY(fcsukf STATIC cukf.c)

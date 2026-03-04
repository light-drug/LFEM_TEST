# 以下是服务器上Makefile文件的实现
# 编译器设置

all: main example

CC = g++
CCOPT = -fdiagnostics-color=always -O3 -fopenmp
PG = -pg
WARNING = -Wall -Wextra -Werror
CallCC17 = -std=c++17
CallCC20 = -std=c++20
OPENMP = -fopenmp
AddressSanitizer = -fsanitize=address

# 路径设置
# Root = /export/home/panliang/LFEM_TEST/src
Root = src
Self1 = src/basis
Self2 = src/config
Self3 = src/general
Self4 = src/src_test
Self5 = src/mesh
Self6 = src/problems
Self7 = src/operator
Self8 = src/problems2D

SelfExample1 = src/example2D
SelfExample2 = src/example1D

LAPACK_DIR = /export/home/panliang/Package_download/lapack_download/lapack-3.12.0/build
Eigen_DIR = /export/home/panliang/eigen_download/eigen-3.4.0
MPI_Lib_DIR = /export/home/panliang/Package_download/mpi_download/OPENMPI/lib
MPI_Include_DIR = /export/home/panliang/Package_download/mpi_download/OPENMPI/include
MKL_Lib_DIR = /opt/intel/oneapi/mkl/2025.1/lib
MKL_Include_DIR = /opt/intel/oneapi/mkl/2025.1/include
# MKL_Lib_DIR = /export/home/panliang/Package_download/mkl_download/mkl/mkl/2024.2/lib
# MKL_Include_DIR = /export/home/panliang/Package_download/mkl_download/mkl/mkl/2024.2/include

# 包含路径和库路径
INCLUDE_PATH = -I. -I $(Eigen_DIR) \
										-I $(MKL_Include_DIR) \
										-I $(Self1) \
										-I $(Self2) \
										-I $(Self3) \
										-I $(Self4) \
										-I $(Self5) \
										-I $(Self6) \
										-I $(Self7) \
										-I $(Self8)
LIB_PATH = -L.	-L $(MKL_Lib_DIR) \
								-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm -ldl \
								-fopenmp

#####################################################################

# 源文件和目标文件
LIB_SRC_FILES = \
            $(wildcard $(Self1)/*.cpp) \
            $(wildcard $(Self2)/*.cpp) \
            $(wildcard $(Self3)/*.cpp) \
            $(wildcard $(Self4)/*.cpp) \
            $(wildcard $(Self5)/*.cpp) \
            $(wildcard $(Self6)/*.cpp) \
            $(wildcard $(Self7)/*.cpp) \
            $(wildcard $(Self8)/*.cpp)

LIB_OBJ_FILES = $(patsubst %.cpp, build/%.o, $(LIB_SRC_FILES))
  # 自动生成对应的 .o 文件列表

# examples（每个cpp一个可执行文件）
EXAMPLE_SRC_FILES = $(wildcard $(SelfExample1)/*.cpp) \
                    $(wildcard $(SelfExample2)/*.cpp)
EXAMPLE_BINS = $(patsubst $(SelfExample1)/%.cpp, bin/%, $(wildcard $(SelfExample1)/*.cpp)) \
               $(patsubst $(SelfExample2)/%.cpp, bin/%, $(wildcard $(SelfExample2)/*.cpp))

# main 自己的源
MAIN_SRC = main.cpp
MAIN_OBJ = $(patsubst %.cpp, build/%.o, $(MAIN_SRC))

# 编译可执行文件
main: $(MAIN_OBJ) $(LIB_OBJ_FILES)
	$(CC) $(CCOPT) $(CallCC20) $(MAIN_OBJ) $(LIB_OBJ_FILES) $(INCLUDE_PATH) $(LIB_PATH) -o main

bin/%: build/$(SelfExample1)/%.o $(LIB_OBJ_FILES)
	@mkdir -p $(dir $@)
	$(CC) $(CCOPT) $(CallCC20) $< $(LIB_OBJ_FILES) $(INCLUDE_PATH) $(LIB_PATH) -o $@

bin/%: build/$(SelfExample2)/%.o $(LIB_OBJ_FILES)
	@mkdir -p $(dir $@)
	$(CC) $(CCOPT) $(CallCC20) $< $(LIB_OBJ_FILES) $(INCLUDE_PATH) $(LIB_PATH) -o $@

# 通用规则
build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CC) $(CCOPT) $(CallCC20) $(INCLUDE_PATH) -c $< -o $@


# 清除生成的文件
clean:
	rm -rf build bin
	rm -f main *.o *.out fort.* main.exe pch.hpp.gch

clean_o:
	find . -type f -name "*.o" -exec rm -f {} \;

clc:
	rm -rf build bin
	find . -type f -name "*.o" -exec rm -f {} \;
	rm -f main *.o fort.* main.exe pch.hpp.gch

example: $(EXAMPLE_BINS)

EXAMPLE_OBJ_FILES = $(patsubst %.cpp, build/%.o, $(EXAMPLE_SRC_FILES))
.SECONDARY: $(EXAMPLE_OBJ_FILES)

.PHONY: clean


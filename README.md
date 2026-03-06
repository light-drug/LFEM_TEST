# LFEM_TEST (Light's FEM)

A lightweight tensor-mesh solver for partial differential equations (work in progress).

This is a simple tensor-mesh based high-order FD/FV/DG framework for hyperbolic equations.

## Dependencies

The following libraries are required:

Eigen: https://gitlab.com/libeigen/eigen/-/releases/3.4.1  
OpenMP: https://www.openmp.org  
Intel MKL: https://www.intel.cn/content/www/cn/zh/developer/articles/guide/intel-math-kernel-library-intel-mkl-2019-install-guide.html  

## Compilation

The code has been tested on the computer **c0196** at Xiamen University.

You can compile the project using the provided `Makefile`.  
Before compiling, please modify the paths of the following variables in the Makefile:

- `Eigen_DIR`
- `MPI_Lib_DIR`
- `MPI_Include_DIR`
- `MKL_Lib_DIR`
- `MKL_Include_DIR`

### Build

Run

```bash
make


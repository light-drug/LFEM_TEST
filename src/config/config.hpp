#ifndef QUEST_CONFIG_HPP
#define QUEST_CONFIG_HPP

#define QUEST_USE_DOUBLE 
#define QUEST_USE_INT
#define QUEST_USE_OPENMP
#define QUEST_USE_MKL
#define QUEST_USE_EIGEN
#define QUEST_SPACE_DIM 1
#define QUEST_MESH_DIM 1

#include <iostream>  // Add this header
#include <string>    // Add this header
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>


#ifdef QUEST_USE_SINGLE
using real_t = float;
#elif defined QUEST_USE_DOUBLE
using real_t = double;
#endif

using creal_t = std::complex<real_t>;

#ifdef QUEST_USE_INT
using QUEST_Int = int;
#elif defined QUEST_USE_SIZET
using QUEST_Int = size_t;
#endif

#ifdef QUEST_USE_OPENMP
#include <omp.h>
#endif

#ifdef QUEST_USE_EIGEN 
  using IntMatrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
  using IntVector = Eigen::Vector<int, Eigen::Dynamic>;
  using Matrix = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Vector<real_t, Eigen::Dynamic>;
  using DiagnalMatrix = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic>;
  using SparseMatrix = Eigen::SparseMatrix<real_t>;

  using cMatrix = Eigen::Matrix<creal_t, Eigen::Dynamic, Eigen::Dynamic>;
  using cVector = Eigen::Vector<creal_t, Eigen::Dynamic>;
  using cDiagnalMatrix = Eigen::DiagonalMatrix<creal_t, Eigen::Dynamic>;
  using cSparseMatrix = Eigen::SparseMatrix<creal_t>;
#endif 

#ifdef QUEST_USE_MKL
#include <mkl.h>
#include <Eigen/PardisoSupport>
#endif 

#ifdef QUEST_USE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#endif 


#if defined(_MSC_VER) && defined(QUEST_SHARED_BUILD)
#ifdef QUEST_EXPORTS
#define QUEST_EXPORT __declspec(dllexport)
#else
#define QUEST_EXPORT __declspec(dllimport)
#endif
#else
#define QUEST_EXPORT
#endif


#define PAUSE() do { \
  std::cout << "Press Enter to continue... " << std::endl; \
  getchar(); \
} while (0)


class Color {
public:
  static constexpr const char* red   = "\033[31m";
  static constexpr const char* green = "\033[32m";
  static constexpr const char* yellow= "\033[33m";
  static constexpr const char* blue  = "\033[34m";
  static constexpr const char* reset = "\033[0m";
};


#endif  // QUEST_CONFIG_HPP

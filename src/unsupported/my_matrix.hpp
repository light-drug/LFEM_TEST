#ifndef LFEM_MY_MATRIX_HPP
#define LFEM_MY_MATRIX_HPP

#include "config.hpp"
#include "error.hpp"

namespace LFEM {

class my_matrix 
{
private:

  real_t* data_;
  LFEM_Int rows_;
  LFEM_Int cols_;
  LFEM_Int size_;

public:
  
  my_matrix(LFEM_Int& rows, LFEM_Int& cols)
          : rows_(rows), cols_(cols), size_(rows * cols)  {
    data_ = new T[size_]();
  };

  my_matrix(LFEM_Int& rows, LFEM_Int& cols, real_t* data)
          : rows_(rows), cols_(cols), size_(rows * cols) {
    data_ = new T[size_];
    std::copy(data, data + size_, data_);
  };

  my_matrix(const my_matrix& other)
    : rows_(other.rows_), cols_(other.cols_), size_(other.size_) {
    data_ = new T[size_];
    std::copy(other.data_, other.data_ + size_, data_);
  };

  void reshape(LFEM_Int rows, LFEM_Int cols) {
    LFEM_VERIFY(size_ == (rows * cols), " Matrix reshape Error, size != (rows * cols)");
    rows_ = rows;
    cols_ = cols;
  };

  T& operator()(LFEM_Int i, LFEM_Int j) {
    LFEM_VERIFY(i >= 0 && i < rows_ && j >= 0 && j < cols_, "Matrix index out of bounds");
    return data_[i * cols_ + j];
  };

  const T& operator()(LFEM_Int i, LFEM_Int j) const {
    LFEM_VERIFY(i >= 0 && i < rows_ && j >= 0 && j < cols_, "index out of bounds");
    return data_[i * cols_ + j];
  };

  // this: M, a * M
  void LeftMulti(const T& a, const my_matrix* result);

  // this: M, M * a
  void RightMulti(const T& a, const my_matrix* result);

  // this: M, A * M
  void LeftMulti(const my_matrix& A, const my_matrix* result);

  // this: M, M * A
  void RightMulti(const my_matrix& A, const my_matrix* result);

  // this: M, a + A = A + a
  void Add(const T& a, const my_matrix* result);

  // this: M, M + A = A + M
  void Add(const T& A, const my_matrix* result);

  // this: M, A^T * M
  void LeftTMulti(const my_matrix& A, const my_matrix* result);

  // this: M, M * A^T
  void RightTMulti(const my_matrix& A, const my_matrix* result);

  // this: M, A^T * M^T
  void LeftTMultiT(const my_matrix& A, const my_matrix* result);

  // this: M, M^T * A^T
  void RightTMultiT(const my_matrix& A, const my_matrix* result);

  const LFEM_Int& rows() const { return rows_; };
  const LFEM_Int& cols() const { return cols_; };
  const LFEM_Int& size() const { return size_; };

  ~my_matrix() {delete[] data_; };
};

} // namespace LFEM

// LFEM_FESPACE_HPP

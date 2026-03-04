/*
 * Copyright (C) 2025 Xiamen University
 *
 * @Author: Liang Pan
 * @Date:   2025-04-12
 * @Last Modified by: Liang Pan
 * @Last Modified time: 2025-04-12
 */

#include "basis.hpp"

namespace QUEST {
  
namespace Internal {
  
real_t basisphi(int j,real_t x,int d) {
  real_t result;
  switch (j) {
  case 0:
      switch (d)
      { case 0: result = 1.e0; break;
      default: result = 0.e0; break; }
      break;
  case 1:
      switch (d)
      { case 0: result = x; break;
      case 1: result = 1.e0; break;
      default: result = 0.e0; break; }
      // result  = result * 2;
      break;
  case 2:
      switch (d)
      { case 0: result = 3.e0 * std::pow(x,2) - 1.e0 / 4.e0; break;
      case 1: result = 6.e0 * x; break;
      case 2: result = 6.e0; break;
      default: result = 0.e0; break; }
      // result = result * 4;
      break;
  case 3:
      switch (d) 
      { case 0: result = 20.e0 * std::pow(x,3) - 3.e0 * x; break;
      case 1: result = 60.e0 * std::pow(x,2) - 3.e0; break;
      case 2: result = 120.e0 * x; break;
      case 3: result = 120.e0; break;
      default:result = 0.e0; break; }
      break;
  case 4:
      switch (d) 
      { case 0: result = 70.e0 * std::pow(x,4) - 15.e0 * std::pow(x,2) + 3.e0 / 8.e0;break;
      case 1: result = 280.e0 * std::pow(x,3) - 30.e0 * x;break;
      case 2: result = 840.e0 * std::pow(x,2) - 30.e0;break;
      case 3: result = 1680.e0 * x;break;
      case 4: result = 1680.e0; break;
      default: result = 0.e0; break; }
      break;
  }
  return result;
};
  
} // namespace Internal

BasisFunction1D::BasisFunction1D(const BasisType& basistype) : basistype_(basistype) {
  using namespace QUEST;
  k1D_ = static_cast<int>(basistype_) - static_cast<int>(BasisType::P0) + 1;
  if (static_cast<int>(basistype_) < static_cast<int>(BasisType::P0) ||
      static_cast<int>(basistype_) > static_cast<int>(BasisType::P4)) {
    QUEST_ERROR(" Unsupported BasisType !!! ");
  }
};

BasisFunction1D& BasisFunction1D::dx(const int& dx_order) {
  dx_order_ = dx_order;
  return *this;
};

void BasisFunction1D::Map(const real_t& x, Vector* result) {
  result->resize(k1D_); result->setZero();
  if (static_cast<int>(basistype_) >= static_cast<int>(BasisType::P0) ||
      static_cast<int>(basistype_) <= static_cast<int>(BasisType::P4)) {
    for (int i = 0; i < k1D_; ++i) {
      (*result)(i) = QUEST::Internal::basisphi(i, x, dx_order_);
    }
  } else {
    QUEST_ERROR(" Unsupported BasisType !!! ");
  }
};

void BasisFunction1D::Map(const Vector& x, Matrix* result) {
  int size_ = x.size();
  result->resize(k1D_, size_); result->setZero();
  if (static_cast<int>(basistype_) >= static_cast<int>(BasisType::P0) ||
      static_cast<int>(basistype_) <= static_cast<int>(BasisType::P4)) {
    for (int i = 0; i < k1D_; i++) {
      for (int qx = 0; qx < size_; qx++) {
        (*result)(i, qx) = QUEST::Internal::basisphi(i, x(qx), dx_order_);
      }      
    }
  } else {
    QUEST_ERROR(" Unsupported BasisType !!! ");
  }
};

const int& BasisFunction1D::getpolydim() const {
  return k1D_;
};

const int& BasisFunction1D::getk1D() const {
  return k1D_;
};

BasisFunction2D::BasisFunction2D(const BasisType& basistype) : basistype_(basistype) {
  using namespace QUEST;
  if (static_cast<int>(basistype_) < static_cast<int>(BasisType::Q0))
  {
    k1D_ = static_cast<int>(basistype_) - static_cast<int>(BasisType::P0) + 1;
    polydim_ = k1D_ * (k1D_ + 1) / 2;
  } 
  else if (static_cast<int>(basistype_) >= static_cast<int>(BasisType::Q0))
  {
    k1D_ = static_cast<int>(basistype_) - static_cast<int>(BasisType::Q0) + 1;
    polydim_ = k1D_ * k1D_;
  }
  
};

BasisFunction2D& BasisFunction2D::dx(const int& dx_order) {
  dx_order_ = dx_order;
  return *this;
};

BasisFunction2D& BasisFunction2D::dy(const int& dy_order) {
  dy_order_ = dy_order;
  return *this;
};

void BasisFunction2D::Map(const real_t& x,
        const real_t& y,
        Vector* result)
{
  result->resize(polydim_); result->setZero();
  if (static_cast<int>(basistype_) < static_cast<int>(BasisType::Q0))
  {
    int index = 0;
    for (int j = 0; j < k1D_; j++){
      for (int i = 0; i < (k1D_ - j); i++){
        (*result)(index) = Internal::basisphi(i, x, dx_order_)
                        * Internal::basisphi(j, y, dy_order_);
        index++;
      } 
    }
  } else if (static_cast<int>(basistype_) >= static_cast<int>(BasisType::Q0))
  {
    int index = 0;
    for (int j = 0; j < k1D_; j++){
      for (int i = 0; i < k1D_; i++){
        (*result)(index) = Internal::basisphi(i, x, dx_order_)
                        * Internal::basisphi(j, y, dy_order_);
        index++;
      } 
    }
  }
};

void BasisFunction2D::Map(const Matrix& Coor,
                        Matrix* result)
{
  const int size = Coor.cols();
  result->resize(polydim_, size); result->setZero();
  if (static_cast<int>(basistype_) < static_cast<int>(BasisType::Q0))
  {
    int index = 0;
    for (int j = 0; j < k1D_; j++){
      for (int i = 0; i < (k1D_ - j); i++){
        for (int qx = 0; qx < size; qx++) {
          (*result)(index, qx) = Internal::basisphi(i, Coor(0, qx), dx_order_)
                                  * Internal::basisphi(j, Coor(1, qx), dy_order_);
        }   
        index++;
      } 
    }
  } else if (static_cast<int>(basistype_) >= static_cast<int>(BasisType::Q0))
  {
    int index = 0;
    for (int j = 0; j < k1D_; j++){
      for (int i = 0; i < k1D_; i++){
        for (int qx = 0; qx < size; qx++) {
          (*result)(index, qx) = Internal::basisphi(i, Coor(0, qx), dx_order_)
                                  * Internal::basisphi(j, Coor(1, qx), dy_order_);
        }   
        index++;
      } 
    }
  }
};

const int& BasisFunction2D::getpolydim() const
{
  return polydim_;
};

const int& BasisFunction2D::getk1D() const
{
  return k1D_;
};

} // namespace QUEST



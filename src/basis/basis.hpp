/*
 * Copyright (C) 2025 Xiamen University
 *
 * @Author: Liang Pan
 * @Date:   2025-04-12
 * @Last Modified by: Liang Pan
 * @Last Modified time: 2025-04-12
 */

#ifndef QUEST_BASIS_HPP
#define QUEST_BASIS_HPP

#include "config.hpp"
#include "error.hpp"

namespace QUEST {
  
enum class BasisType {
  P0 = 0, 
  P1 = 1, 
  P2 = 2, 
  P3 = 3, 
  P4 = 4, 
  Q0 = 15,
  Q1 = 16, 
  Q2 = 17, 
  Q3 = 18, 
  Q4 = 19,
};

namespace Internal {

real_t basisphi(int j,real_t x,int d);
  
} // namespace Internal

class BasisFunction1D {
  public:
  
    ~BasisFunction1D() = default;

    BasisFunction1D(const BasisType& basistype);

    BasisFunction1D& dx(const int& dx_order);

    void Map(const real_t& x, Vector* result);

    void Map(const Vector& x, Matrix* result);

    const int& getpolydim() const;

    const int& getk1D() const;
  
  private:
  
    int dx_order_ = 0;
    BasisType basistype_;
    int k1D_;
    
};

class BasisFunction2D 
{
  public:
  
    ~BasisFunction2D() = default;

    BasisFunction2D(const BasisType& basistype);

    BasisFunction2D& dx(const int& dx_order);

    BasisFunction2D& dy(const int& dy_order);

    void Map(const real_t& x,
            const real_t& y,
            Vector* result);

    void Map(const Matrix& Coor,
            Matrix* result);

    const int& getpolydim() const;

    const int& getk1D() const;
  
  private:
  
    int dx_order_ = 0;
    int dy_order_ = 0;
    BasisType basistype_;
    int polydim_;
    int k1D_;
};

} // namespace LFE 

#endif // QUEST_BASIS_HPP 

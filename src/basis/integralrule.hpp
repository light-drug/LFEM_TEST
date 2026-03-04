#ifndef QUEST_INTEGRALRULE_HPP
#define QUEST_INTEGRALRULE_HPP

#include "config.hpp"
#include "error.hpp"

#include <string_view>

namespace QUEST {

enum class QuadratureType {
  kGaussLegendre = 0,
  kGaussLobatto = 1,
  kGaussChebyshev = 2,
};

namespace Internal {

  // ************* 参考单元 [-1/2, 1/2] ********************* //
  void initializeGaussLegendre(const int n, Vector* x, Vector* w);

  // ************** 参考单元 [-1/2, 1/2] ******************** //
  void initializeGaussLobatto(const int n, Vector* x, Vector* w);

  // ************** 参考单元 [-1/2, 1/2], w(x) != 1 ********* // 
  void initializeGaussChebyshev(const int n, Vector* x, Vector* w);

} // namespace internal

  // void initializeQuadradure(const int& n,
  //   const QuadratureType& quatype,
  //   Vector* q1_,
  //   Vector* w1_);

  void initializeQuadradure(const int& n,
    const QuadratureType& quatype,
    Vector* q1_,
    Vector* w1_);     
    
  std::string name(const QuadratureType& type);

} // namespace QUEST

#endif  // INTEGRALRULE_HPP

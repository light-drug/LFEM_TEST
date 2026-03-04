#ifndef QUEST_HYPERBOLIC_HPP
#define QUEST_HYPERBOLIC_HPP

#include "config.hpp"
#include "basis.hpp"
#include "Tensormesh.hpp"
#include "fespace.hpp"
#include "vector_overload.hpp"

// Compute the following step
//  ∫_T F(u):∇v,   -∫_e F̂(u)⋅[[v]]
// 其中F是一个非线性项或者是一个显性项
// 主要用作显式计算

namespace QUEST
{

class HyperbolicFormBase
{
protected:

  const fespace1D* fe_;

  int num_equations_;

  int dim_;

public: 

  HyperbolicFormBase(const fespace1D* fe);

  const fespace1D* getfespace() const;

  virtual void init_modal(const std::vector<Matrix>& u_modal);

  virtual void ComputeFu(const std::vector<Matrix>& u_nodal, const int& dimid, std::vector<Matrix>& Fu) const = 0;

  void ComputeFuTest(const std::vector<Matrix>& u_nodal, std::vector<Matrix>& Temp) const = 0; 

  virtual void ComputeFlux(const std::vector<real_t>& u_int, 
                          const std::vector<real_t>& u_ext,
                          const Vector& normal, 
                          std::vector<real_t>& flux) const = 0;

  virtual void ComputeFluxAll(const std::vector<Matrix>& u_nodal, 
                          std::vector<Matrix>& flux_int,
                          std::vector<Matrix>& flux_ext) const = 0;

  void ComputeFluxTest(const std::vector<Matrix>& flux_int, const std::vector<Matrix>& flux_ext, 
                      std::vector<Matrix>& Temp) const = 0;
};


class UpwindFluxPeriodicForm: public HyperbolicFormBase
{
private:
  const Vector Coeff_;

  
public:

  UpwindFluxPeriodicForm(const fespace1D* fe, const Vector Coeff_);

  void ComputeFu(const std::vector<Matrix>& u_nodal, const int& dimid, std::vector<Matrix>& Fu) const = 0;

  void ComputeFlux(const std::vector<real_t>& u_int, 
                          const std::vector<real_t>& u_ext,
                          const Vector& normal, 
                          std::vector<real_t>& flux) const = 0;

  void ComputeFluxAll(const std::vector<Matrix>& u_nodal, 
                          std::vector<Matrix>& flux_int,
                          std::vector<Matrix>& flux_ext) const = 0;

  ~UpwindFluxPeriodicForm() = default;

};


class DownwindFluxPeriodicForm: public HyperbolicFormBase
{
public:

  DownwindFluxForm(/* args */);

  ~DownwindFluxForm();

};


  
} // namespace QUEST




#endif // QUEST_HYPERBOLIC_HPP

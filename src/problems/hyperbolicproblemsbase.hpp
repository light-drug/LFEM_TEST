#ifndef QUEST_PROBLEMSBASE_HPP
#define QUEST_PROBLEMSBASE_HPP

#include "fespace.hpp"
#include "config.hpp"
#include <vector>

#include "RK_table.hpp"

namespace QUEST
{

class HyperbolicProblems1DBase
{
public:

  HyperbolicProblems1DBase(const TensorMesh1D* mesh1D,
                          const fespace1D* fe, 
                          const EX_TVDRK* rk_table);

  virtual ~HyperbolicProblems1DBase() = default;

  virtual void init() = 0;
  virtual void setdt(real_t* dt) = 0;
  virtual real_t max_speed_compute(const std::vector<Matrix>& u_modal) = 0;
  virtual void updateAll(const real_t& Trun, const real_t& dt) = 0;
  virtual void fu_compute(const std::vector<Matrix>& u_modal,
                          std::vector<Matrix>* fu_nodal) = 0;
  virtual void Lu_compute(const std::vector<Matrix>& u_modal,
                          const real_t& Trun, const real_t& dt,
                          std::vector<Matrix>* Lu) = 0;
  virtual void fluxint_compute(const std::vector<Matrix>& u_modal,
                                std::vector<Matrix>* flux_int) = 0;
  virtual void fluxext_compute(const std::vector<Matrix>& u_modal,
                                const std::vector<Matrix>& Dirichlet,
                                std::vector<Matrix>* flux_ext) = 0;
  virtual void numerical_flux(const Vector& u_L, const Vector& u_R,
                              const real_t& normal, 
                              const real_t& max_speed,
                              Vector* flux) = 0;

  virtual void setcfl(const real_t& cfl);

  const std::vector<Matrix>& getumodal() const;

protected: 
  
  const TensorMesh1D* mesh1D_;
  const fespace1D* fe_;
  const EX_TVDRK* rk_table_;

  real_t cfl_;

  int num_equations_;

  std::vector<Matrix> u_modal_;
  std::vector<std::vector<Matrix>> u_modal_stages_;
  std::vector<std::vector<Matrix>> Lu_stages_;
  std::vector<Matrix> flux_int_;
  std::vector<Matrix> flux_ext_;

};

class HyperbolicProblems2DBase
{
public:

  HyperbolicProblems2DBase(const TensorMesh2D* mesh2D,
                           const fespace2D* fe,
                           const EX_TVDRK* rk_table);

  virtual ~HyperbolicProblems2DBase() = default;

  virtual void init() = 0;
  virtual void setdt(real_t* dt) = 0;
  virtual real_t max_speed_compute(const std::vector<Matrix>& u_modal) = 0;
  virtual void updateAll(const real_t& Trun, const real_t& dt) = 0;
  virtual void flux_compute(const std::vector<Matrix>& u_modal,
                            std::vector<Matrix>* fu_nodal,
                            std::vector<Matrix>* fv_nodal) = 0;
  virtual void Lu_compute(const std::vector<Matrix>& u_modal,
                          const real_t& Trun,
                          const real_t& dt,
                          std::vector<Matrix>* Lu) = 0;
  virtual void fluxint_compute(const std::vector<Matrix>& u_modal,
                               std::vector<Matrix>* flux_int) = 0;
  virtual void fluxext_compute(const std::vector<Matrix>& u_modal,
                               const std::vector<Matrix>& Dirichlet,
                               std::vector<Matrix>* flux_ext) = 0;
  virtual void numerical_flux(const Vector& u_L, const Vector& u_R,
                              const Eigen::Vector2d& normal,
                              const real_t& max_speed,
                              Vector* flux) = 0;

  virtual void setcfl(const real_t& cfl);

  const std::vector<Matrix>& getumodal() const;

protected:

  const TensorMesh2D* mesh2D_;
  const fespace2D* fe_;
  const EX_TVDRK* rk_table_;

  real_t cfl_;

  int num_equations_;

  std::vector<Matrix> u_modal_;
  std::vector<std::vector<Matrix>> u_modal_stages_;
  std::vector<std::vector<Matrix>> Lu_stages_;
  std::vector<Matrix> flux_int_;
  std::vector<Matrix> flux_ext_;

};

// class IMEXDriftDiffusion1D: public HyperbolicProblems1DBase
// {
// public:

//   IMEXDriftDiffusion1D(/* args */);

//   ~IMEXDriftDiffusion1D();
// };

} // namespace QUEST



#endif // QUEST_PROBLEMSBASE_HPP

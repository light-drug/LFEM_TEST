#ifndef QUEST_EULER1D_ACCTEST_HPP
#define QUEST_EULER1D_ACCTEST_HPP

#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "hyperbolicproblemsbase.hpp"
#include "config.hpp"
#include "error.hpp"

namespace QUEST
{

class EulerPerodicalBase: public HyperbolicProblems1DBase
{
protected: 
  
  static real_t cfl_;
  static real_t gamma_;
  static real_t pi_;

  Matrix pre_nodal_;
  Matrix vx_nodal_;
  real_t max_speed_;
  std::vector<Matrix> fk_;

public:
  EulerPerodicalBase(const fespace1D* fe, const int& t_order);

  ~EulerPerodicalBase() = default;

  static Matrix rho_init(const Matrix& x);
  static Matrix vx_init(const Matrix& x);
  static Matrix rhovx_init(const Matrix& x); 
  static Matrix pre_init(const Matrix& x);
  static Matrix erg_init(const Matrix& x);

  void init() override;
  void setdt(real_t* dt) override;
  void time_stepping(real_t& dt) override;
  void update() override;
  virtual void fk_compute();
  virtual void flux_compute(const real_t& rho_L, const real_t& vx_L, const real_t& pre_L, const real_t& E_L,
                  const real_t& rho_R, const real_t& vx_R, const real_t& pre_R, const real_t& E_R,
                  const real_t& normal, 
                  Vector* flux);
  virtual void fluxint_compute();
  virtual void fluxext_compute();

};

} // namespace QUEST

#endif  // QUEST_EULER1D_ACCTEST_HPP

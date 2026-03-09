#ifndef QUEST_EULER1D_ACCTEST_HPP
#define QUEST_EULER1D_ACCTEST_HPP

#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "hyperbolicproblemsbase.hpp"
#include "config.hpp"
#include "error.hpp"

namespace QUEST
{

class Euler1D_DG_TVDRK
 : public HyperbolicProblems1DBase
{
protected: 
  
  real_t cfl_;
  real_t gamma_;
  real_t pi_;

public:
  Euler1D_DG_TVDRK(const TensorMesh1D* mesh1D,
                  const fespace1D* fe, 
                  const EX_TVDRK* rk_table);
  // u_modal_: rho, rhovx, E
  // fu: rhovx, rho*vx^2 + pre, (E + pre) * vx

  ~Euler1D_DG_TVDRK() override = default;

  virtual Matrix rho_init(const Matrix& x);
  virtual real_t rho_init(const real_t& x);
  virtual Matrix vx_init(const Matrix& x);
  virtual real_t vx_init(const real_t& x);
  virtual Matrix pre_init(const Matrix& x);
  virtual real_t pre_init(const real_t& x);

  virtual real_t rho_bc(const real_t& x, const real_t& t);
  virtual real_t vx_bc(const real_t& x, const real_t& t);
  virtual real_t pre_bc(const real_t& x, const real_t& t);

  void init() override;
  void setdt(real_t* dt) override;
  real_t max_speed_compute(const std::vector<Matrix>& u_modal) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void fu_compute(const std::vector<Matrix>& u_modal,
                  std::vector<Matrix>* fu_nodal) override;
  void Lu_compute(const std::vector<Matrix>& u_modal,
                  const real_t& Trun, 
                  const real_t& dt,
                  std::vector<Matrix>* Lu) override;
  void fluxint_compute(const std::vector<Matrix>& u_modal,
                      std::vector<Matrix>* flux_int) override;
  void fluxext_compute(const std::vector<Matrix>& u_modal,
                      const std::vector<Matrix>& Dirichlet,
                      std::vector<Matrix>* flux_ext) override;
  void numerical_flux(const Vector& u_L,
                      const Vector& u_R,
                      const real_t& normal,
                      const real_t& max_speed,
                      Vector* flux) override;

  virtual void setgamma(const real_t& gamma);
  virtual void setlimiter(const bool enable_limiter, const real_t theta = 1.5e0);

protected:
  void tvd_limiter_characteristic(std::vector<Matrix>* u_modal);
  real_t minmod3(const real_t a, const real_t b, const real_t c) const;

  bool use_limiter_ = false;
  real_t limiter_theta_ = 1.5e0;

};

class Euler1D_DG_TVDRK_period 
  : public Euler1D_DG_TVDRK
{
public:
  Euler1D_DG_TVDRK_period(const TensorMesh1D* mesh1D,
                          const fespace1D* fe, 
                          const EX_TVDRK* rk_table);
  ~Euler1D_DG_TVDRK_period() override = default;

   void Lu_compute(const std::vector<Matrix>& u_modal,
                  const real_t& Trun, 
                  const real_t& dt,
                  std::vector<Matrix>* Lu) override;
  void fluxext_compute(const std::vector<Matrix>& u_modal,
                      const std::vector<Matrix>& Dirichlet,
                      std::vector<Matrix>* flux_ext) override;

  real_t rho_bc(const real_t& x, const real_t& t) override;
  real_t vx_bc(const real_t& x, const real_t& t) override;
  real_t pre_bc(const real_t& x, const real_t& t) override;
  
  virtual Matrix rho_real(const Matrix& x, const real_t& t);
  virtual real_t rho_real(const real_t& x, const real_t& t);
  virtual Matrix vx_real(const Matrix& x, const real_t& t);
  virtual real_t vx_real(const real_t& x, const real_t& t);
  virtual Matrix pre_real(const Matrix& x, const real_t& t);
  virtual real_t pre_real(const real_t& x, const real_t& t);

};

} // namespace QUEST

#endif  // QUEST_EULER1D_ACCTEST_HPP

#ifndef QUEST_KIETICLD_DG_IMEX_EX_HPP
#define QUEST_KIETICLD_DG_IMEX_EX_HPP

#include <iostream>
#include <vector>

#include "config.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "RK_table.hpp"
#include "fespace.hpp"
#include "Tensormesh.hpp"


namespace QUEST
{

// 这个是Riemann problem for telegraph equation
// example 5 in Schur by Peng
class Kinetic1D_LD_DG_IMEX_EX
{
protected:
  const KineticTensorMesh1D* mesh1D_;
  const fespace1D* fe_;
  const IMEX_RK* rk_table_;

  real_t eps_;
  real_t eps2_;
  real_t pi_;
  real_t Chy_;
  real_t Cdif_;
  real_t CR_;
  real_t beta1_; // - 0.e0 for central flux, +-0.5e0 for alternative flux
  real_t sigmas_;
  int NTH_;
  int stages_;

  struct model_data_
  {
    Matrix rho;
    std::vector<Matrix> g;
  };
  
  model_data_ kinetic_modal_;
  std::vector<model_data_> kinetic_modal_stages_;
  std::vector<Matrix> ah_;
  std::vector<std::vector<Matrix>> bh_;
  std::vector<Matrix> dh_;
  std::vector<std::vector<Matrix>> sh_;
  // std::vector<Matrix> rho_stages_;
  // std::vector<std::vector<Matrix>> g_stages_;


public:
  Kinetic1D_LD_DG_IMEX_EX(const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table);
  virtual ~Kinetic1D_LD_DG_IMEX_EX() = default;

  virtual void seteps(const real_t& knu);
  virtual void setChy(const real_t& Chy);
  virtual void setCdif(const real_t& Cdif);
  virtual void setCR(const real_t& CR); // penalty term for LDG boundary flux
  virtual void setbeta1(const real_t& beta1);
  virtual void setsigmas(const real_t& sigmas);
  virtual void setNTH(const int& NTH);

  virtual void init();
  virtual void setdt(real_t* dt);
  // virtual void updaterho(const real_t& Trun, const real_t& dt);
  // virtual void updateg(const real_t& Trun, const real_t& dt);
  // virtual void updateAll(const real_t& Trun, const real_t& dt);
  
  virtual Matrix rho_init(const Matrix& x);
  virtual real_t rho_init(const real_t& x);
  virtual Matrix f_init(const Matrix& x, const real_t& v);
  virtual real_t f_init(const real_t& x, const real_t& v);
  virtual Matrix g_init(const Matrix& x, const real_t& v);
  virtual real_t g_init(const real_t& x, const real_t& v);
  // virtual real_t f_bc(const real_t& x, const real_t& v, const real_t& t,
  //                     const Matrix& rho_modal, const std::vector<Matrix>& g_modal);
  // virtual real_t g_bc(const real_t& x, const real_t& v, const real_t& t,
  //                     const Matrix& rho_modal, const std::vector<Matrix>& g_modal);
  virtual real_t fL_bc(const int& j, const real_t& t,
                      const model_data_& modal);
  virtual real_t fR_bc(const int& j, const real_t& t,
                      const model_data_& modal);
  virtual real_t gL_bc(const int& j, const real_t& t,
                      const model_data_& modal, const real_t& rho_L);
  virtual real_t gR_bc(const int& j, const real_t& t,
                      const model_data_& modal, const real_t& rho_R);
  virtual real_t rho_numericalbc(const real_t& x, const real_t& t,
                      const model_data_& modal);
  virtual void updateAll(const real_t& Trun, const real_t& dt);

  virtual const Matrix& getrho_modal() const;
  virtual const std::vector<Matrix>& getg_modal() const;

  virtual void fluxint_compute(const Matrix& modal, 
                          const real_t& beta,
                          Matrix* flux_int);
  virtual void fluxext_compute(const Matrix& modal, 
                          const real_t& beta,
                          const Matrix& Dirichlet,
                          Matrix* flux_ext);

  virtual void fluxint_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          Matrix* flux_int);
  virtual void fluxext_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          const Matrix& Dirichlet,
                          Matrix* flux_ext);

  virtual void ah_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          const Matrix& boundary_flux,
                          Matrix* ah);
  virtual void ah_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux);

  virtual void bh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          const std::vector<Matrix>& boundary_flux,
                          std::vector<Matrix>* bh);
  virtual void bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux);

  virtual void dh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          const Matrix& boundary_flux,
                          Matrix* dh);
  virtual void dh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          Matrix* boundary_flux);

  virtual void sh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* sh);
};


class Kinetic1D_LD_DG_IMEX_EX_IsentropicBC : 
  public virtual Kinetic1D_LD_DG_IMEX_EX
{
public:
  Kinetic1D_LD_DG_IMEX_EX_IsentropicBC(const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table);
  ~Kinetic1D_LD_DG_IMEX_EX_IsentropicBC() override = default;

  Matrix rho_init(const Matrix& x) override;
  real_t rho_init(const real_t& x) override;
  Matrix f_init(const Matrix& x, const real_t& v) override;
  real_t f_init(const real_t& x, const real_t& v) override;
  Matrix g_init(const Matrix& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;

  real_t fL_bc(const int& j, const real_t& t,
                      const model_data_& modal) override;
  real_t fR_bc(const int& j, const real_t& t,
                      const model_data_& modal) override;

  void init() override;
};


class Kinetic1D_LD_DG_IMEX_EX_twovel_period : 
  public virtual Kinetic1D_LD_DG_IMEX_EX
{
public:
  Kinetic1D_LD_DG_IMEX_EX_twovel_period(const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table);
  ~Kinetic1D_LD_DG_IMEX_EX_twovel_period() override = default;

  Matrix rho_init(const Matrix& x) override;
  real_t rho_init(const real_t& x) override;
  Matrix f_init(const Matrix& x, const real_t& v) override;
  real_t f_init(const real_t& x, const real_t& v) override;
  Matrix g_init(const Matrix& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;

  virtual Matrix rho_real(const Matrix& x, const real_t& t);
  virtual real_t rho_real(const real_t& x, const real_t& t);
  virtual Matrix f_real(const Matrix& x, const real_t& v, const real_t& t);
  virtual real_t f_real(const real_t& x, const real_t& v, const real_t& t);
  virtual Matrix g_real(const Matrix& x, const real_t& v, const real_t& t);
  virtual real_t g_real(const real_t& x, const real_t& v, const real_t& t);

  void init() override;
  void updateAll(const real_t& Trun, const real_t& dt) override;

  void ah_extflux_compute(const model_data_& modal, const real_t& Trun, 
                          Matrix* boundary_flux) override;
  void bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux) override;
  void dh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux) override;

  virtual void getrho_real_modal(const real_t& Tstop, Matrix* rho_real_nodal);
  virtual void getrho_real_nodal(const real_t& Tstop, Matrix* rho_real_nodal);
  virtual void getg_real_nodal(const real_t& Tstop, 
                      std::vector<Matrix>* g_real_nodal);

protected:
  real_t r_;
};

} // namespace QUEST



#endif // QUEST_KIETICLD_DG_IMEX_EX_HPP

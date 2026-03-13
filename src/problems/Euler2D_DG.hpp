#ifndef QUEST_EULER2D_DG_HPP
#define QUEST_EULER2D_DG_HPP

#include "hyperbolicproblemsbase.hpp"

namespace QUEST
{

class Euler2D_DG_TVDRK : public HyperbolicProblems2DBase
{
public:
  Euler2D_DG_TVDRK(const TensorMesh2D* mesh2D,
                   const fespace2D* fe,
                   const EX_TVDRK* rk_table);
  ~Euler2D_DG_TVDRK() override = default;

  virtual Matrix rho_init(const Matrix& x, const Matrix& y);
  virtual Matrix vx_init(const Matrix& x, const Matrix& y);
  virtual Matrix vy_init(const Matrix& x, const Matrix& y);
  virtual Matrix pre_init(const Matrix& x, const Matrix& y);

  virtual real_t rho_bc(const real_t& x, const real_t& y, const real_t& t);
  virtual real_t vx_bc(const real_t& x, const real_t& y, const real_t& t);
  virtual real_t vy_bc(const real_t& x, const real_t& y, const real_t& t);
  virtual real_t pre_bc(const real_t& x, const real_t& y, const real_t& t);

  void init() override;
  void setdt(real_t* dt) override;
  real_t max_speed_compute(const std::vector<Matrix>& u_modal) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void flux_compute(const std::vector<Matrix>& u_modal,
                    std::vector<Matrix>* fu_nodal,
                    std::vector<Matrix>* fv_nodal) override;
  void Lu_compute(const std::vector<Matrix>& u_modal,
                  const real_t& Trun,
                  const real_t& dt,
                  std::vector<Matrix>* Lu) override;
  void fluxint_compute(const std::vector<Matrix>& u_modal,
                       std::vector<Matrix>* flux_int) override;
  void fluxext_compute(const std::vector<Matrix>& u_modal,
                       const std::vector<Matrix>& Dirichlet,
                       std::vector<Matrix>* flux_ext) override;
  void numerical_flux(const Vector& u_L, const Vector& u_R,
                      const Eigen::Vector2d& normal,
                      const real_t& max_speed,
                      Vector* flux) override;

  virtual void setgamma(const real_t& gamma);

protected:
  real_t gamma_;
  real_t cfl_;
  real_t pi_;
};

class Euler2D_DG_TVDRK_period : public Euler2D_DG_TVDRK
{
public:
  Euler2D_DG_TVDRK_period(const TensorMesh2D* mesh2D,
                          const fespace2D* fe,
                          const EX_TVDRK* rk_table);
  ~Euler2D_DG_TVDRK_period() override = default;

  void Lu_compute(const std::vector<Matrix>& u_modal,
                  const real_t& Trun,
                  const real_t& dt,
                  std::vector<Matrix>* Lu) override;
  void fluxext_compute(const std::vector<Matrix>& u_modal,
                       const std::vector<Matrix>& Dirichlet,
                       std::vector<Matrix>* flux_ext) override;

  real_t rho_bc(const real_t& x, const real_t& y, const real_t& t) override;
  real_t vx_bc(const real_t& x, const real_t& y, const real_t& t) override;
  real_t vy_bc(const real_t& x, const real_t& y, const real_t& t) override;
  real_t pre_bc(const real_t& x, const real_t& y, const real_t& t) override;

  virtual Matrix rho_real(const Matrix& x, const Matrix& y, const real_t& t);
  virtual real_t rho_real(const real_t& x, const real_t& y, const real_t& t);
  virtual Matrix vx_real(const Matrix& x, const Matrix& y, const real_t& t);
  virtual real_t vx_real(const real_t& x, const real_t& y, const real_t& t);
  virtual Matrix vy_real(const Matrix& x, const Matrix& y, const real_t& t);
  virtual real_t vy_real(const real_t& x, const real_t& y, const real_t& t);
  virtual Matrix pre_real(const Matrix& x, const Matrix& y, const real_t& t);
  virtual real_t pre_real(const real_t& x, const real_t& y, const real_t& t);
};

} // namespace QUEST

#endif // QUEST_EULER2D_DG_HPP

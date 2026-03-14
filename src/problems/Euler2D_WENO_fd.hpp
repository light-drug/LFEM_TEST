#ifndef QUEST_EULER2D_WENO_FD_HPP
#define QUEST_EULER2D_WENO_FD_HPP

#include <vector>

#include "FDmesh.hpp"
#include "RK_table.hpp"
#include "config.hpp"

namespace QUEST
{

class Euler2D_WENO_FD
{
public:
  Euler2D_WENO_FD(const FDmesh2D* mesh2D,
                  const EX_TVDRK* rk_table,
                  const int x_order);
  virtual ~Euler2D_WENO_FD() = default;

  virtual void setgamma(const real_t& gamma);
  virtual void setCFL(const real_t& cfl);

  virtual void init();
  virtual void setdt(real_t* dt) const;
  virtual void updateAll(const real_t& Trun, const real_t& dt);
  virtual void Lu_compute(const std::vector<Matrix>& u,
                          const real_t& Trun,
                          std::vector<Matrix>* Lu);

  virtual real_t rho_init(const real_t& x, const real_t& y) const;
  virtual real_t vx_init(const real_t& x, const real_t& y) const;
  virtual real_t vy_init(const real_t& x, const real_t& y) const;
  virtual real_t pre_init(const real_t& x, const real_t& y) const;
  virtual Matrix rho_init(const Matrix& x, const Matrix& y) const;
  virtual Matrix vx_init(const Matrix& x, const Matrix& y) const;
  virtual Matrix vy_init(const Matrix& x, const Matrix& y) const;
  virtual Matrix pre_init(const Matrix& x, const Matrix& y) const;

  virtual real_t rho_bc(const real_t& x, const real_t& y, const real_t& t) const;
  virtual real_t vx_bc(const real_t& x, const real_t& y, const real_t& t) const;
  virtual real_t vy_bc(const real_t& x, const real_t& y, const real_t& t) const;
  virtual real_t pre_bc(const real_t& x, const real_t& y, const real_t& t) const;

  virtual Matrix rho_real(const Matrix& x, const Matrix& y, const real_t& t) const;

  const std::vector<Matrix>& getu() const;
  const Matrix& getrho() const;

protected:
  const FDmesh2D* mesh2D_;
  const EX_TVDRK* rk_table_;
  int x_order_;
  real_t gamma_;
  real_t CFL_;
  int num_equations_;

  // u_[k] are all (xDiv, yDiv), values at cell centers.
  std::vector<Matrix> u_;
  std::vector<std::vector<Matrix>> u_stages_;
  std::vector<std::vector<Matrix>> Lu_stages_;

  virtual real_t max_wave_speed(const std::vector<Matrix>& u) const;
  virtual Vector weno3_left_biased(const std::vector<Vector>& unei) const;
  virtual Vector weno5_left_biased(const std::vector<Vector>& unei) const;
  virtual Vector weno3_right_biased(const std::vector<Vector>& unei) const;
  virtual Vector weno5_right_biased(const std::vector<Vector>& unei) const;
  virtual Vector flux_x_eval(const Vector& u) const;
  virtual Vector flux_y_eval(const Vector& u) const;

  virtual Vector get_state(const std::vector<Matrix>& u,
                           const int i,
                           const int j,
                           const real_t& Trun,
                           const bool periodic) const;
};

class Euler2D_WENO_FD_period : public Euler2D_WENO_FD
{
public:
  Euler2D_WENO_FD_period(const FDmesh2D* mesh2D,
                         const EX_TVDRK* rk_table,
                         const int x_order);
  ~Euler2D_WENO_FD_period() override = default;

  real_t rho_init(const real_t& x, const real_t& y) const override;
  real_t rho_bc(const real_t& x, const real_t& y, const real_t& t) const override;
  Matrix rho_real(const Matrix& x, const Matrix& y, const real_t& t) const override;
  void Lu_compute(const std::vector<Matrix>& u,
                  const real_t& Trun,
                  std::vector<Matrix>* Lu) override;
};

} // namespace QUEST

#endif // QUEST_EULER2D_WENO_FD_HPP

#ifndef QUEST_EULER1D_WENO_FD_HPP
#define QUEST_EULER1D_WENO_FD_HPP

#include <vector>

#include "FDmesh.hpp"
#include "RK_table.hpp"
#include "config.hpp"

namespace QUEST
{

class Euler1D_WENO_FD
{
protected:
  const FDmesh* mesh1D_;
  const EX_TVDRK* rk_table_;
  real_t a_;
  int x_order_;

  real_t pi_;
  real_t CFL_;

  std::vector<Vector> u_; // rho, rhovx, E 守恒变量
  std::vector<std::vector<Vector>> u_stages_;
  std::vector<std::vector<Vector>> Lu_stages_;

  virtual std::vector<Vector> extend_with_ghost(const std::vector<Vector>& u, const real_t& Trun) const;
  virtual void extend_left_ghost(const int ghost_id, const real_t& Trun,
                                std::vector<Vector>* u);
  virtual void extend_right_ghost(const int ghost_id, const real_t& Trun,
                                  std::vector<Vector>* u) const;

  virtual real_t weno3_left_biased(const Vector& ue, const int iface) const;
  virtual real_t weno5_left_biased(const Vector& ue, const int iface) const;

public:
  Euler1D_WENO_FD(const FDmesh* mesh1D,
                  const EX_TVDRK* rk_table,
                  const real_t a,
                  const int x_order);
  virtual ~Euler1D_WENO_FD() = default;

  virtual void setCFL(const real_t& cfl);
  virtual const std::vector<Vector>& getu() const;
  virtual const Vector& getrho() const;
  virtual const Vector& getrhovx() const;
  virtual const Vector& getE() const;
  virtual Vector getpre(const std::vector<Vector>& u) const;

  virtual void init();
  virtual void setdt(real_t* dt) const;
  virtual void Lu_compute(const std::vector<Vector>& u, const real_t& Trun, 
                          std::vector<Vector>* Lu);
  virtual void updateAll(const real_t& Trun, const real_t& dt);

  virtual real_t rho_init(const real_t& x) const;
  virtual Vector rho_init(const Vector& x) const;
  virtual real_t rho_bc(const real_t& x, const real_t& t) const;
  virtual Vector rho_bc(const Vector& x, const real_t& t) const;

  virtual real_t vx_init(const real_t& x) const;
  virtual Vector vx_init(const Vector& x) const;
  virtual real_t vx_bc(const real_t& x, const real_t& t) const;
  virtual Vector vx_bc(const Vector& x, const real_t& t) const;

  virtual real_t pre_init(const real_t& x) const;
  virtual Vector pre_init(const Vector& x) const;
  virtual real_t pre_bc(const real_t& x, const real_t& t) const;
  virtual Vector pre_bc(const Vector& x, const real_t& t) const;
};

class Euler1D_WENO_FD_period : public Euler1D_WENO_FD
{
protected:
  const FDmesh_period* mesh1D_period_;

  std::vector<Vector> extend_with_ghost(const std::vector<Vector>& u, 
                                        const real_t& Trun) const override;
  real_t left_ghost_value(const int ghost_id, const real_t& Trun,
                          std::vector<Vector>* u) const override;
  real_t right_ghost_value(const int ghost_id, const real_t& Trun,
                            std::vector<Vector>* u) const override;

public:
  Euler1D_WENO_FD_period(const FDmesh_period* mesh1D,
                          const EX_TVDRK* rk_table,
                          const real_t a,
                          const int x_order);
  ~Euler1D_WENO_FD_period() override = default;

  real_t rho_init(const real_t& x) const override;
  Vector rho_init(const Vector& x) const override;
  real_t rho_bc(const real_t& x, const real_t& t) const override;
  Vector rho_bc(const Vector& x, const real_t& t) const override;

  real_t vx_init(const real_t& x) const override;
  Vector vx_init(const Vector& x) const override;
  real_t vx_bc(const real_t& x, const real_t& t) const override;
  Vector vx_bc(const Vector& x, const real_t& t) const override;

  real_t pre_init(const real_t& x) const override;
  Vector pre_init(const Vector& x) const override;
  real_t pre_bc(const real_t& x, const real_t& t) const override;
  Vector pre_bc(const Vector& x, const real_t& t) const override;
};

} // namespace QUEST

#endif // QUEST_EULER1D_WENO_FD_HPP

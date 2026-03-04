#ifndef QUEST_LINEARADEVECTION_WENO_FD_HPP
#define QUEST_LINEARADEVECTION_WENO_FD_HPP

#include <vector>

#include "FDmesh.hpp"
#include "RK_table.hpp"
#include "config.hpp"

namespace QUEST
{

class LinearAdvection_WENO_FD
{
protected:
  const FDmesh* mesh1D_;
  const EX_TVDRK* rk_table_;
  real_t a_;
  int x_order_;

  real_t pi_;
  real_t CFL_;

  Vector u_;
  std::vector<Vector> u_stages_;
  std::vector<Vector> Lu_stages_;

  virtual Vector extend_with_ghost(const Vector& u, const real_t& Trun) const;
  virtual real_t left_ghost_value(const int ghost_id, const real_t& Trun) const;
  virtual real_t right_ghost_value(const Vector& u, const int ghost_id, const real_t& Trun) const;

  virtual real_t weno3_left_biased(const Vector& ue, const int iface) const;
  virtual real_t weno5_left_biased(const Vector& ue, const int iface) const;

public:
  LinearAdvection_WENO_FD(const FDmesh* mesh1D,
                          const EX_TVDRK* rk_table,
                          const real_t a,
                          const int x_order);
  virtual ~LinearAdvection_WENO_FD() = default;

  virtual void setCFL(const real_t& cfl);
  virtual const Vector& getu() const;

  virtual void init();
  virtual void setdt(real_t* dt) const;
  virtual void Lu_compute(const Vector& u, const real_t& Trun, Vector* Lu);
  virtual void updateAll(const real_t& Trun, const real_t& dt);

  virtual real_t u_init(const real_t& x) const;
  virtual Vector u_init(const Vector& x) const;
  virtual real_t u_bc(const real_t& x, const real_t& t) const;
  virtual Vector u_bc(const Vector& x, const real_t& t) const;

  virtual real_t u_real(const real_t& x, const real_t& t) const;
  virtual Vector u_real(const Vector& x, const real_t& t) const;
};

class LinearAdvection_WENO_FD_period : public LinearAdvection_WENO_FD
{
protected:
  const FDmesh_period* mesh1D_period_;

  Vector extend_with_ghost(const Vector& u, const real_t& Trun) const override;
  real_t left_ghost_value(const int ghost_id, const real_t& Trun) const override;
  real_t right_ghost_value(const Vector& u, const int ghost_id, const real_t& Trun) const override;

public:
  LinearAdvection_WENO_FD_period(const FDmesh_period* mesh1D,
                                 const EX_TVDRK* rk_table,
                                 const real_t a,
                                 const int x_order);
  ~LinearAdvection_WENO_FD_period() override = default;

  real_t u_init(const real_t& x) const override;
  Vector u_init(const Vector& x) const override;
  real_t u_real(const real_t& x, const real_t& t) const override;
  Vector u_real(const Vector& x, const real_t& t) const override;
};

} // namespace QUEST

#endif // QUEST_LINEARADEVECTION_WENO_FD_HPP

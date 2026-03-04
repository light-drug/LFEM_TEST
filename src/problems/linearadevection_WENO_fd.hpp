#ifndef QUEST_LINEARADEVECTION_WENO_FD_HPP
#define QUEST_LINEARADEVECTION_WENO_FD_HPP

#include <iostream>
#include <vector>
#include "config.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "RK_table.hpp"
#include "FDmesh.hpp"

namespace QUEST
{

class LinearAdvection_WENO_FD
{
protected:

  const FDmesh* mesh1D_;
  const EX_TVDRK* rk_table_;
  const real_t a_;
  
  real_t pi_;
  real_t CFL_;
  
  Vector u_;
  std::vector<Vector> u_stages_;
  Vector Lu_;

public: 
  LinearAdvection_WENO_FD(const FDmesh* mesh1D,
                          const EX_TVDRK* rk_table,
                          const real_t a);
  virtual ~LinearAdvection_WENO_FD() = default;

  virtual void setCFL(const real_t& cfl);

  virtual void init();
  virtual void setdt(real_t* dt);
  virtual void Lu_compute(const Vector& u, const real_t& Trun, Vector* Lu);
  virtual void updateAll(const real_t& Trun, const real_t& dt);
  virtual real_t u_init(const real_t& x);
  virtual Vector u_init(const Vector& x);
  virtual real_t u_bc(const real_t& x, const real_t& t);
  virtual real_t u_bc(const Vector& x, const real_t& t);

};

class LinearAdvection_WENO_FD_period
 : public virtual LinearAdvection_WENO_FD
{
protected:
  const FDmesh_period* mesh1D_;
  const EX_TVDRK* rk_table_;
  const real_t a_;

public: 
  LinearAdvection_WENO_FD(const FDmesh_period* mesh1D,
                          const EX_TVDRK* rk_table,
                          const real_t a);
  virtual ~LinearAdvection_WENO_FD() = default;

  void init() override;
  void setdt(real_t* dt) override;
  void Lu_compute(const Vector& u, const real_t& Trun, Vector* Lu) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  real_t u_init(const real_t& x) override;
  Vector u_init(const Vector& x) override;
  real_t u_real(const real_t& x, const real_t& t) override;
  Vector u_real(const Vector& x, const real_t& t) override;

};

};

#endif // QUEST_LINEARADEVECTION_WENO_FD_HPP

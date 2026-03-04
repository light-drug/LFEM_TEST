#ifndef QUEST_DRIFTDIFFUSION1D_SL_LOMAC_HPP
#define QUEST_DRIFTDIFFUSION1D_SL_LOMAC_HPP

#include "config.hpp"
#include "FDmesh.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "Poissonsolver1DFD.hpp"
#include "KineticDriftDiffusion1D_SL.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

namespace QUEST
{
  
enum class TimeSteppingType {
  Exponential,
  Rktype,
};

class KineticDriftDiffusion1D_SL_lomac
  : public virtual KineticDriftDiffusion1D_SL
{
public:
  KineticDriftDiffusion1D_SL_lomac(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticDriftDiffusion1D_SL_lomac() override = default;

  virtual void setTimeSteppingType(const TimeSteppingType& time_stepping);
  virtual void setsigmas(const real_t& sigmas);
  
  void init(const Solver1DTypeFD& soltype) override;
  void character_tracing(const real_t& Trun, const real_t& dt) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf) override;
  
  void SolveS(const real_t& Trun, const Vector& bf) override;
  void update_E(const real_t& Trun) override;
  void update_f(const real_t& Trun, const real_t& dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) override;
  
  const Vector& getphi();

  // virtual void character_tracing_half(const real_t& Trun, const real_t& dt);
protected:
  TimeSteppingType time_stepping_;

  // std::vector<Vector> xback_half_;
  // std::vector<Vector> vback_half_;
  real_t phi_L_;
  real_t phi_R_;
  real_t sigmas_;
};

class KineticDriftDiffusion1D_SL_lomac_period
  : public virtual KineticDriftDiffusion1D_SL_period,
    public virtual KineticDriftDiffusion1D_SL_lomac
{
public:
  KineticDriftDiffusion1D_SL_lomac_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticDriftDiffusion1D_SL_lomac_period() override = default;

  void init(const Solver1DTypeFD& soltype) override;
  void character_tracing(const real_t& Trun, const real_t& dt) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                  const real_t& dt,
                  Vector* bf) override;
  void SolveS(const real_t& Trun, const Vector& bf) override;

  void update_E(const real_t& Trun) override;
  void update_f(const real_t& Trun, const real_t& dt) override;
  void interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) override;

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
  real_t rho_d(const real_t& x) override;
  Vector rho_d(const Vector& x) override;
  real_t source(const real_t& x, const real_t& t) override;
  Vector source(const Vector& x, const real_t& t) override;
  real_t fsource(const real_t& x, const real_t& v, const real_t& t) override;
  Vector fsource(const Vector& x, const real_t& v, const real_t& t) override;
  real_t fsource_dx(const real_t& x, const real_t& t) override;
  Vector fsource_dx(const Vector& x, const real_t& t) override;

  real_t f_init(const real_t& x, const real_t& v) override;
  Vector f_init(const Vector& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;
  Vector g_init(const Vector& x, const real_t& v) override;

  real_t rho_real(const real_t& x, const real_t& t) override;
  Vector rho_real(const Vector& x, const real_t& t) override;
  Vector E_real(const Vector& x, const real_t& t) override;
  real_t g_real(const real_t& x, const real_t& v, const real_t& t) override;
  Vector g_real(const Vector& x, const real_t& v, const real_t& t) override;
};

class KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma: 
  public virtual KineticDriftDiffusion1D_SL_period_Different_Gamma,
  public virtual KineticDriftDiffusion1D_SL_lomac_period
{
public:
  KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma(const KineticFDmesh_period* mesh1D,
                          Poissonsolver1DFD_period* poisol,
                          const int& x_order,
                          const int& t_order);
  ~KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma() override = default;

  void init(const Solver1DTypeFD& soltype) override;
  void character_tracing(const real_t& Trun, const real_t& dt) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                  const real_t& dt,
                  Vector* bf) override;
  void SolveS(const real_t& Trun, const Vector& bf) override;

  void update_E(const real_t& Trun) override;
  void update_f(const real_t& Trun, const real_t& dt) override;
  void interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) override;

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
  real_t rho_d(const real_t& x) override;
  Vector rho_d(const Vector& x) override;
  real_t source(const real_t& x, const real_t& t) override;
  Vector source(const Vector& x, const real_t& t) override;
  real_t fsource(const real_t& x, const real_t& v, const real_t& t) override;
  Vector fsource(const Vector& x, const real_t& v, const real_t& t) override;
  real_t fsource_dx(const real_t& x, const real_t& t) override;
  Vector fsource_dx(const Vector& x, const real_t& t) override;

  real_t f_init(const real_t& x, const real_t& v) override;
  Vector f_init(const Vector& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;
  Vector g_init(const Vector& x, const real_t& v) override;

  real_t rho_real(const real_t& x, const real_t& t) override;
  Vector rho_real(const Vector& x, const real_t& t) override;
  Vector E_real(const Vector& x, const real_t& t) override;
  real_t g_real(const real_t& x, const real_t& v, const real_t& t) override;
  Vector g_real(const Vector& x, const real_t& v, const real_t& t) override;
};

class KineticDriftDiffusion1D_SL_lomac_PNjunction:
  public virtual KineticDriftDiffusion1D_SL_lomac
{
public:
  KineticDriftDiffusion1D_SL_lomac_PNjunction(const KineticFDmesh* mesh1D,
                          Poissonsolver1DFD* poisol,
                          const int& x_order,
                          const int& t_order);
  
  void init(const Solver1DTypeFD& soltype) override;
  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
  real_t rho_numericalbc(const real_t& x, const real_t& t) override;
  real_t phi_bc(const real_t& x, const real_t & t) override;
  real_t rho_d(const real_t& x) override;
  Vector rho_d(const Vector& x) override;
  real_t source(const real_t& x, const real_t& t) override;
  Vector source(const Vector& x, const real_t& t) override;
  real_t fsource(const real_t& x, const real_t& v, const real_t& t) override;
  Vector fsource(const Vector& x, const real_t& v, const real_t& t) override;
  real_t fsource_dx(const real_t& x, const real_t& t) override;
  Vector fsource_dx(const Vector& x, const real_t& t) override;

  real_t f_init(const real_t& x, const real_t& v) override;
  Vector f_init(const Vector& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;
  Vector g_init(const Vector& x, const real_t& v) override;
  real_t fL_bc(const real_t& v, const real_t& t) override;
  real_t fL_bc(const int& v_index, const real_t& t) override;
  real_t fR_bc(const real_t& v, const real_t& t) override;
  real_t fR_bc(const int& v_index, const real_t& t) override;

};

}; // namespace QUEST

#endif // QUEST_DRIFTDIFFUSION1D_SL_LOMAC_HPP

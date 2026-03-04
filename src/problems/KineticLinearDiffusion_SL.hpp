#ifndef QUEST_KINETICLINEARDIFFUSION1D_SL_HPP
#define QUEST_KINETICLINEARDIFFUSION1D_SL_HPP

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

class KineticLinearDiffusion1D_SL :
  public KineticDriftDiffusion1D_SL
{ 
protected:
  real_t sigmas_;
public:
  KineticLinearDiffusion1D_SL(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticLinearDiffusion1D_SL() override = default;

  void init(const Solver1DTypeFD& soltype) override;
  void character_tracing(const real_t& Trun, const real_t& dt) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                  const real_t& dt,
                  Vector* bf) override;
  void SolveS(const real_t& Trun, const Vector& bf) override;

  void update_f(const real_t& Trun, const real_t& dt) override;
  void interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) override;

  virtual void setsigmas(const real_t& sigmas);

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
  real_t rho_numericalbc(const real_t& x, const real_t& t) override;
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

class KineticLinearDiffusion1D_SL_period: 
  public KineticDriftDiffusion1D_SL_period
{
protected:

#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> ldlt_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> ldlt_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;

public:
  KineticLinearDiffusion1D_SL_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticLinearDiffusion1D_SL_period() override = default;

  void character_tracing(const real_t& Trun, const real_t& dt) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                  const real_t& dt,
                  Vector* bf) override;
  void SolveS(const real_t& Trun, const Vector& bf) override;

  void update_f(const real_t& Trun, const real_t& dt) override;
  void interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) override;

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
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
  real_t f_real(const real_t& x, const real_t& v, const real_t& t) override;
  Vector f_real(const Vector& x, const real_t& v, const real_t& t) override;
  real_t g_real(const real_t& x, const real_t& v, const real_t& t) override;
  Vector g_real(const Vector& x, const real_t& v, const real_t& t) override;

};

class KineticLinearDiffusion1D_SL_period_twovel: 
  public KineticLinearDiffusion1D_SL_period
{
protected:
  std::vector<Vector> rho_inter_;
#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> ldlt_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> ldlt_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;

public:
  KineticLinearDiffusion1D_SL_period_twovel(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticLinearDiffusion1D_SL_period_twovel() override = default;

  void init(const Solver1DTypeFD& soltype) override;
  void setdt(real_t* dt) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;

  void update_f(const real_t& Trun, const real_t& dt) override;

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
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
  real_t f_real(const real_t& x, const real_t& v, const real_t& t) override;
  Vector f_real(const Vector& x, const real_t& v, const real_t& t) override;
  real_t g_real(const real_t& x, const real_t& v, const real_t& t) override;
  Vector g_real(const Vector& x, const real_t& v, const real_t& t) override;

};

class KineticLinearDiffusion1D_SL_period_twovel_v2: 
  public KineticLinearDiffusion1D_SL_period_twovel
{
protected:

  Eigen::GMRES<SparseMatrix> gmresv1_;
  Eigen::GMRES<SparseMatrix> gmresv2_;

  SparseMatrix Sv1_;
  SparseMatrix Sv2_;
  Vector bfv1_;
  Vector bfv2_;

public:
  KineticLinearDiffusion1D_SL_period_twovel_v2(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order);

  ~KineticLinearDiffusion1D_SL_period_twovel_v2() override = default;

  void update_f(const real_t& Trun, const real_t& dt) override;
  // virtual void generateSv();
};

} // namespace QUEST

#endif  // QUEST_KINETICLINEARDIFFUSION1D_SL_HPP

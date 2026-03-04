#ifndef QUEST_KINETICDRIFTDIFFUSION1D_SL_HPP
#define QUEST_KINETICDRIFTDIFFUSION1D_SL_HPP

#include "config.hpp"
#include "FDmesh.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "Poissonsolver1DFD.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

namespace QUEST
{

enum class Solver1DTypeFD {
  CG,
  PCG,
  LDLT,
  LU,
  GMRES,
  BICG,
};

class KineticDriftDiffusion1D_SL
{
protected:

  const KineticFDmesh* mesh1D_;
  Poissonsolver1DFD* poisol_;
  const int& x_order_;
  const int& t_order_;

  std::vector<Vector> f_;
  std::vector<Vector> finter_;
  std::vector<Vector> rhointer_;
  std::vector<Vector> xback_;
  std::vector<Vector> vback_;
  std::vector<Vector> Fx_;
  Vector SLflux_;
  SparseMatrix S_;
  Solver1DTypeFD soltype_;
  Vector Source_;
  Vector Source_v_;
  std::vector<Vector> Source_f_;

#ifndef QUEST_USE_MKL
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif 
  Eigen::BiCGSTAB<SparseMatrix> bicg_;
  Eigen::GMRES<SparseMatrix> gmres_;

  Vector rhon_;
  Vector rho_;
  Vector rhod_;
  Vector phi_;
  Vector E_;

  real_t theta_;
  real_t cfl_;
  real_t pi_;
  real_t gamma_;
  real_t eps_;
  real_t sparsetol_;
  real_t iterationtol_;
  real_t D_;
  int NTH_;

  real_t rho_L_;
  real_t rho_R_;

  int iter_;

public:
  KineticDriftDiffusion1D_SL(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  virtual ~KineticDriftDiffusion1D_SL() = default;

  virtual void settheta(const real_t& theta);
  virtual void setcfl(const real_t& cfl);
  virtual void setpi(const real_t& pi);
  virtual void setgamma(const real_t& gamma);
  virtual void seteps(const real_t& eps);
  virtual void setsparsetol(const real_t& sparsetol);
  virtual void setiterationtol(const real_t& iterationtol);
  virtual void setD(const real_t& D);
  virtual void setNTH(const int& NTH);

  virtual const real_t& getcfl() const;
  virtual const int& getNTH() const;
  
  virtual void init(const Solver1DTypeFD& soltype);
  virtual void character_tracing(const real_t& Trun, const real_t& dt);
  virtual void setdt(real_t* dt);
  virtual void update_rho(const real_t& Trun, const real_t& dt);
  virtual void generateS(const real_t& Trun, const real_t& dt);
  virtual void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf);
  virtual void SolveS(const real_t& Trun, const Vector& bf);
  virtual void update_E(const real_t& Trun);
  virtual void update_f(const real_t& Trun, const real_t& dt);
  virtual void updateAll(const real_t& Trun, const real_t& dt);
  virtual void interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp);
  
  virtual real_t rho_init(const real_t& x);
  virtual Vector rho_init(const Vector& x);
  virtual real_t rho_numericalbc(const real_t& x, const real_t& t);
  virtual real_t phi_bc(const real_t& x, const real_t & t);
  virtual real_t rho_d(const real_t& x);
  virtual Vector rho_d(const Vector& x);
  virtual real_t source(const real_t& x, const real_t& t);
  virtual Vector source(const Vector& x, const real_t& t);
  virtual real_t fsource(const real_t& x, const real_t& v, const real_t& t);
  virtual Vector fsource(const Vector& x, const real_t& v, const real_t& t);
  virtual real_t fsource_dx(const real_t& x, const real_t& t);
  virtual Vector fsource_dx(const Vector& x, const real_t& t);

  virtual real_t f_init(const real_t& x, const real_t& v);
  virtual Vector f_init(const Vector& x, const real_t& v);
  virtual real_t g_init(const real_t& x, const real_t& v);
  virtual Vector g_init(const Vector& x, const real_t& v);
  virtual real_t fL_bc(const real_t& v, const real_t& t);
  virtual real_t fL_bc(const int& v_index, const real_t& t);
  virtual real_t fR_bc(const real_t& v, const real_t& t);
  virtual real_t fR_bc(const int& v_index, const real_t& t);
  virtual real_t Maxwell(const real_t& v);

  virtual const Vector& getrho() const;
  virtual const std::vector<Vector>& getf() const;
  virtual const Vector& getE() const;
  virtual const Vector getrhod() const;
  virtual const int& getiter() const;
  virtual Vector getrhoinit();
  virtual std::vector<Vector> getfinit();

};

class KineticDriftDiffusion1D_SL_period: 
  public virtual KineticDriftDiffusion1D_SL
{
protected:
  Vector rho_final_;
  Vector f_final_;
public:
  KineticDriftDiffusion1D_SL_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticDriftDiffusion1D_SL_period() override = default;
  
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
  
  virtual real_t rho_real(const real_t& x, const real_t& t);
  virtual Vector rho_real(const Vector& x, const real_t& t);
  virtual Vector E_real(const Vector& x, const real_t& t);
  virtual void setrhofinal(const real_t& t);
  virtual const Vector& getrhofinal() const;
  virtual Vector getEfinal(const real_t& t);
  virtual std::vector<Vector> getffinal(const real_t& t);
  virtual real_t f_real(const real_t& x, const real_t& v, const real_t& t);
  virtual Vector f_real(const Vector& x, const real_t& v, const real_t& t);
  virtual real_t g_real(const real_t& x, const real_t& v, const real_t& t);
  virtual Vector g_real(const Vector& x, const real_t& v, const real_t& t);

};


// 对不同gamma算例的精度测试, 算法只对\varepsilon AP，并不对\gamma AP
class KineticDriftDiffusion1D_SL_period_Different_Gamma: 
  public virtual KineticDriftDiffusion1D_SL_period
{
public:

  KineticDriftDiffusion1D_SL_period_Different_Gamma(const KineticFDmesh_period* mesh1D,
                          Poissonsolver1DFD_period* poisol,
                          const int& x_order,
                          const int& t_order);
  ~KineticDriftDiffusion1D_SL_period_Different_Gamma() override = default;

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

class KineticDriftDiffusion1D_SL_PNjunction
  : public KineticDriftDiffusion1D_SL
{
public: 
  KineticDriftDiffusion1D_SL_PNjunction(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  virtual ~KineticDriftDiffusion1D_SL_PNjunction() = default;

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


class KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative:
  public KineticDriftDiffusion1D_SL_period_Different_Gamma
{
public:
  KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative(
                          const KineticFDmesh_period* mesh1D,
                          Poissonsolver1DFD_period* poisol,
                          const int& x_order,
                          const int& t_order);
  ~KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative() override = default;

  void update_rho(const real_t& Trun, const real_t& dt) override;
  void update_f(const real_t& Trun, const real_t& dt) override;
};

} // namespace QUEST

#endif // QUEST_KINETICDRIFTDIFFUSION1D_SL_HPP 

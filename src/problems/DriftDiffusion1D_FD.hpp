#ifndef QUEST_DRIFTDIFFUSION1D_HPP
#define QUEST_DRIFTDIFFUSION1D_HPP

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

/* ************************************
  \rho_t + (\phi_x * rho)_x = rho_xx
  finite difference method for drift-diffusion equation
  implicit treatment for diffusion term, implicit treatment for drift term
****************************************** */

class DriftDiffusion1D_FD 
{
protected:
  const FDmesh* mesh1D_;
  Poissonsolver1DFD_Base* poisol_;
  const int& x_order_;
  const int& t_order_;

  SparseMatrix S_;
  Solver1DTypeFD soltype_;
  Vector Source_;

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

  real_t cfl_;
  real_t pi_;
  real_t gamma_;
  real_t sparsetol_;
  real_t iterationtol_;

  real_t rho_L_;
  real_t rho_R_;

  int iter_;

public:
  DriftDiffusion1D_FD(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order);
  virtual ~DriftDiffusion1D_FD() = default;

  virtual void setcfl(const real_t& cfl);
  virtual void setpi(const real_t& pi);
  virtual void setgamma(const real_t& gamma);
  virtual void setsparsetol(const real_t& sparsetol);
  virtual void setiterationtol(const real_t& iterationtol);

  virtual const real_t& getcfl() const;
  virtual const real_t& getgamma() const;

  virtual void init(const Solver1DTypeFD& soltype);
  virtual void setdt(real_t* dt);
  virtual void update_rho(const real_t& Trun, const real_t& dt);
  virtual void generateS(const real_t& Trun, const real_t& dt);
  virtual void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf);
  virtual void SolveS(const real_t& Trun, const Vector& bf);
  virtual void update_E(const real_t& Trun, const real_t& dt);
  virtual void updateAll(const real_t& Trun, const real_t& dt);

  virtual real_t rho_init(const real_t& x);
  virtual Vector rho_init(const Vector& x);
  virtual real_t phi_bc(const real_t& x, const real_t & t);
  virtual real_t rho_bc(const real_t& x, const real_t & t);
  virtual real_t rho_d(const real_t& x);
  virtual Vector rho_d(const Vector& x);
  virtual real_t source(const real_t& x, const real_t& t);
  virtual Vector source(const Vector& x, const real_t& t);

  virtual const Vector& getrho() const;
  virtual const Vector& getE() const;
  virtual const Vector& getphi() const;
  virtual const Vector& getrhod() const;
  virtual Vector getrhoinit();

};

class DriftDiffusion1D_FD_PNjunction 
  : public virtual DriftDiffusion1D_FD 
{
public:
  DriftDiffusion1D_FD_PNjunction(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order);
  ~DriftDiffusion1D_FD_PNjunction() override = default;

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
  real_t phi_bc(const real_t& x, const real_t & t) override;
  real_t rho_bc(const real_t& x, const real_t & t) override;
  real_t rho_d(const real_t& x) override;
  Vector rho_d(const Vector& x) override;
  real_t source(const real_t& x, const real_t& t) override;
  Vector source(const Vector& x, const real_t& t) override;
};

// 这个是迭代法中+罚项法
// (A^{(l-1)} + \mu/\gamma) * \rho^{(l)} = \rho^n + \mu/\gamma * \rho^{(l-1)}
class DriftDiffusion1D_FD_penaltyIter 
  : public virtual DriftDiffusion1D_FD 
{
protected:
  real_t mu_;
public:
  DriftDiffusion1D_FD_penaltyIter(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order);
  ~DriftDiffusion1D_FD_penaltyIter() override = default;

  // void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf) override;

  virtual void setmu(const real_t& mu);
  virtual const real_t& getmu() const;

};

// 待定 还没弄清楚
// 这个是隐式格式+罚项法
// \rho_t + (\phi_x * (\rho - \rho_d))_x + (\rho_d\phi_{x})_x = rho_xx
// \gamma * \phi_{xx}  = \rho - \rho_d
// S * \rho^{n+1} + B * \phi^{n+1} = bf
// \rho^{n+1}     + Poi *  \phi^{n+1} = \rho_d
class DriftDiffusion1D_FD_penaltyImplicit 
  : public virtual DriftDiffusion1D_FD 
{
protected:
  SparseMatrix B_;
public:
  DriftDiffusion1D_FD_penaltyImplicit(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order);
  ~DriftDiffusion1D_FD_penaltyImplicit() override = default;

  void init(const Solver1DTypeFD& soltype) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf) override;
  void update_E(const real_t& Trun, const real_t& dt) override;
  virtual void generateB(const real_t& Trun, const real_t& dt);
};

class DriftDiffusion1D_FD_PNjunction_penaltyImplicit 
  : public DriftDiffusion1D_FD_penaltyImplicit,
    public DriftDiffusion1D_FD_PNjunction
{
public:
  DriftDiffusion1D_FD_PNjunction_penaltyImplicit(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order);
  ~DriftDiffusion1D_FD_PNjunction_penaltyImplicit() override = default;

  using DriftDiffusion1D_FD_penaltyImplicit::init;
  using DriftDiffusion1D_FD_penaltyImplicit::update_rho;
  using DriftDiffusion1D_FD_penaltyImplicit::generateS;
  using DriftDiffusion1D_FD_penaltyImplicit::generatebf;
  using DriftDiffusion1D_FD_penaltyImplicit::update_E;

  using DriftDiffusion1D_FD_PNjunction::rho_init;
  using DriftDiffusion1D_FD_PNjunction::rho_d;
  using DriftDiffusion1D_FD_PNjunction::phi_bc;
  using DriftDiffusion1D_FD_PNjunction::rho_bc;
  using DriftDiffusion1D_FD_PNjunction::source;
};

// class DriftDiffusion1D_FD_APQuasineutral : public DriftDiffusion1D_FD 
// {
// public:
//   DriftDiffusion1D_FD_APQuasineutral(const FDmesh* mesh1D,
//                         Poissonsolver1DFD_Base* poisol,
//                         const int& x_order,
//                         const int& t_order);
//   ~DriftDiffusion1D_FD_APQuasineutral() override = default;

// };

class DriftDiffusion1D_FD_period : public DriftDiffusion1D_FD 
{
protected:
  Vector rho_final_;

public:
  DriftDiffusion1D_FD_period(const FDmesh_period* mesh1D,
                        Poissonsolver1DFD_Base_period* poisol,
                        const int& x_order,
                        const int& t_order);
  ~DriftDiffusion1D_FD_period() override = default;
  
  // void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                  const real_t& dt,
                  Vector* bf) override;
  // void SolveS(const real_t& Trun, const Vector& bf) override;
  void update_E(const real_t& Trun, const real_t& dt) override;

  real_t rho_init(const real_t& x) override;
  Vector rho_init(const Vector& x) override;
  real_t rho_d(const real_t& x) override;
  Vector rho_d(const Vector& x) override;
  real_t source(const real_t& x, const real_t& t) override;
  Vector source(const Vector& x, const real_t& t) override;

  virtual real_t rho_real(const real_t& x, const real_t& t);
  virtual Vector rho_real(const Vector& x, const real_t& t);
  virtual Vector E_real(const Vector& x, const real_t& t);
  virtual void setrhofinal(const real_t& t);
  virtual const Vector& getrhofinal() const;
  virtual Vector getEfinal(const real_t& t);
  virtual Vector getphifinal(const real_t& t);

};




} // namespace QUEST

#endif // QUEST_DRIFTDIFFUSION1D_HPP

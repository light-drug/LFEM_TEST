#ifndef QUEST_APDRIFTDIFFUSION1D_SL_HPP
#define QUEST_APDRIFTDIFFUSION1D_SL_HPP

#include "config.hpp"
#include "FDmesh.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "Poissonsolver1DFD.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "KineticDriftDiffusion1D_SL.hpp"

namespace QUEST
{

// ************************ // 
// This class preserves both quasi-neutral limit and diffusion limit.
// ************************ //
class APDriftDiffusion_SL : public KineticDriftDiffusion1D_SL
{
protected:
  real_t quasipoi_tol_;
  SparseMatrix quasipoi_;

  Solver1DTypeFD soltype_poi_;
#ifndef QUEST_USE_MKL
  Eigen::SparseLU<SparseMatrix> lupoi_;
#else 
  Eigen::PardisoLU<SparseMatrix> lupoi_;
#endif 
  Eigen::GMRES<SparseMatrix> gmrespoi_;
public:
  
  APDriftDiffusion_SL(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  ~APDriftDiffusion_SL() override = default;
  
  void setquasipoi_tol(const real_t& quasipoi_tol);
  // 这个地方将飘移项的离散方式重写, 此时S的组装方式发生改变
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;

  // rho更新矩阵右端项发生改变
  void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf) override;
  
  // 电场更新方程变为重写后的泊松方程
  virtual void update_E(const real_t& Trun, const real_t& dt);

  // 确定求解quasi-neutral的泊松方程求解器
  virtual void setsoltype_poi(const Solver1DTypeFD& soltype_poi);

  // 由于电场更新不再仅仅需要 update_E
  // virtual void update_E_final(const real_t& Trun, const real_t& dt);
};

class APDriftDiffusion1D_SL_PNjunction : public APDriftDiffusion_SL
{
public:
  APDriftDiffusion1D_SL_PNjunction(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  ~APDriftDiffusion1D_SL_PNjunction() override = default;

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
  
  real_t fL_bc(const real_t& v, const real_t& t) override;
  real_t fL_bc(const int& v_index, const real_t& t) override;
  real_t fR_bc(const real_t& v, const real_t& t) override;
  real_t fR_bc(const int& v_index, const real_t& t) override;

  real_t phi_bc(const real_t& x, const real_t & t) override;
};


class APDriftDiffusion1D_SL_period : public KineticDriftDiffusion1D_SL_period
{
protected:
  real_t quasipoi_tol_;
  
  SparseMatrix quasipoi_;
  Eigen::GMRES<SparseMatrix> gmrespoi_;

public:
  APDriftDiffusion1D_SL_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order);
  ~APDriftDiffusion1D_SL_period() override = default;
  
  void setquasipoi_tol(const real_t& quasipoi_tol);
  // 这个地方将飘移项的离散方式重写, 此时S的组装方式发生改变
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;

  // rho更新矩阵右端项发生改变
  void generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf) override;
  
  // 电场更新方程变为重写后的泊松方程
  virtual void update_E(const real_t& Trun, const real_t& dt);

  // 由于电场更新不再仅仅需要 update_E
  virtual void update_E_final(const real_t& Tstop);

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

}; // namespace QUEST

#endif // QUEST_APDRIFTDIFFUSION1D_SL_HPP

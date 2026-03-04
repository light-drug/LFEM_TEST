#ifndef QUEST_KIETICLD_DG_IMEX_SCHUR_HPP
#define QUEST_KIETICLD_DG_IMEX_SCHUR_HPP

#include <iostream>
#include <vector>

#include "config.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "RK_table.hpp"
#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "KineticLD_DG_IMEX_EX.hpp"

namespace QUEST
{

enum class Solver1DType {
  CG = 0,
  PCG = 1,
  LDLT = 2,
  LU = 3,
  GMRES = 4,
  BICG = 5,
};

class KineticLD_DG_IMEX_IM_Schur : 
  public virtual Kinetic1D_LD_DG_IMEX_EX
{
protected:

  SparseMatrix M_;
  SparseMatrix Minv_;
  SparseMatrix Da_;
  SparseMatrix Db_;
  SparseMatrix Theta_;
  SparseMatrix Thetainv_;
  // std::vector<SparseMatrix> Theta_ref_;
  // std::vector<SparseMatrix> Thetainv_ref_;
  SparseMatrix temp_;
  SparseMatrix Da_Theta_;
  SparseMatrix H_;

  real_t dt_ref_;
  Solver1DType schur_solver_type_;
  real_t schur_tol_;

#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> chol_;
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> chol_;
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;
  Eigen::GMRES<SparseMatrix> gmres_;

public:
  KineticLD_DG_IMEX_IM_Schur(const KineticTensorMesh1D* mesh1D,
                              const fespace1D* fe,
                              const IMEX_RK* rk_table,
                              const Solver1DType& schur_solver_type);
  virtual ~KineticLD_DG_IMEX_IM_Schur() = default;

  virtual void setsparsetol(const real_t& sparsetol);

  virtual void D_compute(const real_t& be, SparseMatrix* D);
  virtual void Da_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Da,
                            std::vector<SparseMatrix>* Da_ext);
  

  // virtual void Db_compute(SparseMatrix* Db);
  virtual void Db_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Db,
                            std::vector<SparseMatrix>* Db_ext);

  virtual void M_compute(SparseMatrix* M);
  virtual void Mbc_compute(SparseMatrix* Mbc);
  virtual void Minv_compute(SparseMatrix* Minv);
  
  virtual void SolveRho(const real_t& Trun,
                        const real_t& dt,
                        const real_t& a,
                        model_data_* modal);
  virtual void SolveRho_DichletBoundary(const real_t& Trun,
                        const real_t& dt,
                        const real_t& a,
                        const Matrix& rhoboundary_flux,
                        const Matrix& vgboundary_flux,
                        model_data_* modal);
  virtual void SolveG(const real_t& Trun,
                      const real_t& dt,
                      const real_t& a,
                      model_data_* modal);
  
  void init() override;
  void setdt(real_t* dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
};

class KineticLD_DG_IMEX_IM_Schur_IsentropicBC : 
  public virtual KineticLD_DG_IMEX_IM_Schur
{
public:
  KineticLD_DG_IMEX_IM_Schur_IsentropicBC(const KineticTensorMesh1D* mesh1D,
                                          const fespace1D* fe,
                                          const IMEX_RK* rk_table,
                                          const Solver1DType& schur_solver_type);
  ~KineticLD_DG_IMEX_IM_Schur_IsentropicBC() override = default;

  Matrix rho_init(const Matrix& x) override;
  real_t rho_init(const real_t& x) override;
  Matrix f_init(const Matrix& x, const real_t& v) override;
  real_t f_init(const real_t& x, const real_t& v) override;
  Matrix g_init(const Matrix& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;

  real_t fL_bc(const int& j, const real_t& t,
                      const model_data_& modal) override;
  real_t fR_bc(const int& j, const real_t& t,
                      const model_data_& modal) override;

  void init() override;
};

class KineticLD_DG_IMEX_IM_Schur_twovel_period : 
  public KineticLD_DG_IMEX_IM_Schur,
  public Kinetic1D_LD_DG_IMEX_EX_twovel_period
{
public:
  KineticLD_DG_IMEX_IM_Schur_twovel_period(const KineticTensorMesh1D* mesh1D,
                              const fespace1D* fe,
                              const IMEX_RK* rk_table,
                              const Solver1DType& schur_solver_type);
  ~KineticLD_DG_IMEX_IM_Schur_twovel_period() override = default;

  void init() override;
  Matrix rho_init(const Matrix& x) override;
  real_t rho_init(const real_t& x) override;
  Matrix f_init(const Matrix& x, const real_t& v) override;
  real_t f_init(const real_t& x, const real_t& v) override;
  Matrix g_init(const Matrix& x, const real_t& v) override;
  real_t g_init(const real_t& x, const real_t& v) override;

  Matrix rho_real(const Matrix& x, const real_t& t) override;
  real_t rho_real(const real_t& x, const real_t& t) override;
  Matrix f_real(const Matrix& x, const real_t& v, const real_t& t) override;
  real_t f_real(const real_t& x, const real_t& v, const real_t& t) override;
  Matrix g_real(const Matrix& x, const real_t& v, const real_t& t) override;
  real_t g_real(const real_t& x, const real_t& v, const real_t& t) override;

  void setdt(real_t* dt) override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  void D_compute(const real_t& be, SparseMatrix* D) override;

  void ah_compute(const model_data_& modal, 
                const real_t& Trun, 
                const Matrix& boundary_flux,
                Matrix* ah) override;
  void bh_compute(const model_data_& modal, 
                const real_t& Trun,
                const std::vector<Matrix>& boundary_flux,
                std::vector<Matrix>* bh) override;
  void dh_compute(const model_data_& modal, 
                const real_t& Trun,
                const Matrix& boundary_flux,
                Matrix* dh) override;
  void sh_compute(const model_data_& modal, 
                const real_t& Trun,
                std::vector<Matrix>* sh) override;

  void ah_extflux_compute(const model_data_& modal, const real_t& Trun, 
                          Matrix* boundary_flux) override;
  void bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux) override;
  void dh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux) override;
};
  
} // namespace QUEST

#endif  // QUEST_KIETICLD_DG_IMEX_SCHUR_HPP 

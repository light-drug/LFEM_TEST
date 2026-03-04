#ifndef QUEST_KIETICLD_DG_IMEX_SCHUR_DIRICHLET_HPP
#define QUEST_KIETICLD_DG_IMEX_SCHUR_DIRICHLET_HPP

#include <iostream>
#include <vector>

#include "config.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "RK_table.hpp"
#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "KineticLD_DG_IMEX_EX.hpp"
#include "KineticLD_DG_IMEX_Schur.hpp"

namespace QUEST
{

class KineticLD_DG_IMEX_IM_Schur_New :
  public virtual KineticLD_DG_IMEX_IM_Schur
{
protected:
  SparseMatrix Mbc_;
  std::vector<SparseMatrix> Da_ext_;
  std::vector<SparseMatrix> Db_ext_;
public:
  KineticLD_DG_IMEX_IM_Schur_New(const KineticTensorMesh1D* mesh1D,
                                const fespace1D* fe,
                                const IMEX_RK* rk_table,
                                const Solver1DType& schur_solver_type);
  ~KineticLD_DG_IMEX_IM_Schur_New() override = default;

  void init() override;
  void updateAll(const real_t& Trun, const real_t& dt) override;
  
  virtual real_t fL_explicit_bc(const int& j, const real_t& t);
  virtual real_t fR_explicit_bc(const int& j, const real_t& t);
  
  virtual void dh_explicit_extflux_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun,
                          Matrix* boundary_flux); 
  void SolveRho_DichletBoundary(const real_t& Trun,
                                const real_t& dt,
                                const real_t& a,
                                const Matrix& rhoboundary_flux,
                                const Matrix& vgboundary_flux,
                                model_data_* modal) override;
  void Da_ext_treat(const real_t& Trun,
                    const real_t& dt, 
                    const SparseMatrix& Da,
                    std::vector<SparseMatrix>* Da_ext) override;
  
  void Db_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Db,
                            std::vector<SparseMatrix>* Db_ext) override;
  
  void Mbc_compute(SparseMatrix* Mbc) override;
                     
};

class KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC : 
  public virtual KineticLD_DG_IMEX_IM_Schur_New
{
public:
  KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC(const KineticTensorMesh1D* mesh1D,
                                          const fespace1D* fe,
                                          const IMEX_RK* rk_table,
                                          const Solver1DType& schur_solver_type);
  ~KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC() override = default;

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
  real_t fL_explicit_bc(const int& j, const real_t& t);
  real_t fR_explicit_bc(const int& j, const real_t& t);

  void init() override;
};

}; // namespace QUEST

#endif // QUEST_KIETICLD_DG_IMEX_SCHUR_DIRICHLET_HPP
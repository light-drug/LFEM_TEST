#ifndef QUEST_KIETICDriftD_DG_IMEX_SCHUR_DIRICHLET_HPP
#define QUEST_KIETICDriftD_DG_IMEX_SCHUR_DIRICHLET_HPP

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
#include "KineticLD_DG_IMEX_Schur_Dirichlet.hpp"
#include "poissonsolver1D.hpp"

namespace QUEST
{

class KineticLinearD_DG_IMEX_IM_Schur_withMaxwell
 : public KineticDriftD_DG_IMEX_IM_Schur
{
protected:

public:
  KineticLinearD_DG_IMEX_IM_Schur_withMaxwell(const KineticTensorMesh1D* mesh1D,
                                const fespace1D* fe,
                                const IMEX_RK* rk_table,
                                const PoissonSolver1D* poisson_solver,
                                const Solver1DType& schur_solver_type);
  virtual ~KineticLinearD_DG_IMEX_IM_Schur_withMaxwell() = default;

  void Eh_compute(const Matrix& E_modal, SparseMatrix* ME) override;

  void SolveRho_DichletBoundary(const real_t& Trun,
                                const real_t& dt,
                                const real_t& a,
                                const Matrix& dh_bc,
                                const Matrix& ah_bc,
                                const Matrix& rho_prev,
                                model_data_* modal) override;

  void ch_compute(const model_data_& modal, 
                  const real_t& Trun,
                  std::vector<Matrix>* ch) override;
  virtual void WENO3_VelocityReconstruct(const std::vector<Matrix>& g,
                                        const std::vector<Matrix>& g_nodal,
                                        std::vector<Matrix>* flux_upwind_nodal,
                                        std::vector<Matrix>* flux_downwind_nodal);                

  virtual void dh_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const Matrix& dh_bc,
                          Matrix* dh);
  virtual void dh_bc_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun,
                          Matrix* dh_bc);

  virtual void mh_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const Matrix& E,
                          Matrix* mh);

  virtual void sh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* sh);


  virtual Matrix rho_init(const Matrix& x);
  virtual real_t rho_init(const real_t& x);
  virtual Matrix f_init(const Matrix& x, const real_t& v);
  virtual real_t f_init(const real_t& x, const real_t& v);
  virtual Matrix g_init(const Matrix& x, const real_t& v);
  virtual real_t g_init(const real_t& x, const real_t& v);
  virtual Matrix source(const Matrix& x, const real_t& t);
  virtual real_t source(const real_t& x, const real_t& t);
  virtual Matrix fsource(const Matrix& x, const real_t& v, const real_t& t);
  virtual real_t fsource(const real_t& x, const real_t& v, const real_t& t);
  virtual Matrix rho_d(const Matrix& x);
  virtual real_t rho_d(const real_t& x);

  virtual real_t fL_bc(const int& j, const real_t& t,
                      const model_data_& modal);
  virtual real_t fR_bc(const int& j, const real_t& t,
                      const model_data_& modal);
  virtual real_t fL_bc(const int& j, const real_t& t,
                      const Matrix& rho, const std::vector<Matrix>& g);
  virtual real_t fR_bc(const int& j, const real_t& t,
                      const Matrix& rho, const std::vector<Matrix>& g);
  virtual real_t fL_explicit_bc(const int& j, const real_t& t);
  virtual real_t fR_explicit_bc(const int& j, const real_t& t);
  virtual real_t gL_bc(const int& j, const real_t& t,
                      const model_data_& modal, const real_t& rho_L);
  virtual real_t gR_bc(const int& j, const real_t& t,
                      const model_data_& modal, const real_t& rho_R);
  virtual real_t phi_bc(const real_t& x, const real_t& t);
  virtual real_t rho_numericalbc(const real_t& x, const real_t& t,
                      const Matrix& rho, const std::vector<Matrix>& g);
  virtual real_t Maxwell(const real_t& v);

  virtual void init();
  virtual void setdt(real_t* dt);
  virtual void updateAll(const real_t& Trun, const real_t& dt);
  virtual void SolveG(const real_t& Trun, const real_t& dt, const real_t& a, model_data_* modal);
};

// class KineticDriftD_DG_IMEX_IM_Schur_New :
//   public virtual KineticDriftD_DG_IMEX_IM_Schur
// {
// public:
//   KineticDriftD_DG_IMEX_IM_Schur_New(const KineticTensorMesh1D* mesh1D,
//                                 const fespace1D* fe,
//                                 const IMEX_RK* rk_table,
//                                 const PoissonSolver1D* poisson_solver,
//                                 const Solver1DType& schur_solver_type);
//   ~KineticDriftD_DG_IMEX_IM_Schur_New() = default;

//   void bh_compute(const model_data_& modal, 
//                   const real_t& Trun,
//                   const std::vector<Matrix>& boundary_flux,
//                   std::vector<Matrix>* bh);

//   void ch_compute(const model_data_& modal, 
//                   const real_t& Trun,
//                   std::vector<Matrix>* ch);
//   void WENO3_VelocityReconstruct(const std::vector<Matrix>& g,
//                                 std::vector<Matrix>* flux_upwind_modal,
//                                 std::vector<Matrix>* flux_downwind_modal);   
// };


}; // namespace QUEST

#endif // QUEST_KIETICLD_DG_IMEX_SCHUR_DIRICHLET_HPP

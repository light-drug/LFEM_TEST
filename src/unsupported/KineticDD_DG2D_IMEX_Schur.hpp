#ifndef QUEST_KIETICDD_DG2D_IMEX_SCHUR_DIRICHLET_HPP
#define QUEST_KIETICDD_DG2D_IMEX_SCHUR_DIRICHLET_HPP

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

class KineticDD_DG2d_IMEX_IM_Schur 
{
protected:
  const TensorMesh1D* mesh1D_;
  const TensorMesh2D* kinetic_mesh2D_;
  const fespace1D* fe1D_;
  const fespace2D* kinetic_fe2D_;
  const IMEX_RK* rk_table_;
  const PoissonSolver1D* poisson_solver_;

  real_t eps_;
  real_t eps2_;
  real_t pi_;
  real_t Chy_;
  real_t Cdif_;
  real_t CR_;
  real_t beta1_; // - 0.e0 for central flux, +-0.5e0 for alternative flux
  real_t sigmas_;
  real_t theta_;
  real_t schur_tol_;
  real_t iter_tol_;
  real_t gamma_; // Debye length
  int NTH_;
  int stages_;
  real_t Maxwell_sum_;

  SparseMatrix M_;
  SparseMatrix Mbc_;
  SparseMatrix Minv_;
  SparseMatrix Da_;
  SparseMatrix Db_;
  SparseMatrix ME_;
  SparseMatrix Theta_;
  SparseMatrix Thetainv_;
  // std::vector<SparseMatrix> Theta_ref_;
  // std::vector<SparseMatrix> Thetainv_ref_;
  SparseMatrix temp_;
  SparseMatrix H_;
  Solver1DType schur_solver_type_;

#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> chol_;
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> chol_;
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;
  Eigen::GMRES<SparseMatrix> gmres_;

  struct model_data_
  {
    Matrix rho;
    Matrix g;
    Matrix phi;
    Matrix E;
  };

  Matrix rho_d_modal_;
  Matrix rho_d_nodal_;
  model_data_ kinetic_modal_;
  std::vector<model_data_> kinetic_modal_stages_;
  std::vector<Matrix> ah_;
  std::vector<std::vector<Matrix>> ch_;
  std::vector<std::vector<Matrix>> bh_;
  std::vector<Matrix> dh_;
  std::vector<Matrix> mh_;
  std::vector<Matrix> sh_;
  std::vector<Matrix> sourceh_;
  std::vector<Matrix> source_sumh_;
  Matrix zero_modal_mat_;
  Matrix zero_nodal_mat_;

  std::vector<std::vector<Matrix>> mixedflux_u_v_;

public:
  KineticDD_DG2d_IMEX_IM_Schur(const TensorMesh1D* mesh1D,
                                const TensorMesh2D* kinetic_mesh2D,
                                const fespace1D* fe1D,
                                const fespace2D* kinetic_fe2D,
                                const IMEX_RK* rk_table,
                                const PoissonSolver1D* poisson_solver,
                                const Solver1DType& schur_solver_type);
  virtual ~KineticDD_DG2d_IMEX_IM_Schur() = default;

  virtual void seteps(const real_t& knu);
  virtual void setCR(const real_t& CR); // penalty term for LDG boundary flux
  virtual void setbeta1(const real_t& beta1);
  virtual void setsigmas(const real_t& sigmas);
  virtual void setgamma(const real_t& gamma);
  virtual void settheta(const real_t& theta);
  virtual void setNTH(const int& NTH);
  virtual void setsparsetol(const real_t& schur_tol);
  virtual void setitertol(const real_t& iter_tol);

  virtual const Matrix& getrho_modal() const;
  virtual const Matrix& getg_modal() const;
  virtual const Matrix& getrho_d_modal() const;

  virtual void D_compute(const real_t& be, SparseMatrix* D);
  virtual void Da_compute(const real_t& Trun,
                            const real_t& dt,
                            SparseMatrix* Da);

  virtual void Db_compute(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Db,
                            SparseMatrix* Db);

  virtual void M_compute(SparseMatrix* M);
  virtual void Mbc_compute(SparseMatrix* Mbc);
  virtual void Minv_compute(SparseMatrix* Minv);
  virtual void Eh_compute(const Matrix& E_modal, SparseMatrix* ME);
  
  // virtual void FourierMatrix(const real_t& xita, std::string& OutfilePath);
  // virtual void FourierDMatrix_compute(const real_t& be, const real_t& xita, cMatrix* D);
  // virtual void FourierUMatrix_compute(const cMatrix& D_NegativeWind,
  //                                     const cMatrix& D_PositiveWind,
  //                                     const real_t& xita, 
  //                                     cMatrix* U);

  virtual void fluxint_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          Matrix* flux_int);
  virtual void fluxext_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          const Matrix& Dirichlet,
                          Matrix* flux_ext);

  virtual void ah_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          const Matrix& ah_bc,
                          Matrix* ah);
  virtual void ah_bc_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun, 
                          Matrix* ah_bc);
  virtual void SolveRho_DichletBoundary(const real_t& Trun,
                                        const real_t& dt,
                                        const real_t& a,
                                        const Matrix& dh_bc,
                                        const Matrix& ah_bc,
                                        const Matrix& rho_prev,
                                        model_data_* modal);

  virtual void bh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          const std::vector<Matrix>& boundary_flux,
                          std::vector<Matrix>* bh);
  virtual void bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux);
  virtual void ch_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* ch);                

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

  virtual void sourceh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          Matrix* source_sumh,
                          std::vector<Matrix>* sourceh);
  virtual void MixedFlux_compute();


  virtual Matrix rho_init(const Matrix& x);
  virtual real_t rho_init(const real_t& x);
  virtual Matrix f_init(const Matrix& x, const Matrix& v);
  virtual real_t f_init(const real_t& x, const real_t& v);
  virtual Matrix g_init(const Matrix& x, const Matrix& v);
  virtual real_t g_init(const real_t& x, const real_t& v);
  virtual Matrix source(const Matrix& x, const Matrix& t);
  virtual real_t source(const real_t& x, const real_t& t);
  virtual Matrix fsource(const Matrix& x, const Matrix& v, const real_t& t);
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
  virtual Matrix Maxwell(const Matrix& v);
  virtual real_t Maxwell(const real_t& v);

  virtual void init();
  virtual void setMaxwell();
  virtual void setdt(real_t* dt);
  virtual void updateAll(const real_t& Trun, const real_t& dt);
  virtual void SolveG(const real_t& Trun, const real_t& dt, const real_t& a, model_data_* modal);
};


// class KineticDriftD_DG_IMEX_IM_Schur_period
//  : public KineticDriftD_DG_IMEX_IM_Schur
// {
// public:
//   KineticDriftD_DG_IMEX_IM_Schur_period(const KineticTensorMesh1D* mesh1D,
//                                       const fespace1D* fe,
//                                       const IMEX_RK* rk_table,
//                                       const PoissonSolver1D_period* poisson_solver,
//                                       const Solver1DType& schur_solver_type);
//   ~KineticDriftD_DG_IMEX_IM_Schur_period() override = default;

//   void init() override;
//   void setdt(real_t* dt) override;
//   void updateAll(const real_t& Trun, const real_t& dt) override;
//   void D_compute(const real_t& be, SparseMatrix* D) override;
//   void Da_ext_treat(const real_t& Trun,
//                     const real_t& dt, 
//                     const SparseMatrix& Da,
//                     std::vector<SparseMatrix>* Da_ext) override;

//   void Db_ext_treat(const real_t& Trun,
//                     const real_t& dt, 
//                     const SparseMatrix& Db,
//                     SparseMatrix* Db_ext) override;  

//   void ah_compute(const model_data_& modal, 
//                   const real_t& Trun, 
//                   const Matrix& ah_bc,
//                   Matrix* ah) override;
//   void ah_bc_compute(const Matrix& rho, 
//                     const std::vector<Matrix>& g,
//                     const real_t& Trun, 
//                     Matrix* ah_bc) override;

//   void dh_compute(const Matrix& rho, 
//                   const std::vector<Matrix>& g,
//                   const Matrix& dh_bc,
//                   Matrix* dh) override;
//   void dh_bc_compute(const Matrix& rho, 
//                     const std::vector<Matrix>& g,
//                     const real_t& Trun,
//                     Matrix* dh_bc) override;
  
//   void Mbc_compute(SparseMatrix* Mbc) override;

//   void mh_compute(const Matrix& rho, 
//                   const std::vector<Matrix>& g,
//                   const Matrix& E,
//                   Matrix* mh) override;

//   void bh_compute(const model_data_& modal, 
//                   const real_t& Trun,
//                   const std::vector<Matrix>& boundary_flux,
//                   std::vector<Matrix>* bh) override;
//   void bh_extflux_compute(const model_data_& modal, 
//                           const real_t& Trun,
//                           std::vector<Matrix>* boundary_flux) override;
  
//   void FourierMatrix(const real_t& xita, std::string& OutfilePath) override;
//   void FourierDMatrix_compute(const real_t& be, const real_t& xita, cMatrix* D) override;
//   void FourierUMatrix_compute(const cMatrix& D_NegativeWind,
//                               const cMatrix& D_PositiveWind,
//                               const real_t& xita, cMatrix* U) override;

//   void SolveRho_DichletBoundary(const real_t& Trun,
//                                 const real_t& dt,
//                                 const real_t& a,
//                                 const Matrix& dh_bc,
//                                 const Matrix& ah_bc,
//                                 const Matrix& rho_prev,
//                                 model_data_* modal) override;
  
//   virtual void SolveRho_PeriodBoundary(const real_t& Trun,
//                                       const real_t& dt,
//                                       const real_t& a,
//                                       const Matrix& rho_prev,
//                                       model_data_* modal);

//   Matrix rho_init(const Matrix& x) override;
//   real_t rho_init(const real_t& x) override;
//   Matrix f_init(const Matrix& x, const real_t& v) override;
//   real_t f_init(const real_t& x, const real_t& v) override;
//   Matrix g_init(const Matrix& x, const real_t& v) override;
//   real_t g_init(const real_t& x, const real_t& v) override;
//   Matrix source(const Matrix& x, const real_t& t) override;
//   real_t source(const real_t& x, const real_t& t) override;
//   Matrix fsource(const Matrix& x, const real_t& v, const real_t& t) override;
//   real_t fsource(const real_t& x, const real_t& v, const real_t& t) override;
//   Matrix rho_d(const Matrix& x) override;
//   real_t rho_d(const real_t& x) override;

//   virtual Matrix rho_real(const Matrix& x, const real_t& t);
//   virtual real_t rho_real(const real_t& x, const real_t& t);
//   virtual Matrix f_real(const Matrix& x, const real_t& v, const real_t& t);
//   virtual real_t f_real(const real_t& x, const real_t& v, const real_t& t);
//   virtual Matrix g_real(const Matrix& x, const real_t& v, const real_t& t);
//   virtual real_t g_real(const real_t& x, const real_t& v, const real_t& t);

//   real_t fL_bc(const int& j, const real_t& t,
//                       const model_data_& modal) override;
//   real_t fR_bc(const int& j, const real_t& t,
//                       const model_data_& modal) override;
//   real_t fL_bc(const int& j, const real_t& t,
//                       const Matrix& rho, const std::vector<Matrix>& g) override;
//   real_t fR_bc(const int& j, const real_t& t,
//                       const Matrix& rho, const std::vector<Matrix>& g) override;
//   real_t fL_explicit_bc(const int& j, const real_t& t) override;
//   real_t fR_explicit_bc(const int& j, const real_t& t) override;
//   real_t gL_bc(const int& j, const real_t& t,
//                       const model_data_& modal, const real_t& rho_L) override;
//   real_t gR_bc(const int& j, const real_t& t,
//                       const model_data_& modal, const real_t& rho_R) override;
//   real_t phi_bc(const real_t& x, const real_t& t) override;
//   real_t rho_numericalbc(const real_t& x, const real_t& t,
//                       const Matrix& rho, const std::vector<Matrix>& g) override;

//   virtual void getrho_real_modal(const real_t& Tstop, Matrix* rho_real_nodal);
//   virtual void getrho_real_nodal(const real_t& Tstop, Matrix* rho_real_nodal);
//   virtual void getg_real_nodal(const real_t& Tstop, 
//                       std::vector<Matrix>* g_real_nodal);
// };

// class KineticLinearD_DG_IMEX_IM_Schur_period
//   : public KineticDriftD_DG_IMEX_IM_Schur_period
// {
// public:
//   KineticLinearD_DG_IMEX_IM_Schur_period(const KineticTensorMesh1D* mesh1D,
//                                       const fespace1D* fe,
//                                       const IMEX_RK* rk_table,
//                                       const PoissonSolver1D_period* poisson_solver,
//                                       const Solver1DType& schur_solver_type);
//   ~KineticLinearD_DG_IMEX_IM_Schur_period() override = default;

//   void setitertol(const real_t& iter_tol) override;

//   void init() override;
//   void updateAll(const real_t& Trun, const real_t& dt) override;
//   void setdt(real_t* dt) override;

//   void SolveRho_PeriodBoundary(const real_t& Trun,
//                                 const real_t& dt,
//                                 const real_t& a,
//                                 const Matrix& rho_prev,
//                                 model_data_* modal) override;

//   void Eh_compute(const Matrix& E_modal, SparseMatrix* ME) override;

//   Matrix rho_init(const Matrix& x) override;
//   real_t rho_init(const real_t& x) override;
//   Matrix f_init(const Matrix& x, const real_t& v) override;
//   real_t f_init(const real_t& x, const real_t& v) override;
//   Matrix g_init(const Matrix& x, const real_t& v) override;
//   real_t g_init(const real_t& x, const real_t& v) override;
//   Matrix source(const Matrix& x, const real_t& t) override;
//   real_t source(const real_t& x, const real_t& t) override;
//   Matrix fsource(const Matrix& x, const real_t& v, const real_t& t) override;
//   real_t fsource(const real_t& x, const real_t& v, const real_t& t) override;
//   Matrix rho_d(const Matrix& x) override;
//   real_t rho_d(const real_t& x) override;

//   Matrix rho_real(const Matrix& x, const real_t& t) override;
//   real_t rho_real(const real_t& x, const real_t& t) override;
//   Matrix f_real(const Matrix& x, const real_t& v, const real_t& t) override;
//   real_t f_real(const real_t& x, const real_t& v, const real_t& t) override;
//   Matrix g_real(const Matrix& x, const real_t& v, const real_t& t) override;
//   real_t g_real(const real_t& x, const real_t& v, const real_t& t) override;
// };

}; // namespace QUEST

#endif // QUEST_KIETICDD_DG2D_IMEX_SCHUR_DIRICHLET_HPP

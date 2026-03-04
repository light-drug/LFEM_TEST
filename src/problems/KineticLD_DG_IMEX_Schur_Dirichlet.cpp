#include "KineticLD_DG_IMEX_Schur_Dirichlet.hpp"

namespace QUEST
{

KineticLD_DG_IMEX_IM_Schur_New::KineticLD_DG_IMEX_IM_Schur_New(const KineticTensorMesh1D* mesh1D,
                              const fespace1D* fe,
                              const IMEX_RK* rk_table,
                              const Solver1DType& schur_solver_type)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table),
    KineticLD_DG_IMEX_IM_Schur(mesh1D, fe, rk_table, schur_solver_type) {};

real_t KineticLD_DG_IMEX_IM_Schur_New::fL_explicit_bc(const int& j, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(j) >= 0) {
    f = 2.e0;
  } else {
    QUEST_ERROR(" The velocity should be inflow (v >= 0) in fL_explicit_bc ! ");
  }
  return f;
};

real_t KineticLD_DG_IMEX_IM_Schur_New::fR_explicit_bc(const int& j, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(j) <= 0) {
    f = 1.e0;
  } else {
    QUEST_ERROR(" The velocity should be inflow (v <= 0) in fR_explicit_bc ! ");
  }
  return f;
};

void KineticLD_DG_IMEX_IM_Schur_New::init() 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweight = mesh1D_->getvweights();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  const int& NTdofs = fe_->getNTdofs();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();
  // QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  // QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");

  ah_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_); 

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  };
  kinetic_modal_stages_.resize(stages_);
  D_compute(beta1_, &Da_);
  
  Da_ = - Da_;
  Da_ext_treat(0.e0, 0.e0, Da_, &Da_ext_);
  D_compute(- beta1_, &Db_); // 交替通量
  PAUSE();
  Db_ext_treat(0.e0, 0.e0, Db_, &Db_ext_);
  PAUSE();
  M_compute(&M_);
  Mbc_compute(&Mbc_);
  Minv_compute(&Minv_);
  temp_.resize(NTdofs, NTdofs);
  for (int j = 0; j < Nv; j++)
  {
    temp_ += vweight(j) * (Da_ * V(j) + Da_ext_[j])  * Minv_ * (Db_ * V(j) + Db_ext_[j]);
  };
  // Da_Theta_ = Da_ * Minv_;
};

void KineticLD_DG_IMEX_IM_Schur_New::SolveRho_DichletBoundary(const real_t& Trun,
                        const real_t& dt,
                        const real_t& a,
                        const Matrix& rhoboundary_flux,
                        const Matrix& vgboundary_flux,
                        model_data_* modal)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  real_t JacobiDet1 = JacobiDet(0);
  // ******************************** //
  // Theta_compute(Trun, dt, a, &Theta_);
  // Thetainv_compute(Trun, dt, a, &Thetainv_);
  real_t parainv = 1.e0 / (eps2_ + a * sigmas_ * dt);
  // Matrix b0_hat = Matrix::Zero(polydim, ncell);
  
  Matrix b0_bc = modal->rho;
  Matrix b0_temp;

  fe_->Assemble_Flux_bc(vgboundary_flux, &b0_temp);
  b0_temp = b0_temp * JacobiDet1;
  b0_bc = b0_bc + (dt * a) * b0_temp;  // 入流部分边界条件的处理

  b0_temp.setZero();
  fe_->Assemble_Flux_bc(rhoboundary_flux, &b0_temp); // 微观方程中rho的边界显式处理
  b0_temp = b0_temp * JacobiDet1;

  Vector b0_hat = Vector::Zero(NTdofs);
  Vector tmp = Vector::Zero(NTdofs);
  Matrix stag = Matrix::Zero(polydim, ncell);
  real_t vj;
  SparseMatrix Coef = a * dt * Minv_ * parainv;
  for (int j = 0; j < Nv; j++)
  { 
    vj = V(j);
    stag = modal->g[j] + (dt * a * V(j)) * b0_temp; // g的边界处理
    Eigen::Map<const Vector> gj(stag.data(), NTdofs);
    tmp = Coef * (vweights(j) * gj);
    b0_hat += vj * (Da_ * tmp) + (Da_ext_[j] * tmp);
  };
  Eigen::Map<Vector> rho_old_vector(b0_bc.data(), NTdofs, 1);

  Vector b0 = rho_old_vector - b0_hat; // Schur补右端项
  H_ = M_ + (a * dt) * Mbc_ + (dt * a * dt * a * parainv) * temp_; // Schur补矩阵
  
  switch (schur_solver_type_)
  {
  case Solver1DType::LU:
    // std::cout << "  Solve the Schur complement by LDLT ......\n";
    lu_.analyzePattern(H_);
    lu_.compute(H_);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case Solver1DType::GMRES:
    // std::cout << "  Solve the Schur complement by GMRES ......\n";
    gmres_.setTolerance(schur_tol_);
    gmres_.set_restart(30);
    gmres_.setMaxIterations(NTdofs);
    gmres_.compute(H_);
    // std::cout << "  iterations = " << gmres_.iterations() << std::endl;
    // std::cout << "  error = " << gmres_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;

  default:
    QUEST_ERROR(" The Schur complement Solver is not been implemented ! ");
    break;
  }
  Vector rhotemp;
  switch (schur_solver_type_)
  {
  case Solver1DType::LU:
    rhotemp = lu_.solve(b0);
    break;
  
  case Solver1DType::GMRES:
    rhotemp = gmres_.solve(b0);
    if (gmres_.info() != Eigen::Success) {
      double rel_error = gmres_.error();
      int iters = gmres_.iterations();
      QUEST_ERROR("GMRES failed to converge! "
              << " iterations = " << iters
              << ", relative residual = " << rel_error);
    }
    break;

  default:
    QUEST_ERROR(" The Schur complement Solver is not been implemented ! ");
    break;
  }
  Eigen::Map<Matrix> rho(rhotemp.data(), polydim, ncell);
  modal->rho = rho;
};

void KineticLD_DG_IMEX_IM_Schur_New::dh_explicit_extflux_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun,
                          Matrix* boundary_flux)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const int& ncell = mesh1D_->getncell();
  // *** 
  boundary_flux->resize(1, 2);
  boundary_flux->setZero();

  // Peng Schur中提到的边界取法
  real_t rho_L = 0.e0;
  real_t f;
  Vector temp;
  for (int j = 0; j < Nv; j++) {
    if (V(j) >= 0) {
      f = fL_explicit_bc(j, Trun);
    } else {
      f = (rho.col(0) + eps_ * g[j].col(0)).dot(boundary_u[1]);
    }
    rho_L += Vweights(j) * f;
  }
  real_t rho_R = 0.e0;
  for (int j = 0; j < Nv; j++) {
    if (V(j) <= 0) {
      f = fR_explicit_bc(j, Trun);
    } else {
      f = (rho.col(ncell - 1) + eps_ * g[j].col(ncell - 1)).dot(boundary_u[0]);
    }
    rho_R += Vweights(j) * f;
  }

  (*boundary_flux)(0, 0) = rho_L;
  (*boundary_flux)(0, 1) = rho_R;
}; 

void KineticLD_DG_IMEX_IM_Schur_New::Da_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Da,
                            std::vector<SparseMatrix>* Da_ext)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dv_u = fe_->getdv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = mesh1D_->getJacobiDet();
  const Vector& Jx = mesh1D_->getJx();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const std::vector<real_t>& intnormals = 
    mesh1D_->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    mesh1D_->getintboundaryneighbors();
  const int& intboundaryNum = 
    mesh1D_->getintboundaryNum();
  const std::vector<real_t>& extnormals = 
    mesh1D_->getextboundarynormal();
  const std::vector<int>& extNei =
    mesh1D_->getextboundaryneighbors();
  const int& extboundaryNum =
    mesh1D_->getextboundaryNum();
  // ******************************** //
  Da_ext->resize(Nv);
  
  // SparseMatrix da_bc(NTdofs, NTdofs); 
  
  int estimatedNonZeros = extboundaryNum * polydim * polydim;
  
#pragma omp parallel num_threads(NTH_), default(shared)
  {
    int test_cell_Index, trial_cell_Index;
    int test_basis_index, trial_basis_index;
    Vector test_qua_value, trial_qua_value;
    real_t test_normal, trial_normal;
    int alpha, beta;
    real_t ini_value;
    real_t vj;
    int j, i;
    std::vector<Eigen::Triplet<real_t>> tripletList_bc;
#pragma omp for schedule(static)
    for (j = 0; j < Nv; j++)
    {
      vj = V(j);
      tripletList_bc.clear();
      tripletList_bc.reserve(estimatedNonZeros);
      (*Da_ext)[j].resize(NTdofs, NTdofs);
      // da_bc.setZero();
      for (i = 0; i < extboundaryNum; i++) 
      {
        test_cell_Index = extNei[i];
        test_qua_value = boundary_u[1-i];
        test_normal = extnormals[i];

        trial_cell_Index = extNei[i];
        trial_qua_value = boundary_u[1-i];
        trial_normal = extnormals[i];

        for (test_basis_index = 0; test_basis_index < polydim; test_basis_index++) 
        {
          alpha = Tm(test_basis_index, test_cell_Index);
          for (trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) 
          {
            beta = Tm(trial_basis_index, trial_cell_Index);
            // std::cout << " vj = " << vj << std::endl;
            // std::cout << " test_normal = " << test_normal << std::endl;
            // std::cout << " real_t(vj * test_normal >= 0.e0 = " << real_t(vj * test_normal >= 0.e0) << std::endl;
            ini_value = vj * trial_qua_value(trial_basis_index)
                        * test_qua_value(test_basis_index) * test_normal 
                        - eps_ * CR_ * trial_qua_value(trial_basis_index) 
                        * test_qua_value(test_basis_index) * real_t(vj * test_normal > 0.e0);
            tripletList_bc.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
      // da_bc.setFromTriplets(tripletList_bc.begin(), tripletList_bc.end());
      (*Da_ext)[j].setFromTriplets(tripletList_bc.begin(), tripletList_bc.end());
      // (*Da_ext)[j] = vj * Da_ + da_bc;
    };
  }
};

void KineticLD_DG_IMEX_IM_Schur_New::Db_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Db,
                            std::vector<SparseMatrix>* Db_ext)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dv_u = fe_->getdv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = mesh1D_->getJacobiDet();
  const Vector& Jx = mesh1D_->getJx();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const std::vector<real_t>& intnormals = 
    mesh1D_->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    mesh1D_->getintboundaryneighbors();
  const int& intboundaryNum = 
    mesh1D_->getintboundaryNum();
  const std::vector<real_t>& extnormals = 
    mesh1D_->getextboundarynormal();
  const std::vector<int>& extNei =
    mesh1D_->getextboundaryneighbors();
  const int& extboundaryNum =
    mesh1D_->getextboundaryNum();
  // ******************************** //
  Db_ext->resize(Nv);
  
  // SparseMatrix db_bc(NTdofs, NTdofs); 
  
  int estimatedNonZeros = extboundaryNum * polydim * polydim;

#pragma omp parallel num_threads(NTH_), default(shared)
  {
    int test_cell_Index, trial_cell_Index;
    int test_basis_index, trial_basis_index;
    Vector test_qua_value, trial_qua_value;
    real_t test_normal, trial_normal;
    int alpha, beta;
    real_t ini_value;
    real_t vj;
    int j, i;
    std::vector<Eigen::Triplet<real_t>> tripletList_bc;
#pragma omp for schedule(static)
    for (j = 0; j < Nv; j++)
    {
      vj = V(j);
      tripletList_bc.clear();
      tripletList_bc.reserve(estimatedNonZeros);
      (*Db_ext)[j].resize(NTdofs, NTdofs);
      // db_bc.setZero();
      for (i = 0; i < extboundaryNum; i++) 
      {
        test_cell_Index = extNei[i];
        test_qua_value = boundary_u[1-i];
        test_normal = extnormals[i];

        trial_cell_Index = extNei[i];
        trial_qua_value = boundary_u[1-i];
        trial_normal = extnormals[i];

        for (test_basis_index = 0; test_basis_index < polydim; test_basis_index++) 
        {
          alpha = Tm(test_basis_index, test_cell_Index);
          for (trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) 
          {
            beta = Tm(trial_basis_index, trial_cell_Index);
            ini_value = - 0.5e0 * vj * trial_qua_value(trial_basis_index)
                * test_qua_value(test_basis_index) * test_normal;
            tripletList_bc.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
      // db_bc.setFromTriplets(tripletList_bc.begin(), tripletList_bc.end());
      (*Db_ext)[j].setFromTriplets(tripletList_bc.begin(), tripletList_bc.end());
      // (*Db_ext)[j] = vj * Db_ + db_bc;
    }
  }
};

void KineticLD_DG_IMEX_IM_Schur_New::Mbc_compute(SparseMatrix* Mbc)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = mesh1D_->getJacobiDet();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& extnormals = 
    mesh1D_->getextboundarynormal();
  const std::vector<int>& extNei = 
    mesh1D_->getextboundaryneighbors();
  const int& extboundaryNum = 
    mesh1D_->getextboundaryNum();
  // ******************************** //

  Mbc->resize(NTdofs, NTdofs);
  std::vector<Eigen::Triplet<real_t>> tripletList_Mbc;
  int estimatedNonZeros = extboundaryNum * polydim * polydim;
  tripletList_Mbc.reserve(estimatedNonZeros);

  int test_cell_Index, trial_cell_Index;
  Vector test_qua_value, trial_qua_value;
  real_t test_normal, trial_normal;
  int alpha, beta;
  real_t ini_value;
  for (int i = 0; i < extboundaryNum; i++) {
    test_cell_Index = extNei[i];
    test_qua_value = boundary_u[1-i];
    test_normal = extnormals[i];

    trial_cell_Index = extNei[i];
    trial_qua_value = boundary_u[1-i];
    trial_normal = extnormals[i];

    for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
      alpha = Tm(test_basis_index, test_cell_Index);
      for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
        beta = Tm(trial_basis_index, trial_cell_Index);
        ini_value =  (0.5e0 * CR_ * trial_qua_value(trial_basis_index)) 
            * test_qua_value(test_basis_index) * test_normal * trial_normal;
        tripletList_Mbc.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };

  Mbc->setFromTriplets(tripletList_Mbc.begin(), tripletList_Mbc.end());
};

void KineticLD_DG_IMEX_IM_Schur_New::updateAll(const real_t& Trun, 
                                          const real_t& dt)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Matrix& Vweights = mesh1D_->getvweights();
  const Matrix& Ae = rk_table_->getA();
  const Matrix& Ai = rk_table_->getAi();
  const Vector& be = rk_table_->getb();
  const Vector& bi = rk_table_->getbi();
  const Vector& ce = rk_table_->getc();
  const Vector& ci = rk_table_->getci();
  const DiagnalMatrix& M_ref = fe_->getv_u_diag();
  const DiagnalMatrix& Minv_ref = fe_->getv_u_diaginv();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const int& NTdofs = fe_->getNTdofs();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const int& extboundarynum = mesh1D_->getextboundaryNum();
  // *** 
  Matrix rhoboundary_flux = Matrix::Zero(1,extboundarynum);
  Matrix vgboundary_flux = Matrix::Zero(1,extboundarynum);
  std::vector<Matrix> gboundary_flux;
  const Matrix Matrixrho = JacobiDet(0) * M_ref * kinetic_modal_.rho;
  const std::vector<Matrix> Matrixg = eps2_ * JacobiDet(0) * M_ref * kinetic_modal_.g;
  real_t rho_L = 0.e0, rho_R = 0.e0;
  real_t f;
  Vector temp;

  for (int s = 0; s < stages_; s++)
  {   
    kinetic_modal_stages_[s].rho = Matrixrho;
    kinetic_modal_stages_[s].g = Matrixg;
    
    rhoboundary_flux.setZero();
    vgboundary_flux.setZero();
    const std::vector<Matrix>* gg_ptr = nullptr;
    if (s == 0) { gg_ptr = &kinetic_modal_.g; } 
    else { gg_ptr = &kinetic_modal_stages_[s-1].g; };
    
    rho_L = 0.e0;  rho_R = 0.e0;
    for (int j = 0; j < Nv; j++) {
      if (V(j) >= 0) {
        f = fL_explicit_bc(j, Trun);
      } else {
        f = (eps_ * gg_ptr->at(j).col(0)).dot(boundary_u[1]);
      };
      rho_L += Vweights(j) * f;
    }
    for (int j = 0; j < Nv; j++) {
      if (V(j) <= 0) {
        f = fR_explicit_bc(j, Trun);
      } else {
        f = (eps_ * gg_ptr->at(j).col(ncell - 1)).dot(boundary_u[0]);
      }
      rho_R += Vweights(j) * f;
    }
    rhoboundary_flux(0, 0) = - rho_L;
    rhoboundary_flux(0, 1) = - rho_R;

    rho_L = 0.e0;  rho_R = 0.e0;
    for (int j = 0; j < Nv; j++) {
      if (V(j) >= 0) {
        f = fL_explicit_bc(j, Trun);
        rho_L += Vweights(j) * f;
      }
    }
    for (int j = 0; j < Nv; j++) {
      if (V(j) <= 0) {
        f = fR_explicit_bc(j, Trun);
        rho_R += Vweights(j) * f;
      }
    }
    vgboundary_flux(0, 0) = - CR_ * rho_L;
    vgboundary_flux(0, 1) = CR_ * rho_R;

    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - (dt * Ai(s, i)) * ah_[i];
      #pragma omp parallel num_threads(NTH_), default(shared)
      { 
        #pragma omp for schedule(static)
        for (int j = 0; j < Nv; j++)
        {
          kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                        - (dt * Ae(s, i) * eps_)  * bh_[i][j]
                                        + (dt * Ai(s, i) * V(j)) * dh_[i]
                                        - (dt * sigmas_ * Ai(s, i)) * sh_[i][j];
        };
      };
    };

    SolveRho_DichletBoundary(Trun, dt, Ai(s, s),
                              rhoboundary_flux, vgboundary_flux,
                              &(kinetic_modal_stages_[s]));
    
    if (s == 0)
    {
      dh_explicit_extflux_compute(kinetic_modal_stages_[s].rho, 
                                  kinetic_modal_.g,
                                  Trun + ci(s) * dt, &rhoboundary_flux);
    } else 
    {
      dh_explicit_extflux_compute(kinetic_modal_stages_[s].rho, 
                                  kinetic_modal_stages_[s-1].g,
                                  Trun + ci(s) * dt, &rhoboundary_flux);
    }
    dh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, rhoboundary_flux, &(dh_[s]));
    dh_[s] = JacobiDet(0) * dh_[s];

    for (int j = 0; j < Nv; j++)
    {
      kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j] 
                                    + (dt * Ai(s, s) * V(j)) * dh_[s];
    };
    SolveG(Trun, dt, Ai(s, s), &(kinetic_modal_stages_[s]));
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));
    sh_[s] = JacobiDet(0) * sh_[s];

    ah_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &vgboundary_flux);
    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));
    ah_[s] = JacobiDet(0) * ah_[s];
    
    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
    bh_[s] = JacobiDet(0) * bh_[s];
  };
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
};

KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC(
                      const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table,
                      const Solver1DType& schur_solver_type)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table),
    KineticLD_DG_IMEX_IM_Schur(mesh1D, fe, rk_table, schur_solver_type),
    KineticLD_DG_IMEX_IM_Schur_New(mesh1D, fe, rk_table, schur_solver_type) {};

void KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::init() 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweight = mesh1D_->getvweights();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  const int& NTdofs = fe_->getNTdofs();
  // *** 
  QUEST_VERIFY(std::abs(v1 + 1.e0) < 1.e-12, " v1 must be equal to - 1.e0 !");
  QUEST_VERIFY(std::abs(v2 - 1.e0) < 1.e-12, " v2 must be equal to 1.e0 !");
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();
  // QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  // QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");

  ah_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_); 

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  };
  kinetic_modal_stages_.resize(stages_);
  
  D_compute(beta1_, &Da_);
  Da_ = - Da_;
  Da_ext_treat(0.e0, 0.e0, Da_, &Da_ext_);
  D_compute(- beta1_, &Db_); // 交替通量
  Db_ext_treat(0.e0, 0.e0, Db_, &Db_ext_);

  M_compute(&M_);
  Mbc_compute(&Mbc_);
  Minv_compute(&Minv_);
  temp_.resize(NTdofs, NTdofs);
  for (int j = 0; j < Nv; j++)
  {
    temp_ += vweight(j) * (Da_ * V(j) + Da_ext_[j])  * Minv_ * (Db_ * V(j) + Db_ext_[j]);
  };
  // Da_Theta_ = Da_ * Minv_;
};

Matrix KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::rho_init(const Matrix& x) 
{
  Matrix rho = x;
  rho.setZero();
  return rho;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::rho_init(const real_t& x) 
{
  return 0.e0;
};

Matrix KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::f_init(const Matrix& x, const real_t& v) 
{
  Matrix f = x;
  f.setZero();
  return f;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::f_init(const real_t& x, const real_t& v) 
{
  return 0.e0;
};

Matrix KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::g_init(const Matrix& x, const real_t& v) 
{
  Matrix g = Matrix::Zero(x.rows(), x.cols());
  return g;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::g_init(const real_t& x, const real_t& v) 
{
  real_t g = 0.e0;
  return g;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::fL_bc(const int& j,  const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  real_t f;
  if (V(j) >= 0) {
    f = 1.e0;
  } else {
    Vector temp = modal.rho.col(0) + eps_ * modal.g[j].col(0);
    f = temp.dot(boundary_u[1]);
  }
  return f;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::fR_bc(const int& j,  const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  real_t f;
  if (V(j) <= 0) 
  {
    f = 0.e0;
  } else {
    Vector temp = modal.rho.col(ncell - 1) + eps_ * modal.g[j].col(ncell - 1);
    f = temp.dot(boundary_u[0]);
  }
  return f;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::fL_explicit_bc(const int& j, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(j) >= 0) {
    f = 1.e0;
  } else {
    QUEST_ERROR(" The velocity should be inflow (v >= 0) in fL_explicit_bc ! ");
  }
  return f;
};

real_t KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC::fR_explicit_bc(const int& j, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(j) <= 0) {
    f = 0.e0;
  } else {
    QUEST_ERROR(" The velocity should be inflow (v <= 0) in fR_explicit_bc ! ");
  }
  return f;
};


} // namespace QUEST 
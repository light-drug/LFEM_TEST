#include "KineticLD_DG_IMEX_Schur.hpp"


namespace QUEST
{

KineticLD_DG_IMEX_IM_Schur::KineticLD_DG_IMEX_IM_Schur(const KineticTensorMesh1D* mesh1D,
                              const fespace1D* fe,
                              const IMEX_RK* rk_table,
                              const Solver1DType& schur_solver_type)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table),
    schur_solver_type_(schur_solver_type) {};

void KineticLD_DG_IMEX_IM_Schur::setsparsetol(const real_t& schur_tol)
{
  schur_tol_ = schur_tol;
};

void KineticLD_DG_IMEX_IM_Schur::setdt(real_t* dt)
{
  // ** 传入相关变量
  const int& polydim = fe_->getbasis()->getpolydim();
  const real_t& hx = mesh1D_->gethx();
  // *** 
  switch (polydim)
  {
  case 1:
    if (eps_ <= 0.5e0 * sigmas_ * hx)
    {
      *dt = 0.75e0 * hx;
    } else {
      *dt = std::min(0.75e0 * hx, eps2_ * hx / (eps_ - 0.5e0 * sigmas_ * hx));
    }
    break;

  case 2:
    if (eps_ <= 0.025e0 * sigmas_ * hx)
    {
      *dt = 0.75e0 * hx;
    } else {
      *dt = std::min(0.75e0 * hx, eps2_ * hx / (eps_ - 0.025e0 * sigmas_ * hx) / std::sqrt(10));
    }
    break;
  
  case 3:
    if (eps_ <= 0.05e0 * sigmas_ * hx)
    {
      *dt = 0.75e0 * hx;
    } else {
      *dt = std::min(0.75e0 * hx, eps2_ * hx * 0.1e0 / (eps_ - 0.05e0 * sigmas_ * hx));
    }
    break;
    
  default:
    break;
  }
};

void KineticLD_DG_IMEX_IM_Schur::init() 
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
  D_compute(- beta1_, &Db_); // 交替通量

  M_compute(&M_);
  Minv_compute(&Minv_);
  real_t dv2_sum = 0.e0;
  for (int j = 0; j < Nv; j++)
  {
    dv2_sum += vweight(j) * V(j) * V(j);
  };

  temp_ = dv2_sum * Da_ * Minv_ * Db_;
  Da_Theta_ = Da_ * Minv_;
};

void KineticLD_DG_IMEX_IM_Schur::D_compute(const real_t& be, SparseMatrix* D)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dv_u = fe_->getdv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intnormals = 
    fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    fe_->getmesh1D()->getintboundaryneighbors();
  const int& intboundaryNum = 
    fe_->getmesh1D()->getintboundaryNum();
  // ******************************** //
  D->resize(NTdofs, NTdofs);
  int estimatedNonZeros = ncell * polydim * polydim + 
          intboundaryNum * polydim * polydim * 4;
  std::vector<Eigen::Triplet<real_t>> tripletList_D;
  tripletList_D.reserve(estimatedNonZeros);
  
  Matrix Bref;
  Bref = dv_u;
  for (int i = 0; i < ncell; i++) {
    int alpha_start = i * polydim;
    for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
      for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jx(i);
        tripletList_D.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };

  int test_cell_Index, trial_cell_Index;
  Vector test_qua_value, trial_qua_value;
  real_t test_normal, trial_normal;
  int alpha, beta;
  real_t ini_value;
  for (int i = 0; i < intboundaryNum; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      test_cell_Index = intNei[i](test_cell);
      test_qua_value = boundary_u[test_cell];
      test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        trial_cell_Index = intNei[i](trial_cell);
        trial_qua_value = boundary_u[trial_cell];
        trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
          alpha = Tm(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
            beta = Tm(trial_basis_index, trial_cell_Index);
            ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
                + trial_qua_value(trial_basis_index) * be * trial_normal) 
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = - ini_value;   
            // 这是一维代码不需要积分
            tripletList_D.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };
  
  D->setFromTriplets(tripletList_D.begin(), tripletList_D.end());
};

void KineticLD_DG_IMEX_IM_Schur::Da_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Da,
                            std::vector<SparseMatrix>* Da_ext)
{
  QUEST_ERROR(" Da_ext_treat() is not implmented in this class, only exist in the Schur_New one !");
};

void KineticLD_DG_IMEX_IM_Schur::Db_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Db,
                            std::vector<SparseMatrix>* Db_ext)
{
  QUEST_ERROR(" Db_ext_treat() is not implmented in this class, only exist in the Schur_New one !");
};


void KineticLD_DG_IMEX_IM_Schur::M_compute(SparseMatrix* M)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const IntMatrix& Tm = fe_->getTm();
  // ******************************** //

  M->resize(NTdofs, NTdofs);
  std::vector<Eigen::Triplet<real_t>> tripletList_M;
  int estimatedNonZeros = NTdofs;
  tripletList_M.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell; i++) {
    int alpha_start = i * polydim;
    for (int basis_index = 0; basis_index < polydim; basis_index++) {
      real_t ini_value = v_u(basis_index,basis_index) * JacobiDet(i);
      int alpha = alpha_start + basis_index;
      tripletList_M.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  M->setFromTriplets(tripletList_M.begin(), tripletList_M.end());
};

void KineticLD_DG_IMEX_IM_Schur::Mbc_compute(SparseMatrix* Mbc)
{
  QUEST_ERROR(" Mbc_compute() is not implmented in this class, only exist in the Schur_New one !");
};

void KineticLD_DG_IMEX_IM_Schur::Minv_compute(SparseMatrix* Minv)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  // ******************************** //

  Minv->resize(NTdofs, NTdofs);
  std::vector<Eigen::Triplet<real_t>> tripletList_Minv;
  int estimatedNonZeros = NTdofs;
  tripletList_Minv.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell; i++) {
    int alpha_start = i * polydim;
    for (int basis_index = 0; basis_index < polydim; basis_index++) {
      real_t ini_value = v_u(basis_index,basis_index) * JacobiDet(i);
      ini_value = 1.e0 / ini_value;
      int alpha = alpha_start + basis_index;
      tripletList_Minv.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  Minv->setFromTriplets(tripletList_Minv.begin(), tripletList_Minv.end());
};

void KineticLD_DG_IMEX_IM_Schur::updateAll(const real_t& Trun, 
                                          const real_t& dt)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Matrix& vweights = mesh1D_->getvweights();
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
  // *** 
  Matrix rhoboundary_flux, vgboundary_flux;
  std::vector<Matrix> gboundary_flux;
  const Matrix Matrixrho = JacobiDet(0) * M_ref * kinetic_modal_.rho;
  const std::vector<Matrix> Matrixg = eps2_ * JacobiDet(0) * M_ref * kinetic_modal_.g;

  for (int s = 0; s < stages_; s++)
  { 
    kinetic_modal_stages_[s].rho = Matrixrho;
    kinetic_modal_stages_[s].g = Matrixg;
    if (s == 0)
    {
      dh_extflux_compute(kinetic_modal_, Trun, &rhoboundary_flux);
      ah_extflux_compute(kinetic_modal_, Trun, &vgboundary_flux);
    } else 
    {
      dh_extflux_compute(kinetic_modal_stages_[s-1], Trun + ci(s-1) * dt, &rhoboundary_flux);
      ah_extflux_compute(kinetic_modal_stages_[s-1], Trun + ci(s-1) * dt, &vgboundary_flux);
    }
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - (dt * Ai(s, i)) * ah_[i];
      #pragma omp parallel num_threads(NTH_), default(shared)
      { 
        real_t vj, Mj;
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

    dh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, rhoboundary_flux, &(dh_[s]));
    dh_[s] = JacobiDet(0) * dh_[s];

    #pragma omp parallel num_threads(NTH_), default(shared)
    { 
    #pragma omp for schedule(static)
      for (int j = 0; j < Nv; j++)
      {
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j] 
                                      + (dt * Ai(s, s) * V(j)) * dh_[s];
      };
    };
    SolveG(Trun, dt, Ai(s, s), &(kinetic_modal_stages_[s]));
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));
    sh_[s] = JacobiDet(0) * sh_[s];

    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));
    ah_[s] = JacobiDet(0) * ah_[s];
    
    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
    bh_[s] = JacobiDet(0) * bh_[s];
  };
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
};

void KineticLD_DG_IMEX_IM_Schur::SolveRho(const real_t& Trun,
                        const real_t& dt,
                        const real_t& a,
                        model_data_* modal)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // ******************************** //
  // Theta_compute(Trun, dt, a, &Theta_);
  // Thetainv_compute(Trun, dt, a, &Thetainv_);
  real_t parainv = 1.e0 / (eps2_ + a * sigmas_ * dt);
  Matrix b0_hat = Matrix::Zero(polydim, ncell);
  for (int j = 0; j < Nv; j++)
  {
    b0_hat += a * dt * vweights(j) * V(j) * modal->g[j];
  };
  Eigen::Map<Vector> b0(b0_hat.data(), NTdofs, 1);
  Eigen::Map<Vector> rho_old_vector(modal->rho.data(), NTdofs, 1);

  b0 = rho_old_vector - Da_Theta_ * parainv * b0;
  H_ = M_ + (dt * a * dt * a * parainv) * temp_;
  
  switch (schur_solver_type_)
  {
  case Solver1DType::CG:
    // std::cout << "  Solve the Schur complement by CG ......\n";
    cg_.setTolerance(schur_tol_);
    cg_.compute(H_);
    // std::cout << "  iterations = " << cg_.iterations() << std::endl;
    // std::cout << "  error = " << cg_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case Solver1DType::LDLT:
    // std::cout << "  Solve the Schur complement by LDLT ......\n";
    chol_.analyzePattern(H_);
    chol_.compute(H_);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case Solver1DType::GMRES:
    // std::cout << "  Solve the Schur complement by GMRES ......\n";
    gmres_.setTolerance(schur_tol_);
    gmres_.set_restart(200);
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
  case Solver1DType::CG:
    rhotemp = cg_.solve(b0);
    if (cg_.info() != Eigen::Success) {
      QUEST_ERROR(" CG failed to converge! ");
    }
    break;
  
  case Solver1DType::LDLT:
    rhotemp = chol_.solve(b0);
    break;
  
  case Solver1DType::GMRES:
    rhotemp = gmres_.solve(b0);
    if (gmres_.info() != Eigen::Success) {
      QUEST_ERROR(" GMRES failed to converge! ");
    }
    break;

  default:
    QUEST_ERROR(" The Schur complement Solver is not been implemented ! ");
    break;
  }
  Eigen::Map<Matrix> rho(rhotemp.data(), polydim, ncell);
  modal->rho = rho;
};

void KineticLD_DG_IMEX_IM_Schur::SolveRho_DichletBoundary(const real_t& Trun,
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
  Matrix b0_hat = Matrix::Zero(polydim, ncell);
  Matrix b0_bc = modal->rho;
  Matrix b0_temp;

  fe_->Assemble_Flux_bc(vgboundary_flux, &b0_temp);
  b0_temp = b0_temp * JacobiDet1;
  b0_bc = b0_bc - (dt * a) * b0_temp;  // 边界处显式处理

  b0_temp.setZero();
  fe_->Assemble_Flux_bc(rhoboundary_flux, &b0_temp); // 微观方程中rho的边界显式处理
  // std::cout << " rhoboundary_flux = \n" << rhoboundary_flux << std::endl;
  // std::cout << " JacobiDet1 = " << JacobiDet1 << std::endl;
  // PAUSE();
  b0_temp = b0_temp * JacobiDet1;
  for (int j = 0; j < Nv; j++)
  {
    b0_hat += (a * dt * vweights(j) * V(j)) * (modal->g[j] - (dt * a * V(j)) * b0_temp);
  };
  Eigen::Map<Vector> b0(b0_hat.data(), NTdofs, 1);
  Eigen::Map<Vector> rho_old_vector(b0_bc.data(), NTdofs, 1);

  b0 = rho_old_vector - Da_Theta_ * parainv * b0; // Schur补右端项
  H_ = M_ + (dt * a * dt * a * parainv) * temp_; // Schur补矩阵
  
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
    gmres_.set_restart(100);
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

void KineticLD_DG_IMEX_IM_Schur::FourierMatrix(const real_t& xita, std::string& OutfilePath,
                                              cMatrix* M_fourier, 
                                              cMatrix* D_neg_fourier,
                                              cMatrix* D_plus_fourier,
                                              cMatrix* U_fourier)
{
  // ********* 传入相关变量 ********** //
  const Matrix& v_u = fe_->getv_u();
  const Matrix& dv_u = fe_->getdv_u();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& polydim = fe_->getbasis()->getpolydim();
  // ******************************** //

  cMatrix D_NegativeWind, D_PositiveWind;
  FourierDMatrix_compute(0.5e0, xita, &D_PositiveWind);
  D_PositiveWind = - D_PositiveWind;
  FourierDMatrix_compute(- 0.5e0, xita, &D_NegativeWind);
  D_NegativeWind = - D_NegativeWind;

  *M_fourier = v_u;
  D_plus_fourier->resize(polydim, Nv * polydim);
  D_neg_fourier->resize(Nv * polydim, polydim);
  for (int i = 0; i < Nv; i++)
  {
    D_plus_fourier->block(0, i * polydim, polydim, polydim) = vweights(i) * V(i) * D_PositiveWind;
    D_neg_fourier->block(i * polydim, 0, polydim, polydim) = V(i) * D_NegativeWind;
  };
  FourierUMatrix_compute(D_NegativeWind, D_PositiveWind, xita, U_fourier);

};

void KineticLD_DG_IMEX_IM_Schur::FourierDMatrix_compute(const real_t& be, const real_t& xita, cMatrix* D) 
{
  // ********* 传入相关变量 ********** //
  const int& polydim = fe_->getbasis()->getpolydim();
  const Matrix& dv_u = fe_->getdv_u();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // ******************************** //
  D->resize(polydim, polydim);
  
  *D = - dv_u;
  creal_t ini_value;
  for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) 
  {
    for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) 
    {
      ini_value = (0.5e0 + be) * boundary_u[0](trial_basis_index) 
                  * boundary_u[1](test_basis_index) 
                  * std::exp(creal_t(0.e0, -1.e0 * xita))
                + (0.5e0 - be) * boundary_u[1](trial_basis_index) 
                  * boundary_u[1](test_basis_index)
                - (0.5e0 + be) * boundary_u[0](trial_basis_index) 
                  * boundary_u[0](test_basis_index)
                - (0.5e0 - be) * boundary_u[1](trial_basis_index) 
                  * boundary_u[0](test_basis_index) 
                  * std::exp(creal_t(0.e0, xita));
      ini_value = - ini_value;
      (*D)(test_basis_index, trial_basis_index) += ini_value;
    };
  };
};

void KineticLD_DG_IMEX_IM_Schur::FourierUMatrix_compute(
                              const cMatrix& D_NegativeWind,
                              const cMatrix& D_PositiveWind,
                              const real_t& xita, 
                              cMatrix* U)
{
  // ********* 传入相关变量 ********** //
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& polydim = fe_->getbasis()->getpolydim();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dv_u = fe_->getdv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  // ******************************** //
  U->resize(Nv * polydim, Nv * polydim);
  U->setZero();
  real_t test_v, trial_v;
  for (int test_i = 0; test_i < Nv; test_i++)
  {
    test_v = V(test_i);
    for (int trial_i = 0; trial_i < Nv; trial_i++)
    {
      trial_v = V(trial_i);
      
      for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++)
      {
        for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++)
        {
          int alpha = test_i * polydim + test_basis_index;
          int beta = trial_i * polydim + trial_basis_index;

          if(test_i == trial_i)
          {
            (*U)(alpha, beta) = (test_v <= 0) ? 
                test_v * D_NegativeWind(test_basis_index, trial_basis_index) :
                test_v * D_PositiveWind(test_basis_index, trial_basis_index);
          }
          (*U)(alpha, beta) -= (test_v <= 0) ? 
                vweights(trial_i) * trial_v * D_NegativeWind(test_basis_index, trial_basis_index) :
                vweights(trial_i) * trial_v * D_PositiveWind(test_basis_index, trial_basis_index);
        }
      }
    }
  }
};

void KineticLD_DG_IMEX_IM_Schur::SolveG(const real_t& Trun,
                      const real_t& dt,
                      const real_t& a,
                      model_data_* modal)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const IntMatrix& Tm = fe_->getTm();
  const DiagnalMatrix& M_ref = fe_->getv_u_diag();
  const DiagnalMatrix& Minv_ref = fe_->getv_u_diaginv();
  const int& Nv = mesh1D_->getNv();
  // ******************************** //
  real_t parainv = 1.e0 / (eps2_ + a * sigmas_ * dt) / JacobiDet(0);
  modal->g = (parainv * Minv_ref) * modal->g;
};

KineticLD_DG_IMEX_IM_Schur_IsentropicBC::KineticLD_DG_IMEX_IM_Schur_IsentropicBC(
                      const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table,
                      const Solver1DType& schur_solver_type)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table),
    KineticLD_DG_IMEX_IM_Schur(mesh1D, fe, rk_table, schur_solver_type) {};

void KineticLD_DG_IMEX_IM_Schur_IsentropicBC::init() 
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
  D_compute(- beta1_, &Db_); // 交替通量

  M_compute(&M_);
  Minv_compute(&Minv_);
  real_t dv2_sum = 0.e0;
  for (int j = 0; j < Nv; j++)
  {
    dv2_sum += vweight(j) * V(j) * V(j);
  };

  temp_ = dv2_sum * Da_ * Minv_ * Db_;
  Da_Theta_ = Da_ * Minv_;
};

Matrix KineticLD_DG_IMEX_IM_Schur_IsentropicBC::rho_init(const Matrix& x) 
{
  Matrix rho = x;
  rho.setZero();
  return rho;
};

real_t KineticLD_DG_IMEX_IM_Schur_IsentropicBC::rho_init(const real_t& x) 
{
  return 0.e0;
};

Matrix KineticLD_DG_IMEX_IM_Schur_IsentropicBC::f_init(const Matrix& x, const real_t& v) 
{
  Matrix f = x;
  f.setZero();
  return f;
};

real_t KineticLD_DG_IMEX_IM_Schur_IsentropicBC::f_init(const real_t& x, const real_t& v) 
{
  return 0.e0;
};

Matrix KineticLD_DG_IMEX_IM_Schur_IsentropicBC::g_init(const Matrix& x, const real_t& v) 
{
  Matrix g = Matrix::Zero(x.rows(), x.cols());
  return g;
};

real_t KineticLD_DG_IMEX_IM_Schur_IsentropicBC::g_init(const real_t& x, const real_t& v) 
{
  real_t g = 0.e0;
  return g;
};

real_t KineticLD_DG_IMEX_IM_Schur_IsentropicBC::fL_bc(const int& j,  const real_t& t,
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

real_t KineticLD_DG_IMEX_IM_Schur_IsentropicBC::fR_bc(const int& j,  const real_t& t,
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

KineticLD_DG_IMEX_IM_Schur_twovel_period::KineticLD_DG_IMEX_IM_Schur_twovel_period(
                              const KineticTensorMesh1D* mesh1D,
                              const fespace1D* fe,
                              const IMEX_RK* rk_table,
                              const Solver1DType& schur_solver_type)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table),
    KineticLD_DG_IMEX_IM_Schur(mesh1D, fe, rk_table, schur_solver_type),
    Kinetic1D_LD_DG_IMEX_EX_twovel_period(mesh1D, fe, rk_table) {};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::init() 
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
  QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");
  QUEST_VERIFY(Nv == 2, " Nv must equal to 2 !");
  r_ = 2.e0 / (1.e0 + std::sqrt(1.e0 - 4.e0 * eps2_));

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
  D_compute(- beta1_, &Db_); // 交替通量

  M_compute(&M_);
  Minv_compute(&Minv_);
  real_t dv2_sum = 0.e0;
  for (int j = 0; j < Nv; j++)
  {
    dv2_sum += vweight(j) * V(j) * V(j);
  };

  temp_ = dv2_sum * Da_ * Minv_ * Db_;
  Da_Theta_ = Da_ * Minv_;
};

Matrix KineticLD_DG_IMEX_IM_Schur_twovel_period::
  rho_init(const Matrix& x)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_init(x);
};

real_t KineticLD_DG_IMEX_IM_Schur_twovel_period::
  rho_init(const real_t& x)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_init(x);
};

Matrix KineticLD_DG_IMEX_IM_Schur_twovel_period::
  f_init(const Matrix& x, const real_t& v)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_init(x, v);
};

real_t KineticLD_DG_IMEX_IM_Schur_twovel_period::
  f_init(const real_t& x, const real_t& v)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_init(x, v);
};

Matrix KineticLD_DG_IMEX_IM_Schur_twovel_period::
  g_init(const Matrix& x, const real_t& v)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_init(x, v);
};

real_t KineticLD_DG_IMEX_IM_Schur_twovel_period::
  g_init(const real_t& x, const real_t& v)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_init(x, v);
};

Matrix KineticLD_DG_IMEX_IM_Schur_twovel_period::
  rho_real(const Matrix& x, const real_t& t) 
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_real(x, t);
};

real_t KineticLD_DG_IMEX_IM_Schur_twovel_period::
  rho_real(const real_t& x, const real_t& t)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_real(x, t);
};

Matrix KineticLD_DG_IMEX_IM_Schur_twovel_period::
  f_real(const Matrix& x, const real_t& v, const real_t& t)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_real(x, v, t);
};

real_t KineticLD_DG_IMEX_IM_Schur_twovel_period::
  f_real(const real_t& x, const real_t& v, const real_t& t)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_real(x, v, t);
};

Matrix KineticLD_DG_IMEX_IM_Schur_twovel_period::
  g_real(const Matrix& x, const real_t& v, const real_t& t)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_real(x, v, t);
};

real_t KineticLD_DG_IMEX_IM_Schur_twovel_period::
  g_real(const real_t& x, const real_t& v, const real_t& t)
{
  return Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_real(x, v, t);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  setdt(real_t* dt)
{
  KineticLD_DG_IMEX_IM_Schur::setdt(dt);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  updateAll(const real_t& Trun, const real_t& dt)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Matrix& vweights = mesh1D_->getvweights();
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
  // *** 
  Matrix rhoboundary_flux, vgboundary_flux;
  std::vector<Matrix> gboundary_flux;
  const Matrix Matrixrho = JacobiDet(0) * M_ref * kinetic_modal_.rho;
  const std::vector<Matrix> Matrixg = eps2_ * JacobiDet(0) * M_ref * kinetic_modal_.g;
  PAUSE();
  for (int s = 0; s < stages_; s++)
  { 
    kinetic_modal_stages_[s].rho = Matrixrho;
    kinetic_modal_stages_[s].g = Matrixg;
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - (dt * Ai(s, i)) * ah_[i];
      for (int j = 0; j < Nv; j++)
      {
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                      - (dt * Ae(s, i) * eps_)  * bh_[i][j]
                                      + (dt * Ai(s, i) * V(j)) * dh_[i]
                                      - (dt * sigmas_ * Ai(s, i)) * sh_[i][j];
      };
    };
    
    SolveRho(Trun, dt, Ai(s, s), &(kinetic_modal_stages_[s]));

    dh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &rhoboundary_flux);
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

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  ah_compute(const model_data_& modal, 
                const real_t& Trun, 
                const Matrix& boundary_flux,
                Matrix* ah)
{
  Kinetic1D_LD_DG_IMEX_EX::ah_compute(modal, Trun, boundary_flux, ah);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  bh_compute(const model_data_& modal, 
              const real_t& Trun,
              const std::vector<Matrix>& boundary_flux,
              std::vector<Matrix>* bh)
{
  Kinetic1D_LD_DG_IMEX_EX::bh_compute(modal, Trun, boundary_flux, bh);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  dh_compute(const model_data_& modal, 
              const real_t& Trun,
              const Matrix& boundary_flux,
              Matrix* dh)
{
  Kinetic1D_LD_DG_IMEX_EX::dh_compute(modal, Trun, boundary_flux, dh);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  sh_compute(const model_data_& modal, 
              const real_t& Trun,
              std::vector<Matrix>* sh)
{
  Kinetic1D_LD_DG_IMEX_EX::sh_compute(modal, Trun, sh);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  ah_extflux_compute(const model_data_& modal, 
                    const real_t& Trun, 
                    Matrix* boundary_flux)
{
  Kinetic1D_LD_DG_IMEX_EX_twovel_period::ah_extflux_compute(modal, Trun, boundary_flux);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  bh_extflux_compute(const model_data_& modal, 
                    const real_t& Trun, 
                    std::vector<Matrix>* boundary_flux)
{
  Kinetic1D_LD_DG_IMEX_EX_twovel_period::bh_extflux_compute(modal, Trun, boundary_flux);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  dh_extflux_compute(const model_data_& modal, 
                    const real_t& Trun, 
                    Matrix* boundary_flux)
{
  Kinetic1D_LD_DG_IMEX_EX_twovel_period::dh_extflux_compute(modal, Trun, boundary_flux);
};

void KineticLD_DG_IMEX_IM_Schur_twovel_period::
  D_compute(const real_t& be, SparseMatrix* D)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dv_u = fe_->getdv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intnormals = 
    fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    fe_->getmesh1D()->getintboundaryneighbors();
  const int& intboundaryNum = 
    fe_->getmesh1D()->getintboundaryNum();
  const std::vector<real_t>& extnormals = 
    fe_->getmesh1D()->getextboundarynormal();
  const std::vector<Eigen::Vector2i>& extNei_period =
    fe_->getmesh1D()->getextboundaryneighbors_period();
  const int& extboundaryNum =
    fe_->getmesh1D()->getextboundaryNum();
  // ******************************** //
  D->resize(NTdofs, NTdofs);
  int estimatedNonZeros = ncell * polydim * polydim + 
          (intboundaryNum + 1) * polydim * polydim * 4;
  std::vector<Eigen::Triplet<real_t>> tripletList_D;
  tripletList_D.reserve(estimatedNonZeros);
  
  Matrix Bref;
  Bref = dv_u;
  for (int i = 0; i < ncell; i++) {
    int alpha_start = i * polydim;
    for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
      for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jx(i);
        tripletList_D.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };

  int test_cell_Index, trial_cell_Index;
  Vector test_qua_value, trial_qua_value;
  real_t test_normal, trial_normal;
  int alpha, beta;
  real_t ini_value;
  for (int i = 0; i < intboundaryNum; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      test_cell_Index = intNei[i](test_cell);
      test_qua_value = boundary_u[test_cell];
      test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        trial_cell_Index = intNei[i](trial_cell);
        trial_qua_value = boundary_u[trial_cell];
        trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
          alpha = Tm(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
            beta = Tm(trial_basis_index, trial_cell_Index);
            ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
                + trial_qua_value(trial_basis_index) * be * trial_normal) 
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = - ini_value;   
            // 这是一维代码不需要积分
            tripletList_D.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int test_cell = 0; test_cell < 2; test_cell++) {
    test_cell_Index = extNei_period[0](test_cell);
    test_qua_value = boundary_u[1 - test_cell];
    test_normal = extnormals[0] * std::pow(-1.e0,test_cell);
    for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
      trial_cell_Index = extNei_period[0](trial_cell);
      trial_qua_value = boundary_u[1 - trial_cell];
      trial_normal = extnormals[0] * std::pow(-1.e0,trial_cell);
      for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
        alpha = Tm(test_basis_index, test_cell_Index);
        for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
          beta = Tm(trial_basis_index, trial_cell_Index);
          ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
              + trial_qua_value(trial_basis_index) * be * trial_normal) 
              * test_qua_value(test_basis_index) * test_normal;
          ini_value = - ini_value;   
          // 这是一维代码不需要积分
          tripletList_D.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
        };
      };
    };
  };
  
  D->setFromTriplets(tripletList_D.begin(), tripletList_D.end());
};

} // namespace QUEST
  
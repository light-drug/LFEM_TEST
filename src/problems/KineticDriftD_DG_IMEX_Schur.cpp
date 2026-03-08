#include "KineticDriftD_DG_IMEX_Schur.hpp"

namespace QUEST
{

void KineticDriftD_DG_IMEX_IM_Schur::seteps(const real_t& knu)
{
  eps_ = knu;
  eps2_ = knu * knu;
};

void KineticDriftD_DG_IMEX_IM_Schur::setCR(const real_t& CR)
{
  CR_ = CR;
};

void KineticDriftD_DG_IMEX_IM_Schur::setbeta1(const real_t& beta1)
{
  beta1_ = beta1;
};

void KineticDriftD_DG_IMEX_IM_Schur::setsigmas(const real_t& sigmas)
{
  sigmas_ = sigmas;
};

void KineticDriftD_DG_IMEX_IM_Schur::setsparsetol(const real_t& schur_tol)
{
  schur_tol_ = schur_tol;
};

void KineticDriftD_DG_IMEX_IM_Schur::setitertol(const real_t& iter_tol)
{
  iter_tol_ = iter_tol;
};

void KineticDriftD_DG_IMEX_IM_Schur::setNTH(const int& NTH)
{
  NTH_ = NTH;
};

void KineticDriftD_DG_IMEX_IM_Schur::settheta(const real_t& theta)
{
  theta_ = theta;
};

void KineticDriftD_DG_IMEX_IM_Schur::setgamma(const real_t& gamma)
{
  gamma_ = gamma;
};

const Matrix& KineticDriftD_DG_IMEX_IM_Schur::getrho_modal() const 
{
  return kinetic_modal_.rho;
};

const std::vector<Matrix>& KineticDriftD_DG_IMEX_IM_Schur::getg_modal() const 
{
  return kinetic_modal_.g;
};

const Matrix& KineticDriftD_DG_IMEX_IM_Schur::getrho_d_modal() const 
{
  return rho_d_modal_;
};


KineticDriftD_DG_IMEX_IM_Schur::KineticDriftD_DG_IMEX_IM_Schur(const KineticTensorMesh1D* mesh1D,
                                const fespace1D* fe,
                                const IMEX_RK* rk_table,
                                const PoissonSolver1D* poisson_solver,
                                const Solver1DType& schur_solver_type)
  : mesh1D_(mesh1D), fe_(fe), rk_table_(rk_table), poisson_solver_(poisson_solver),
   schur_solver_type_(schur_solver_type) {};

void KineticDriftD_DG_IMEX_IM_Schur::init()
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& numqua = fe_->getnumqua(); 
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweight = mesh1D_->getvweights();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  const int& NTdofs = fe_->getNTdofs();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();
  setMaxwell();
  std::cout << " Maxwell sum = " << Maxwell_sum_ << std::endl;
  // QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  // QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");
  ah_.resize(stages_);
  ch_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_); 
  mh_.resize(stages_);
  source_sumh_.resize(stages_);
  sourceh_.resize(stages_);

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_d(x); }, &rho_d_modal_);
  fe_->Interpolate_Initial(
      [this](const Matrix& x) { return this->rho_d(x); }, &rho_d_nodal_);
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  };
  kinetic_modal_stages_.resize(stages_);
  zero_modal_mat_ = Matrix::Zero(polydim, ncell);
  zero_nodal_mat_ = Matrix::Zero(numqua, ncell);

  D_compute(beta1_, &Da_);
  Da_ = - Da_;
  Da_ext_treat(0.e0, 0.e0, Da_, &Da_ext_);
  D_compute(- beta1_, &Db_); // 交替通量
  Db_ext_treat(0.e0, 0.e0, Db_, &Db_ext_);

  M_compute(&M_);
  Mbc_compute(&Mbc_);
  Minv_compute(&Minv_);
  temp_.resize(NTdofs, NTdofs);
  real_t Mj;
  for (int j = 0; j < Nv; j++)
  {
    Mj = Maxwell(V(j));
    temp_ += vweight(j) * Mj * (Da_ * V(j) + Da_ext_[j])  * (Minv_ * (Db_ * V(j) + V(j) * Db_ext_));
  };

  Matrix phi_Dirichlet(1, 2);
  phi_Dirichlet(0, 0) = phi_bc(x1, 0.e0);
  phi_Dirichlet(0, 1) = phi_bc(x2, 0.e0);
  Matrix rho_nodal;
  fe_->modal_to_nodal1D(kinetic_modal_.rho, &rho_nodal);
  Matrix poi_lhs = (rho_d_nodal_ - rho_nodal) / gamma_;  // 泊松方程右端项
  poisson_solver_->solveall(phi_Dirichlet, poi_lhs, &(kinetic_modal_.E), &(kinetic_modal_.phi)); 
};

void KineticDriftD_DG_IMEX_IM_Schur::setMaxwell()
{
  // ** 传入相关变量
  const int& Nv = mesh1D_->getNv();
  const Vector& vweight = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t vj, Mj;
  Maxwell_sum_ = 0.e0;
  for (int j = 0; j < Nv; j++)
  {
    vj = V(j);
    Mj = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) * std::exp( - vj * vj / ( 2 * theta_));
    Maxwell_sum_ += vweight(j) * Mj;
  }
};

void KineticDriftD_DG_IMEX_IM_Schur::setdt(real_t* dt)
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
  *dt = gamma_ * (*dt);
};

void KineticDriftD_DG_IMEX_IM_Schur::D_compute(const real_t& be, SparseMatrix* D)
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

void KineticDriftD_DG_IMEX_IM_Schur::Da_ext_treat(const real_t& Trun,
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
  };
};

void KineticDriftD_DG_IMEX_IM_Schur::Db_ext_treat(const real_t& Trun,
                            const real_t& dt, 
                            const SparseMatrix& Db,
                            SparseMatrix* Db_ext)
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

  // SparseMatrix db_bc(NTdofs, NTdofs); 
  std::vector<Eigen::Triplet<real_t>> tripletList_bc;
  int estimatedNonZeros = extboundaryNum * polydim * polydim;
  tripletList_bc.reserve(estimatedNonZeros);

  int test_cell_Index, trial_cell_Index;
  int test_basis_index, trial_basis_index;
  Vector test_qua_value, trial_qua_value;
  real_t test_normal, trial_normal;
  int alpha, beta;
  real_t ini_value;
  int i;
  Db_ext->resize(NTdofs, NTdofs);
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
        ini_value = - 0.5e0 * trial_qua_value(trial_basis_index)
            * test_qua_value(test_basis_index) * test_normal;
        tripletList_bc.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };
  Db_ext->setFromTriplets(tripletList_bc.begin(), tripletList_bc.end());
};

void KineticDriftD_DG_IMEX_IM_Schur::M_compute(SparseMatrix* M)
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

void KineticDriftD_DG_IMEX_IM_Schur::Mbc_compute(SparseMatrix* Mbc)
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

void KineticDriftD_DG_IMEX_IM_Schur::Minv_compute(SparseMatrix* Minv)
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

void KineticDriftD_DG_IMEX_IM_Schur::Eh_compute(const Matrix& E_modal, SparseMatrix* ME)
{
  // ********* 传入相关变量 ********** //
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const Matrix& test_ref = fe_->gettest_ref();
  const Matrix& test_ref_T = fe_->gettest_ref_T();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  // ******************************** //
  Matrix E_nodal;
  fe_->modal_to_nodal1D(E_modal, &E_nodal);
  ME->resize(NTdofs, NTdofs);
  std::vector<Eigen::Triplet<real_t>> tripletList_ME;
  int estimatedNonZeros = ncell * polydim * polydim;
  tripletList_ME.reserve(estimatedNonZeros);
  
  Matrix local_ME;
  real_t ini_value;
  int alpha, beta;
  for (int i = 0; i < ncell; i++) {
    int alpha_start = i * polydim;
    local_ME = test_ref.array().rowwise() * E_nodal.col(i).transpose().array() * JacobiDet(i);
    local_ME = local_ME * wqua_diag * test_ref_T;
    for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) 
    {
      for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) 
      {
        ini_value = local_ME(test_basis_index, trial_basis_index);
        alpha = alpha_start + test_basis_index;
        beta = alpha_start + trial_basis_index;
        tripletList_ME.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      }
    };
  };
  ME->setFromTriplets(tripletList_ME.begin(), tripletList_ME.end());
};

void KineticDriftD_DG_IMEX_IM_Schur::FourierMatrix(const real_t& xita, std::string& OutfilePath)
{

};

void KineticDriftD_DG_IMEX_IM_Schur::FourierDMatrix_compute(const real_t& be, const real_t& xita, cMatrix* D)
{

};

void KineticDriftD_DG_IMEX_IM_Schur::FourierUMatrix_compute(const cMatrix& D_NegativeWind,
                                    const cMatrix& D_PositiveWind,
                                    const real_t& xita, 
                                    cMatrix* U)
{

};

void KineticDriftD_DG_IMEX_IM_Schur::fluxint_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          Matrix* flux_int)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& intboundarynum = fe_->getmesh1D()->getintboundaryNum();
  // const int& extboundarynum = fe_->getmesh1D()->getextboundarynum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  // const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  // const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  
  *flux_int = Matrix::Zero(1, intboundarynum);
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal;
  for (int i = 0; i < intboundarynum; i++)
  {
    Cellindex0 = IntBNei[i](0);
    Cellindex1 = IntBNei[i](1);
    val0 = a * modal.col(Cellindex0).dot(boundary_u[0]);
    val1 = a * modal.col(Cellindex1).dot(boundary_u[1]);
    normal = IntBNormal[i];
    (*flux_int)(0, i) = (a * normal > 0) ? val0 : val1;
  };
};

void KineticDriftD_DG_IMEX_IM_Schur::fluxext_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          const Matrix& Dirichlet,
                          Matrix* flux_ext)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  // const int& intboundarynum = fe_->getmesh1D()->getintboundarynum();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  // const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  // const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  
  *flux_ext = Matrix::Zero(1, extboundarynum);
  int Cellindex0, Cellindex1;
  Vector temp0, temp1;
  real_t val0, val1;
  real_t normal;
  for (int i = 0; i < extboundarynum; i++)
  {
    Cellindex0 = ExtBNei[i];
    temp0 = modal.col(Cellindex0).transpose();
    val0 = a * temp0.dot(boundary_u[1 - i]);
    val1 = a * Dirichlet(0, i);
    normal = ExtBNormal[i];
    (*flux_ext)(0, i) = (a * normal > 0) ? val0 : val1;
  };
};

void KineticDriftD_DG_IMEX_IM_Schur::ah_compute(const model_data_& modal, 
        const real_t& Trun,
        const Matrix& ah_bc,
        Matrix* ah)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  // Matrix vg_modal = Matrix::Zero(polydim, ncell);
  // for (int j = 0; j < Nv; j++) {
  //   vg_modal += vweights(j) * V(j) * modal.g[j];
  // };
  // Matrix vg_nodal;
  // fe_->modal_to_nodal1D(vg_modal, &vg_nodal);
  // Matrix rhs;
  // fe_->Assemble_F(vg_nodal, 1, &rhs);
  // Matrix flux_int;
  // fluxint_compute(vg_modal, beta1_, &flux_int);
  // // fluxext_compute(vg_modal, beta1_, boundary_flux, &flux_ext);  // 交错通量 beta
  // Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  // fe_->Assemble_Flux(flux_int, &flux_rhs);
  // *ah = - rhs + flux_rhs;

  real_t vj;
  ah->resize(polydim, ncell);
  ah->setZero();
  Eigen::Map<Vector> ah_vec(ah->data(), NTdofs);
  Eigen::Map<const Vector> rho_vec(modal.rho.data(), NTdofs);
  ah_vec.setZero();
  for (int j = 0; j < Nv; j++)
  { 
    vj = V(j);
    Eigen::Map<const Vector> gj(modal.g[j].data(), NTdofs);
    ah_vec += (vj * (Da_ * gj) + (Da_ext_[j] * gj)) * vweights(j);
  };
  ah_vec += Mbc_ * rho_vec;
  *ah -= ah_bc;
};

void KineticDriftD_DG_IMEX_IM_Schur::ah_bc_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun, 
                          Matrix* ah_bc)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  // *** 
  Matrix vgboundary_flux(1, extboundarynum);
  real_t rho_L = 0.e0;  
  real_t rho_R = 0.e0;
  real_t f;
  for (int j = 0; j < Nv; j++) {
    if (V(j) >= 0) {
      f = fL_explicit_bc(j, Trun);
      rho_L += vweights(j) * f;
    }
  }
  for (int j = 0; j < Nv; j++) {
    if (V(j) <= 0) {
      f = fR_explicit_bc(j, Trun);
      rho_R += vweights(j) * f;
    }
  }
  vgboundary_flux(0, 0) = - CR_ * rho_L;
  vgboundary_flux(0, 1) = CR_ * rho_R;
  fe_->Assemble_Flux_bc(vgboundary_flux, ah_bc);
  *ah_bc *= JacobiDet1;
};

void KineticDriftD_DG_IMEX_IM_Schur::bh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          const std::vector<Matrix>& boundary_flux,
                          std::vector<Matrix>* bh) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 
  std::vector<Matrix> g_nodal;
  fe_->modal_to_nodal1D(modal.g, &g_nodal);
  Matrix rhs;
  Matrix flux_int, flux_ext;
  Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  Matrix dh = Matrix::Zero(polydim, ncell);
  Matrix temp;
  bh->resize(Nv);
  for (int j = 0; j < Nv; j++) {
    temp = V(j) * g_nodal[j];
    fe_->Assemble_F(temp, 1, &rhs);
    fluxint_upwind_compute(modal.g[j], V(j), &flux_int);
    // fluxext_upwind_compute(modal.g[j], V(j), boundary_flux[j], &flux_ext);
    fe_->Assemble_Flux(flux_int, boundary_flux[j], &flux_rhs);
    (*bh)[j] = (- rhs + flux_rhs);
    dh = dh + vweights(j) * (*bh)[j];  // <vg_x>
  };
  // (I - \Pi)(vg_x)
  // \Pi(vg_x) = <vg_x> M
  real_t Mj;
  for (int j = 0; j < Nv; j++)
  {
    Mj = Maxwell(V(j));
    (*bh)[j] = (*bh)[j] - dh * Mj;
    (*bh)[j] *= JacobiDet1;
  };

};

void KineticDriftD_DG_IMEX_IM_Schur::bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  boundary_flux->resize(Nv);
  real_t rho_L = rho_numericalbc(x1, Trun, modal.rho, modal.g);
  real_t rho_R = rho_numericalbc(x2, Trun, modal.rho, modal.g);

  int Cellindex0;
  real_t temp0, temp1;
  for (int j = 0; j < Nv; j++) 
  { 
    boundary_flux->at(j).resize(1, 2);
    boundary_flux->at(j).setZero();
    if (V(j) >= 0)
    {
      Cellindex0 = ExtBNei[1];
      temp0 = V(j) * modal.g[j].col(Cellindex0).dot(boundary_u[0]);
      boundary_flux->at(j)(0, 0) = V(j) * gL_bc(j, Trun, modal, rho_L);
      boundary_flux->at(j)(0, 1) = temp0; 
    } else if (V(j) < 0)
    {
      Cellindex0 = ExtBNei[0];
      temp0 = V(j) * modal.g[j].col(Cellindex0).dot(boundary_u[1]);
      boundary_flux->at(j)(0, 0) = temp0;
      boundary_flux->at(j)(0, 1) = V(j) * gR_bc(j, Trun, modal, rho_R);
    };
  };
};

void KineticDriftD_DG_IMEX_IM_Schur::dh_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const Matrix& dh_bc,
                          Matrix* dh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  // Matrix rho_nodal;
  // fe_->modal_to_nodal1D(modal.rho, &rho_nodal);
  // Matrix rhs;
  // fe_->Assemble_F(rho_nodal, 1, &rhs);
  // Matrix flux_int;
  // fluxint_compute(modal.rho, - beta1_, &flux_int);
  // Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  // fe_->Assemble_Flux_int(flux_int, &flux_rhs);
  
  // *dh = rhs - flux_rhs;
  dh->resize(polydim, ncell);
  dh->setZero();
  Eigen::Map<Vector> dh_vec(dh->data(), NTdofs);

  Eigen::Map<const Vector> rho_vec(rho.data(), NTdofs);
  dh_vec = Db_ * rho_vec + Db_ext_ * rho_vec;
  
  *dh += dh_bc;
};

void KineticDriftD_DG_IMEX_IM_Schur::dh_bc_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun, 
                          Matrix* dh_bc)
{ 
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const Vector& V = mesh1D_->getV();
  const int& ncell = mesh1D_->getncell();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  // *** 
  real_t rho_L = 0.e0, rho_R = 0.e0;
  real_t f;
  Matrix rhoboundary_flux(1, extboundarynum);
  for (int j = 0; j < Nv; j++) {
    if (V(j) >= 0) {
      f = fL_explicit_bc(j, Trun);
    } else {
      f = (eps_ * g[j].col(0)).dot(boundary_u[1]);
    };
    rho_L += Vweights(j) * f;
  }
  for (int j = 0; j < Nv; j++) {
    if (V(j) <= 0) {
      f = fR_explicit_bc(j, Trun);
    } else {
      f = (eps_ * g[j].col(ncell - 1)).dot(boundary_u[0]);
    }
    rho_R += Vweights(j) * f;
  }
  rhoboundary_flux(0, 0) = - rho_L;
  rhoboundary_flux(0, 1) = - rho_R;

  fe_->Assemble_Flux_bc(rhoboundary_flux, dh_bc); // 微观方程中rho的边界显式处理
  *dh_bc *= JacobiDet1;
};

void KineticDriftD_DG_IMEX_IM_Schur::mh_compute(const Matrix& rho, 
                const std::vector<Matrix>& g,
                const Matrix& E,
                Matrix* mh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const Vector& V = mesh1D_->getV();
  const int& ncell = mesh1D_->getncell();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 
  Matrix E_nodal;
  fe_->modal_to_nodal1D(E, &E_nodal);
  
  Matrix rho_nodal;
  fe_->modal_to_nodal1D(rho, &rho_nodal);
  E_nodal = E_nodal.array() * rho_nodal.array();
  
  fe_->Assemble_F(E_nodal, 0, mh);
  *mh *= JacobiDet1;
};

void KineticDriftD_DG_IMEX_IM_Schur::sh_compute(const model_data_& modal,  
                          const real_t& Trun, 
                          std::vector<Matrix>* sh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 
  
  std::vector<Matrix> g_nodal;
  fe_->modal_to_nodal1D(modal.g, &g_nodal);
  std::vector<Matrix> rhs;
  fe_->Assemble_F(g_nodal, 0, &rhs);
  *sh = rhs * JacobiDet1;
};

void KineticDriftD_DG_IMEX_IM_Schur::ch_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* ch)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& numqua = fe_->getnumqua(); 
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const real_t& hv = mesh1D_->gethv();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 
  real_t hvinv = 1.e0 / hv;
  real_t hvinv2 = 1.e0 / (2.e0 * hv);
  int wind;
  Matrix E_nodal;
  fe_->modal_to_nodal1D(modal.E, &E_nodal);
  std::vector<Matrix> g_nodal;
  fe_->modal_to_nodal1D(modal.g, &g_nodal);
  Matrix g_diff_modal = Matrix::Zero(polydim, ncell);
  Matrix g_diff_nodal = Matrix::Zero(numqua, ncell);
  ch->resize(Nv);
  
  switch(polydim)
  {
  case 1:
  {
    const Matrix* gj_1;
    for (int j = 0; j < Nv; j++) 
    {
      g_diff_nodal.setZero();
      for (int i = 0; i < ncell; i++)
      {
        wind = modal.E(0, i) >= 0 ? -1 : 1; 
        if (j-wind >= 0 && j-wind < Nv) { gj_1 = &(g_nodal[j-wind]); }
        else { gj_1 = &zero_nodal_mat_; };
        
        g_diff_nodal.col(i) = (g_nodal[j].col(i) - gj_1->col(i)) * hvinv * real_t(wind);
      }
      g_diff_nodal = - 1.e0 * g_diff_nodal.array() * E_nodal.array();
      fe_->Assemble_F(g_diff_nodal, 0, &((*ch)[j]));
      (*ch)[j] *= JacobiDet1;
    };
    break;
  }
  case 2:
  {
    const Matrix* gj_1; const Matrix* gj_2;
    for (int j = 0; j < Nv; j++) 
    {
      g_diff_nodal.setZero();
      for (int i = 0; i < ncell; i++)
      {
        wind = modal.E(0, i) >= 0 ? -1 : 1;

        if (j-wind >= 0 && j-wind < Nv) { gj_1 = &(g_nodal[j-wind]); }
        else { gj_1 = &zero_nodal_mat_; };
        if (j-2*wind >= 0 && j-2*wind < Nv) { gj_2 = &(g_nodal[j-2*wind]); }
        else { gj_2 = &zero_nodal_mat_; };

        g_diff_nodal.col(i) = (3.e0 * g_nodal[j].col(i) 
                              - 4.e0 * gj_1->col(i)
                              + gj_2->col(i)) * hvinv2 * real_t(wind);
      }
      g_diff_nodal = - 1.e0 * g_diff_nodal.array() * E_nodal.array();
      fe_->Assemble_F(g_diff_nodal, 0, &((*ch)[j]));
      (*ch)[j] *= JacobiDet1;
    };
    break;
  }
  case 3:
  {
    std::vector<Matrix> flux_upwind_nodal, flux_downwind_nodal;
    WENO3_VelocityReconstruct(modal.g, g_nodal, &flux_upwind_nodal, &flux_downwind_nodal);
    for (int j = 0; j < Nv; j++) 
    {
      g_diff_nodal.setZero();
      for (int i = 0; i < ncell; i++)
      { 
        for (int q = 0; q < numqua; q++)
        {
          wind = E_nodal(q, i) >= 0 ? -1 : 1;
          switch (wind)
          {
            case -1:
              g_diff_nodal(q, i) = (flux_downwind_nodal[j+1](q, i) - 
                                    flux_downwind_nodal[j](q, i)) * hvinv;
              
              break;
            case 1:
              g_diff_nodal(q, i) = (flux_upwind_nodal[j+1](q, i) - 
                                    flux_upwind_nodal[j](q, i)) * hvinv;
              break;
          };
        };
      }
      g_diff_nodal = - 1.e0 * g_diff_nodal.array() * E_nodal.array();
      fe_->Assemble_F(g_diff_nodal, 0, &((*ch)[j]));
      (*ch)[j] *= JacobiDet1;
    };
    break;
  }
  default:
    QUEST_ERROR(" high order velocity discretization is not implemented yet ! ");
    break;
  }
  
};

void KineticDriftD_DG_IMEX_IM_Schur::WENO3_VelocityReconstruct(
                                      const std::vector<Matrix>& g,
                                      const std::vector<Matrix>& g_nodal,
                                      std::vector<Matrix>* flux_upwind_nodal,
                                      std::vector<Matrix>* flux_downwind_nodal)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const real_t& hv = mesh1D_->gethv();
  const int& numqua = fe_->getnumqua();
  // *** 
  flux_upwind_nodal->resize(Nv + 1);
  flux_downwind_nodal->resize(Nv + 1);
  flux_upwind_nodal->at(0) = zero_nodal_mat_;
  flux_downwind_nodal->at(0) = zero_nodal_mat_;
  flux_upwind_nodal->at(Nv) = zero_nodal_mat_;
  flux_downwind_nodal->at(Nv) = zero_nodal_mat_;

  Vector beta0, beta1;
  const real_t d0 = 2.e0 / 3.e0;
  const real_t d1 = 1.e0 / 3.e0;
  Vector alpha0, alpha1, alpha_sum;
  const real_t eps_weno = 1.e-6;
  const Matrix* g2;
  const Matrix* g1;
  const Matrix* g0;
  for (int j = 1; j < Nv; j++)
  {
    flux_upwind_nodal->at(j) = zero_nodal_mat_;
    flux_downwind_nodal->at(j) = zero_nodal_mat_;
    if (j-2 >= 0 && j-2 < Nv) { g2 = &(g_nodal[j-2]); }
    else { g2 = &zero_nodal_mat_; };
    if (j-1 >= 0 && j-1 < Nv) { g1 = &(g_nodal[j-1]); }
    else { g1 = &zero_nodal_mat_; };
    g0 = &(g_nodal[j]);
    Matrix& flux_up   = flux_upwind_nodal->at(j);
    Matrix& flux_down = flux_downwind_nodal->at(j-1);
    for(int i = 0; i < ncell; i++)
    {
      Eigen::Ref<const Vector> g2i = g2->col(i);
      Eigen::Ref<const Vector> g1i = g1->col(i);
      Eigen::Ref<const Vector> g0i = g0->col(i);
      beta0 = (g2i - g1i).array().square();
      beta1 = (g1i - g0i).array().square();
      alpha0 = (d0 / (eps_weno + beta0.array()).square()).matrix();
      alpha1 = (d1 / (eps_weno + beta1.array()).square()).matrix();
      alpha_sum = alpha0 + alpha1;
      alpha0 = alpha0.array() / alpha_sum.array();
      alpha1 = alpha1.array() / alpha_sum.array();
      if (j == 1)
      {
        flux_up.col(i) = alpha0.array() * (0.5e0 * g1i.array() + 0.5e0 * g0i.array())
                        + alpha1.array() * (- 0.5e0 * g2i.array() + 1.5e0 * g1i.array());
      } else
      {
        flux_up.col(i) = alpha0.array() * (0.5e0 * g1i.array() + 0.5e0 * g0i.array())
                        + alpha1.array() * (- 0.5e0 * g2i.array() + 1.5e0 * g1i.array());
        flux_down.col(i) = alpha0.array() * (1.5e0 * g1i.array() - 0.5e0 * g0i.array())
                          + alpha1.array() * (0.5e0 * g2i.array() + 0.5e0 * g1i.array());
      };  
    };
  }  
};

void KineticDriftD_DG_IMEX_IM_Schur::sourceh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          Matrix* source_sumh,
                          std::vector<Matrix>* sourceh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& numqua = fe_->getnumqua(); 
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const real_t& hv = mesh1D_->gethv();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 
  std::vector<Matrix> source_nodal;
  source_nodal.resize(Nv);
  sourceh->resize(Nv);
  Matrix source_temp_nodal;
  source_temp_nodal = zero_nodal_mat_;
  
  
  for (int j = 0; j < Nv; j++)
  {
    real_t vj = V(j);
    fe_->Interpolate_Final(
      [this, vj](const Matrix& x, const real_t& t) { return this->fsource(x, vj, t); },
      Trun,
      &(source_nodal.at(j)));
    source_temp_nodal += vweights(j) * source_nodal[j];
  };
  fe_->Assemble_F(source_temp_nodal, 0, source_sumh);
  *source_sumh *= JacobiDet1;
  // std::cout << " source_sumh = " << *source_sumh << std::endl;
  // PAUSE();
  #pragma omp parallel num_threads(NTH_), default(shared)
  {
    real_t vj, Mj;
  #pragma omp for schedule(static)
    for (int j = 0; j < Nv; j++)
    {
      vj = V(j);
      Mj = Maxwell(vj);
      source_nodal[j] = source_nodal[j] - source_temp_nodal * Mj;
      fe_->Assemble_F(source_nodal[j], 0, &((*sourceh)[j]));
      (*sourceh)[j] *= JacobiDet1;
    };
  };
  
};

void KineticDriftD_DG_IMEX_IM_Schur::updateAll(const real_t& Trun, 
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
  Matrix ah_bc = Matrix::Zero(1,extboundarynum);
  Matrix dh_bc = Matrix::Zero(1,extboundarynum);
  Matrix rhoboundary_flux, vgboundary_flux;
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
    kinetic_modal_stages_[s].E = kinetic_modal_.E;
    kinetic_modal_stages_[s].phi = kinetic_modal_.phi;
    
    rhoboundary_flux.setZero();
    vgboundary_flux.setZero();
    const std::vector<Matrix>* gg_ptr = nullptr;
    if (s == 0) { gg_ptr = &kinetic_modal_.g; } 
    else { gg_ptr = &kinetic_modal_stages_[s-1].g; };

    sourceh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(source_sumh_[s]), &(sourceh_[s]));
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - (dt * Ai(s, i)) * ah_[i]
                                    + (dt * Ae(s, i)) * source_sumh_[i];
      
      #pragma omp parallel num_threads(NTH_), default(shared)
      { 
      real_t vj, Mj;
      #pragma omp for schedule(static)
        for (int j = 0; j < Nv; j++)
        {
          vj = V(j);
          Mj = Maxwell(vj);
          kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                        - (dt * Ae(s, i) * eps_) * bh_[i][j]
                                        - (dt * Ae(s, i) * eps_) * ch_[i][j]
                                        + (dt * Ai(s, i) * vj * Mj) * dh_[i]
                                        - (dt * Ai(s, i) * vj * Mj) / theta_ * mh_[i]
                                        - (dt * sigmas_ * Ai(s, i)) * sh_[i][j]
                                        + (dt * Ae(s, i) * eps2_) * sourceh_[i][j];
        };
      };
    };
    // kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
    //                               + (dt * Ai(s, s)) * source_sumh_[s];
    // #pragma omp parallel num_threads(NTH_), default(shared)
    // {
    // real_t vj;
    // #pragma omp for schedule(static)
    //   for (int j = 0; j < Nv; j++)
    //   {
    //     vj = V(j);
    //     kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
    //                                   + (dt * Ai(s, s) * eps2_) * sourceh_[s][j];
    //   };
    // };
    if (s == 0)
    {
      dh_bc_compute(kinetic_modal_stages_[s].rho, 
                    kinetic_modal_.g,
                    Trun + ci(s) * dt, 
                    &rhoboundary_flux);
    } else 
    {
      dh_bc_compute(kinetic_modal_stages_[s].rho, 
                    kinetic_modal_stages_[s-1].g,
                    Trun + ci(s) * dt, 
                    &rhoboundary_flux);
    }
    ah_bc_compute(kinetic_modal_stages_[s].rho,
                  kinetic_modal_stages_[s].g,
                  Trun + ci(s) * dt, &vgboundary_flux);

    SolveRho_DichletBoundary(Trun + ci(s) * dt,  dt,  Ai(s, s),
                            rhoboundary_flux, 
                            vgboundary_flux,
                            kinetic_modal_.rho,
                            &(kinetic_modal_stages_[s]));
    
    
    dh_compute(kinetic_modal_stages_[s].rho,
              kinetic_modal_stages_[s].g,
              rhoboundary_flux, &(dh_[s]));
    
    mh_compute(kinetic_modal_stages_[s].rho,
                kinetic_modal_stages_[s].g, 
                kinetic_modal_stages_[s].E, 
                &(mh_[s]));
    #pragma omp parallel num_threads(NTH_), default(shared)
    { 
    real_t vj, Mj;
    #pragma omp for schedule(static)
      for (int j = 0; j < Nv; j++)
      {
        vj = V(j);
        Mj = Maxwell(vj);
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j] 
                                      + (dt * Ai(s, s) * vj * Mj) * dh_[s]
                                      - (dt * Ai(s, s) * vj * Mj) / theta_ * mh_[s];
      };
    };
    
    SolveG(Trun, dt, Ai(s, s), &(kinetic_modal_stages_[s]));
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));

    ch_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(ch_[s]));
    
    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));
    
    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
  };
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
};

void KineticDriftD_DG_IMEX_IM_Schur::SolveRho_DichletBoundary(const real_t& Trun,
                                        const real_t& dt,
                                        const real_t& a,
                                        const Matrix& dh_bc,
                                        const Matrix& ah_bc,
                                        const Matrix& rho_prev,
                                        model_data_* modal)
{
  // ********* 传入相关变量 ********** //
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
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
  real_t parainv = 1.e0 / (eps2_ + a * sigmas_ * dt);
  Matrix b0_bc = modal->rho; //宏观方程入流边界条件的处理
  b0_bc = b0_bc + (dt * a) * ah_bc;

  Matrix stag;
  Vector tmp;
  Vector b0_hat = Vector::Zero(NTdofs);
  real_t vj, Mj;
  SparseMatrix Coef = a * dt * Minv_ * parainv;
  for (int j = 0; j < Nv; j++)
  {
    vj = V(j);
    Mj = Maxwell(V(j));
    stag = modal->g[j] + (dt * a * vj * Mj) * dh_bc;  // 微观方程中边界条件的处理
    Eigen::Map<const Vector> gj(stag.data(), NTdofs);
    tmp = Coef * (vweights(j) * gj);
    b0_hat += vj * (Da_ * tmp) + (Da_ext_[j] * tmp);
  };
  
  Eigen::Map<Vector> rho_old_vector(b0_bc.data(), NTdofs);
  Vector b0 = rho_old_vector - b0_hat; // Schur补的右端项
  real_t ppp = dt * a * dt * a * parainv;
  H_ = M_ + (a * dt) * Mbc_ + ppp * temp_;

  H_ += 0.5e0 / gamma_ * M_;
 
  int iter = 0;
  real_t rho_error = 10;
  Matrix phi_Dirichlet(1, 2);
  Matrix poi_lhs, rhotemp_nodal;
  Matrix rho_prev_nodal;
  Eigen::Map<const Vector> rho_prev_vec(rho_prev.data(), NTdofs);
  fe_->modal_to_nodal1D(rho_prev, &rho_prev_nodal);
  Vector rhotemp;
  phi_Dirichlet(0, 0) = phi_bc(x1, Trun);
  phi_Dirichlet(0, 1) = phi_bc(x2, Trun);
  while (rho_error > iter_tol_)
  {
    Eh_compute(modal->E, &ME_);
    Coef = Minv_ * ME_;
    H_ -= ppp * Da_ * Coef;
    for (int j = 0; j < Nv; j++)
    { 
      vj = V(j);
      Mj = Maxwell(vj);
      H_ -= (ppp * vweights(j) * Mj / theta_) * Da_ext_[j] * vj * Coef;
    }
    b0 += 0.5e0 / gamma_ * M_ * rho_prev_vec;
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
      gmres_.set_restart(300);
      gmres_.setMaxIterations(NTdofs);
      gmres_.compute(H_);
      // std::cout << "  iterations = " << gmres_.iterations() << std::endl;
      // std::cout << "  The end of Solving ! " << std::endl;
      break;

    default:
      QUEST_ERROR(" The Schur complement Solver is not been implemented ! ");
      break;
    }
    
    switch (schur_solver_type_)
    {
    case Solver1DType::LU:
      rhotemp = lu_.solve(b0);
      break;
    case Solver1DType::GMRES:
      rhotemp = gmres_.solve(b0);
      std::cout << "  error = " << gmres_.error()      << std::endl;
      if (gmres_.info() != Eigen::Success) {
        real_t rel_error = gmres_.error();
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
    Eigen::Map<Matrix> rhotemp_modal(rhotemp.data(), polydim, ncell);
    fe_->modal_to_nodal1D(rhotemp_modal, &rhotemp_nodal);
    fe_->computerrorL2(rhotemp_nodal, rho_prev_nodal, &rho_error); // 计算相对残差
    rho_prev_nodal = rhotemp_nodal;

    poi_lhs = (rho_d_nodal_ - rhotemp_nodal) / gamma_;  // 泊松方程右端项
    // 求解泊松方程并且进入下一次循环
    poisson_solver_->solveall(phi_Dirichlet, poi_lhs, &(modal->E), &(modal->phi)); 
    iter++;
    std::cout << "Schur iteration: " << iter
              << ", rho_error = " << rho_error << std::endl;
  }
  Eigen::Map<Matrix> rhotemp_modal(rhotemp.data(), polydim, ncell);
  modal->rho = rhotemp_modal;
};

void KineticDriftD_DG_IMEX_IM_Schur::SolveG(const real_t& Trun,
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
  DiagnalMatrix Mtemp = parainv * Minv_ref;
  #pragma omp parallel num_threads(NTH_), default(shared)
  {
  #pragma omp for schedule(static)
    for (int j = 0; j < Nv; j++)
    {
      modal->g[j] = Mtemp * modal->g[j];
    }
  }
};

real_t KineticDriftD_DG_IMEX_IM_Schur::rho_init(const real_t& x) 
{ 
  return 1.e0;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur::rho_init(const Matrix& x) 
{
  Matrix rho = x;
  rho.setConstant(1.e0);
  return rho;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::rho_d(const real_t& x) 
{
  return 1.e0 - (1.e0 - 0.001)/2.e0 * (std::tanh((x - 0.3e0)/0.02e0) - std::tanh((x - 0.7e0)/0.02e0));
};

Matrix KineticDriftD_DG_IMEX_IM_Schur::rho_d(const Matrix& x) 
{
  Matrix x1 = (x.array() - 0.3) / 0.02;
  Matrix x2 = (x.array() - 0.7) / 0.02;
  Matrix rhod = 1.0 - (1.0 - 0.001) / 2.0 * (x1.array().tanh() - x2.array().tanh());
  return rhod;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::source(const real_t& x, const real_t& t) 
{
  return 0.e0;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur::source(const Matrix& x, const real_t& t) 
{
  return Matrix::Zero(x.rows(), x.cols());
};

real_t KineticDriftD_DG_IMEX_IM_Schur::fsource(const real_t& x, const real_t& v, const real_t& t) 
{
  return 0.e0;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur::fsource(const Matrix& x, const real_t& v, const real_t& t) 
{
  return Matrix::Zero(x.rows(), x.cols());
};

real_t KineticDriftD_DG_IMEX_IM_Schur::f_init(const real_t& x, const real_t& v) 
{
  return Maxwell(v);
};

Matrix KineticDriftD_DG_IMEX_IM_Schur::f_init(const Matrix& x, const real_t& v) 
{
  Matrix f = x;
  f.setConstant(1.e0);
  f = f * Maxwell(v);
  return f;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::g_init(const real_t& x, const real_t& v) 
{
  return 0.e0;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur::g_init(const Matrix& x, const real_t& v) 
{
  Matrix g = x;
  g.setZero();
  return g;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::fL_bc(const int& j,  const real_t& t,
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

real_t KineticDriftD_DG_IMEX_IM_Schur::fR_bc(const int& j,  const real_t& t,
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
    f = 1.e0;
  } else {
    Vector temp = modal.rho.col(ncell - 1) + eps_ * modal.g[j].col(ncell - 1);
    f = temp.dot(boundary_u[0]);
  }
  return f;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::fL_bc(const int& j,  const real_t& t,
      const Matrix& rho, const std::vector<Matrix>& g) 
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
    Vector temp = rho.col(0) + eps_ * g[j].col(0);
    f = temp.dot(boundary_u[1]);
  }
  return f;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::fR_bc(const int& j,  const real_t& t,
      const Matrix& rho, const std::vector<Matrix>& g) 
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
    f = 1.e0;
  } else {
    Vector temp = rho.col(ncell - 1) + eps_ * g[j].col(ncell - 1);
    f = temp.dot(boundary_u[0]);
  }
  return f;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::Maxwell(const real_t& v) {
  real_t M = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) * std::exp( - v * v / ( 2 * theta_));
  M = M / Maxwell_sum_;
  return M;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::fL_explicit_bc(const int& j, const real_t& t)
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

real_t KineticDriftD_DG_IMEX_IM_Schur::fR_explicit_bc(const int& j, const real_t& t)
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

real_t KineticDriftD_DG_IMEX_IM_Schur::gL_bc(const int& j,  const real_t& t,
      const model_data_& modal, const real_t& rho_L) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t g;
  g = (fL_bc(j, t, modal) - rho_L) / eps_;
  return g;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::gR_bc(const int& j,  const real_t& t,
      const model_data_& modal, const real_t& rho_R) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t g;
  g = (fR_bc(j, t, modal) - rho_R) / eps_;
  return g;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::phi_bc(const real_t& x, const real_t& t)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t phi;
  if (std::abs(x - x1) < 1.e-10) {
    return 0.e0;
  } else if (std::abs(x - x2) < 1.e-10) {
    return 5.e0;
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return phi;
};

real_t KineticDriftD_DG_IMEX_IM_Schur::rho_numericalbc(
      const real_t& x, const real_t& t,
      const Matrix& rho, const std::vector<Matrix>& g) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  // *** 
  real_t rhobc = 0.e0;
  if (std::abs(x - x1) < 1.e-11) {
    for (int j = 0; j < Nv; j++) {
      rhobc += Vweights(j) * fL_bc(j, t, rho, g);
    }
  } else if (std::abs(x - x2) < 1.e-11) {
    for (int j = 0; j < Nv; j++) {
      rhobc += Vweights(j) * fR_bc(j, t, rho, g);
    }
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return rhobc;
};


KineticDriftD_DG_IMEX_IM_Schur_period::KineticDriftD_DG_IMEX_IM_Schur_period(
                  const KineticTensorMesh1D* mesh1D,
                  const fespace1D* fe,
                  const IMEX_RK* rk_table,
                  const PoissonSolver1D_period* poisson_solver,
                  const Solver1DType& schur_solver_type)
  : KineticDriftD_DG_IMEX_IM_Schur(mesh1D, fe, rk_table, poisson_solver, schur_solver_type) {};

void KineticDriftD_DG_IMEX_IM_Schur_period::init()
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& numqua = fe_->getnumqua(); 
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweight = mesh1D_->getvweights();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  const int& NTdofs = fe_->getNTdofs();
  // *** 
  QUEST_VERIFY(fe_->getmesh1D()->IsPeriodBoundary(), "The mesh is not periodical !");
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();
  setMaxwell();
  std::cout << " Maxwell sum = " << Maxwell_sum_ << std::endl;
  // QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  // QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");
  ah_.resize(stages_);
  ch_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_); 
  mh_.resize(stages_);
  source_sumh_.resize(stages_);
  sourceh_.resize(stages_);

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_d(x); }, &rho_d_modal_);
  fe_->Interpolate_Initial(
      [this](const Matrix& x) { return this->rho_d(x); }, &rho_d_nodal_);
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  };
  kinetic_modal_stages_.resize(stages_);
  zero_modal_mat_ = Matrix::Zero(polydim, ncell);
  zero_nodal_mat_ = Matrix::Zero(numqua, ncell);

  D_compute(beta1_, &Da_);
  Da_ = - Da_;
  // Da_ext_treat(0.e0, 0.e0, Da_, &Da_ext_);
  D_compute(- beta1_, &Db_); // 交替通量
  // Db_ext_treat(0.e0, 0.e0, Db_, &Db_ext_);

  M_compute(&M_);
  // Mbc_compute(&Mbc_);
  Minv_compute(&Minv_);
  temp_.resize(NTdofs, NTdofs);
  // real_t Mj;
  // for (int j = 0; j < Nv; j++)
  // {
  //   Mj = Maxwell(V(j));
  //   temp_ += vweight(j) * Mj * (Da_ * V(j))  * (Minv_ * (Db_ * V(j)));
  // };
  temp_ = Da_ * Minv_ * Db_;
  temp_ = temp_ * theta_;

  Matrix phi_Dirichlet(1, 2);
  phi_Dirichlet.setZero();
  Matrix rho_nodal;
  fe_->modal_to_nodal1D(kinetic_modal_.rho, &rho_nodal);
  Matrix poi_lhs = (rho_d_nodal_ - rho_nodal) / gamma_;  // 泊松方程右端项
  poisson_solver_->solveall(phi_Dirichlet, poi_lhs, &(kinetic_modal_.E), &(kinetic_modal_.phi)); 
};

void KineticDriftD_DG_IMEX_IM_Schur_period::setdt(real_t* dt)
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
  *dt = gamma_ * (*dt);
};

void KineticDriftD_DG_IMEX_IM_Schur_period::D_compute(const real_t& be, SparseMatrix* D)
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

void KineticDriftD_DG_IMEX_IM_Schur_period::Da_ext_treat(const real_t& Trun,
                    const real_t& dt, 
                    const SparseMatrix& Da,
                    std::vector<SparseMatrix>* Da_ext) 
{
  QUEST_ERROR("KineticDriftD_DG_IMEX_IM_Schur_period do not need Da_ext_treat() !");
};

void KineticDriftD_DG_IMEX_IM_Schur_period::Db_ext_treat(const real_t& Trun,
                    const real_t& dt, 
                    const SparseMatrix& Db,
                    SparseMatrix* Db_ext) 
{
  QUEST_ERROR("KineticDriftD_DG_IMEX_IM_Schur_period do not need Db_ext_treat() !");
};

void KineticDriftD_DG_IMEX_IM_Schur_period::ah_compute(const model_data_& modal, 
        const real_t& Trun,
        const Matrix& ah_bc,
        Matrix* ah)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  
  real_t vj;
  Matrix vg_modal = zero_modal_mat_;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    vg_modal += vweights(j) * V(j) * modal.g[j];
  };

  ah->resize(polydim, ncell);
  ah->setZero();
  Eigen::Map<Vector> ah_vec(ah->data(), NTdofs);
  Eigen::Map<const Vector> vgj(vg_modal.data(), NTdofs);
  ah_vec = Da_ * vgj;
};

void KineticDriftD_DG_IMEX_IM_Schur_period::ah_bc_compute(const Matrix& rho, 
                          const std::vector<Matrix>& g,
                          const real_t& Trun, 
                          Matrix* ah_bc)
{
  QUEST_ERROR("KineticDriftD_DG_IMEX_IM_Schur_period do not need ah_bc_compute() !");
};

void KineticDriftD_DG_IMEX_IM_Schur_period::dh_compute(const Matrix& rho, 
                                                      const std::vector<Matrix>& g,
                                                      const Matrix& dh_bc,
                                                      Matrix* dh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 

  dh->resize(polydim, ncell);
  dh->setZero();
  Eigen::Map<Vector> dh_vec(dh->data(), NTdofs);

  Eigen::Map<const Vector> rho_vec(rho.data(), NTdofs);
  dh_vec = Db_ * rho_vec;
};

void KineticDriftD_DG_IMEX_IM_Schur_period::dh_bc_compute(const Matrix& rho, 
                  const std::vector<Matrix>& g,
                  const real_t& Trun,
                  Matrix* dh_bc) 
{
  QUEST_ERROR("KineticDriftD_DG_IMEX_IM_Schur_period do not need dh_bc_compute() !");
};

void KineticDriftD_DG_IMEX_IM_Schur_period::Mbc_compute(SparseMatrix* Mbc) 
{
  QUEST_ERROR("KineticDriftD_DG_IMEX_IM_Schur_period do not need dh_bc_compute() !");
};
    
void KineticDriftD_DG_IMEX_IM_Schur_period::mh_compute(const Matrix& rho, 
                  const std::vector<Matrix>& g,
                  const Matrix& E,
                  Matrix* mh) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const Vector& V = mesh1D_->getV();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& NTdofs = fe_->getNTdofs();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 

  mh->resize(polydim, ncell);
  mh->setZero();
  Eigen::Map<Vector> mh_vec(mh->data(), NTdofs);

  Eigen::Map<const Vector> rho_vec(rho.data(), NTdofs);
  mh_vec = ME_ * rho_vec;
  // std::cout << " mh = " << *mh << std::endl;

  // Matrix E_nodal;
  // fe_->modal_to_nodal1D(E, &E_nodal);
  
  // Matrix rho_nodal;
  // fe_->modal_to_nodal1D(rho, &rho_nodal);
  // E_nodal = E_nodal.array() * rho_nodal.array();
  
  // Matrix mmmm;
  // fe_->Assemble_F(E_nodal, 0, &mmmm);
  // mmmm *= JacobiDet1;
  // mmmm = mmmm - *mh;
  // std::cout << " distance = " << mmmm.norm() << std::endl;
  // PAUSE();
};

void KineticDriftD_DG_IMEX_IM_Schur_period::bh_compute(const model_data_& modal, 
                  const real_t& Trun,
                  const std::vector<Matrix>& boundary_flux,
                  std::vector<Matrix>* bh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  real_t JacobiDet1 = JacobiDet(0);
  // *** 
  std::vector<Matrix> g_nodal;
  fe_->modal_to_nodal1D(modal.g, &g_nodal);
  Matrix rhs;
  Matrix flux_int, flux_ext;
  Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  Matrix dh = Matrix::Zero(polydim, ncell);
  Matrix temp;
  bh->resize(Nv);
  for (int j = 0; j < Nv; j++) {
    temp = V(j) * g_nodal[j];
    fe_->Assemble_F(temp, 1, &rhs);
    fluxint_upwind_compute(modal.g[j], V(j), &flux_int);
    // fluxext_upwind_compute(modal.g[j], V(j), boundary_flux[j], &flux_ext);
    fe_->Assemble_Flux(flux_int, boundary_flux[j], &flux_rhs);
    (*bh)[j] = (- rhs + flux_rhs);
    dh = dh + vweights(j) * (*bh)[j];  // <vg_x>
  };
  // (I - \Pi)(vg_x)
  // \Pi(vg_x) = <vg_x> M
  real_t Mj;
  for (int j = 0; j < Nv; j++)
  {
    Mj = Maxwell(V(j));
    bh->at(j) = (bh->at(j) - dh * Mj) * JacobiDet1;
  };
  // Matrix ttt = Matrix::Zero(polydim, ncell);
  // for (int j = 0; j < Nv; j++)
  // {
  //   ttt += bh->at(j) * vweights(j);
  // };
  // std::cout << " ttt = " << ttt << std::endl;
  // PAUSE();
  
};

void KineticDriftD_DG_IMEX_IM_Schur_period::bh_extflux_compute(const model_data_& modal, 
                      const real_t& Trun,
                      std::vector<Matrix>* boundary_flux) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  const std::vector<Eigen::Vector2i>& ExtBNei_period =  mesh1D_->getextboundaryneighbors_period();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  boundary_flux->resize(Nv);
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal = ExtBNormal[0];
  real_t v;
  Cellindex0 = ExtBNei_period[0](0);
  Cellindex1 = ExtBNei_period[0](1);
  for (int j = 0; j < Nv; j++) 
  { 
    v = V(j);
    boundary_flux->at(j).setZero();
    boundary_flux->at(j).resize(1, extboundarynum);
    val0 = v * modal.g[j].col(Cellindex0).dot(boundary_u[1]);
    val1 = v * modal.g[j].col(Cellindex1).dot(boundary_u[0]);
    boundary_flux->at(j)(0, 0) = (v * normal > 0) ? val0 : val1;
    // std::cout << " boundary_flux->at(j)(0, 0) = " << boundary_flux->at(j)(0, 0) << std::endl;
    boundary_flux->at(j)(0, 1) = boundary_flux->at(j)(0, 0);
  };
};

void KineticDriftD_DG_IMEX_IM_Schur_period::FourierMatrix(const real_t& xita, std::string& OutfilePath)
{
  cMatrix D_NegativeWind, D_PositiveWind, U;
  FourierDMatrix_compute(0.5e0, xita, &D_PositiveWind);
  D_PositiveWind = - D_PositiveWind;
  FourierDMatrix_compute(- 0.5e0, xita, &D_NegativeWind);
  D_NegativeWind = - D_NegativeWind;
  FourierUMatrix_compute(D_NegativeWind, D_PositiveWind, xita, &U);
};

void KineticDriftD_DG_IMEX_IM_Schur_period::FourierDMatrix_compute(const real_t& be, const real_t& xita, cMatrix* D) 
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
  // D->resize(NTdofs, NTdofs);
  // int estimatedNonZeros = ncell * polydim * polydim + 
  //         (intboundaryNum + 1) * polydim * polydim * 4;
  // std::vector<Eigen::Triplet<creal_t>> tripletList_D;
  // tripletList_D.reserve(estimatedNonZeros);
  
  // Matrix Bref;
  // Bref = dv_u;
  // for (int i = 0; i < ncell; i++) {
  //   int alpha_start = i * polydim;
  //   for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
  //     for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
  //       real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jx(i);
  //       tripletList_D.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
  //                       alpha_start + trial_basis_index, creal_t(ini_value, 0.e0)));
  //     };
  //   };
  // };

  // int test_cell_Index, trial_cell_Index;
  // Vector test_qua_value, trial_qua_value;
  // real_t test_normal, trial_normal;
  // int alpha, beta;
  // real_t ini_value_temp;
  // creal_t ini_value;
  // for (int i = 0; i < intboundaryNum; i++) {
  //   for (int test_cell = 0; test_cell < 2; test_cell++) {
  //     test_cell_Index = intNei[i](test_cell);
  //     test_qua_value = boundary_u[test_cell];
  //     test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
  //     for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
  //       trial_cell_Index = intNei[i](trial_cell);
  //       trial_qua_value = boundary_u[trial_cell];
  //       trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);
  //       for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
  //         alpha = Tm(test_basis_index, test_cell_Index);
  //         for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
  //           beta = Tm(trial_basis_index, trial_cell_Index);
  //           ini_value_temp =  (0.5e0 * trial_qua_value(trial_basis_index) 
  //               + trial_qua_value(trial_basis_index) * be * trial_normal) 
  //               * test_qua_value(test_basis_index) * test_normal;
  //           ini_value = - ini_value_temp * std::exp(creal_t(0.e0, real_t(trial_cell_Index - test_cell_Index) * xita));   
  //           // 这是一维代码不需要积分
  //           tripletList_D.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
  //         };
  //       };
  //     };
  //   };
  // };

  // for (int test_cell = 0; test_cell < 2; test_cell++) {
  //   test_cell_Index = extNei_period[0](test_cell);
  //   test_qua_value = boundary_u[1 - test_cell];
  //   test_normal = extnormals[0] * std::pow(-1.e0,test_cell);
  //   for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
  //     trial_cell_Index = extNei_period[0](trial_cell);
  //     trial_qua_value = boundary_u[1 - trial_cell];
  //     trial_normal = extnormals[0] * std::pow(-1.e0,trial_cell);
  //     for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
  //       alpha = Tm(test_basis_index, test_cell_Index);
  //       for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
  //         beta = Tm(trial_basis_index, trial_cell_Index);
  //         ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
  //             + trial_qua_value(trial_basis_index) * be * trial_normal) 
  //             * test_qua_value(test_basis_index) * test_normal;
  //         ini_value = - ini_value_temp * std::exp(creal_t(0.e0, real_t(trial_cell_Index - test_cell_Index) * xita));   
  //         tripletList_D.push_back(Eigen::Triplet<creal_t>(alpha, beta, ini_value));
  //       };
  //     };
  //   };
  // };
  
  // D->setFromTriplets(tripletList_D.begin(), tripletList_D.end());

  // Matrix Bref;
  // Bref = dv_u;
  
  // int test_cell_Index, trial_cell_Index;
  // Vector test_qua_value, trial_qua_value;
  // real_t test_normal, trial_normal;
  // int alpha, beta;
  // real_t ini_value_temp;
  // creal_t ini_value;
  // for (int test_cell = 0; test_cell < 2; test_cell++) {
  //   test_cell_Index = test_cell;
  //   test_qua_value = boundary_u[1 - test_cell];
  //   test_normal = extnormals[0] * std::pow(- 1.e0, test_cell);
  //   for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
  //     trial_cell_Index = trial_cell;
  //     trial_qua_value = boundary_u[1 - trial_cell];
  //     trial_normal = extnormals[0] * std::pow(- 1.e0,trial_cell);
  //     for (int test_basis_index = 0; test_basis_index < polydim; test_basis_index++) {
  //       alpha = test_cell_Index;
  //       for (int trial_basis_index = 0; trial_basis_index < polydim; trial_basis_index++) {
  //         beta = trial_cell_Index;
  //         ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
  //             + trial_qua_value(trial_basis_index) * be * trial_normal) 
  //             * test_qua_value(test_basis_index) * test_normal;
  //         ini_value = - ini_value_temp * std::exp(creal_t(0.e0, real_t(trial_cell_Index - test_cell_Index) * xita));   
  //         tripletList_D.push_back(Eigen::Triplet<creal_t>(alpha, beta, ini_value));
  //       };
  //     };
  //   };
  // };

};

void KineticDriftD_DG_IMEX_IM_Schur_period::FourierUMatrix_compute(
                              const cMatrix& D_NegativeWind,
                              const cMatrix& D_PositiveWind,
                              const real_t& xita, 
                              cMatrix* U)
{
  // ********* 传入相关变量 ********** //
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
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
  // U->resize(Nv * NTdofs, Nv * NTdofs);
  // int estimatedNonZeros = ncell * polydim * polydim + 
  //         (intboundaryNum + 1) * polydim * polydim * 4;
  // std::vector<Eigen::Triplet<creal_t>> tripletList_U;
  // tripletList_U.reserve(estimatedNonZeros * Nv * Nv);

  // real_t vtest, vtrial;
  // real_t Mtest, Mtrial;
  // creal_t ini_value;
  // cSparseMatrix* D_point;
  // for (int test_val = 0; test_val < Nv; test_val++)
  // {
  //   for (int trial_val = 0; trial_val < Nv; trial_val++)
  //   {
  //     vtest = V(test_val); 
  //     vtrial = V(trial_val);
  //     Mtest = Maxwell(vtest);
  //     Mtrial = Maxwell(vtrial);
  //     if (vtrial > 0)
  //     {
  //       D_point = &D_PositiveWind;
  //     } else if (vtrial < 0)
  //     {
  //       D_point = &D_NegativeWind;
  //     }
  //     if (test_val == trial_val)
  //     {
  //       for (int k = 0; k < D_point->outerSize(); ++k) 
  //       {
  //         for (cSparseMatrix::InnerIterator it(&D_point, k); it; ++it) {
  //           alpha = test_val * Nv + it.row();
  //           beta = trial_val * Nv + it.col();
  //           ini_value = (vtest - vweights(trial_val) * vtrial * Mtest) * it.value();
  //           tripletList_U.push_back(Eigen::Triplet<creal_t>(alpha, beta, ini_value));
  //         }
  //       }
  //     }
  //     for (int k = 0; k < D_point->outerSize(); ++k) 
  //     {
  //       for (cSparseMatrix::InnerIterator it(&D_point, k); it; ++it) {
  //         alpha = test_val * Nv + it.row();
  //         beta = trial_val * Nv + it.col();
  //         ini_value = - vweights(trial_val) * vtrial * Mtest * it.value()
  //         tripletList_U.push_back(Eigen::Triplet<creal_t>(alpha, beta, ini_value));
  //       }
  //     }
  //   }
  // }
  // U->setFromTriplets(tripletList_U.begin(), tripletList_U.end());
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::rho_init(const Matrix& x)
{
  Matrix rho = rho_d(x).array() + gamma_ * x.array().sin();
  return rho;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::rho_init(const real_t& x)
{
  return rho_d(x) + gamma_ * std::sin(x);
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::f_init(const Matrix& x, const real_t& v)
{
  Matrix f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::f_init(const real_t& x, const real_t& v)
{
  real_t f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::g_init(const Matrix& x, const real_t& v)
{
  real_t M = Maxwell(v);
  Matrix g = - v * M * (x.array().cos() + 1.e0 / theta_ * x.array().cos() * x.array().sin());
  g = g.array() * gamma_ - v * M * x.array().cos() / theta_;
  return g;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::g_init(const real_t& x, const real_t& v)
{
  real_t M = Maxwell(v);
  real_t g = - v * M * (std::cos(x) + 1.e0 / theta_ * std::cos(x) * std::sin(x));
  g = g * gamma_ - v * M * std::cos(x) / theta_;
  return g;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::source(const Matrix& x, const real_t& t)
{
  Matrix S = 2.e0 * x;
  S = - std::exp(- 2.e0 * theta_ * t) * S.array().cos() * gamma_;
  S = S.array() + std::exp(- theta_ * t) * x.array().sin();
  return S;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::source(const real_t& x, const real_t& t)
{
  real_t S = - std::exp(- 2.e0 * theta_  * t) * std::cos(2.e0 * x);
  S = S * gamma_;
  S = S + std::exp(- theta_ * t) * std::sin(x);
  return S;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::fsource(const Matrix& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  Matrix S;
  Matrix sinx = x.array().sin();
  Matrix cosx = x.array().cos();
  Matrix sin2x = (2.e0 * x).array().sin();
  Matrix cos2x = (2.e0 * x).array().cos();
  Matrix E = - std::exp(- theta_ * t) * cosx;
  S = std::exp(- theta_ * t) * M * (eps_ * v * theta_ * cosx.array()
                                - theta_ * sinx.array() + v2 * sinx.array()) +
      std::exp(- 2 * theta_ * t) * M * (eps_ * v * theta_ * sin2x.array()/theta_ 
                                - 1.e0 / theta_ * v2 * cos2x.array() + 
                                cosx.array() * cosx.array() * (1.e0 - v2 / theta_)) +
      std::exp(- 3 * theta_ * t) * M * (cosx.array() * cosx.array()) * sinx.array() * 
                            (1.e0 / theta_  - v2 / (theta_ * theta_));
  S = S * gamma_ + v2 * M / theta_ * std::exp(- theta_ * t) * sinx;
  S = S.array() + E.array().square() * (M / theta_ - v2 * M / std::pow(theta_, 2)) + eps_ * v * M / theta_ * (- theta_ * E.array());
  return S;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fsource(const real_t& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  real_t S;
  real_t E = - std::exp(- theta_ * t) * std::cos(x);
  S = std::exp(- theta_ * t) * M * (eps_ * v * theta_ * std::cos(x) - 
                                theta_ * sin(x) + v2 * sin(x)) +
      std::exp(- 2 * theta_ * t) * M * (eps_ * v * theta_ * std::sin(2*x)/theta_ 
                                - v2 / theta_ * std::cos(2*x) + 
                                std::cos(x) * std::cos(x) * (1.e0 - v2 / theta_)) +
      std::exp(- 3 * theta_ * t) * M * (std::cos(x) * std::cos(x)) * std::sin(x) * 
                            (1.e0 / theta_  - v2 / (theta_ * theta_));
  S = S * gamma_ + v2 * M / theta_ * std::exp(- theta_ * t) * std::sin(x);
  S = S + std::pow(E, 2) * (M / theta_ - v2 * M / std::pow(theta_, 2)) + eps_ * v * M / theta_ * (- theta_ * E);
  return S;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::rho_d(const Matrix& x)
{
  Matrix rhod = x;
  rhod.setConstant(1.e0);
  return rhod;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::rho_d(const real_t& x)
{
  return 1.e0;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::rho_real(const Matrix& x, const real_t& t)
{
  Matrix rho = rho_d(x).array() + gamma_ * std::exp(- theta_ * t) * x.array().sin();
  return rho;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::rho_real(const real_t& x, const real_t& t)
{
  real_t rho = rho_d(x) + gamma_ * std::exp(- theta_ * t) * std::sin(x);
  return rho;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::f_real(const Matrix& x, const real_t& v, const real_t& t)
{
  Matrix f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::f_real(const real_t& x, const real_t& v, const real_t& t)
{
  real_t f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

Matrix KineticDriftD_DG_IMEX_IM_Schur_period::g_real(const Matrix& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  Matrix g =  - v * M * (std::exp(- theta_ * t) * x.array().cos() + 
            std::exp(- 2.e0 * theta_ * t) / theta_ * x.array().cos() * x.array().sin());
  g = g.array() * gamma_ - std::exp(- theta_ * t) * v * M * x.array().cos() / theta_;
  return g;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::g_real(const real_t& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  real_t g =  - v * M * (std::exp(- theta_ * t) * std::cos(x) + 
            std::exp(- 2.e0 * theta_ * t) / theta_ * std::cos(x) * std::sin(x));
  g = g * gamma_ - std::exp(- theta_ * t) * v * M * std::cos(x) / theta_;
  return g;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fL_bc(const int& j, const real_t& t,
                    const model_data_& modal)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need fL_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fR_bc(const int& j, const real_t& t,
                    const model_data_& modal)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need fR_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fL_bc(const int& j, const real_t& t,
                    const Matrix& rho, const std::vector<Matrix>& g)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need fL_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fR_bc(const int& j, const real_t& t,
                    const Matrix& rho, const std::vector<Matrix>& g)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need fR_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fL_explicit_bc(const int& j, const real_t& t)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need fL_explicit_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::fR_explicit_bc(const int& j, const real_t& t)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need fR_explicit_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::gL_bc(const int& j, const real_t& t,
                    const model_data_& modal, const real_t& rho_L)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need gL_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::gR_bc(const int& j, const real_t& t,
                    const model_data_& modal, const real_t& rho_R)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need gR_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::phi_bc(const real_t& x, const real_t& t)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need phi_bc() ! ");
  return 0.e0;
};

real_t KineticDriftD_DG_IMEX_IM_Schur_period::rho_numericalbc(const real_t& x, const real_t& t,
                    const Matrix& rho, const std::vector<Matrix>& g)
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need rho_numericalbc() ! ");
  return 0.e0;
};

void KineticDriftD_DG_IMEX_IM_Schur_period::SolveRho_DichletBoundary(
                                                const real_t& Trun,
                                                const real_t& dt,
                                                const real_t& a,
                                                const Matrix& dh_bc,
                                                const Matrix& ah_bc,
                                                const Matrix& rho_prev,
                                                model_data_* modal) 
{
  QUEST_ERROR(" KineticDriftD_DG_IMEX_IM_Schur_period do not need SolveRho_DichletBoundary() ! ");
};

void KineticDriftD_DG_IMEX_IM_Schur_period::SolveRho_PeriodBoundary(const real_t& Trun,
                              const real_t& dt,
                              const real_t& a,
                              const Matrix& rho_prev,
                              model_data_* modal)
{
  // ********* 传入相关变量 ********** //
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
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
  real_t parainv = 1.e0 / (eps2_ + a * sigmas_ * dt);
  Matrix b0_bc = modal->rho; // 周期边界不需要处理边界条件

  Matrix stag;
  Vector tmp;
  Vector b0_hat = Vector::Zero(NTdofs);
  real_t vj, Mj;
  SparseMatrix Coef = a * dt * Minv_ * parainv;
  SparseMatrix H_temp;
  
  for (int j = 0; j < Nv; j++)
  {
    vj = V(j);
    Mj = Maxwell(V(j));
    Eigen::Map<const Vector> gj(modal->g[j].data(), NTdofs); // 周期边界不需要处理边界条件
    tmp = Coef * (vweights(j) * gj);
    b0_hat += vj * (Da_ * tmp);
  };
  Eigen::Map<Vector> rho_old_vector(b0_bc.data(), NTdofs);
  Vector b0_temp = rho_old_vector - b0_hat; // Schur补的右端项
  real_t ppp = dt * a * dt * a * parainv;
  H_temp = M_ + ppp * temp_;

  // H_ += 0.5e0 / gamma_ * M_; // 这个是看要不要加罚项
 
  int iter = 0;
  real_t rho_error = 10;
  Matrix phi_Dirichlet(1, 2);
  phi_Dirichlet.setZero(); 
  Matrix poi_lhs, rhotemp_nodal;
  Matrix rho_prev_nodal;
  Eigen::Map<const Vector> rho_prev_vec(rho_prev.data(), NTdofs);
  fe_->modal_to_nodal1D(rho_prev, &rho_prev_nodal);
  Vector rhotemp;
  Vector b0 = b0_temp;
  while (rho_error > iter_tol_)
  {
    Eh_compute(modal->E, &ME_);
    Coef = Minv_ * ME_;
    H_ = H_temp - ppp * Da_ * Coef;
    // b0 = b0_temp + 0.5e0 / gamma_ * M_ * rho_prev_vec; // 添加罚项对应的右端项的变化
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
      gmres_.set_restart(300);
      gmres_.setMaxIterations(NTdofs);
      gmres_.compute(H_);
      // std::cout << "  iterations = " << gmres_.iterations() << std::endl;
      // std::cout << "  The end of Solving ! " << std::endl;
      break;

    default:
      QUEST_ERROR(" The Schur complement Solver is not been implemented ! ");
      break;
    }
    
    switch (schur_solver_type_)
    {
    case Solver1DType::LU:
      rhotemp = lu_.solve(b0);
      break;
    case Solver1DType::GMRES:
      rhotemp = gmres_.solve(b0);
      // std::cout << "  error = " << gmres_.error() << std::endl;
      // std::cout << "  iterations = " << gmres_.iterations() << std::endl;
      if (gmres_.info() != Eigen::Success) {
        real_t rel_error = gmres_.error();
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
    Eigen::Map<Matrix> rhotemp_modal(rhotemp.data(), polydim, ncell);
    fe_->modal_to_nodal1D(rhotemp_modal, &rhotemp_nodal);
    fe_->computerrorL2(rhotemp_nodal, rho_prev_nodal, &rho_error); // 计算相对残差
    rho_prev_nodal = rhotemp_nodal;

    poi_lhs = (rho_d_nodal_ - rhotemp_nodal) / gamma_;  // 泊松方程右端项
    // 求解泊松方程并且进入下一次循环
    poisson_solver_->solveall(phi_Dirichlet, poi_lhs, &(modal->E), &(modal->phi)); 
    iter++;
    std::cout << "Schur iteration: " << iter  << ", rho_error = " << rho_error << std::endl;
  }
  Eigen::Map<Matrix> rhotemp_modal(rhotemp.data(), polydim, ncell);
  modal->rho = rhotemp_modal;
};

void KineticDriftD_DG_IMEX_IM_Schur_period::updateAll(const real_t& Trun, 
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
  Matrix ah_bc = Matrix::Zero(1,extboundarynum);
  Matrix dh_bc = Matrix::Zero(1,extboundarynum);
  Matrix rhoboundary_flux, vgboundary_flux;
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
    kinetic_modal_stages_[s].E = kinetic_modal_.E;
    kinetic_modal_stages_[s].phi = kinetic_modal_.phi;
    
    rhoboundary_flux.setZero();
    vgboundary_flux.setZero();
    
    sourceh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(source_sumh_[s]), &(sourceh_[s]));
    
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - (dt * Ai(s, i)) * ah_[i]
                                    + (dt * Ai(s, i)) * source_sumh_[i];
      #pragma omp parallel num_threads(NTH_), default(shared)
      { 
      real_t vj, Mj;
      #pragma omp for schedule(static)
        for (int j = 0; j < Nv; j++)
        {
          vj = V(j);
          Mj = Maxwell(vj);
          kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                        - (dt * Ae(s, i) * eps_) * bh_[i][j]
                                        - (dt * Ae(s, i) * eps_) * ch_[i][j]
                                        + (dt * Ai(s, i) * vj * Mj) * dh_[i]
                                        - (dt * Ai(s, i) * vj * Mj) / theta_ * mh_[i]
                                        - (dt * sigmas_ * Ai(s, i)) * sh_[i][j]
                                        + (dt * Ai(s, i) * eps2_) * sourceh_[i][j];
        };
      };
    };    
   
    kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                  + (dt * Ai(s, s)) * source_sumh_[s];
    #pragma omp parallel num_threads(NTH_), default(shared)
    {
    #pragma omp for schedule(static)
      for (int j = 0; j < Nv; j++)
      {
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                      + (dt * Ai(s, s) * eps2_) * sourceh_[s][j];
      };
    };
    
    SolveRho_PeriodBoundary(Trun + ci(s) * dt,  dt,  Ai(s, s),
                            kinetic_modal_.rho,
                            &(kinetic_modal_stages_[s]));
    
    dh_compute(kinetic_modal_stages_[s].rho,
              kinetic_modal_stages_[s].g,
              rhoboundary_flux, &(dh_[s]));
    
    Eh_compute(kinetic_modal_stages_[s].E, &ME_);
    mh_compute(kinetic_modal_stages_[s].rho,
                kinetic_modal_stages_[s].g, 
                kinetic_modal_stages_[s].E, 
                &(mh_[s]));
    #pragma omp parallel num_threads(NTH_), default(shared)
    {
      real_t vj, Mj;
    #pragma omp for schedule(static)
      for (int j = 0; j < Nv; j++)
      {
        vj = V(j);
        Mj = Maxwell(vj);
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                      + (dt * Ai(s, s) * vj * Mj) * dh_[s]
                                      - (dt * Ai(s, s) * vj * Mj) / theta_ * mh_[s];
      };
    }
    
    SolveG(Trun, dt, Ai(s, s), &(kinetic_modal_stages_[s]));
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));

    ch_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(ch_[s]));
    
    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));
    
    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
  };
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
  kinetic_modal_.E = kinetic_modal_stages_[stages_ - 1].E;
  kinetic_modal_.phi = kinetic_modal_stages_[stages_ - 1].phi;
};

void KineticDriftD_DG_IMEX_IM_Schur_period::getrho_real_modal(const real_t& Tstop, Matrix* rho_real_modal)
{
  fe_->Project_Final(
      [this](const Matrix& x, const real_t& t) { return this->rho_real(x, t); }, 
      Tstop, 
      rho_real_modal);
};

void KineticDriftD_DG_IMEX_IM_Schur_period::getrho_real_nodal(const real_t& Tstop, Matrix* rho_real_nodal)
{
  fe_->Interpolate_Final(
      [this](const Matrix& x, const real_t& t) { return this->rho_real(x, t); },
      Tstop, 
      rho_real_nodal);
};

void KineticDriftD_DG_IMEX_IM_Schur_period::getg_real_nodal(const real_t& Tstop, 
                    std::vector<Matrix>* g_real_nodal)
{
  // ** 传入相关变量
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t vj;
  g_real_nodal->resize(Nv);
  for (int j = 0; j < Nv; j++)
  {
    vj = V(j);
    fe_->Interpolate_Final(
      [this, vj](const Matrix& x, const real_t& t) { return this->g_real(x, vj, t); },
      Tstop,
      &(g_real_nodal->at(j)));
  };
};

KineticLinearD_DG_IMEX_IM_Schur_period::KineticLinearD_DG_IMEX_IM_Schur_period(
                  const KineticTensorMesh1D* mesh1D,
                  const fespace1D* fe,
                  const IMEX_RK* rk_table,
                  const PoissonSolver1D_period* poisson_solver,
                  const Solver1DType& schur_solver_type)
  : KineticDriftD_DG_IMEX_IM_Schur_period(mesh1D, fe, rk_table, poisson_solver, schur_solver_type) {};

void KineticLinearD_DG_IMEX_IM_Schur_period::init()
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& numqua = fe_->getnumqua(); 
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweight = mesh1D_->getvweights();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  const int& NTdofs = fe_->getNTdofs();
  // *** 
  QUEST_VERIFY(fe_->getmesh1D()->IsPeriodBoundary(), "The mesh is not periodical !");
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();
  setMaxwell();
  std::cout << " Maxwell sum = " << Maxwell_sum_ << std::endl;
  // QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  // QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");
  ah_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_); 
  source_sumh_.resize(stages_);
  sourceh_.resize(stages_);

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_d(x); }, &rho_d_modal_);
  fe_->Interpolate_Initial(
      [this](const Matrix& x) { return this->rho_d(x); }, &rho_d_nodal_);
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  };
  kinetic_modal_stages_.resize(stages_);
  zero_modal_mat_ = Matrix::Zero(polydim, ncell);
  zero_nodal_mat_ = Matrix::Zero(numqua, ncell);

  D_compute(beta1_, &Da_);
  Da_ = - Da_;
  // Da_ext_treat(0.e0, 0.e0, Da_, &Da_ext_);
  D_compute(- beta1_, &Db_); // 交替通量
  // Db_ext_treat(0.e0, 0.e0, Db_, &Db_ext_);

  M_compute(&M_);
  // Mbc_compute(&Mbc_);
  Minv_compute(&Minv_);
  temp_.resize(NTdofs, NTdofs);
  // real_t Mj;
  // for (int j = 0; j < Nv; j++)
  // {
  //   Mj = Maxwell(V(j));
  //   temp_ += vweight(j) * Mj * (Da_ * V(j))  * (Minv_ * (Db_ * V(j)));
  // };
  temp_ = Da_ * Minv_ * Db_;
  temp_ = temp_ * theta_;

};

void KineticLinearD_DG_IMEX_IM_Schur_period::setdt(real_t* dt)
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
  *dt = gamma_ * (*dt) / 5;
  // std::cout << " dt = " << *dt << std::endl;
  // PAUSE();
};

void KineticLinearD_DG_IMEX_IM_Schur_period::updateAll(const real_t& Trun, 
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
  Matrix ah_bc = Matrix::Zero(1,extboundarynum);
  Matrix dh_bc = Matrix::Zero(1,extboundarynum);
  Matrix rhoboundary_flux, vgboundary_flux;
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

    sourceh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(source_sumh_[s]), &(sourceh_[s]));
    
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - (dt * Ai(s, i)) * ah_[i]
                                    + (dt * Ai(s, i)) * source_sumh_[i];
      #pragma omp parallel num_threads(NTH_), default(shared)
      { 
      real_t vj, Mj;
      #pragma omp for schedule(static)
        for (int j = 0; j < Nv; j++)
        {
          vj = V(j);
          Mj = Maxwell(vj);
          kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                        - (dt * Ae(s, i) * eps_) * bh_[i][j]
                                        + (dt * Ai(s, i) * vj * Mj) * dh_[i]
                                        - (dt * sigmas_ * Ai(s, i)) * sh_[i][j]
                                        + (dt * Ai(s, i) * eps2_) * sourceh_[i][j];
        };
      };
    };    
   
    kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho + dt * Ai(s, s) * source_sumh_[s];
    #pragma omp parallel num_threads(NTH_), default(shared)
    {
    real_t vj;
    #pragma omp for schedule(static)
      for (int j = 0; j < Nv; j++)
      {
        vj = V(j);
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j]
                                      + (dt * Ai(s, s) * eps2_) * sourceh_[s][j];
      };
    };
    
    SolveRho_PeriodBoundary(Trun + ci(s) * dt,  dt,  Ai(s, s),
                            kinetic_modal_.rho,
                            &(kinetic_modal_stages_[s]));
    
    dh_compute(kinetic_modal_stages_[s].rho,
              kinetic_modal_stages_[s].g,
              rhoboundary_flux, &(dh_[s]));

    #pragma omp parallel num_threads(NTH_), default(shared)
    {
      real_t vj, Mj;
    #pragma omp for schedule(static)
      for (int j = 0; j < Nv; j++)
      {
        vj = V(j);
        Mj = Maxwell(vj);
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j] 
                                      + (dt * Ai(s, s) * vj * Mj) * dh_[s];
      };
    }
    
    SolveG(Trun, dt, Ai(s, s), &(kinetic_modal_stages_[s]));
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));
    
    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));
    
    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
  };
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
};

void KineticLinearD_DG_IMEX_IM_Schur_period::SolveRho_PeriodBoundary(const real_t& Trun,
                                const real_t& dt,
                                const real_t& a,
                                const Matrix& rho_prev,
                                model_data_* modal)
{
// ********* 传入相关变量 ********** //
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
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
  real_t parainv = 1.e0 / (eps2_ + a * sigmas_ * dt);
  Matrix b0_bc = modal->rho; // 周期边界不需要处理边界条件

  Matrix stag;
  Vector tmp;
  Vector b0_hat = Vector::Zero(NTdofs);
  real_t vj, Mj;
  SparseMatrix Coef = a * dt * Minv_ * parainv;
  
  for (int j = 0; j < Nv; j++)
  {
    vj = V(j);
    Mj = Maxwell(V(j));
    Eigen::Map<const Vector> gj(modal->g[j].data(), NTdofs); // 周期边界不需要处理边界条件
    tmp = Coef * (vweights(j) * gj);
    b0_hat += vj * (Da_ * tmp);
  };
  Eigen::Map<Vector> rho_old_vector(b0_bc.data(), NTdofs);
  Vector b0 = rho_old_vector - b0_hat; // Schur补的右端项
  real_t ppp = dt * a * dt * a * parainv;
  H_ = M_ + ppp * temp_;

  int iter = 0;
  real_t rho_error = 10;
  Matrix rhotemp_nodal;
  Vector rhotemp;

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
    gmres_.set_restart(300);
    gmres_.setMaxIterations(NTdofs);
    gmres_.compute(H_);
    // std::cout << "  iterations = " << gmres_.iterations() << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;

  default:
    QUEST_ERROR(" The Schur complement Solver is not been implemented ! ");
    break;
  }
  
  switch (schur_solver_type_)
  {
  case Solver1DType::LU:
    rhotemp = lu_.solve(b0);
    break;
  case Solver1DType::GMRES:
    rhotemp = gmres_.solve(b0);
    // std::cout << "  error = " << gmres_.error() << std::endl;
    // std::cout << "  iterations = " << gmres_.iterations() << std::endl;
    if (gmres_.info() != Eigen::Success) {
      real_t rel_error = gmres_.error();
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

  Eigen::Map<Matrix> rhotemp_modal(rhotemp.data(), polydim, ncell);
  modal->rho = rhotemp_modal;
};

void KineticLinearD_DG_IMEX_IM_Schur_period::setitertol(const real_t& iter_tol)
{
  QUEST_ERROR(" KineticLinearD_DG_IMEX_IM_Schur_period do not need setitertol() ! ");
};

void KineticLinearD_DG_IMEX_IM_Schur_period::Eh_compute(const Matrix& E_modal, SparseMatrix* ME) 
{
  QUEST_ERROR(" KineticLinearD_DG_IMEX_IM_Schur_period do not need Eh_compute() ! ");
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::rho_init(const Matrix& x)
{
  Matrix rho = rho_d(x).array() + gamma_ * x.array().sin();
  return rho;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::rho_init(const real_t& x)
{
  return rho_d(x) + gamma_ * std::sin(x);
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::f_init(const Matrix& x, const real_t& v)
{
  Matrix f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::f_init(const real_t& x, const real_t& v)
{
  real_t f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::g_init(const Matrix& x, const real_t& v)
{
  real_t M = Maxwell(v);
  Matrix g = - v * M * x.array().cos();
  g = g * gamma_;
  return g;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::g_init(const real_t& x, const real_t& v)
{
  real_t M = Maxwell(v);
  real_t g = - v * M * std::cos(x);
  g = g * gamma_;
  return g;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::source(const Matrix& x, const real_t& t)
{
  QUEST_ERROR(" KineticLinearD_DG_IMEX_IM_Schur_period do not need source() ! ");
  Matrix S = x;
  S.setZero();
  return S;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::source(const real_t& x, const real_t& t)
{
  QUEST_ERROR(" KineticLinearD_DG_IMEX_IM_Schur_period do not need source() ! ");
  return 0.e0;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::fsource(const Matrix& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  Matrix S;
  Matrix sinx = x.array().sin();
  Matrix cosx = x.array().cos();
  S = std::exp(- theta_ * t) * M * (eps_ * v * theta_ * cosx.array()
                                - theta_ * sinx.array() + v2 * sinx.array());
  S = S * gamma_;
  return S;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::fsource(const real_t& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  real_t S;
  S = std::exp(- theta_ * t) * M * (eps_ * v * theta_ * std::cos(x) - 
                                theta_ * sin(x) + v2 * sin(x));
  S = S * gamma_;
  return S;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::rho_d(const Matrix& x)
{
  Matrix rhod = x;
  rhod.setConstant(1.e0);
  return rhod;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::rho_d(const real_t& x)
{
  return 1.e0;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::rho_real(const Matrix& x, const real_t& t)
{
  Matrix rho = rho_d(x).array() + gamma_ * std::exp(- theta_ * t) * x.array().sin();
  return rho;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::rho_real(const real_t& x, const real_t& t)
{
  real_t rho = rho_d(x) + gamma_ * std::exp(- theta_ * t) * std::sin(x);
  return rho;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::f_real(const Matrix& x, const real_t& v, const real_t& t)
{
  Matrix f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::f_real(const real_t& x, const real_t& v, const real_t& t)
{
  real_t f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

Matrix KineticLinearD_DG_IMEX_IM_Schur_period::g_real(const Matrix& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  Matrix g =  - v * M * (std::exp(- theta_ * t) * x.array().cos());
  g = g.array() * gamma_;
  return g;
};

real_t KineticLinearD_DG_IMEX_IM_Schur_period::g_real(const real_t& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  real_t g =  - v * M * (std::exp(- theta_ * t) * std::cos(x));
  g = g * gamma_;
  return g;
};

}; // namespace QUEST

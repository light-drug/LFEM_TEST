#include "KineticDD_DG2D_IMEX_Schur.hpp"

namespace QUEST
{

void KineticDD_DG2d_IMEX_IM_Schur::seteps(const real_t& knu)
{
  eps_ = knu;
  eps2_ = knu * knu;
};

void KineticDD_DG2d_IMEX_IM_Schur::setCR(const real_t& CR)
{
  CR_ = CR;
};

void KineticDD_DG2d_IMEX_IM_Schur::setbeta1(const real_t& beta1)
{
  beta1_ = beta1;
};

void KineticDD_DG2d_IMEX_IM_Schur::setsigmas(const real_t& sigmas)
{
  sigmas_ = sigmas;
};

void KineticDD_DG2d_IMEX_IM_Schur::setsparsetol(const real_t& schur_tol)
{
  schur_tol_ = schur_tol;
};

void KineticDD_DG2d_IMEX_IM_Schur::setitertol(const real_t& iter_tol)
{
  iter_tol_ = iter_tol;
};

void KineticDD_DG2d_IMEX_IM_Schur::setNTH(const int& NTH)
{
  NTH_ = NTH;
};

void KineticDD_DG2d_IMEX_IM_Schur::settheta(const real_t& theta)
{
  theta_ = theta;
};

void KineticDD_DG2d_IMEX_IM_Schur::setgamma(const real_t& gamma)
{
  gamma_ = gamma;
};

const Matrix& KineticDD_DG2d_IMEX_IM_Schur::getrho_modal() const 
{
  return kinetic_modal_.rho;
};

const std::vector<Matrix>& KineticDD_DG2d_IMEX_IM_Schur::getg_modal() const 
{
  return kinetic_modal_.g;
};

const Matrix& KineticDD_DG2d_IMEX_IM_Schur::getrho_d_modal() const 
{
  return rho_d_modal_;
};


KineticDD_DG2d_IMEX_IM_Schur::KineticDD_DG2d_IMEX_IM_Schur(const TensorMesh1D* mesh1D,
                                const TensorMesh2D* kinetic_mesh2D,
                                const fespace1D* fe1D,
                                const fespace2D* kinetic_fe2D,
                                const IMEX_RK* rk_table,
                                const PoissonSolver1D* poisson_solver,
                                const Solver1DType& schur_solver_type)
  : mesh1D_(mesh1D),  kinetic_mesh2D_(kinetic_mesh2D), 
  fe1D_(fe1D), kinetic_fe2D_(kinetic_fe2D), 
  rk_table_(rk_table), poisson_solver_(poisson_solver),
  schur_solver_type_(schur_solver_type) {};

void KineticDD_DG2d_IMEX_IM_Schur::MixedFlux_compute()
{
  // ** 传入相关变量
  const real_t& x1 = kinetic_mesh2D_->getx1();
  const real_t& x2 = kinetic_mesh2D_->getx2();
  const real_t& v1 = kinetic_mesh2D_->gety1();
  const real_t& v2 = kinetic_mesh2D_->gety2();
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& test_ref_rho = fe1D_->gettest_ref();
  const Matrix& test_ref_dx_rho = fe1D_->gettest_dx_ref();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const int& Nv = kinetic_mesh2D_->getyDiv();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  const std::vector<Vector>& boundary_u_rho_temp = fe1D_->getboundary_u();
  const std::vector<Matrix>& boundary_u_g = kinetic_fe2D_->getboundary_u();
  const DiagnalMatrix& wqua_diag1d_ref = 
    kinetic_fe2D_->getwqua_diag1d();
  // *** 

  boundary_u_rho_.resize(2);
  Matrix ones(1, numqua_rho);
  for (int i = 0; i < boundary_u_rho_temp.size(); i++)
  {
    boundary_u_rho_[i] = boundary_u_rho_temp[i] * ones;
  }
  mixedflux_u_v_.resize(2);
  for (int i = 0; i < 2; i++)
  {
    mixedflux_u_v_[i].resize(4);
    for (int j = 0; j < 4; j++)
    {
      mixedflux_u_v_[i][j] = boundary_u_rho[i] * wqua_diag1d_ref * boundary_u_g[j].transpose();
    };
  };

  test_ref1D_.resize(polydim_rho, numqua_g);
  test_ref1D_dx_.resize(polydim_rho, numqua_g);
  int index = 0;
  for (int qy = 0; qy < numqua_rho; qy++)
  {
    for (int qx = 0; qx < numqua_rho; qx++)
    {
      for (int i = 0; i < ploydim_rho; i++)
      {
        test_ref1D_(i, index) = test_ref_rho(i, qx);
        test_ref1D_dx_(i, index) = test_ref_dx_rho(i, qx);
      }
      index++;
    };
  };
  test_ref1D_T_ = test_ref1D_.transpose();
  test_ref1D_dx_T_ = test_ref1D_dx_.transpose();
};

Vector KineticDD_DG2d_IMEX_IM_Schur::vector_nodal1Dto2D(const Eigen::Ref<const Vector>& vec)
{
  // ** 传入相关变量
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const int& Nv = kinetic_mesh2D_->getyDiv();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  // ***
  Vector temp(numqua_g);
  int index = 0;
  for (int qy = 0; qy < numqua_rho; qy++)
  {
    for (int qx = 0; qx < numqua_rho; qx++)
    {
      temp(index) = vec(qx);
      index++;
    };
  };
};

void KineticDD_DG2d_IMEX_IM_Schur::init()
{
  // ** 传入相关变量
  const real_t& x1 = kinetic_mesh2D_->getx1();
  const real_t& x2 = kinetic_mesh2D_->getx2();
  const real_t& v1 = kinetic_mesh2D_->gety1();
  const real_t& v2 = kinetic_mesh2D_->gety2();
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const int& Nv = kinetic_mesh2D_->getyDiv();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
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
  kinetic_fe2D_->Interpolate_Initial(
              [this](const Matrix& x, const Matrix& v) { return this->V_value(x, v); }, &V_nodal_);
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
  Da_compute(beta1_, &Da_);
  Db_compute(- beta1_, &Db_); // 交替通量

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

void KineticDD_DG2d_IMEX_IM_Schur::setMaxwell()
{
  // ** 传入相关变量
  const std::vector<Matrix>& Pb = kinetic_fe2D_->getPb();
  const int& Nv = kinetic_fe2D_->getyDiv();
  // *** 
  Matrix M_Pb = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) *
                  (- Pb[1].array() * Pb[1].array() / ( 2 * theta_)).exp();
  Matrix M_modal;
  kinetic_fe2D_.nodal_to_modal2D(M_Pb, &M_modal);
  Maxwell_sum_ = 0.e0;
  for (int j = 0; j < Nv; j++)
  {
    Maxwell_sum_ += M_modal(0, mesh2D_->CellIndex(0, j));
  };
};

void KineticDD_DG2d_IMEX_IM_Schur::setdt(real_t* dt)
{
  // ** 传入相关变量
  const int& polydim = fe1D_->getbasis()->getpolydim();
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

void KineticDD_DG2d_IMEX_IM_Schur::Da_compute(const real_t& beta,
                                              SparseMatrix* Da)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& Tm_rho = fe1D_->getTm();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const std::vector<Matrix>& boundary_u_g = kinetic_fe2D_->getboundary_u();
  const Matrix& Tm_g = kinetic_fe2D_->getTm();
  const DiagnalMatrix& wqua_diag1d_ref = 
    kinetic_fe2D_->getwqua_diag1d();
  const DiagnalMatrix& wqua_diag_ref = 
    kinetic_fe2D_->getwqua_diag();
  const Matrix& test_ref2D = kinetic_mesh2D_->gettest_ref();
  const int& Nx = kinetic_mesh2D_->getxDiv();
  const int& Nv = kinetic_mesh2D_->getyDiv();
  const real_t& hx = kinetic_mesh2D_->gethx();
  const real_t& hv = kinetic_mesh2D_->gethy();

  const Matrix& dv_u = fe1D_->getdv_u();
  const Matrix& v_u = fe1D_->getv_u();
  const Vector& JacobiDet_rho = mesh1D_->getJacobiDet();
  const Vector& JacobiDet_g = mesh2D_->getJacobiDet();
  const Vector& Jx = mesh1D_->getJx();

  const std::vector<Eigen::Vector2d>& intnormals = 
    kinetic_mesh2D_->getintboundarynormal();
  const std::vector<real_t>& intlength = 
    kinetic_mesh2D_->getintboundarylength();
  const std::vector<Eigen::Vector2i>& intNei = 
    kinetic_mesh2D_->getintboundaryneighbors();
  const int& intboundaryNum = 
    kinetic_mesh2D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex =
    kinetic_mesh2D_->getintboundarytypeindex();
  const std::vector<real_t>&  extlength = 
    kinetic_mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& extnormals = 
    kinetic_mesh2D_->getextboundarynormal();
  const std::vector<int>& extNei =
    kinetic_mesh2D_->getextboundaryneighbors();
  const int& extboundaryNum =
    kinetic_mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBTypeIndex =
    kinetic_mesh2D_->getextboundarytypeindex();
  const std::vector<Eigen::Vector2d>& extboundarycenter = 
    kinetic_mesh2D_->getextboundarycenter();
  const Matrix& CellCenter = 
    kinetic_mesh2D_->getCellCenter();
  const std::vector<Matrix>& CoorBdrRef2d = 
    kinetic_fe2D_->getCoorBdrRef();
  // ******************************** //
  
  std::vector<Eigen::Triplet<real_t>> tripletList_Da;
  int estimatedNonZeros = ncell_g * polydim_rho * ploydim_g
                        + intboundaryNum * polydim_rho * ploydim_g * 4
                        + extboundaryNum * polydim_rho * polydim_g * 4;
  tripletList_Da.reserve(estimatedNonZeros);
  Da->resize(NTdofs_rho, NTdofs_g);
  Vector V_vec;
  V_vec.resize(numqua_g);
  for (int i = 0; i < ncell_g; i++) 
  {
    int x_cell_index = mesh2D_->xCellIndex(i);
    V_vec = V_nodal_.col(i);
    Bref = test_ref2D.array().rowwise() * V_vec.array();
    Bref = test_ref1D_dx_ * wqua_diag_ref * Bref.transpose();
    for (int test_basis_index = 0; test_basis_index < polydim_rho; test_basis_index++) 
    {
      int alpha = Tm_rho(test_basis_index, x_cell_index);
      for (int trial_basis_index = 0; trial_basis_index < polydim_g; trial_basis_index++) 
      {
        int beta = Tm_g(trial_basis_index, i);
        real_t ini_value = - Bref(test_basis_index, trial_basis_index) * JacobiDet_g(0) * Jx(0);
        tripletList_Da.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };

  int test_cell_Index, trial_cell_Index;
  int test_x_cell_index;
  Matrix test_qua_value;
  Matrix trial_qua_value;
  Matrix temp;
  real_t test_normal, trial_normal;
  real_t normal;
  real_t len;
  int alpha, beta;
  real_t ini_value;
  V_vec.resize(numqua_rho);
  for (int i = 0; i < intboundaryNum; i++) 
  {
    normal = intnormals[i](0);
    len = intlength[i];
    if (std::abs(normal) > 0.5e0)
    {
      for (int q = 0; q < numqua_rho; q++)
      {
        V_vec(q) = CellCenter(1, intNei[i](0)) 
                  + CoorBdrRef2d[IntBTypeIndex[i](0)](1, q) * hv;
      }
      for (int test_cell = 0; test_cell < 2; test_cell++) {
        test_cell_Index = intNei[i](test_cell);
        test_x_cell_index = kinetic_mesh2D_->xCellIndex(test_cell_index);
        test_qua_value = boundary_u_rho_[test_cell];
        test_normal = intnormals[i](0) * std::pow(-1.e0, test_cell);

        for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
          trial_cell_Index = intNei[i](trial_cell);
          trial_qua_value = boundary_u_g[IntBTypeIndex[i](trial_cell)];
          trial_normal = intnormals[i](0) * std::pow(-1.e0, trial_cell);  
      
          temp = test_qua_value * V_vec.asDiagonal() 
                  * wqua_diag1d_ref * trial_qua_value.transpose() * len;
          for (int test_basis_index = 0; test_basis_index < polydim_rho; test_basis_index++) 
          {
            alpha = Tm_rho(test_basis_index, test_x_cell_index);
            for (int trial_basis_index = 0; trial_basis_index < polydim_g; trial_basis_index++) 
            {
              beta = Tm_g(trial_basis_index, trial_cell_Index);
              ini_value =  (0.5e0 + be * trial_normal) 
                          * temp(test_basis_index, trial_basis_index) * test_normal;
              // 这是一维代码不需要积分
              if (std::abs(ini_value >= 1.e-14) > 1.e-14)
              {
                tripletList_Da.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
              }
            };
          };
        };
      };
    };
  };
  
  Matrix temp_penalty;
  for (int i = 0; i < extboundaryNum; i++) 
  {
    normal = extnormals[i](0);
    len = extlength[i];
    if (std::abs(normal) > 0.5e0)
    {
      for (int q = 0; q < numqua_rho; q++)
      {
        V_vec(q) = CellCenter(1, extNei[i]) 
                  + CoorBdrRef2d[ExtBTypeIndex[i]](1, q) * hv;
      }
      test_cell_Index = ExtNei[i];
      test_x_cell_index = kinetic_mesh2D_->xCellIndex(test_cell_index);
      if (test_x_cell_index == 0)
      {
        test_qua_value = boundary_u_rho_[1];
      } else if (test_x_cell_index == Nx - 1)
      {
        test_qua_value = boundary_u_rho_[0];
      }
      test_normal = normal;

      trial_cell_Index = ExtNei[i];
      trial_qua_value = boundary_u_g[ExtBTypeIndex[i]];
      trial_normal = normal;  
  
      temp = test_qua_value * V_vec.asDiagonal() 
              * wqua_diag1d_ref * trial_qua_value.transpose() * len;

      V_vec = ((V_vec.array() * normal) > 0).cast<double>();
      temp_penalty = test_qua_value * V_vec.asDiagonal() 
                    * wqua_diag1d_ref * trial_qua_value.transpose() * len;
      for (int test_basis_index = 0; test_basis_index < polydim_rho; test_basis_index++) 
      {
        alpha = Tm_rho(test_basis_index, test_x_cell_index);
        for (int trial_basis_index = 0; trial_basis_index < polydim_g; trial_basis_index++) 
        {
          beta = Tm_g(trial_basis_index, trial_cell_Index);
          ini_value = temp(test_basis_index, trial_basis_index) * test_normal;
          ini_value -= temp_penalty(test_basis_index, trial_basis_index);
          // 这是一维代码不需要积分
          if (std::abs(ini_value >= 1.e-14) > 1.e-14)
          {
            tripletList_Da.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          }
        };
      };
    };
  };
  
  Da->setFromTriplets(tripletList_Da.begin(), tripletList_Da.end());
  // (*Da_ext)[j] = vj * Da_ + da_bc;
};

void KineticDD_DG2d_IMEX_IM_Schur::Db_compute(const real_t& beta,
                                                SparseMatrix* Db)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& Tm_rho = fe1D_->getTm();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const std::vector<Matrix>& boundary_u_g = kinetic_fe2D_->getboundary_u();
  const Matrix& Tm_g = kinetic_fe2D_->getTm();
  const DiagnalMatrix& wqua_diag1d_ref = 
    kinetic_fe2D_->getwqua_diag1d();
  const DiagnalMatrix& wqua_diag_ref = 
    kinetic_fe2D_->getwqua_diag();
  const Matrix& test_ref2D = kinetic_mesh2D_->gettest_ref();
  const Matrix& test_ref2D_dx = kinetic_mesh2D_->gettest_ref_dx();
  const int& Nx = kinetic_mesh2D_->getxDiv();
  const int& Nv = kinetic_mesh2D_->getyDiv();
  const real_t& hx = kinetic_mesh2D_->gethx();
  const real_t& hv = kinetic_mesh2D_->gethy();

  const Matrix& dv_u = fe1D_->getdv_u();
  const Matrix& v_u = fe1D_->getv_u();
  const Vector& JacobiDet_rho = mesh1D_->getJacobiDet();
  const Vector& JacobiDet_g = mesh2D_->getJacobiDet();
  const Vector& Jx = mesh1D_->getJx();

  const std::vector<Eigen::Vector2d>& intnormals = 
    kinetic_mesh2D_->getintboundarynormal();
  const std::vector<real_t>& intlength = 
    kinetic_mesh2D_->getintboundarylength();
  const std::vector<Eigen::Vector2i>& intNei = 
    kinetic_mesh2D_->getintboundaryneighbors();
  const int& intboundaryNum = 
    kinetic_mesh2D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex =
    kinetic_mesh2D_->getintboundarytypeindex();
  const std::vector<real_t>&  extlength = 
    kinetic_mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& extnormals = 
    kinetic_mesh2D_->getextboundarynormal();
  const std::vector<int>& extNei =
    kinetic_mesh2D_->getextboundaryneighbors();
  const int& extboundaryNum =
    kinetic_mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBTypeIndex =
    kinetic_mesh2D_->getextboundarytypeindex();
  const std::vector<Eigen::Vector2d>& extboundarycenter = 
    kinetic_mesh2D_->getextboundarycenter();
  const Matrix& CellCenter = 
    kinetic_mesh2D_->getCellCenter();
  const std::vector<Matrix>& CoorBdrRef2d = 
    kinetic_fe2D_->getCoorBdrRef();
  // ******************************** //
  
  std::vector<Eigen::Triplet<real_t>> tripletList_Db;
  int estimatedNonZeros = ncell_g * polydim_rho * ploydim_g
                        + intboundaryNum * polydim_rho * ploydim_g * 4
                        + extboundaryNum * polydim_rho * polydim_g * 4;
  tripletList_Db.reserve(estimatedNonZeros);
  Db->resize(NTdofs_g, NTdofs_rho);
  Vector V_vec;
  V_vec.resize(numqua_g);
  for (int i = 0; i < ncell_g; i++) 
  {
    int x_cell_index = mesh2D_->xCellIndex(i);
    V_vec = V_nodal_.col(i);
    V_vec = V_vec.array() * Maxwell_nodal_.col(i).array();
    Bref = test_ref1D_.array().colwise() * V_vec.array();
    Bref = test_ref2D_dx * wqua_diag_ref * Bref.transpose();
    for (int test_basis_index = 0; test_basis_index < polydim_g; test_basis_index++) 
    {
      int alpha = Tm_g(test_basis_index, i);
      for (int trial_basis_index = 0; trial_basis_index < polydim_rho; trial_basis_index++) 
      {
        int beta = Tm_rho(trial_basis_index, x_cell_index);
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet_g(0) * Jx(0);
        tripletList_Db.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };

  int test_cell_Index, trial_cell_Index;
  int trial_x_cell_index;
  Matrix test_qua_value;
  Matrix trial_qua_value;
  Matrix temp;
  real_t test_normal, trial_normal;
  real_t normal;
  real_t len;
  int alpha, beta;
  real_t ini_value;
  V_vec.resize(numqua_rho);
  for (int i = 0; i < intboundaryNum; i++) 
  {
    normal = intnormals[i](0);
    len = intlength[i];
    if (std::abs(normal) > 0.5e0)
    {
      for (int q = 0; q < numqua_rho; q++)
      {
        V_vec(q) = CellCenter(1, intNei[i](0)) 
                  + CoorBdrRef2d[IntBTypeIndex[i](0)](1, q) * hv;
      }
      V_vec = V_vec.array() * Maxwell(V_vec).array();
      for (int test_cell = 0; test_cell < 2; test_cell++) {
        test_cell_Index = intNei[i](test_cell);
        test_qua_value = boundary_u_g[IntBTypeIndex[i](test_cell)];
        test_normal = intnormals[i](0) * std::pow(-1.e0, test_cell);

        for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
          trial_cell_Index = intNei[i](trial_cell);
          trial_x_cell_index = kinetic_mesh2D_->xCellIndex(trial_cell_index);
          trial_qua_value = boundary_u_rho_[trial_cell];
          trial_normal = intnormals[i](0) * std::pow(-1.e0, trial_cell);  
      
          temp = test_qua_value * V_vec.asDiagonal() 
                  * wqua_diag1d_ref * trial_qua_value.transpose() * len;
          for (int test_basis_index = 0; test_basis_index < polydim_g; test_basis_index++) 
          {
            alpha = Tm_g(test_basis_index, test_cell_Index);
            for (int trial_basis_index = 0; trial_basis_index < polydim_rho; trial_basis_index++) 
            {
              beta = Tm_rho(trial_basis_index, trial_x_cell_index);
              ini_value =  (0.5e0 + be * trial_normal) 
                          * temp(test_basis_index, trial_basis_index) * test_normal;
              // 这是一维代码不需要积分
              if (std::abs(ini_value >= 1.e-14) > 1.e-14)
              {
                tripletList_Db.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
              }
            };
          };
        };
      };
    };
  };
  
  for (int i = 0; i < extboundaryNum; i++) 
  {
    normal = extnormals[i](0);
    len = extlength[i];
    if (std::abs(normal) > 0.5e0)
    {
      for (int q = 0; q < numqua_rho; q++)
      {
        V_vec(q) = CellCenter(1, extNei[i]) 
                  + CoorBdrRef2d[ExtBTypeIndex[i]](1, q) * hv;
      }
      V_vec = V_vec.array() * Maxwell(V_vec).array();
      test_cell_Index = ExtNei[i];
      test_qua_value = boundary_u_g[ExtBTypeIndex[i]];
      test_normal = normal;

      trial_cell_Index = ExtNei[i];
      trial_x_cell_index = kinetic_mesh2D_->xCellIndex(trial_cell_index);
      if (test_x_cell_index == 0)
      {
        trial_qua_value = boundary_u_rho_[1];
      } else if (test_x_cell_index == Nx - 1)
      {
        trial_qua_value = boundary_u_rho_[0];
      }
      trial_normal = normal;  
  
      temp = test_qua_value * V_vec.asDiagonal() 
              * wqua_diag1d_ref * trial_qua_value.transpose() * len;
      
      for (int test_basis_index = 0; test_basis_index < polydim_g; test_basis_index++) 
      {
        alpha = Tm_g(test_basis_index, test_cell_index);
        for (int trial_basis_index = 0; trial_basis_index < polydim_rho; trial_basis_index++) 
        {
          beta = Tm_rho(trial_basis_index, trial_x_cell_index);
          ini_value = - 0.5e0 * temp(test_basis_index, trial_basis_index) * test_normal;
          // 这是一维代码不需要积分
          if (std::abs(ini_value >= 1.e-14) > 1.e-14)
          {
            tripletList_Db.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          }
        };
      };
    };
  };
  
  Db->setFromTriplets(tripletList_Db.begin(), tripletList_Db.end());
  // (*Da_ext)[j] = vj * Da_ + da_bc;
};

void KineticDD_DG2d_IMEX_IM_Schur::Mrho_compute(SparseMatrix* M)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& v_u = fe1D_->getv_u();
  const Vector& JacobiDet_rho = mesh1D_->getJacobiDet();
  const IntMatrix& Tm_rho = fe_->getTm();
  // ******************************** //

  M->resize(NTdofs_rho, NTdofs_rho);
  std::vector<Eigen::Triplet<real_t>> tripletList_M;
  int estimatedNonZeros = NTdofs_rho;
  tripletList_M.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_rho; i++) 
  {
    for (int basis_index = 0; basis_index < polydim_rho; basis_index++) 
    {
      real_t ini_value = v_u(basis_index,basis_index) * JacobiDet_rho(i);
      int alpha = Tm_rho(basis_index, i);
      tripletList_M.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  M->setFromTriplets(tripletList_M.begin(), tripletList_M.end());
};

void KineticDD_DG2d_IMEX_IM_Schur::Mbc_compute(SparseMatrix* Mbc)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& v_u = fe1D_->getv_u();
  const Vector& JacobiDet = mesh1D_->getJacobiDet();
  const std::vector<Vector>& boundary_u = fe1D_->getboundary_u();
  const IntMatrix& Tm_rho = fe1D_->getTm();
  const std::vector<real_t>& extnormals = 
    mesh1D_->getextboundarynormal();
  const std::vector<int>& extNei = 
    mesh1D_->getextboundaryneighbors();
  const int& extboundaryNum = 
    mesh1D_->getextboundaryNum();
  // ******************************** //

  Mbc->resize(NTdofs_rho, NTdofs_rho);
  std::vector<Eigen::Triplet<real_t>> tripletList_Mbc;
  int estimatedNonZeros = extboundaryNum * polydim * polydim;
  tripletList_Mbc.reserve(estimatedNonZeros);

  int test_cell_Index, trial_cell_Index;
  Vector test_qua_value, trial_qua_value;
  real_t test_normal, trial_normal;
  int alpha, beta;
  real_t ini_value;
  for (int i = 0; i < extboundaryNum; i++) 
  {
    test_cell_Index = extNei[i];
    test_qua_value = boundary_u[1-i];
    test_normal = extnormals[i];

    trial_cell_Index = extNei[i];
    trial_qua_value = boundary_u[1-i];
    trial_normal = extnormals[i];

    for (int test_basis_index = 0; test_basis_index < polydim_rho; test_basis_index++) 
    {
      alpha = Tm_rho(test_basis_index, test_cell_Index);
      for (int trial_basis_index = 0; trial_basis_index < polydim_rho; trial_basis_index++) 
      {
        beta = Tm_rho(trial_basis_index, trial_cell_Index);
        ini_value =  (0.5e0 * CR_ * trial_qua_value(trial_basis_index)) 
            * test_qua_value(test_basis_index) * test_normal * trial_normal;
        tripletList_Mbc.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };

  Mbc->setFromTriplets(tripletList_Mbc.begin(), tripletList_Mbc.end());
};

void KineticDD_DG2d_IMEX_IM_Schur::Mg_compute(SparseMatrix* M)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const Matrix& v_u = kinetic_fe2D_->getv_u();
  const Vector& JacobiDet_g = kinetic_mesh2D_->getJacobiDet();
  const IntMatrix& Tm_g = kinetic_fe2D_->getTm();
  // ******************************** //

  M->resize(NTdofs_g, NTdofs_g);
  std::vector<Eigen::Triplet<real_t>> tripletList_M;
  int estimatedNonZeros = NTdofs_g;
  tripletList_M.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_g; i++) 
  {
    for (int basis_index = 0; basis_index < polydim_g; basis_index++) 
    {
      real_t ini_value = v_u(basis_index, basis_index) * JacobiDet_g(i);
      int alpha = Tm_g(basis_index, i);
      tripletList_M.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  M->setFromTriplets(tripletList_M.begin(), tripletList_M.end());
};

void KineticDD_DG2d_IMEX_IM_Schur::Mginv_compute(SparseMatrix* Minv)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const Matrix& v_u = kinetic_fe2D_->getv_u();
  const Vector& JacobiDet_g = kinetic_mesh2D_->getJacobiDet();
  const IntMatrix& Tm_g = kinetic_fe2D_->getTm();
  // ******************************** //

  Minv->resize(NTdofs_g, NTdofs_g);
  std::vector<Eigen::Triplet<real_t>> tripletList_M;
  int estimatedNonZeros = NTdofs_g;
  tripletList_Minv.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_g; i++) 
  {
    for (int basis_index = 0; basis_index < polydim_g; basis_index++) 
    {
      real_t ini_value = v_u(basis_index, basis_index) * JacobiDet_g(i);
      ini_value = 1.e0 / ini_value;
      int alpha = Tm_g(basis_index, i);
      tripletList_Minv.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  Minv->setFromTriplets(tripletList_Minv.begin(), tripletList_Minv.end());
};

void KineticDD_DG2d_IMEX_IM_Schur::Eh_compute(const Matrix& E_modal, SparseMatrix* ME)
{
  // ********* 传入相关变量 ********** //
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& Tm_rho = fe1D_->getTm();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  const std::vector<Matrix>& boundary_u_g = kinetic_fe2D_->getboundary_u();
  const Matrix& Tm_g = kinetic_fe2D_->getTm();
  const DiagnalMatrix& wqua_diag1d_ref = 
    kinetic_fe2D_->getwqua_diag1d();
  const DiagnalMatrix& wqua_diag_ref = 
    kinetic_fe2D_->getwqua_diag();
  const Matrix& test_ref2D = kinetic_mesh2D_->gettest_ref();
  const Matrix& test_ref2D_T = kinetic_mesh2D_->gettest_ref_T();
  const Matrix& test_ref2D_dx = kinetic_mesh2D_->gettest_ref_dx();
  const int& Nx = kinetic_mesh2D_->getxDiv();
  const int& Nv = kinetic_mesh2D_->getyDiv();
  const real_t& hx = kinetic_mesh2D_->gethx();
  const real_t& hv = kinetic_mesh2D_->gethy();

  const Matrix& dv_u = fe1D_->getdv_u();
  const Matrix& v_u = fe1D_->getv_u();
  const Vector& JacobiDet_rho = mesh1D_->getJacobiDet();
  const Vector& JacobiDet_g = mesh2D_->getJacobiDet();
  const Vector& Jx = mesh1D_->getJx();
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
  Vector V_vec;
  V_vec.resize(numqua_g);
  Vector E_nodal_local;
  for (int i = 0; i < ncell_g; i++) 
  {
    int x_cell_index = mesh2D_->xCellIndex(i);
    E_nodal_local = vector_nodal1Dto2D(E_nodal.col(x_cell_index));
    local_ME = test_ref2D.array().rowwise() * E_nodal_local.transpose().array();
    V_vec = V_nodal_.col(i);
    V_vec = V_vec.array() * Maxwell_nodal_.col(i).array();
    local_ME = local_ME.array().rowwise() * V_vec.transpose().array();
    local_ME = local_ME * wqua_diag_ref * test_ref1D_T_;
    for (int test_basis_index = 0; test_basis_index < polydim_g; test_basis_index++) 
    {
      int alpha = Tm_g(test_basis_index, i);
      for (int trial_basis_index = 0; trial_basis_index < polydim_rho; trial_basis_index++) 
      {
        int beta = Tm_rho(trial_basis_index, x_cell_index);
        real_t ini_value = local_ME(test_basis_index, trial_basis_index) * JacobiDet_g(0);
        tripletList_ME.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };
  ME->setFromTriplets(tripletList_ME.begin(), tripletList_ME.end());
};

void KineticDD_DG2d_IMEX_IM_Schur::FourierMatrix(const real_t& xita, std::string& OutfilePath)
{

};

void KineticDD_DG2d_IMEX_IM_Schur::FourierDMatrix_compute(const real_t& be, const real_t& xita, cMatrix* D)
{

};

void KineticDD_DG2d_IMEX_IM_Schur::FourierUMatrix_compute(const cMatrix& D_NegativeWind,
                                    const cMatrix& D_PositiveWind,
                                    const real_t& xita, 
                                    cMatrix* U)
{

};

void KineticDD_DG2d_IMEX_IM_Schur::fluxint_upwind_compute(const Matrix& modal, 
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

void KineticDD_DG2d_IMEX_IM_Schur::fluxext_upwind_compute(const Matrix& modal, 
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

void KineticDD_DG2d_IMEX_IM_Schur::ah_compute(const model_data_& modal, 
                                              Matrix* ah)
{
  // ** 传入相关变量
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const Matrix& test_ref_rho = fe1D_->gettest_ref();
  const Matrix& test_ref_dx_rho = fe1D_->gettest_dx_ref();
  const int& ncell_g = kinetic_mesh2D_->getncell();
  const int& polydim_g = kinetic_fe2D_->getbasis()->getpolydim();
  const int& NTdofs_g = kinetic_fe2D_->getNTdofs();
  const int& numqua_g = kinetic_fe2D_->getNTdofs();
  // *** 

  ah->resize(polydim_rho, ncell_rho);
  ah->setZero();
  Eigen::Map<Vector> ah_vec(ah->data(), NTdofs_rho);
  Eigen::Map<Vector> g_vec(modal.g, NTdofs_g);
  ah_vec.setZero();

  ah_vec = Da_ * g_vec;
};

void KineticDD_DG2d_IMEX_IM_Schur::ah_bc_compute(const model_data_& modal, 
                                                  Matrix* ah_bc)
{
  // ** 传入相关变量
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  // *** 
  ah_bc->resize(polydim_rho, ncell_rho);
  ah_bc->setZero();
  Eigen::Map<Vector> ah_bc_vec(ah_bc->data(), NTdofs_rho);
  Eigen::Map<Vector> rho_vec(modal.rho, NTdofs_rho);
  ah_bc_vec = Mbc_ * rho_vec;
};

void KineticDD_DG2d_IMEX_IM_Schur::ah_bc_inflow_compute(const model_data_& modal,
                                                        const real_t& Trun,
                                                        Matrix* ah_bc_inflow)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell_rho = mesh1D_->getncell();
  const int& polydim_rho = fe1D_->getbasis()->getpolydim();
  const int& NTdofs_rho = fe1D_->getNTdofs();
  const int& numqua_rho = fe1D_->getnumqua();
  const std::vector<real_t>&  extlength = 
    kinetic_mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& extnormals = 
    kinetic_mesh2D_->getextboundarynormal();
  const std::vector<int>& extNei =
    kinetic_mesh2D_->getextboundaryneighbors();
  const int& extboundaryNum =
    kinetic_mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBTypeIndex =
    kinetic_mesh2D_->getextboundarytypeindex();
  const std::vector<Eigen::Vector2d>& extboundarycenter = 
    kinetic_mesh2D_->getextboundarycenter();
  const Matrix& CellCenter = 
    kinetic_mesh2D_->getCellCenter();
  const std::vector<Matrix>& CoorBdrRef2d = 
    kinetic_fe2D_->getCoorBdrRef();
  // *** 
  ah_bc_inflow->resize(polydim_rho, ncell_rho);
  ah_bc_inflow->setZero();
  Eigen::Map<Vector> ah_bc_vec(ah_bc->data(), NTdofs_rho);

  Vector V_vec(numqua_rho), f_vec(numqua_rho);
  real_t normal;
  Vector temp(polydim_rho);
  real_t xbc;
  for (int i = 0; i < extboundaryNum; i++) 
  {
    normal = extnormals[i](0);
    len = extlength[i];
    if (std::abs(normal) > 0.5e0)
    {
      test_cell_Index = ExtNei[i];
      test_x_cell_index = kinetic_mesh2D_->xCellIndex(test_cell_index);
      if (test_x_cell_index == 0)
      {
        test_qua_value = boundary_u_rho_[1];
        xbc = x1;
      } else if (test_x_cell_index == Nx - 1)
      {
        test_qua_value = boundary_u_rho_[0];
        xbc = x2;
      }
      test_normal = normal;

      for (int q = 0; q < numqua_rho; q++)
      {
        V_vec(q) = CellCenter(1, extNei[i]) 
                  + CoorBdrRef2d[ExtBTypeIndex[i]](1, q) * hv;
        f_vec(q) = fin_bc(xbc, V_vec(q), Trun);
      }
      V_vec = ((V_vec.array() * normal) <= 0).cast<double>();
      temp = test_qua_value * wqua_diag1d_ref * f_vec * len;

      for (int test_basis_index = 0; test_basis_index < polydim_rho; test_basis_index++) 
      {
        alpha = Tm_rho(test_basis_index, test_x_cell_index);
        ini_value = temp(test_basis_index) * CR_;
      };
    };
  };
};

void KineticDD_DG2d_IMEX_IM_Schur::bh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          const std::vector<Matrix>& boundary_flux,
                          Matrix* bh) 
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

void KineticDD_DG2d_IMEX_IM_Schur::bh_extflux_compute(const model_data_& modal, 
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

void KineticDD_DG2d_IMEX_IM_Schur::dh_compute(const Matrix& rho, 
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

void KineticDD_DG2d_IMEX_IM_Schur::dh_bc_compute(const Matrix& rho, 
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

void KineticDD_DG2d_IMEX_IM_Schur::mh_compute(const Matrix& rho, 
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

void KineticDD_DG2d_IMEX_IM_Schur::sh_compute(const model_data_& modal,  
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

void KineticDD_DG2d_IMEX_IM_Schur::ch_compute(const model_data_& modal, 
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
        wind = modal.E(0, i) >= 0 ? -1 : 1;
        switch (wind)
        {
          case -1:
            g_diff_nodal.col(i) = (flux_downwind_nodal[j+1].col(i) - 
                                  flux_downwind_nodal[j].col(i)) * hvinv;
            
            break;
          case 1:
            g_diff_nodal.col(i) = (flux_upwind_nodal[j+1].col(i) - 
                                  flux_upwind_nodal[j].col(i)) * hvinv;
            break;
        }
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

void KineticDD_DG2d_IMEX_IM_Schur::WENO3_VelocityReconstruct(
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
  flux_upwind_nodal->at(0) = Matrix::Zero(numqua, ncell);
  flux_downwind_nodal->at(0) = Matrix::Zero(numqua, ncell);
  flux_upwind_nodal->at(Nv) = Matrix::Zero(numqua, ncell);
  flux_downwind_nodal->at(Nv) = Matrix::Zero(numqua, ncell);

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

void KineticDD_DG2d_IMEX_IM_Schur::sourceh_compute(const model_data_& modal, 
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

void KineticDD_DG2d_IMEX_IM_Schur::updateAll(const real_t& Trun, 
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

void KineticDD_DG2d_IMEX_IM_Schur::SolveRho_DichletBoundary(const real_t& Trun,
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

void KineticDD_DG2d_IMEX_IM_Schur::SolveG(const real_t& Trun,
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

real_t KineticDD_DG2d_IMEX_IM_Schur::rho_init(const real_t& x) 
{ 
  return 1.e0;
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::rho_init(const Matrix& x) 
{
  Matrix rho = x;
  rho.setConstant(1.e0);
  return rho;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::rho_d(const real_t& x) 
{
  return 1.e0 - (1.e0 - 0.001)/2.e0 * (std::tanh((x - 0.3e0)/0.02e0) - std::tanh((x - 0.7e0)/0.02e0));
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::rho_d(const Matrix& x) 
{
  Matrix x1 = (x.array() - 0.3) / 0.02;
  Matrix x2 = (x.array() - 0.7) / 0.02;
  Matrix rhod = 1.0 - (1.0 - 0.001) / 2.0 * (x1.array().tanh() - x2.array().tanh());
  return rhod;
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::V_value(const Matrix& x, const Matrix& v)
{
  return v;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::V_value(const real_t& x, const real_t& v)
{
  return v;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::source(const real_t& x, const real_t& t) 
{
  return 0.e0;
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::source(const Matrix& x, const real_t& t) 
{
  return Matrix::Zero(x.rows(), x.cols());
};

real_t KineticDD_DG2d_IMEX_IM_Schur::fsource(const real_t& x, const real_t& v, const real_t& t) 
{
  return 0.e0;
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::fsource(const Matrix& x, const Matrix& v, const real_t& t) 
{
  return Matrix::Zero(x.rows(), x.cols());
};

real_t KineticDD_DG2d_IMEX_IM_Schur::f_init(const real_t& x, const real_t& v) 
{
  return Maxwell(v);
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::f_init(const Matrix& x, const Matrix& v) 
{
  Matrix f = x;
  f.setConstant(1.e0);
  f = f.array() * Maxwell(v).array();
  return f;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::g_init(const real_t& x, const real_t& v) 
{
  return 0.e0;
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::g_init(const Matrix& x, const Matrix& v) 
{
  Matrix g = x;
  g.setZero();
  return g;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::fin_bc(const real_t& x, 
                                            const real_t& v, 
                                            const real_t& t)
{
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();

  real_t temp;
  if (std::abs(x - x1) < 1.e-13 && v >= 0)
  {
    temp = Maxwell(v);
  } else if (std::abs(x2) < 1.e-13 && v < 0)
  {
    temp = Maxwell(v);
  };
};

Matrix KineticDD_DG2d_IMEX_IM_Schur::Maxwell(const Matrix& v) {
  Matrix M = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) *
                  (- v.array() * v.array() / (2 * theta_)).exp();
  M = M / Maxwell_sum_;
  return M;
};

Vector KineticDD_DG2d_IMEX_IM_Schur::Maxwell(const Vector& v) {
  Vector M = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) *
                  (- v.array() * v.array() / (2 * theta_)).exp();
  M = M / Maxwell_sum_;
  return M;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::Maxwell(const real_t& v) {
  real_t M = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) * std::exp( - v * v / ( 2 * theta_));
  M = M / Maxwell_sum_;
  return M;
};

real_t KineticDD_DG2d_IMEX_IM_Schur::phi_bc(const real_t& x, const real_t& t)
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


}; // namespace QUEST

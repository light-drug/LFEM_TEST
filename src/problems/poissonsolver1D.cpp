#include "poissonsolver1D.hpp"

namespace QUEST
{

PoissonSolver1D::PoissonSolver1D(const fespace1D* fe, const PoissonSolver1DParameter& pa)
  : fe_(fe), pa_(pa) {};

void PoissonSolver1D::generateB() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<real_t>& intnormals_ = 
    fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei_ = 
    fe_->getmesh1D()->getintboundaryneighbors();
  // ******************************** //
  B_.resize(NTdofs_ * dim_, NTdofs_);
  int estimatedNonZeros = ncell_ * polydim_ * polydim_ + 
          intboundaryNum_ * polydim_ * polydim_ * 4;
  std::vector<Eigen::Triplet<real_t>> tripletList_B;
  tripletList_B.reserve(estimatedNonZeros);
  
  Matrix Bref;
  Bref = dv_u_;
  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet_(i) * Jx_(i);
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };
  for (int i = 0; i < intboundaryNum_; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = intNei_[i](test_cell);
      Vector test_qua_value = boundary_u_[test_cell];
      real_t test_normal = intnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = intNei_[i](trial_cell);
        Vector trial_qua_value = boundary_u_[trial_cell];
        real_t trial_normal = intnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
                + trial_qua_value(trial_basis_index) * pa_.C12 * trial_normal) 
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = - ini_value;   // 这是一维代码不需要积分哈哈哈
            tripletList_B.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  B_.setFromTriplets(tripletList_B.begin(), tripletList_B.end());
};

void PoissonSolver1D::generateC() {

  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<real_t>& intnormals_ = 
    fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei_ = 
    fe_->getmesh1D()->getintboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  const std::vector<int>& extNei_ = 
    fe_->getmesh1D()->getextboundaryneighbors();
  // ******************************** //

  C_.resize(NTdofs_ * dim_, NTdofs_);
  int estimatedNonZeros = intboundaryNum_ * polydim_ * polydim_ * 4
    + extboundaryNum_ * polydim_ * polydim_;
  std::vector<Eigen::Triplet<real_t>> tripletList_C;
  tripletList_C.reserve(estimatedNonZeros);

  for (int i = 0; i < intboundaryNum_; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = intNei_[i](test_cell);
      Vector test_qua_value = boundary_u_[test_cell];
      real_t test_normal = intnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = intNei_[i](trial_cell);
        Vector trial_qua_value = boundary_u_[trial_cell];
        real_t trial_normal = intnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value =  (pa_.C11 * trial_qua_value(trial_basis_index)) 
                * test_qua_value(test_basis_index) * test_normal * trial_normal;
            tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_; i++) {

    int test_cell_Index = extNei_[i];
    Vector test_qua_value = boundary_u_[1-i];
    real_t test_normal = extnormals_[i];

    int trial_cell_Index = extNei_[i];
    Vector trial_qua_value = boundary_u_[1-i];
    real_t trial_normal = extnormals_[i];

    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      int alpha = Tm_(test_basis_index, test_cell_Index);
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        int beta = Tm_(trial_basis_index, trial_cell_Index);
        real_t ini_value =  (pa_.C11 * trial_qua_value(trial_basis_index)) 
            * test_qua_value(test_basis_index) * test_normal * trial_normal;
        tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };

  C_.setFromTriplets(tripletList_C.begin(), tripletList_C.end());
};

void PoissonSolver1D::generateK() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  // ******************************** //

  K_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_K;
  int estimatedNonZeros = NTdofs_;
  tripletList_K.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u_(basis_index,basis_index) * JacobiDet_(i);
      int alpha = alpha_start + basis_index;
      tripletList_K.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  K_.setFromTriplets(tripletList_K.begin(), tripletList_K.end());
};

void PoissonSolver1D::generateKinv() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  // ******************************** //

  Kinv_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_K;
  int estimatedNonZeros = NTdofs_;
  tripletList_K.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u_(basis_index,basis_index) * JacobiDet_(i);
      ini_value = 1.e0 / ini_value;
      int alpha = alpha_start + basis_index;
      tripletList_K.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  Kinv_.setFromTriplets(tripletList_K.begin(), tripletList_K.end());
};

void PoissonSolver1D::generateS() {
  BT_ = B_.transpose();
  S_ = C_ + BT_ * Kinv_ * B_;
  std::cout << " Poisson Discrete Matrix S has been set up !! " << std::endl;
}

void PoissonSolver1D::init(const PoissonSolver1DType& poitype,
      const real_t& tol) {

  // ********* 传入相关变量 ********** //
  intboundaryNum_ = fe_->getmesh1D()->getintboundaryNum();
  extboundaryNum_ = fe_->getmesh1D()->getextboundaryNum();
  num_equations_ = fe_->getnum_equations();
  polydim_ = fe_->getbasis()->getpolydim();
  ncell_ = fe_->getmesh1D()->getncell();
  dim_ = fe_->getdim();
  numqua_ = fe_->getnumqua();
  hx_ = fe_->getmesh1D()->gethx();
  // ******************************** //
  generateB();
  generateK();
  generateC();
  generateKinv();
  generateS();
  
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::CG:
    // std::cout << "  Preparation for Conjugate Gradient the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    cg_.compute(S_);
    cg_.setTolerance(tol);
    // std::cout << "  The end of preparation for Conjugate Gradient method ! " << std::endl;
    break;

  case PoissonSolver1DType::PCG:
    QUEST_ERROR(" The PCG method is not been implemented ! ");
    break;
  
  case PoissonSolver1DType::LDLT:
// #ifndef QUEST_USE_MKL
//     std::cout << "  Preparation for LDLT the matrix by Eigen ......\n";
// #else 
//     std::cout << "  Preparation for LDLT the matrix by Pardiso ......\n";
// #endif 
//     std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    chol_.analyzePattern(S_);
    chol_.factorize(S_);
    // std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  TOC;
}

void PoissonSolver1D::solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const {
  All->resize(2);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::CG:
    std::cout << "  Solve the Poisson by CG ......\n";
    xs = cg_.solve(f1);
    std::cout << "  iterations = " << cg_.iterations() << std::endl;
    std::cout << "  eerror = " << cg_.error()      << std::endl;
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DType::LDLT:
    std::cout << "  Solve the Poisson by LDLT ......\n";
    xs = chol_.solve(f1);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;

  Vector f2 = bg - B_ * xs;
  Vector xsd = Kinv_ * f2;
  All->at(0) = xs;
  All->at(1) = xsd;
}

void PoissonSolver1D::solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Matrix>* All) const 
{ 
  All->resize(2);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::CG:
    std::cout << "  Solve the Poisson by CG ......\n";
    xs = cg_.solve(f1);
    std::cout << "  iterations = " << cg_.iterations() << std::endl;
    std::cout << "  eerror = " << cg_.error()      << std::endl;
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DType::LDLT:
    std::cout << "  Solve the Poisson by LDLT ......\n";
    xs = chol_.solve(f1);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;

  Vector f2 = bg - B_ * xs;
  Vector xsd = Kinv_ * f2;
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdtemp(xsd.data(), polydim_, ncell_);
  All->at(0) = xstemp;
  All->at(1) = xsdtemp;
};

// -\Delta \Phi = f_nodal, E = - \nabla \Phi
void PoissonSolver1D::solveall(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Matrix* E,
                Matrix* Phi) const
{
  Vector bg;
  Vector bf;
  Assemble_bg(Dirichlet, &bg);
  Assemble_bf(Dirichlet, f_nodal, &bf);

  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::CG:
    // std::cout << "  Solve the Poisson by CG ......\n";
    xs = cg_.solve(f1);
    // std::cout << "  iterations = " << cg_.iterations() << std::endl;
    // std::cout << "  eerror = " << cg_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DType::LDLT:
    // std::cout << "  Solve the Poisson by LDLT ......\n";
    xs = chol_.solve(f1);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;

  Vector f2 = bg - B_ * xs;
  Vector xsd = Kinv_ * f2;
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdtemp(xsd.data(), polydim_, ncell_);
  *E = - xsdtemp;
  *Phi = xstemp;
};


void PoissonSolver1D::Assemble_bg(const Matrix& Dirichlet,
                Vector* bg) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<int>& extNei_ = 
    fe_->getmesh1D()->getextboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  // ******************************** //

  bg->resize(NTdofs_);
  bg->setZero();

  for (int i = 0; i < extboundaryNum_; i++) {

    int test_cell_Index = extNei_[i];
    Vector test_qua_value = boundary_u_[1-i];
    real_t test_normal = extnormals_[i];

    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      int alpha = Tm_(test_basis_index, test_cell_Index);
      real_t ini_value =  Dirichlet(0, i)
          * test_qua_value(test_basis_index) * test_normal;
      (*bg)(alpha) += ini_value;
    };
  };
};


void PoissonSolver1D::Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<int>& extNei_ = 
    fe_->getmesh1D()->getextboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  // ******************************** //

  bf->resize(NTdofs_);
  bf->setZero();
  Matrix bf_temp;
  fe_->Assemble_F(f_nodal, 0, &bf_temp);
  Eigen::Map<Vector> reshapebf(bf_temp.data(), NTdofs_, 1);
  *bf = reshapebf * JacobiDet_(0);

  for (int i = 0; i < extboundaryNum_; i++) {

    int test_cell_Index = extNei_[i];
    Vector test_qua_value = boundary_u_[1-i];
    real_t test_normal = extnormals_[i];

    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      int alpha = Tm_(test_basis_index, test_cell_Index);
      real_t ini_value =  pa_.C11 * Dirichlet(0, i)
          * test_qua_value(test_basis_index);
      (*bf)(alpha) += ini_value;
    };
  };

};

// period部分
PoissonSolver1D_period::PoissonSolver1D_period(const fespace1D* fe, const PoissonSolver1DParameter& pa)
  : PoissonSolver1D(fe, pa) {};

void PoissonSolver1D_period::init(const PoissonSolver1DType& poitype,
      const real_t& tol) {

  // ********* 传入相关变量 ********** //
  intboundaryNum_ = fe_->getmesh1D()->getintboundaryNum();
  extboundaryNum_ = fe_->getmesh1D()->getextboundaryNum();
  num_equations_ = fe_->getnum_equations();
  polydim_ = fe_->getbasis()->getpolydim();
  ncell_ = fe_->getmesh1D()->getncell();
  dim_ = fe_->getdim();
  numqua_ = fe_->getnumqua();
  hx_ = fe_->getmesh1D()->gethx();
  // ******************************** //
  QUEST_VERIFY(fe_->getmesh1D()->IsPeriodBoundary(), "The mesh is not periodical !");
  generateB();
  generateBT();
  generateK();
  generateC();
  generateKinv();
  generateS();
  
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::LU:
    // std::cout << "  Preparation for Conjugate Gradient the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.compute(S_);
    // std::cout << "  The end of preparation for Conjugate Gradient method ! " << std::endl;
    break;

  case PoissonSolver1DType::GMRES:
    gmres_.setTolerance(tol);
    gmres_.set_restart(1000);
    gmres_.setMaxIterations(ncell_ * polydim_);
    gmres_.compute(S_);
    // std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  TOC;
};


void PoissonSolver1D_period::generateB() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<real_t>& intnormals_ = 
    fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei_ = 
    fe_->getmesh1D()->getintboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  const std::vector<Eigen::Vector2i>& extNei_period_ = 
    fe_->getmesh1D()->getextboundaryneighbors_period();
  // ******************************** //
  B_.resize(NTdofs_ * dim_, NTdofs_);
  int estimatedNonZeros = ncell_ * polydim_ * polydim_ + 
          intboundaryNum_ * polydim_ * polydim_ * 4 + 
          extboundaryNum_ / 2 * polydim_ * polydim_ * 4;
  std::vector<Eigen::Triplet<real_t>> tripletList_B;
  tripletList_B.reserve(estimatedNonZeros);
  
  Matrix Bref;
  Bref = dv_u_;
  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet_(i) * Jx_(i);
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };
  for (int i = 0; i < intboundaryNum_; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = intNei_[i](test_cell);
      const Vector& test_qua_value = boundary_u_[test_cell];
      real_t test_normal = intnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = intNei_[i](trial_cell);
        const Vector& trial_qua_value = boundary_u_[trial_cell];
        real_t trial_normal = intnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
                + trial_qua_value(trial_basis_index) * pa_.C12 * trial_normal) 
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = - ini_value;   // 这是一维代码不需要积分哈哈哈
            tripletList_B.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_ / 2; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = extNei_period_[i](test_cell);
      const Vector& test_qua_value = boundary_u_[1-test_cell];
      real_t test_normal = extnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = extNei_period_[i](trial_cell);
        const Vector& trial_qua_value = boundary_u_[1-trial_cell];
        real_t trial_normal = extnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) 
                + trial_qua_value(trial_basis_index) * pa_.C12 * trial_normal) 
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = - ini_value;   // 这是一维代码不需要积分哈哈哈
            tripletList_B.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  B_.setFromTriplets(tripletList_B.begin(), tripletList_B.end());
};

void PoissonSolver1D_period::generateBT() {
  BT_ = B_.transpose();
};

void PoissonSolver1D_period::generateC() 
{
// ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<real_t>& intnormals_ = 
    fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei_ = 
    fe_->getmesh1D()->getintboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  const std::vector<Eigen::Vector2i>& extNei_period_ = 
    fe_->getmesh1D()->getextboundaryneighbors_period();
  // ******************************** //

  C_.resize(NTdofs_ * dim_, NTdofs_);
  int estimatedNonZeros = intboundaryNum_ * polydim_ * polydim_ * 4
    + extboundaryNum_ / 2 * polydim_ * polydim_ * 4 + ncell_;
  std::vector<Eigen::Triplet<real_t>> tripletList_C;
  tripletList_C.reserve(estimatedNonZeros);

  for (int i = 0; i < intboundaryNum_; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = intNei_[i](test_cell);
      Vector test_qua_value = boundary_u_[test_cell];
      real_t test_normal = intnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = intNei_[i](trial_cell);
        Vector trial_qua_value = boundary_u_[trial_cell];
        real_t trial_normal = intnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value =  (pa_.C11 * trial_qua_value(trial_basis_index)) 
                * test_qua_value(test_basis_index) * test_normal * trial_normal;
            tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_ / 2; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = extNei_period_[i](test_cell);
      const Vector& test_qua_value = boundary_u_[1-test_cell];
      real_t test_normal = extnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = extNei_period_[i](trial_cell);
        const Vector& trial_qua_value = boundary_u_[1-trial_cell];
        real_t trial_normal = extnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value =  (pa_.C11 * trial_qua_value(trial_basis_index)) 
                * test_qua_value(test_basis_index) * test_normal * trial_normal;
            tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  const real_t vol = fe_->getmesh1D()->getCellVol()(0);
  int center = ncell_ / 2;
  int alpha = Tm_(0, center);
  real_t value = vol;
  for (int i = 0; i < ncell_; i++) {
    int beta = Tm_(0, i);
    tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, vol));
  }

  C_.setFromTriplets(tripletList_C.begin(), tripletList_C.end());
};

void PoissonSolver1D_period::generateS() {
  S_ = C_ + BT_ * Kinv_ * B_;
  std::cout << " Poisson Discrete Matrix S has been set up !! " << std::endl;
};

void PoissonSolver1D_period::Assemble_bg(const Matrix& Dirichlet,
                Vector* bg) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<int>& extNei_ = 
    fe_->getmesh1D()->getextboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  // ******************************** //

  bg->resize(NTdofs_);
  bg->setZero();
};

void PoissonSolver1D_period::Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<int>& extNei_ = 
    fe_->getmesh1D()->getextboundaryneighbors();
  const std::vector<real_t>& extnormals_ = 
    fe_->getmesh1D()->getextboundarynormal();
  // ******************************** //

  bf->resize(NTdofs_);
  bf->setZero();
  Matrix bf_temp;
  fe_->Assemble_F(f_nodal, 0, &bf_temp);
  Eigen::Map<Vector> reshapebf(bf_temp.data(), NTdofs_, 1);
  *bf = reshapebf * JacobiDet_(0);
};

void PoissonSolver1D_period::solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const {
  All->resize(2);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::LU:
    std::cout << "  Solve the Poisson by LU ......\n";
    xs = lu_.solve(f1);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DType::GMRES:
    std::cout << "  Solve the Poisson by LDLT ......\n";
    xs = gmres_.solve(f1);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;

  Vector f2 = bg - B_ * xs;
  Vector xsd = Kinv_ * f2;
  All->at(0) = xs;
  All->at(1) = xsd;
}

void PoissonSolver1D_period::solveall(const Vector& bg,
                                      const Vector& bf,
                                      std::vector<Matrix>* All) const 
{ 
  All->resize(2);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::LU:
    std::cout << "  Solve the Poisson by LU ......\n";
    xs = lu_.solve(f1);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DType::GMRES:
    std::cout << "  Solve the Poisson by GMRES ......\n";
    xs = gmres_.solve(f1);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;

  Vector f2 = bg - B_ * xs;
  Vector xsd = Kinv_ * f2;
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdtemp(xsd.data(), polydim_, ncell_);
  All->at(0) = xstemp;
  All->at(1) = xsdtemp;
};

// -\Delta \Phi = f_nodal, E = - \nabla \Phi
void PoissonSolver1D_period::solveall(const Matrix& Dirichlet,
                                      const Matrix& f_nodal,
                                      Matrix* E,
                                      Matrix* Phi) const
{
  Vector bg;
  Vector bf;
  Assemble_bg(Dirichlet, &bg);
  Assemble_bf(Dirichlet, f_nodal, &bf);

  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DType::LU:
    // std::cout << "  Solve the Poisson by LU ......\n";
    xs = lu_.solve(f1);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DType::GMRES:
    // std::cout << "  Solve the Poisson by GMRES ......\n";
    xs = gmres_.solve(f1);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;

  Vector f2 = bg - B_ * xs;
  Vector xsd = Kinv_ * f2;
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdtemp(xsd.data(), polydim_, ncell_);
  *E = - xsdtemp;
  *Phi = xstemp;

};

Poisson_acctest_1D::Poisson_acctest_1D(const fespace1D* fe, const PoissonSolver1DParameter& pa)
  : PoissonSolver1D(fe, pa) {};

Matrix Poisson_acctest_1D::RHS(const Matrix& x) 
{
  // Matrix rhs = - ( x.array().exp() + 2.e0 );
  Matrix rhs = - ( x.array().sin());
  return rhs;
};

real_t Poisson_acctest_1D::u_bc(const real_t& x) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x);
};

Matrix Poisson_acctest_1D::u_real(const Matrix& x) 
{
  // return x.array().exp() + x.array().square();
  return - x.array().sin();
};

real_t Poisson_acctest_1D::u_real(const real_t& x) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x);
};

real_t Poisson_acctest_1D::u_real_dx(const real_t& x) 
{
  // return std::exp(x) + 2.e0 * x;
  return - std::cos(x);
};

Matrix Poisson_acctest_1D::u_real_dx(const Matrix& x) 
{
  // return x.array().exp() + 2.e0 * x.array();
  return - x.array().cos();
};

Poisson_acctest_1D_period::Poisson_acctest_1D_period(const fespace1D* fe, const PoissonSolver1DParameter& pa)
  : PoissonSolver1D_period(fe, pa) {};

Matrix Poisson_acctest_1D_period::RHS(const Matrix& x) 
{
  // Matrix rhs = - ( x.array().exp() + 2.e0 );
  Matrix rhs = - ( x.array().sin());
  return rhs;
};

real_t Poisson_acctest_1D_period::u_bc(const real_t& x) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x);
};

Matrix Poisson_acctest_1D_period::u_real(const Matrix& x) 
{
  // return x.array().exp() + x.array().square();
  return - x.array().sin();
};

real_t Poisson_acctest_1D_period::u_real(const real_t& x) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x);
};

real_t Poisson_acctest_1D_period::u_real_dx(const real_t& x) 
{
  // return std::exp(x) + 2.e0 * x;
  return - std::cos(x);
};

Matrix Poisson_acctest_1D_period::u_real_dx(const Matrix& x) 
{
  // return x.array().exp() + 2.e0 * x.array();
  return - x.array().cos();
};


} // namespace QUEST

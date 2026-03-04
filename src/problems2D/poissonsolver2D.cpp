#include "poissonsolver2D.hpp"

namespace QUEST
{

PoissonSolver2D::PoissonSolver2D(const fespace2D* fe, const PoissonSolver2DParameter& pa)
  : fe_(fe), pa_(pa) {};

void PoissonSolver2D::generateB() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const real_t& hx = fe_->getmesh2D()->gethx();
  const real_t& hy = fe_->getmesh2D()->gethy();
  const Matrix& test_ref = fe_->gettest_ref();
  const Matrix& test_dx_ref = fe_->gettest_dx_ref();
  const Matrix& test_dy_ref = fe_->gettest_dy_ref();
  const Matrix& dxv_u = fe_->getdxv_u();
  const Matrix& dyv_u = fe_->getdyv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh2D()->getJx();
  const Vector& Jy = fe_->getmesh2D()->getJy();
  const std::vector<Matrix>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intlength = 
    fe_->getmesh2D()->getintboundarylength();
  const std::vector<Eigen::Vector2d>& intnormals = 
    fe_->getmesh2D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    fe_->getmesh2D()->getintboundaryneighbors();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = 
    fe_->getmesh2D()->getintboundarytypeindex();
  const std::vector<std::vector<Matrix>>& flux_u_v = 
    fe_->getflux_u_v();
  // ******************************** //
  B_.resize(NTdofs * dim_, NTdofs);
  int estimatedNonZeros = ncell_ * polydim_ * polydim_ + 
          intboundaryNum_ * polydim_ * polydim_ * 4;
  estimatedNonZeros = estimatedNonZeros * dim_;
  std::vector<Eigen::Triplet<real_t>> tripletList_B;
  tripletList_B.reserve(estimatedNonZeros);
  // std::cout << " v_u = \n" << v_u << std::endl;
  // std::cout << " test_ref = \n" << test_ref << std::endl;
  // std::cout << " test_dx_ref = \n" << test_dx_ref << std::endl;
  // std::cout << " test_dy_ref = \n" << test_dy_ref << std::endl;
  // std::cout << " dxv_u = \n" << dxv_u << std::endl;
  // std::cout << " dyv_u = \n" << dyv_u << std::endl;
  // std::cout << " hx = " << hx << std::endl;
  // std::cout << " hy = " << hy << std::endl;
  // std::cout << " JacobiDet = \n" << JacobiDet.transpose() << std::endl;
  // std::cout << " Jx = \n" << Jx.transpose() << std::endl;
  // PAUSE();
  Matrix Bref;
  Bref = dxv_u;
  for (int i = 0; i < ncell_; i++) 
  {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
    {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
      {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jx(i);
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };
  Bref = dyv_u;
  for (int i = 0; i < ncell_; i++) 
  {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
    {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
      {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jy(i);
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index + NTdofs, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };

  for (int d = 0; d < 2; d++)
  {
    Eigen::Vector2d v;
    v.setZero();
    v(d) = 1.e0;
    Matrix temp;
    for (int i = 0; i < intboundaryNum_; i++) 
    {
      real_t len = intlength[i];
      for (int test_cell = 0; test_cell < 2; test_cell++) 
      {
        int test_cell_Index = intNei[i](test_cell);
        // Matrix test_qua_value = boundary_u[IntBTypeIndex[i](test_cell)];
        Eigen::Vector2d test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
        for (int trial_cell = 0; trial_cell < 2; trial_cell++) 
        {
          int trial_cell_Index = intNei[i](trial_cell);
          // Matrix trial_qua_value = boundary_u[IntBTypeIndex[i](trial_cell)];
          Eigen::Vector2d trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);

          temp = flux_u_v[IntBTypeIndex[i](test_cell)][IntBTypeIndex[i](trial_cell)] * len;
          // std::cout << " temp = \n" << flux_u_v[IntBTypeIndex[i](test_cell)][IntBTypeIndex[i](trial_cell)] << std::endl;
          // PAUSE();
          for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
          {
            int alpha = Tm(test_basis_index, test_cell_Index) + d * NTdofs;
            for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
            {
              int beta = Tm(trial_basis_index, trial_cell_Index);
              real_t ini_value =  (0.5e0 + pa_.C12.dot(trial_normal)) 
                  * temp(test_basis_index, trial_basis_index) * v.dot(test_normal);
              ini_value = - ini_value;
              // std::cout << "alpha = " << alpha << std::endl;
              // std::cout << "beta = " << beta << std::endl;
              // std::cout << "length = " << len << std::endl;
              // std::cout << "ini_value = " << ini_value / len << std::endl;
              if ( std::abs(ini_value) > 1.e-14 )
              {
                tripletList_B.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
              };
            };
          };
        };
      };
    };
  };
  
  B_.setFromTriplets(tripletList_B.begin(), tripletList_B.end());
};

void PoissonSolver2D::generateC() {

  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh2D()->getJx();
  const Vector& Jy = fe_->getmesh2D()->getJy();
  const std::vector<Matrix>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intlength = 
    fe_->getmesh2D()->getintboundarylength();
  const std::vector<Eigen::Vector2d>& intnormals = 
    fe_->getmesh2D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    fe_->getmesh2D()->getintboundaryneighbors();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = 
    fe_->getmesh2D()->getintboundarytypeindex();

  const std::vector<real_t>& extlength = 
    fe_->getmesh2D()->getintboundarylength();
  const std::vector<Eigen::Vector2d>& extnormals = 
    fe_->getmesh2D()->getextboundarynormal();
  const std::vector<int>& extNei = 
    fe_->getmesh2D()->getextboundaryneighbors();
  const std::vector<int>& ExtBTypeIndex = 
    fe_->getmesh2D()->getextboundarytypeindex();
  const std::vector<std::vector<Matrix>>& flux_u_v = 
    fe_->getflux_u_v();
  // ******************************** //

  C_.resize(NTdofs, NTdofs);
  int estimatedNonZeros = intboundaryNum_ * polydim_ * polydim_ * 4
    + extboundaryNum_ * polydim_ * polydim_;
  std::vector<Eigen::Triplet<real_t>> tripletList_C;
  tripletList_C.reserve(estimatedNonZeros);

  Matrix temp;
  for (int i = 0; i < intboundaryNum_; i++) 
  {
    real_t len = intlength[i];
    for (int test_cell = 0; test_cell < 2; test_cell++) 
    {
      int test_cell_Index = intNei[i](test_cell);
      Eigen::Vector2d test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) 
      {
        int trial_cell_Index = intNei[i](trial_cell);
        Eigen::Vector2d trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);

        temp = flux_u_v[IntBTypeIndex[i](test_cell)][IntBTypeIndex[i](trial_cell)] * len;
        
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
        {
          int alpha = Tm(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
          {
            int beta = Tm(trial_basis_index, trial_cell_Index);
            real_t ini_value =  pa_.C11 * temp(test_basis_index, trial_basis_index)
                              * test_normal.dot(trial_normal);
            if ( std::abs(ini_value) > 1.e-14 )
            {
              tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
            };
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_; i++) 
  {
    real_t len = extlength[i];

    int test_cell_Index = extNei[i];
    Eigen::Vector2d test_normal = extnormals[i];

    int trial_cell_Index = extNei[i];
    Eigen::Vector2d trial_normal = extnormals[i];

    temp = flux_u_v[ExtBTypeIndex[i]][ExtBTypeIndex[i]] * len;
    // std::cout << " len = " << len << std::endl;
    // std::cout << " temp = \n" << flux_u_v[ExtBTypeIndex[i]][ExtBTypeIndex[i]] << std::endl;
    // PAUSE();
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
    {
      int alpha = Tm(test_basis_index, test_cell_Index);
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
      {
        int beta = Tm(trial_basis_index, trial_cell_Index);
        real_t ini_value =  pa_.C11 * temp(test_basis_index, trial_basis_index);
        if ( std::abs(ini_value) > 1.e-14 )
        {
          tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
        }
      };
    };
  };

  C_.setFromTriplets(tripletList_C.begin(), tripletList_C.end());
};

void PoissonSolver2D::generateK() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const IntMatrix& Tm = fe_->getTm();
  // ******************************** //

  K_.resize(NTdofs * dim_, NTdofs * dim_);
  std::vector<Eigen::Triplet<real_t>> tripletList_K;
  int estimatedNonZeros = NTdofs * dim_;
  tripletList_K.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u(basis_index,basis_index) * JacobiDet(i);
      int alpha = alpha_start + basis_index;
      tripletList_K.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
      alpha += NTdofs;
      tripletList_K.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  K_.setFromTriplets(tripletList_K.begin(), tripletList_K.end());
};

void PoissonSolver2D::generateKinv() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const IntMatrix& Tm = fe_->getTm();
  // ******************************** //

  Kinv_.resize(NTdofs * dim_, NTdofs * dim_);
  std::vector<Eigen::Triplet<real_t>> tripletList_Kinv;
  int estimatedNonZeros = NTdofs * dim_;
  tripletList_Kinv.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u(basis_index,basis_index) * JacobiDet(i);
      ini_value = 1.e0 / ini_value;
      int alpha = alpha_start + basis_index;
      tripletList_Kinv.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
      alpha += NTdofs;
      tripletList_Kinv.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  Kinv_.setFromTriplets(tripletList_Kinv.begin(), tripletList_Kinv.end());
};

void PoissonSolver2D::generateS() {
  BT_ = B_.transpose();
  S_ = C_ + BT_ * Kinv_ * B_;
  std::cout << " Poisson Discrete Matrix S has been set up !! " << std::endl;
}

void PoissonSolver2D::init(const PoissonSolver2DType& poitype,
      const real_t& tol) {

  // ********* 传入相关变量 ********** //
  intboundaryNum_ = fe_->getmesh2D()->getintboundaryNum();
  extboundaryNum_ = fe_->getmesh2D()->getextboundaryNum();
  num_equations_ = fe_->getnum_equations();
  polydim_ = fe_->getbasis()->getpolydim();
  ncell_ = fe_->getmesh2D()->getncell();
  dim_ = fe_->getdim();
  numqua_ = fe_->getnumqua();
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
  case PoissonSolver2DType::CG:
    // std::cout << "  Preparation for Conjugate Gradient the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    cg_.compute(S_);
    cg_.setTolerance(tol);
    // std::cout << "  The end of preparation for Conjugate Gradient method ! " << std::endl;
    break;

  case PoissonSolver2DType::PCG:
    QUEST_ERROR(" The PCG method is not been implemented ! ");
    break;
  
  case PoissonSolver2DType::LDLT:
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

void PoissonSolver2D::solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const 
{
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // ******************************** //

  All->resize(dim_ + 1);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::CG:
    std::cout << "  Solve the Poisson by CG ......\n";
    xs = cg_.solve(f1);
    std::cout << "  iterations = " << cg_.iterations() << std::endl;
    std::cout << "  eerror = " << cg_.error()      << std::endl;
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver2DType::LDLT:
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
  Vector xsdx = xsd.segment(0, NTdofs);
  Vector xsdy = xsd.segment(NTdofs, NTdofs);
  All->at(0) = xs;
  All->at(1) = xsdx;
  All->at(2) = xsdy;
}

void PoissonSolver2D::solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Matrix>* All) const 
{ 
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // ******************************** //

  All->resize(dim_ + 1);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::CG:
    std::cout << "  Solve the Poisson by CG ......\n";
    xs = cg_.solve(f1);
    std::cout << "  iterations = " << cg_.iterations() << std::endl;
    std::cout << "  eerror = " << cg_.error()      << std::endl;
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver2DType::LDLT:
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
  Vector xsdx = xsd.segment(0, NTdofs);
  Vector xsdy = xsd.segment(NTdofs, NTdofs);
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdxtemp(xsdx.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdytemp(xsdy.data(), polydim_, ncell_);
  All->at(0) = xstemp;
  All->at(1) = xsdxtemp;
  All->at(2) = xsdytemp;
};

// -\Delta \Phi = f_nodal, E = - \nabla \Phi
void PoissonSolver2D::solveall(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                std::vector<Matrix>* E,
                Matrix* Phi) const
{
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // ******************************** //

  E->resize(dim_);
  Vector bg;
  Vector bf;
  Assemble_bg(Dirichlet, &bg);
  Assemble_bf(Dirichlet, f_nodal, &bf);

  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::CG:
    // std::cout << "  Solve the Poisson by CG ......\n";
    xs = cg_.solve(f1);
    // std::cout << "  iterations = " << cg_.iterations() << std::endl;
    // std::cout << "  eerror = " << cg_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver2DType::LDLT:
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
  Vector xsdx = xsd.segment(0, NTdofs);
  Vector xsdy = xsd.segment(NTdofs, NTdofs);
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdxtemp(xsdx.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdytemp(xsdy.data(), polydim_, ncell_);
  E->at(0) = - xsdxtemp;
  E->at(1) = - xsdytemp;
  *Phi = xstemp;
};


void PoissonSolver2D::Assemble_bg(const Matrix& Dirichlet,
                Vector* bg) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const DiagnalMatrix& wqua_diag1d = fe_->getwqua_diag1d();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const std::vector<Matrix>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intlength = 
    fe_->getmesh2D()->getintboundarylength();
  const std::vector<real_t>& extlength = 
    fe_->getmesh2D()->getextboundarylength();
  const std::vector<int>& extNei = 
    fe_->getmesh2D()->getextboundaryneighbors();
  const std::vector<Eigen::Vector2d>& extnormals = 
    fe_->getmesh2D()->getextboundarynormal();
  const std::vector<int>& ExtBTypeIndex = 
    fe_->getmesh2D()->getextboundarytypeindex();
  // ******************************** //

  bg->resize(NTdofs * dim_);
  bg->setZero();
  // for (int i = 0; i < intboundaryNum_; i++) 
  // {
  //   std::cout << " length = " << intlength[i] << std::endl;
  // }
  // PAUSE();
  for (int d = 0; d < 2; d++)
  {
    Eigen::Vector2d v;
    v.setZero();
    v(d) = 1.e0;
    for (int i = 0; i < extboundaryNum_; i++) 
    {
      real_t len = extlength[i];

      int test_cell_Index = extNei[i];
      Eigen::Vector2d test_normal = extnormals[i];

      Vector temp = boundary_u[ExtBTypeIndex[i]] * wqua_diag1d * Dirichlet.col(i) * len
                    * v.dot(test_normal);
      for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
      {
        int alpha = Tm(test_basis_index, test_cell_Index) + d * NTdofs;
        real_t ini_value =  temp(test_basis_index);
        (*bg)(alpha) += ini_value;
      };
    };
  }
};


void PoissonSolver2D::Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const DiagnalMatrix& wqua_diag1d = fe_->getwqua_diag1d();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const std::vector<Matrix>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();

  const std::vector<real_t>& extlength = 
    fe_->getmesh2D()->getextboundarylength();
  const std::vector<int>& extNei = 
    fe_->getmesh2D()->getextboundaryneighbors();
  const std::vector<Eigen::Vector2d>& extnormals = 
    fe_->getmesh2D()->getextboundarynormal();
  const std::vector<int>& ExtBTypeIndex = 
    fe_->getmesh2D()->getextboundarytypeindex();
  const std::vector<Eigen::Vector2d> extboundarycenter = 
    fe_->getmesh2D()->getextboundarycenter();
  // ******************************** //

  bf->resize(NTdofs);
  bf->setZero();
  Matrix bf_temp;
  fe_->Assemble_F(f_nodal, 0, 0, &bf_temp);
  Eigen::Map<Vector> reshapebf(bf_temp.data(), NTdofs, 1);
  *bf = reshapebf * JacobiDet(0);

  for (int i = 0; i < extboundaryNum_; i++) 
  {
    real_t len = extlength[i];

    int test_cell_Index = extNei[i];
    Eigen::Vector2d test_normal = extnormals[i];

    Vector temp = pa_.C11 * boundary_u[ExtBTypeIndex[i]] * wqua_diag1d * Dirichlet.col(i) * len;
    // std::cout << " len = " << len << std::endl;
    // std::cout << " cell_index = " << test_cell_Index << std::endl;  
    // std::cout << " normal = " << test_normal.transpose() << std::endl;
    // std::cout << " temp = \n" << temp.transpose() << std::endl;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
    {
      int alpha = Tm(test_basis_index, test_cell_Index);
      real_t ini_value = temp(test_basis_index);
      (*bf)(alpha) += ini_value;
    };
  };
  // PAUSE();
};

// period部分
PoissonSolver2D_period::PoissonSolver2D_period(const fespace2D* fe, const PoissonSolver2DParameter& pa)
  : PoissonSolver2D(fe, pa) {};

void PoissonSolver2D_period::init(const PoissonSolver2DType& poitype,
      const real_t& tol) {

  // ********* 传入相关变量 ********** //
  intboundaryNum_ = fe_->getmesh2D()->getintboundaryNum();
  extboundaryNum_ = fe_->getmesh2D()->getextboundaryNum();
  num_equations_ = fe_->getnum_equations();
  polydim_ = fe_->getbasis()->getpolydim();
  ncell_ = fe_->getmesh2D()->getncell();
  dim_ = fe_->getdim();
  numqua_ = fe_->getnumqua();
  // ******************************** //
  QUEST_VERIFY(fe_->getmesh2D()->IsPeriodBoundary(), "The mesh is not periodical !");
  generateB();
  generateBT();
  generateK();
  generateC();
  generateKinv();
  generateS();
  std::cout << "PoissonSolver2D_period::init begin !" << std::endl;
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::LU:
    // std::cout << "  Preparation for Conjugate Gradient the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.compute(S_);
    // std::cout << "  The end of preparation for Conjugate Gradient method ! " << std::endl;
    break;

  case PoissonSolver2DType::GMRES:
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
  std::cout << "PoissonSolver2D_period::init end !" << std::endl;

};


void PoissonSolver2D_period::generateB() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dxv_u = fe_->getdxv_u();
  const Matrix& dyv_u = fe_->getdyv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh2D()->getJx();
  const Vector& Jy = fe_->getmesh2D()->getJy();
  const std::vector<Matrix>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intlength = 
    fe_->getmesh2D()->getintboundarylength();
  const std::vector<Eigen::Vector2d>& intnormals = 
    fe_->getmesh2D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    fe_->getmesh2D()->getintboundaryneighbors();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = 
    fe_->getmesh2D()->getintboundarytypeindex();
  const std::vector<std::vector<Matrix>>& flux_u_v = 
    fe_->getflux_u_v();

  const std::vector<real_t>& extlength = 
    fe_->getmesh2D()->getextboundarylength();
  const std::vector<Eigen::Vector2d>& extnormals = 
    fe_->getmesh2D()->getextboundarynormal();
  const std::vector<Eigen::Vector2i>& extNei_period = 
    fe_->getmesh2D()->getextboundaryneighbors_period();
  const std::vector<Eigen::Vector2i>& ExtBTypeIndex_period = 
    fe_->getmesh2D()->getextboundarytypeindex_period();
  // ******************************** //
  B_.resize(NTdofs * dim_, NTdofs);
  int estimatedNonZeros = ncell_ * polydim_ * polydim_ + 
          intboundaryNum_ * polydim_ * polydim_ * 4 * dim_ + 
          extboundaryNum_ / 2 * polydim_ * polydim_ * 4;
  estimatedNonZeros = estimatedNonZeros * 2;
  std::vector<Eigen::Triplet<real_t>> tripletList_B;
  tripletList_B.reserve(estimatedNonZeros);
  
  Matrix Bref;
  Bref = dxv_u;
  for (int i = 0; i < ncell_; i++) 
  {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
    {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
      {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jx(i);
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };
  Bref = dyv_u;
  for (int i = 0; i < ncell_; i++) 
  {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
    {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
      {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet(i) * Jy(i);
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index + NTdofs, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };

  for (int d = 0; d < 2; d++)
  {
    Eigen::Vector2d v;
    v.setZero();
    v(d) = 1.e0;
    Matrix temp;
    for (int i = 0; i < intboundaryNum_; i++) 
    {
      real_t len = intlength[i];
      for (int test_cell = 0; test_cell < 2; test_cell++) 
      {
        int test_cell_Index = intNei[i](test_cell);
        // Matrix test_qua_value = boundary_u[IntBTypeIndex[i](test_cell)];
        Eigen::Vector2d test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
        for (int trial_cell = 0; trial_cell < 2; trial_cell++) 
        {
          int trial_cell_Index = intNei[i](trial_cell);
          // Matrix trial_qua_value = boundary_u[IntBTypeIndex[i](trial_cell)];
          Eigen::Vector2d trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);

          temp = flux_u_v[IntBTypeIndex[i](test_cell)][IntBTypeIndex[i](trial_cell)] * len;
          for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
          {
            int alpha = Tm(test_basis_index, test_cell_Index) + d * NTdofs;
            for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
            {
              int beta = Tm(trial_basis_index, trial_cell_Index);
              real_t ini_value =  (0.5e0 + pa_.C12.dot(trial_normal)) 
                  * temp(test_basis_index, trial_basis_index) * v.dot(test_normal);
              ini_value = - ini_value;
              if ( std::abs(ini_value) > 1.e-14 )
              {
                tripletList_B.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
              }
            };
          };
        };
      };
    };
    for (int i = 0; i < extboundaryNum_ / 2; i++) 
    {
      real_t len = extlength[i];
      for (int test_cell = 0; test_cell < 2; test_cell++) 
      {
        int test_cell_Index = extNei_period[i](test_cell);
        Eigen::Vector2d test_normal = extnormals[i] * std::pow(-1.e0,test_cell);
        for (int trial_cell = 0; trial_cell < 2; trial_cell++) 
        {
          int trial_cell_Index = extNei_period[i](trial_cell);
          Eigen::Vector2d trial_normal = extnormals[i] * std::pow(-1.e0,trial_cell);

          temp = flux_u_v[ExtBTypeIndex_period[i](test_cell)][ExtBTypeIndex_period[i](trial_cell)] 
                * len;
          // std::cout << " ExtBTypeIndex_period[i] = " 
          //           << ExtBTypeIndex_period[i](test_cell) <<  ", "
          //           << ExtBTypeIndex_period[i](trial_cell) << std::endl;
          // std::cout << " temp = " << temp.rows() << " * " << temp.cols() << std::endl;
          for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
          {
            int alpha = Tm(test_basis_index, test_cell_Index) + d * NTdofs;
            for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
            {
              int beta = Tm(trial_basis_index, trial_cell_Index);
              real_t ini_value =  (0.5e0 + pa_.C12.dot(trial_normal)) 
                  * temp(test_basis_index, trial_basis_index) * v.dot(test_normal);
              ini_value = - ini_value;
              if ( std::abs(ini_value) > 1.e-14 )
              {
                tripletList_B.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
              };
            };
          };
        };
      };
    };
  };

  B_.setFromTriplets(tripletList_B.begin(), tripletList_B.end());
};

void PoissonSolver2D_period::generateBT() {
  BT_ = B_.transpose();
};

void PoissonSolver2D_period::generateC() 
{
// ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag = fe_->getwqua_diag();
  const Matrix& dxv_u = fe_->getdxv_u();
  const Matrix& dyv_u = fe_->getdyv_u();
  const Matrix& v_u = fe_->getv_u();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  const Vector& Jx = fe_->getmesh2D()->getJx();
  const Vector& Jy = fe_->getmesh2D()->getJy();
  const std::vector<Matrix>& boundary_u = fe_->getboundary_u();
  const IntMatrix& Tm = fe_->getTm();
  const std::vector<real_t>& intlength = 
    fe_->getmesh2D()->getintboundarylength();
  const std::vector<Eigen::Vector2d>& intnormals = 
    fe_->getmesh2D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei = 
    fe_->getmesh2D()->getintboundaryneighbors();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = 
    fe_->getmesh2D()->getintboundarytypeindex();

  const std::vector<real_t>& extlength = 
    fe_->getmesh2D()->getextboundarylength();
  const std::vector<Eigen::Vector2d>& extnormals = 
    fe_->getmesh2D()->getextboundarynormal();
  const std::vector<Eigen::Vector2i>& extNei_period = 
    fe_->getmesh2D()->getextboundaryneighbors_period();
  const std::vector<Eigen::Vector2i>& ExtBTypeIndex_period = 
    fe_->getmesh2D()->getextboundarytypeindex_period();
  const std::vector<std::vector<Matrix>>& flux_u_v = 
    fe_->getflux_u_v();
  // ******************************** //

  C_.resize(NTdofs, NTdofs);
  int estimatedNonZeros = intboundaryNum_ * polydim_ * polydim_ * 4
    + extboundaryNum_ * polydim_ * polydim_;
  std::vector<Eigen::Triplet<real_t>> tripletList_C;
  tripletList_C.reserve(estimatedNonZeros);

  Matrix temp;
  for (int i = 0; i < intboundaryNum_; i++) 
  {
    real_t len = intlength[i];
    for (int test_cell = 0; test_cell < 2; test_cell++) 
    {
      int test_cell_Index = intNei[i](test_cell);
      Eigen::Vector2d test_normal = intnormals[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) 
      {
        int trial_cell_Index = intNei[i](trial_cell);
        Eigen::Vector2d trial_normal = intnormals[i] * std::pow(-1.e0,trial_cell);

        temp = flux_u_v[IntBTypeIndex[i](test_cell)][IntBTypeIndex[i](trial_cell)] * len;
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
        {
          int alpha = Tm(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
          {
            int beta = Tm(trial_basis_index, trial_cell_Index);
            real_t ini_value =  pa_.C11 * temp(test_basis_index, trial_basis_index)
                             * test_normal.dot(trial_normal);
            tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_ / 2; i++) 
  {
    real_t len = extlength[i];
    for (int test_cell = 0; test_cell < 2; test_cell++) 
    {
      int test_cell_Index = extNei_period[i](test_cell);
      Eigen::Vector2d test_normal = extnormals[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) 
      {
        int trial_cell_Index = extNei_period[i](trial_cell);
        Eigen::Vector2d trial_normal = extnormals[i] * std::pow(-1.e0,trial_cell);

        temp = flux_u_v[ExtBTypeIndex_period[i](test_cell)][ExtBTypeIndex_period[i](trial_cell)] 
              * len;
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) 
        {
          int alpha = Tm(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) 
          {
            int beta = Tm(trial_basis_index, trial_cell_Index);
            real_t ini_value =  pa_.C11 * temp(test_basis_index, trial_basis_index)
                             * test_normal.dot(trial_normal);
            tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  const real_t vol = fe_->getmesh2D()->getCellVol()(0);
  int center = ncell_ / 2;
  int alpha = Tm(0, center);
  real_t value = vol;
  for (int i = 0; i < ncell_; i++) {
    int beta = Tm(0, i);
    tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, vol));
  }

  C_.setFromTriplets(tripletList_C.begin(), tripletList_C.end());
};

void PoissonSolver2D_period::generateS() {
  S_ = C_ + BT_ * Kinv_ * B_;
  std::cout << " Poisson Discrete Matrix S has been set up !! " << std::endl;
};

void PoissonSolver2D_period::Assemble_bg(const Matrix& Dirichlet,
                Vector* bg) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  // const Matrix& dv_u_ = fe_->getdv_u();
  // const Matrix& v_u_ = fe_->getv_u();
  // const Vector& JacobiDet_ = fe_->getmesh2D()->getJacobiDet();
  // const Vector& Jx_ = fe_->getmesh2D()->getJx();
  // const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  // const IntMatrix& Tm_ = fe_->getTm();
  // const std::vector<int>& extNei_ = 
  //   fe_->getmesh2D()->getextboundaryneighbors();
  // const std::vector<real_t>& extnormals_ = 
  //   fe_->getmesh2D()->getextboundarynormal();
  // ******************************** //
  bg->resize(NTdofs * dim_);
  bg->setZero();
};

void PoissonSolver2D_period::Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const {
  
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  const Vector& JacobiDet = fe_->getmesh2D()->getJacobiDet();
  // ******************************** //

  bf->resize(NTdofs);
  bf->setZero();
  Matrix bf_temp;
  fe_->Assemble_F(f_nodal, 0, 0, &bf_temp);
  Eigen::Map<Vector> reshapebf(bf_temp.data(), NTdofs, 1);
  *bf = reshapebf * JacobiDet(0);
};

void PoissonSolver2D_period::solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const 
{
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // ******************************** //

  All->resize(dim_ + 1);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::LU:
    // std::cout << "  Solve the Poisson by LU ......\n";
    xs = lu_.solve(f1);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver2DType::GMRES:
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
  Vector xsdx = xsd.segment(0, NTdofs);
  Vector xsdy = xsd.segment(NTdofs, NTdofs);
  All->at(0) = xs;
  All->at(1) = xsdx;
  All->at(2) = xsdy;
}

void PoissonSolver2D_period::solveall(const Vector& bg,
                                      const Vector& bf,
                                      std::vector<Matrix>* All) const 
{ 
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // ******************************** //

  All->resize(dim_ + 1);
  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::LU:
    // std::cout << "  Solve the Poisson by LU ......\n";
    xs = lu_.solve(f1);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver2DType::GMRES:
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
  Vector xsdx = xsd.segment(0, NTdofs);
  Vector xsdy = xsd.segment(NTdofs, NTdofs);
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdxtemp(xsdx.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdytemp(xsdy.data(), polydim_, ncell_);
  All->at(0) = xstemp;
  All->at(1) = xsdxtemp;
  All->at(2) = xsdytemp;
};

// -\Delta \Phi = f_nodal, E = - \nabla \Phi
void PoissonSolver2D_period::solveall(const Matrix& Dirichlet,
                                      const Matrix& f_nodal,
                                      std::vector<Matrix>* E,
                                      Matrix* Phi) const
{
  // ********* 传入相关变量 ********** //
  const int& NTdofs = fe_->getNTdofs();
  // ******************************** //

  E->resize(dim_);
  Vector bg;
  Vector bf;
  Assemble_bg(Dirichlet, &bg);
  Assemble_bf(Dirichlet, f_nodal, &bf);

  Vector f1 = bf + BT_ * Kinv_ * bg;
  Vector xs;
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver2DType::LU:
    // std::cout << "  Solve the Poisson by LU ......\n";
    xs = lu_.solve(f1);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver2DType::GMRES:
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
  Vector xsdx = xsd.segment(0, NTdofs);
  Vector xsdy = xsd.segment(NTdofs, NTdofs);
  Eigen::Map<Matrix> xstemp(xs.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdxtemp(xsdx.data(), polydim_, ncell_);
  Eigen::Map<Matrix> xsdytemp(xsdy.data(), polydim_, ncell_);
  E->at(0) = - xsdxtemp;
  E->at(1) = - xsdytemp;
  *Phi = xstemp;

};

Poisson_acctest_2D::Poisson_acctest_2D(const fespace2D* fe, const PoissonSolver2DParameter& pa)
  : PoissonSolver2D(fe, pa) {};

Matrix Poisson_acctest_2D::RHS(const Matrix& x, const Matrix& y) 
{
  // Matrix rhs = - ( x.array().exp() + 2.e0 );
  Matrix rhs = - ( 2.e0 * x.array().sin() * y.array().cos());
  return rhs;
};

real_t Poisson_acctest_2D::u_bc(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x) * std::cos(y);
};

Matrix Poisson_acctest_2D::u_real(const Matrix& x, const Matrix& y) 
{
  // return x.array().exp() + x.array().square();
  return - x.array().sin() * y.array().cos();
};

real_t Poisson_acctest_2D::u_real(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x) * std::cos(y);
};

real_t Poisson_acctest_2D::u_real_dx(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + 2.e0 * x;
  return - std::cos(x) * std::cos(y);
};

Matrix Poisson_acctest_2D::u_real_dx(const Matrix& x, const Matrix& y) 
{
  // return x.array().exp() + 2.e0 * x.array();
  return - x.array().cos() * y.array().cos();
};

real_t Poisson_acctest_2D::u_real_dy(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + 2.e0 * x;
  return std::sin(x) * std::sin(y);
};

Matrix Poisson_acctest_2D::u_real_dy(const Matrix& x, const Matrix& y) 
{
  return x.array().sin() * y.array().sin();
};

Poisson_acctest_2D_period::Poisson_acctest_2D_period(const fespace2D* fe, const PoissonSolver2DParameter& pa)
  : PoissonSolver2D_period(fe, pa) {};

Matrix Poisson_acctest_2D_period::RHS(const Matrix& x, const Matrix& y) 
{
  // Matrix rhs = - ( x.array().exp() + 2.e0 );
  Matrix rhs = - ( 2.e0 * x.array().sin() * y.array().cos());
  return rhs;
};

real_t Poisson_acctest_2D_period::u_bc(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x) * std::cos(y);
};

Matrix Poisson_acctest_2D_period::u_real(const Matrix& x, const Matrix& y) 
{
  // return x.array().exp() + x.array().square();
  return - x.array().sin() * y.array().cos();
};

real_t Poisson_acctest_2D_period::u_real(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + x * x;
  return - std::sin(x) * std::cos(y);
};

real_t Poisson_acctest_2D_period::u_real_dx(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + 2.e0 * x;
  return - std::cos(x) * std::cos(y);
};

Matrix Poisson_acctest_2D_period::u_real_dx(const Matrix& x, const Matrix& y) 
{
  // return x.array().exp() + 2.e0 * x.array();
  return - x.array().cos() * y.array().cos();
};

real_t Poisson_acctest_2D_period::u_real_dy(const real_t& x, const real_t& y) 
{
  // return std::exp(x) + 2.e0 * x;
  return std::sin(x) * std::sin(y);
};

Matrix Poisson_acctest_2D_period::u_real_dy(const Matrix& x, const Matrix& y) 
{
  return x.array().sin() * y.array().sin();
};

} // namespace QUEST

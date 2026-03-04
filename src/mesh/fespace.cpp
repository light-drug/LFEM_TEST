#include "fespace.hpp"

void QUEST::fespace1D::setNTH(const int& NTH)
{
  NTH_ = NTH;
};

QUEST::fespace1D::fespace1D(const TensorMesh1D* mesh1D, BasisFunction1D* basis,
  const int& qua_order1D, const QuadratureType& quatype)
  : mesh1D_(mesh1D), basis_(basis), qua_order1D_(qua_order1D), quatype_(quatype),
    num_equations_(1) {};

QUEST::fespace1D::fespace1D(const TensorMesh1D* mesh1D, BasisFunction1D* basis,
  const int& qua_order1D, const QuadratureType& quatype, const int num_equations)
  : mesh1D_(mesh1D), basis_(basis), qua_order1D_(qua_order1D), quatype_(quatype),
    num_equations_(num_equations) {};

void QUEST::fespace1D::init() 
{
  generateTm();
  generateDG();
  generatePb();
};

const int& QUEST::fespace1D::getdim() const {
  return dim_;
}

const QUEST::TensorMesh1D* QUEST::fespace1D::getmesh1D() const {
  return mesh1D_;
}

const QUEST::BasisFunction1D* QUEST::fespace1D::getbasis() const {
  return basis_;
}

const int& QUEST::fespace1D::getnum_equations() const {
  return num_equations_;
}

const int& QUEST::fespace1D::getnumqua() const {
  return qua_order1D_;
}

const IntMatrix& QUEST::fespace1D::getTm() const {
  return Tm_;
};

const int& QUEST::fespace1D::getNTdofs() const {
  return NTdofs_;
};

const std::vector<Matrix>& QUEST::fespace1D::getflux_u_v() const {
  return flux_u_v_;
};

const Matrix& QUEST::fespace1D::getdv_u() const {
  return dv_u_;
};

const Matrix& QUEST::fespace1D::getv_u() const {
  return v_u_;
};

const DiagnalMatrix& QUEST::fespace1D::getv_u_diag() const {
  return v_u_diag_;
};

const DiagnalMatrix& QUEST::fespace1D::getv_u_diaginv() const {
  return v_u_diaginv_;
};

const std::vector<Vector>& QUEST::fespace1D::getboundary_u() const {
  return boundary_u_;
};

const Matrix& QUEST::fespace1D::gettest_ref() const {
  return test_ref_;
};

const Matrix& QUEST::fespace1D::gettest_ref_T() const {
  return test_ref_T_;
};

const Matrix& QUEST::fespace1D::gettest_dx_ref() const {
  return test_dx_ref_;
};

const DiagnalMatrix& QUEST::fespace1D::getwqua_diag() const {
  return wqua_diag_;
};

const Matrix& QUEST::fespace1D::getPb() const {
  return Pb_;
};

const int& QUEST::fespace1D::getNTH() const {
  return NTH_;
};

void QUEST::fespace1D::generateTm() {
  // ****** 传入相关变量  ******* //
  const int& ncell_ = mesh1D_->getncell();
  const int& k1D_ = basis_->getpolydim();
  // ************************** //
  NTdofs_ = k1D_ * ncell_;
  
  // ******* Tm 为基函数索引矩阵 ***** //
  int Index = 0;
  Tm_.resize(k1D_, ncell_); Tm_.setZero();
  for (int j = 0; j < ncell_; ++j) {
    for (int i = 0; i < k1D_; ++i) {
      Tm_(i, j) = Index;
      Index++;
    };
  };
}

void QUEST::fespace1D::generateDG() {
  // ****** 传入相关变量  ******* //
  const int& k1D_ = basis_->getpolydim();
  // *************************** //
  flux_u_v_.resize(4);
  boundary_u_.resize(2);

  Vector xqua, wqua;
  real_t xleft = - 1.e0 / 2.e0;
  real_t xright = 1.e0 / 2.e0;
  QUEST::initializeQuadradure(qua_order1D_, quatype_, &xqua, &wqua);
  basis_->dx(0).Map(xqua, &test_ref_);
  basis_->dx(1).Map(xqua, &test_dx_ref_);
  basis_->dx(0).Map(xright, &(boundary_u_[0]));
  basis_->dx(0).Map(xleft, &(boundary_u_[1]));
  wqua_diag_  = DiagnalMatrix(wqua);
  
  test_ref_T_ = test_ref_.transpose();
  dv_u_ = test_dx_ref_ * wqua_diag_ * test_ref_.transpose();
  v_u_ = test_ref_ * wqua_diag_ * test_ref_.transpose();
  v_u_diag_ = DiagnalMatrix(v_u_.diagonal());
  v_u_diaginv_ = v_u_diag_.inverse();
  int index = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      flux_u_v_[index] = boundary_u_[i] * boundary_u_[j].transpose();
      index++;
    };
  };
}

void QUEST::fespace1D::generatePb() {
  // ****** 传入相关变量  ******* //
  const int& ncell_ = mesh1D_->getncell();
  const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
  const real_t& hx_ = mesh1D_->gethx();
  Vector xqua, wqua;  
  QUEST::initializeQuadradure(qua_order1D_, quatype_, &xqua, &wqua);
  // *************************** //
  Pb_.resize(qua_order1D_, ncell_);Pb_.setZero();
  for (int i = 0; i < ncell_; i++) {
    for (int qu = 0; qu < qua_order1D_; qu++) {
      Pb_(qu, i) = Cellcenter_(0, i) + xqua(qu) * hx_;
    }
  }
}

// void QUEST::fespace1D::Project_Initial(Matrix (*func)(const Matrix&), Matrix* M) const {
//   // ****** 传入相关变量  ******* //
//   const int& ncell_ = mesh1D_->getncell();
//   const int& k1D_ = basis_->getpolydim();
//   const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
//   // ************************** //
//   Matrix f_quavalue(qua_order1D_, ncell_);
//   f_quavalue = func(Pb_);
//   nodal_to_modal1D(f_quavalue, M);
// }

// void QUEST::fespace1D::Project_Final(Matrix (*func)(const Matrix&, const real_t&), 
//                             const real_t& Tstop, Matrix* M) const {
//   // ****** 传入相关变量  ******* //
//   const int& ncell_ = mesh1D_->getncell();
//   const int& k1D_ = basis_->getpolydim();
//   const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
//   // ************************** //
//   Matrix f_quavalue(qua_order1D_, ncell_);
//   f_quavalue = func(Pb_, Tstop);
//   nodal_to_modal1D(f_quavalue, M);
// }

// void QUEST::fespace1D::Interpolate_Initial(Matrix (*func)(const Matrix&), Matrix* M) const {
//   // ****** 传入相关变量  ******* //
//   const int& ncell_ = mesh1D_->getncell();
//   const int& k1D_ = basis_->getpolydim();
//   const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
//   // ************************** //
//   M->resize(qua_order1D_, ncell_);
//   *M = func(Pb_);
// }

// void QUEST::fespace1D::Interpolate_Final(Matrix (*func)(const Matrix&, const real_t&), 
//                           const real_t& Tstop, Matrix* M) const {
//   // ****** 传入相关变量  ******* //
//   const int& ncell_ = mesh1D_->getncell();
//   const int& k1D_ = basis_->getpolydim();
//   const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
//   // ************************** //
//   M->resize(qua_order1D_, ncell_);
//   *M = func(Pb_, Tstop);
// }

void QUEST::fespace1D::Assemble_Flux(const std::vector<Matrix>& flux_int, const std::vector<Matrix>& flux_ext, 
                      std::vector<Matrix>* RHS) const {
    // ****** 传入相关变量  ******* //
    const int& ncell_ = mesh1D_->getncell();
    const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
    const real_t& hx_ = mesh1D_->gethx();
    const int& k1D_ = basis_->getk1D();
    const int& intboundaryNum_ = mesh1D_->getintboundaryNum();
    const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
    const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
    const int& extboundaryNum_ = mesh1D_->getextboundaryNum();
    const std::vector<int>& ExtBNei_ =  mesh1D_->getextboundaryneighbors();
    const std::vector<real_t>& ExtBNormal_ = mesh1D_->getextboundarynormal();
    int sys_ = flux_int.size();
    real_t hxinv = 1.e0 / hx_;
    // *************************** //
    Vector From = boundary_u_[0];
    Vector To = boundary_u_[1];
    RHS->resize(sys_);
    for (int j = 0; j < sys_; j++) {
      RHS->at(j).resize(k1D_, ncell_);
      RHS->at(j).setZero();
    };
    std::vector<Matrix> temp0, temp1;
    temp0 = From * flux_int * hxinv;
    temp1 = To * flux_int * hxinv;
// #pragma omp parallel num_threads(NTH_), default(shared)
//   {
    int Cellindex0, Cellindex1;
    int i, test_basis_index, j;
// #pragma omp for schedule(static)
    for (i = 0; i < intboundaryNum_; i++) {
      Cellindex0 = IntBNei[i](0);
      Cellindex1 = IntBNei[i](1);
      for (test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
        for (j = 0; j < sys_; j++) {
          #pragma omp atomic
          (*RHS)[j](test_basis_index, Cellindex0) += temp0[j](test_basis_index, i);
          #pragma omp atomic
          (*RHS)[j](test_basis_index, Cellindex1) -= temp1[j](test_basis_index, i);
        }
      }
    }
  // }
  
  int Cellindex = ExtBNei_[0];
  for (int j = 0; j < sys_; j++) {
    for (int test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
      real_t ini_value = flux_ext[j](0, 0) * To(test_basis_index) * hxinv;
      (*RHS)[j](test_basis_index, Cellindex) -= ini_value;
    };
  };

  Cellindex = ExtBNei_[1];
  for (int j = 0; j < sys_; j++) {
    for (int test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
      real_t ini_value = flux_ext[j](0, 1) * From(test_basis_index) * hxinv;
      (*RHS)[j](test_basis_index, Cellindex) += ini_value;
    };
  };
};

void QUEST::fespace1D::Assemble_Flux(const Matrix& flux_int, const Matrix& flux_ext, 
                      Matrix* RHS) const {
  // ****** 传入相关变量  ******* //
  const int& ncell_ = mesh1D_->getncell();
  const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
  const real_t& hx_ = mesh1D_->gethx();
  const int& k1D_ = basis_->getk1D();
  const int& intboundaryNum_ = mesh1D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  const int& extboundaryNum_ = mesh1D_->getextboundaryNum();
  const std::vector<int>& ExtBNei_ =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal_ = mesh1D_->getextboundarynormal();
  real_t hxinv = 1.e0 / hx_;
  // *************************** //
  Vector From = boundary_u_[0];
  Vector To = boundary_u_[1];
  RHS->resize(k1D_, ncell_);
  RHS->setZero();
  Matrix temp0, temp1;
  temp0 = From * flux_int * hxinv;
  temp1 = To * flux_int * hxinv;
// #pragma omp parallel num_threads(NTH_), default(shared)
//   {
    int Cellindex0, Cellindex1;
    int i, test_basis_index;
// #pragma omp for schedule(static)
    for (i = 0; i < intboundaryNum_; i++) {
      Cellindex0 = IntBNei[i](0);
      Cellindex1 = IntBNei[i](1);
      for (test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
        #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex0) += temp0(test_basis_index, i);
        #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex1) -= temp1(test_basis_index, i);
      }
    }
  // }
  
  int Cellindex = ExtBNei_[0];
  for (int test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
    real_t ini_value = flux_ext(0, 0) * To(test_basis_index) * hxinv;
    (*RHS)(test_basis_index, Cellindex) -= ini_value;
  };

  Cellindex = ExtBNei_[1];
  for (int test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
    real_t ini_value = flux_ext(0, 1) * From(test_basis_index) * hxinv;
    (*RHS)(test_basis_index, Cellindex) += ini_value;
  };
};

void QUEST::fespace1D::Assemble_Flux_int(const Matrix& flux_int, Matrix* RHS) const {
  // ****** 传入相关变量  ******* //
  const int& ncell_ = mesh1D_->getncell();
  const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
  const real_t& hx_ = mesh1D_->gethx();
  const int& k1D_ = basis_->getk1D();
  const int& intboundaryNum_ = mesh1D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  real_t hxinv = 1.e0 / hx_;
  // *************************** //
  Vector From = boundary_u_[0];
  Vector To = boundary_u_[1];
  RHS->resize(k1D_, ncell_);
  RHS->setZero();
  Matrix temp0, temp1;
  temp0 = From * flux_int * hxinv;
  temp1 = To * flux_int * hxinv;
// #pragma omp parallel num_threads(NTH_), default(shared)
//   {
    int Cellindex0, Cellindex1;
    int i, test_basis_index;
// #pragma omp for schedule(static)
    for (i = 0; i < intboundaryNum_; i++) {
      Cellindex0 = IntBNei[i](0);
      Cellindex1 = IntBNei[i](1);
      for (test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
        #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex0) += temp0(test_basis_index, i);
        #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex1) -= temp1(test_basis_index, i);
      }
    }
  // }
};

void QUEST::fespace1D::Assemble_Flux_bc(const Matrix& flux_ext, Matrix* RHS) const {
  // ****** 传入相关变量  ******* //
  const int& ncell_ = mesh1D_->getncell();
  const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
  const real_t& hx_ = mesh1D_->gethx();
  const int& k1D_ = basis_->getk1D();
  const int& intboundaryNum_ = mesh1D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  const int& extboundaryNum_ = mesh1D_->getextboundaryNum();
  const std::vector<int>& ExtBNei_ =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal_ = mesh1D_->getextboundarynormal();
  real_t hxinv = 1.e0 / hx_;
  // *************************** //
  Vector From = boundary_u_[0];
  Vector To = boundary_u_[1];
  RHS->resize(k1D_, ncell_);
  RHS->setZero();

  int Cellindex = ExtBNei_[0];
  for (int test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
    real_t ini_value = flux_ext(0, 0) * To(test_basis_index) * hxinv;
    (*RHS)(test_basis_index, Cellindex) -= ini_value;
  };

  Cellindex = ExtBNei_[1];
  for (int test_basis_index = 0; test_basis_index < k1D_; test_basis_index++) {
    real_t ini_value = flux_ext(0, 1) * From(test_basis_index) * hxinv;
    (*RHS)(test_basis_index, Cellindex) += ini_value;
  };
};

void QUEST::fespace1D::DisplayResult(const Vector& numerical, const Vector& exact,  
    const std::string& title, std::ostream& Outfile) const {
  // ****** 传入相关变量  ******* //
  const int& xDiv = mesh1D_->getxDiv();
  const int& ncell = mesh1D_->getncell();
  const Matrix& Cellcenter = mesh1D_->getCellCenter();
  // *************************** //

  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = \"X\", \"Numerical Solution\", \"Exact Solution\"" << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << xDiv << " F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < ncell; i++) {
    real_t x0 = Cellcenter(0, i);
    Outfile << x0 << " " 
            << numerical(i) << " "
            << exact(i) << " " << std::endl;
  };
};

void QUEST::fespace1D::DisplayResult(const std::vector<Vector>& output, 
  const std::vector<std::string>& outputname,
  const std::string title, std::ostream& Outfile) const {
  // ****** 传入相关变量  ******* //
  const int& xDiv = mesh1D_->getxDiv();
  const int& ncell = mesh1D_->getncell();
  const Matrix& Cellcenter = mesh1D_->getCellCenter();
  int size = output.size();
  // *************************** //

  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = ";
  for (int k = 0; k < size; k++) {
    Outfile << "\"" << outputname[k] <<"\", ";
  }
  Outfile << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << xDiv << " F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < ncell; i++) {
    real_t x0 = Cellcenter(0, i);
    Outfile << x0 << " ";
    for (int k = 0; k < size; k++) {
      Outfile << output[k](i) <<" ";
    };
    Outfile << std::endl;
  };
};

void QUEST::fespace1D::computerrorL1(const Matrix& u_numerical_nodal, 
                          const Matrix& u_exact_nodal, 
                          real_t* error) const {
  Matrix error_nodal = (u_numerical_nodal - u_exact_nodal).cwiseAbs();
  DiagnalMatrix wdiag = wqua_diag_ * mesh1D_->getCellVol()(0);
  error_nodal = wdiag * error_nodal;
  *error = error_nodal.sum();
  // std::cout << u_numerical_nodal << std::endl;
  // std::cout << u_exact_nodal << std::endl;
};

void QUEST::fespace1D::computerrorL2(const Matrix& u_numerical_nodal, 
                          const Matrix& u_exact_nodal, 
                          real_t* error) const {
  Matrix error_nodal = (u_numerical_nodal - u_exact_nodal).cwiseAbs2();
  DiagnalMatrix wdiag = wqua_diag_ * mesh1D_->getCellVol()(0);
  error_nodal = wdiag * error_nodal;
  *error = error_nodal.sum();
  *error = std::sqrt(*error);
};

void QUEST::fespace1D::computerrorLinf(const Matrix& u_numerical_nodal, 
                          const Matrix& u_exact_nodal, 
                          real_t* error) const {
  Matrix error_nodal = (u_numerical_nodal - u_exact_nodal).cwiseAbs();
  *error = error_nodal.maxCoeff();
};

// fespace2D
QUEST::fespace2D::fespace2D(const TensorMesh2D* mesh2D, BasisFunction2D* basis, 
            const int& qua_order1D, const QuadratureType& quatype)
  : mesh2D_(mesh2D), basis_(basis), qua_order1D_(qua_order1D), quatype_(quatype),
    num_equations_(1) {};

QUEST::fespace2D::fespace2D(const TensorMesh2D* mesh2D, BasisFunction2D* basis,
  const int& qua_order1D, const QuadratureType& quatype, const int num_equations)
  : mesh2D_(mesh2D), basis_(basis), qua_order1D_(qua_order1D), quatype_(quatype),
    num_equations_(num_equations) {};

void QUEST::fespace2D::setNTH(const int& NTH)
{
  NTH_ = NTH;
};

const int& QUEST::fespace2D::getdim() const {
  return dim_;
}

const QUEST::TensorMesh2D* QUEST::fespace2D::getmesh2D() const {
  return mesh2D_;
}

const QUEST::BasisFunction2D* QUEST::fespace2D::getbasis() const {
  return basis_;
}

const int& QUEST::fespace2D::getnum_equations() const {
  return num_equations_;
}

const int& QUEST::fespace2D::getnumqua() const {
  return numqua_;
}

const int& QUEST::fespace2D::getnumqua1d() const {
  return qua_order1D_;
}

const IntMatrix& QUEST::fespace2D::getTm() const {
  return Tm_;
};

const int& QUEST::fespace2D::getNTdofs() const {
  return NTdofs_;
};

const std::vector<std::vector<Matrix>>& QUEST::fespace2D::getflux_u_v() const {
  return flux_u_v_;
};

const Matrix& QUEST::fespace2D::getdxv_u() const {
  return dxv_u_;
};

const Matrix& QUEST::fespace2D::getdyv_u() const {
  return dyv_u_;
};

const Matrix& QUEST::fespace2D::getv_u() const {
  return v_u_;
};

const DiagnalMatrix& QUEST::fespace2D::getv_u_diag() const {
  return v_u_diag_;
};

const DiagnalMatrix& QUEST::fespace2D::getv_u_diaginv() const {
  return v_u_diaginv_;
};

const std::vector<Matrix>& QUEST::fespace2D::getboundary_u() const {
  return boundary_u_;
};

const Matrix& QUEST::fespace2D::getCoorRef() const
{
  return CoorRef_;
};

const std::vector<Matrix>& QUEST::fespace2D::getCoorBdrRef() const
{
  return CoorBdrRef_;
};

const Matrix& QUEST::fespace2D::gettest_ref() const {
  return test_ref_;
};

const Matrix& QUEST::fespace2D::gettest_ref_T() const {
  return test_ref_T_;
};

const Matrix& QUEST::fespace2D::gettest_dx_ref() const {
  return test_dx_ref_;
};

const Matrix& QUEST::fespace2D::gettest_dy_ref() const {
  return test_dy_ref_;
};

const DiagnalMatrix& QUEST::fespace2D::getwqua_diag() const {
  return wqua_diag_;
};

const DiagnalMatrix& QUEST::fespace2D::getwqua_diag1d() const {
  return wqua_diag1d_;
};

const std::vector<Matrix>& QUEST::fespace2D::getPb() const {
  return Pb_;
};

const int& QUEST::fespace2D::getNTH() const {
  return NTH_;
};

void QUEST::fespace2D::Assemble_Flux(const std::vector<Matrix>& flux_int, const std::vector<Matrix>& flux_ext, 
                      std::vector<Matrix>* RHS) const
{
  // ****** 传入相关变量  ******* //
  const int& ncell = mesh2D_->getncell();
  const real_t& hx = mesh2D_->gethx();
  const int& k1D = basis_->getk1D();
  const int& ploydim = basis_->getpolydim();
  const int& intboundaryNum = mesh2D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh2D_->getintboundaryneighbors();
  const std::vector<Eigen::Vector2d>& IntBNormal = mesh2D_->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = mesh2D_->getintboundarytypeindex();
  const std::vector<real_t>& IntLength = mesh2D_->getintboundarylength();
  const int& extboundaryNum = mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh2D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtLength = mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& ExtBNormal = mesh2D_->getextboundarynormal();
  const std::vector<int>& ExtBTypeIndex = mesh2D_->getextboundarytypeindex();
  int sys = flux_int.size();
  const Vector& JacobiDet = mesh2D_->getJacobiDet();
  real_t hxinv = 1.e0 / JacobiDet(0);
  // *************************** //
  RHS->resize(sys);
  for (int j = 0; j < sys; j++) {
    RHS->at(j).resize(ploydim, ncell);
    RHS->at(j).setZero();
  };

// #pragma omp parallel num_threads(NTH_), default(shared)
  {
    int Cellindex0, Cellindex1;
    int i, test_basis_index, j;
    real_t length;
    Eigen::Vector2i typeindex;
    Matrix From, To;
// #pragma omp for schedule(static)
    for (i = 0; i < intboundaryNum; i++) 
    {
      Cellindex0 = IntBNei[i](0);
      Cellindex1 = IntBNei[i](1);
      length = IntLength[i];
      typeindex = IntBTypeIndex[i];
      for (j = 0; j < sys; j++) 
      {
        From = boundary_u_[typeindex(0)] * wqua_diag1d_ * flux_int[j].col(i) * length * hxinv;
        To = boundary_u_[typeindex(1)] * wqua_diag1d_ * flux_int[j].col(i) * length * hxinv;
        for (test_basis_index = 0; test_basis_index < ploydim; test_basis_index++) 
        {
          // #pragma omp atomic
          (*RHS)[j](test_basis_index, Cellindex0) += From(test_basis_index);
          // #pragma omp atomic
          (*RHS)[j](test_basis_index, Cellindex1) -= To(test_basis_index);
        }
      }
    }
  }
  {
    int Cellindex;
    int i, test_basis_index, j;
    real_t length;
    int  typeindex;
    Matrix From, To;
    for (i = 0; i < extboundaryNum; i++) 
    {
      Cellindex = ExtBNei[0];
      length = ExtLength[i];
      typeindex = ExtBTypeIndex[i];
      for (j = 0; j < sys; j++) 
      {
        From = boundary_u_[typeindex] * wqua_diag1d_ * flux_ext[j].col(i) * length * hxinv;
        for (test_basis_index = 0; test_basis_index < ploydim; test_basis_index++) 
        {
          // #pragma omp atomic
          (*RHS)[j](test_basis_index, Cellindex) += From(test_basis_index);
        }
      }
    }
  }
};

void QUEST::fespace2D::Assemble_Flux(const Matrix& flux_int, const Matrix& flux_ext, 
                    Matrix* RHS) const
{
  // ****** 传入相关变量  ******* //
  const int& ncell = mesh2D_->getncell();
  const real_t& hx = mesh2D_->gethx();
  const int& k1D = basis_->getk1D();
  const int& ploydim = basis_->getpolydim();
  const int& intboundaryNum = mesh2D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh2D_->getintboundaryneighbors();
  const std::vector<Eigen::Vector2d>& IntBNormal = mesh2D_->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = mesh2D_->getintboundarytypeindex();
  const std::vector<real_t>& IntLength = mesh2D_->getintboundarylength();

  const int& extboundaryNum = mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh2D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtLength = mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& ExtBNormal = mesh2D_->getextboundarynormal();
  const std::vector<int>& ExtBTypeIndex = mesh2D_->getextboundarytypeindex();
  const Vector& JacobiDet = mesh2D_->getJacobiDet();
  real_t hxinv = 1.e0 / JacobiDet(0);
  // *************************** //
  RHS->resize(ploydim, ncell);

// #pragma omp parallel num_threads(NTH_), default(shared)
  {
    int Cellindex0, Cellindex1;
    int i, test_basis_index, j;
    real_t length;
    Eigen::Vector2i typeindex;
    Matrix From, To;
// #pragma omp for schedule(static)
    for (i = 0; i < intboundaryNum; i++) 
    {
      Cellindex0 = IntBNei[i](0);
      Cellindex1 = IntBNei[i](1);
      length = IntLength[i];
      typeindex = IntBTypeIndex[i];

      From = boundary_u_[typeindex(0)] * wqua_diag1d_ * flux_int.col(i) * length * hxinv;
      To = boundary_u_[typeindex(1)] * wqua_diag1d_ * flux_int.col(i) * length * hxinv;
      for (test_basis_index = 0; test_basis_index < ploydim; test_basis_index++) 
      {
        // #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex0) += From(test_basis_index);
        // #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex1) -= To(test_basis_index);
      }
    }
  }
  {
    int Cellindex;
    int i, test_basis_index, j;
    real_t length;
    int  typeindex;
    Matrix From;
    for (i = 0; i < extboundaryNum; i++) 
    {
      Cellindex = ExtBNei[0];
      length = ExtLength[i];
      typeindex = ExtBTypeIndex[i];

      From = boundary_u_[typeindex] * wqua_diag1d_ * flux_ext.col(i) * length * hxinv;
      for (test_basis_index = 0; test_basis_index < ploydim; test_basis_index++) 
      {
        // #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex) += From(test_basis_index);
      }
    }
  }
};

void QUEST::fespace2D::Assemble_Flux_int(const Matrix& flux_int, Matrix* RHS) const
{
  // ****** 传入相关变量  ******* //
  const int& ncell = mesh2D_->getncell();
  const real_t& hx = mesh2D_->gethx();
  const int& k1D = basis_->getk1D();
  const int& ploydim = basis_->getpolydim();
  const int& intboundaryNum = mesh2D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh2D_->getintboundaryneighbors();
  const std::vector<Eigen::Vector2d>& IntBNormal = mesh2D_->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = mesh2D_->getintboundarytypeindex();
  const std::vector<real_t>& IntLength = mesh2D_->getintboundarylength();

  const int& extboundaryNum = mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh2D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtLength = mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& ExtBNormal = mesh2D_->getextboundarynormal();
  const std::vector<int>& ExtBTypeIndex = mesh2D_->getextboundarytypeindex();
  const Vector& JacobiDet = mesh2D_->getJacobiDet();
  real_t hxinv = 1.e0 / JacobiDet(0);
  // *************************** //
  RHS->resize(ploydim, ncell);

// #pragma omp parallel num_threads(NTH_), default(shared)
  {
    int Cellindex0, Cellindex1;
    int i, test_basis_index, j;
    real_t length;
    Eigen::Vector2i typeindex;
    Matrix From, To;
// #pragma omp for schedule(static)
    for (i = 0; i < intboundaryNum; i++) 
    {
      Cellindex0 = IntBNei[i](0);
      Cellindex1 = IntBNei[i](1);
      length = IntLength[i];
      typeindex = IntBTypeIndex[i];

      From = boundary_u_[typeindex(0)] * wqua_diag1d_ * flux_int.col(i) * length * hxinv;
      To = boundary_u_[typeindex(1)] * wqua_diag1d_ * flux_int.col(i) * length * hxinv;
      for (test_basis_index = 0; test_basis_index < ploydim; test_basis_index++) 
      {
        // #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex0) += From(test_basis_index);
        // #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex1) -= To(test_basis_index);
      }
    }
  }
};

void QUEST::fespace2D::Assemble_Flux_bc(const Matrix& flux_ext, Matrix* RHS) const
{
  // ****** 传入相关变量  ******* //
  const int& ncell = mesh2D_->getncell();
  const real_t& hx = mesh2D_->gethx();
  const int& k1D = basis_->getk1D();
  const int& ploydim = basis_->getpolydim();
  const int& intboundaryNum = mesh2D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh2D_->getintboundaryneighbors();
  const std::vector<Eigen::Vector2d>& IntBNormal = mesh2D_->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& IntBTypeIndex = mesh2D_->getintboundarytypeindex();
  const std::vector<real_t>& IntLength = mesh2D_->getintboundarylength();

  const int& extboundaryNum = mesh2D_->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh2D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtLength = mesh2D_->getextboundarylength();
  const std::vector<Eigen::Vector2d>& ExtBNormal = mesh2D_->getextboundarynormal();
  const std::vector<int>& ExtBTypeIndex = mesh2D_->getextboundarytypeindex();
  const Vector& JacobiDet = mesh2D_->getJacobiDet();
  real_t hxinv = 1.e0 / JacobiDet(0);
  // *************************** //
  RHS->resize(ploydim, ncell);

  {
    int Cellindex;
    int i, test_basis_index, j;
    real_t length;
    int  typeindex;
    Matrix From;
    for (i = 0; i < extboundaryNum; i++) 
    {
      Cellindex = ExtBNei[0];
      length = ExtLength[i];
      typeindex = ExtBTypeIndex[i];

      From = boundary_u_[typeindex] * wqua_diag1d_ * flux_ext.col(i) * length * hxinv;
      for (test_basis_index = 0; test_basis_index < ploydim; test_basis_index++) 
      {
        // #pragma omp atomic
        (*RHS)(test_basis_index, Cellindex) += From(test_basis_index);
      }
    }
  }
};

void QUEST::fespace2D::init() 
{
  generateTm();
  generateDG();
  generatePb();
};

void QUEST::fespace2D::generateTm() {
  // ****** 传入相关变量  ******* //
  const int& ncell = mesh2D_->getncell();
  const int& polydim = basis_->getpolydim();
  // ************************** //
  NTdofs_ = polydim * ncell;
  
  // ******* Tm 为基函数索引矩阵 ***** //
  int Index = 0;
  Tm_.resize(polydim, ncell); Tm_.setZero();
  for (int j = 0; j < ncell; ++j) {
    for (int i = 0; i < polydim; ++i) {
      Tm_(i, j) = Index;
      Index++;
    };
  };
}

void QUEST::fespace2D::generateDG() {
  // ****** 传入相关变量  ******* //
  const int& polydim = basis_->getpolydim();
  // *************************** //
  boundary_u_.resize(2);

  Vector xqua, wqua;
  real_t xleft = - 1.e0 / 2.e0;
  real_t xright = 1.e0 / 2.e0;
  QUEST::initializeQuadradure(qua_order1D_, quatype_, &xqua, &wqua);
  numqua_ = qua_order1D_ * qua_order1D_;
  wqua2d_.resize(numqua_);
  CoorRef_.resize(2, numqua_);
  int index = 0;
  for (int qy = 0; qy < qua_order1D_; qy++)
  {
    for (int qx = 0; qx < qua_order1D_; qx++)
    {
      wqua2d_(index) = wqua(qy) * wqua(qx);
      CoorRef_(0, index) = xqua(qx);
      CoorRef_(1, index) = xqua(qy);
      index++;
    }
  }
  basis_->dx(0).dy(0).Map(CoorRef_, &test_ref_);
  basis_->dx(1).dy(0).Map(CoorRef_, &test_dx_ref_);
  basis_->dx(0).dy(1).Map(CoorRef_, &test_dy_ref_);
  int kdim = 2 * dim_;
  CoorBdrRef_.resize(kdim);
  boundary_u_.resize(kdim);
  for (int i = 0; i < kdim; i++)
  {
    CoorBdrRef_[i].resize(2, qua_order1D_);
  }
  for (int q = 0; q < qua_order1D_; q++)
  {
    CoorBdrRef_[0](0, q) = xright;
    CoorBdrRef_[0](1, q) = xqua(q);
    CoorBdrRef_[1](0, q) = xleft;
    CoorBdrRef_[1](1, q) = xqua(q);
    CoorBdrRef_[2](0, q) = xqua(q);
    CoorBdrRef_[2](1, q) = xright;
    CoorBdrRef_[3](0, q) = xqua(q);
    CoorBdrRef_[3](1, q) = xleft;
  }
  for (int i = 0; i < kdim; i++)
  {
    basis_->dx(0).dy(0).Map(CoorBdrRef_[i], &(boundary_u_[i]));
  }
  wqua_diag1d_ = DiagnalMatrix(wqua);
  wqua_diag_ = DiagnalMatrix(wqua2d_);
  
  test_ref_T_ = test_ref_.transpose();
  dxv_u_ = test_dx_ref_ * wqua_diag_ * test_ref_.transpose();
  dyv_u_ = test_dy_ref_ * wqua_diag_ * test_ref_.transpose();
  v_u_ = test_ref_ * wqua_diag_ * test_ref_.transpose();
  v_u_diag_ = DiagnalMatrix(v_u_.diagonal());
  v_u_diaginv_ = v_u_diag_.inverse();

  flux_u_v_.resize(kdim);
  for (int i = 0; i < kdim; i++)
  {
    flux_u_v_[i].resize(kdim);
    for (int j = 0; j < kdim; j++)
    {
      flux_u_v_[i][j] = boundary_u_[i] * wqua_diag1d_ * boundary_u_[j].transpose();
    };
  };
  
};

void QUEST::fespace2D::generatePb() {
  // ****** 传入相关变量  ******* //
  const int& ncell = mesh2D_->getncell();
  const Matrix& Cellcenter = mesh2D_->getCellCenter();
  const real_t& hx = mesh2D_->gethx();
  const real_t& hy = mesh2D_->gethy();
  Vector xqua, wqua;  
  QUEST::initializeQuadradure(qua_order1D_, quatype_, &xqua, &wqua);
  // *************************** //
  Pb_.resize(2);
  Pb_[0].resize(numqua_, ncell);
  Pb_[0].setZero();
  Pb_[1].resize(numqua_, ncell);
  Pb_[1].setZero();
  int index;
  for (int i = 0; i < ncell; i++) 
  {
    index = 0;
    for (int qu = 0; qu < numqua_; qu++) 
    {
      Pb_[0](qu, i) = Cellcenter(0, i) + CoorRef_(0, qu) * hx;
      Pb_[1](qu, i) = Cellcenter(1, i) + CoorRef_(1, qu) * hy;
    }
  }
}

void QUEST::fespace2D::DisplayResult(const Vector& numerical, const Vector& exact,  
    const std::string& title, std::ostream& Outfile) const {
  // ****** 传入相关变量  ******* //
  const int& xDiv = mesh2D_->getxDiv();
  const int& yDiv = mesh2D_->getyDiv();
  const int& ncell = mesh2D_->getncell();
  const Matrix& Cellcenter = mesh2D_->getCellCenter();
  // *************************** //

  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = \"X\", \"Y\" \"Numerical Solution\", \"Exact Solution\"" << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << xDiv << ", J =  " << yDiv << ", F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < ncell; i++) {
    real_t x0 = Cellcenter(0, i);
    real_t y0 = Cellcenter(1, i);
    Outfile << x0 << " " 
            << y0 << " "
            << numerical(i) << " "
            << exact(i) << " " << std::endl;
  };
};

void QUEST::fespace2D::DisplayResult(const std::vector<Vector>& output, 
  const std::vector<std::string>& outputname,
  const std::string title, std::ostream& Outfile) const {
  // ****** 传入相关变量  ******* //
  const int& xDiv = mesh2D_->getxDiv();
  const int& yDiv = mesh2D_->getyDiv();
  const int& ncell = mesh2D_->getncell();
  const Matrix& Cellcenter = mesh2D_->getCellCenter();
  int size = output.size();
  // *************************** //

  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = ";
  for (int k = 0; k < size; k++) {
    Outfile << "\"" << outputname[k] <<"\", ";
  }
  Outfile << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << xDiv 
          << ", J =  " << yDiv << ", F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < ncell; i++) {
    real_t x0 = Cellcenter(0, i);
    real_t y0 = Cellcenter(1, i);
    Outfile << x0 << " " 
            << y0 << " ";
    for (int k = 0; k < size; k++) {
      Outfile << output[k](i) <<" ";
    };
    Outfile << std::endl;
  };
};

void QUEST::fespace2D::computerrorL1(const Matrix& u_numerical_nodal, 
                          const Matrix& u_exact_nodal, 
                          real_t* error) const {
  Matrix error_nodal = (u_numerical_nodal - u_exact_nodal).cwiseAbs();
  DiagnalMatrix wdiag = wqua_diag_ * mesh2D_->getCellVol()(0);
  error_nodal = wdiag * error_nodal;
  *error = error_nodal.sum();
  // std::cout << u_numerical_nodal << std::endl;
  // std::cout << u_exact_nodal << std::endl;
};

void QUEST::fespace2D::computerrorL2(const Matrix& u_numerical_nodal, 
                          const Matrix& u_exact_nodal, 
                          real_t* error) const {
  Matrix error_nodal = (u_numerical_nodal - u_exact_nodal).cwiseAbs2();
  DiagnalMatrix wdiag = wqua_diag_ * mesh2D_->getCellVol()(0);
  error_nodal = wdiag * error_nodal;
  *error = error_nodal.sum();
  *error = std::sqrt(*error);
};

void QUEST::fespace2D::computerrorLinf(const Matrix& u_numerical_nodal, 
                          const Matrix& u_exact_nodal, 
                          real_t* error) const {
  Matrix error_nodal = (u_numerical_nodal - u_exact_nodal).cwiseAbs();
  *error = error_nodal.maxCoeff();
};



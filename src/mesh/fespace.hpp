#ifndef QUEST_FESPACE_HPP
#define QUEST_FESPACE_HPP

#include "Tensormesh.hpp"
#include "basis.hpp"
#include "config.hpp"
#include "error.hpp"
#include "integralrule.hpp"
#include "vector_overload.hpp"

#include <iomanip>
#include <fstream>
#include <iostream>

namespace QUEST {

class fespace1D 
{
public:
  // **** 标量有限元函数空间 ***** // 
  fespace1D(const TensorMesh1D* mesh1D, BasisFunction1D* basis, 
            const int& qua_order1D, const QuadratureType& quatype);

  // **** 向量有限元函数空间 ***** // 
  fespace1D(const TensorMesh1D* mesh1D, BasisFunction1D* basis, 
            const int& qua_order1D, const QuadratureType& quatype, 
            const int num_equations);
  virtual ~fespace1D() = default;

  int NTH_;

  virtual void setNTH(const int& NTH);
  
  virtual void init();
  const int& getdim() const;
  const TensorMesh1D* getmesh1D() const;
  const BasisFunction1D* getbasis() const;
  const int& getnum_equations() const;
  const int& getnumqua() const;
  const IntMatrix& getTm() const;
  const int& getNTdofs() const;
  const std::vector<Matrix>& getflux_u_v() const;
  const Matrix& getdv_u() const;
  const Matrix& getv_u() const;
  const DiagnalMatrix& getv_u_diag() const;
  const DiagnalMatrix& getv_u_diaginv() const;
  const std::vector<Vector>& getboundary_u() const;
  const Matrix& gettest_ref() const;
  const Matrix& gettest_ref_T() const;
  const Matrix& gettest_dx_ref() const;
  const DiagnalMatrix& getwqua_diag() const;
  const Matrix& getPb() const;
  const int& getNTH() const;
  
  template<typename Func>
  void Project_Initial(Func&& func, Matrix* M) const
  {
    // ****** 传入相关变量  ******* //
    const int& ncell_ = mesh1D_->getncell();
    const int& k1D_ = basis_->getpolydim();
    const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
    // ************************** //
    Matrix f_quavalue(qua_order1D_, ncell_);
    f_quavalue = std::forward<Func>(func)(Pb_);
    nodal_to_modal1D(f_quavalue, M);
  };

  template<typename Func>
  void Project_Final(Func&& func, const real_t& Tstop, Matrix* M) const
  {
    // ****** 传入相关变量  ******* //
    const int& ncell_ = mesh1D_->getncell();
    const int& k1D_ = basis_->getpolydim();
    const Matrix& Cellcenter_ = mesh1D_->getCellCenter();
    // ************************** //
    Matrix f_quavalue(qua_order1D_, ncell_);
    f_quavalue = std::forward<Func>(func)(Pb_, Tstop);
    nodal_to_modal1D(f_quavalue, M);
  };

  template<typename Func>
  void Interpolate_Initial(Func&& func, Matrix* M) const
  {
    *M = std::forward<Func>(func)(Pb_);
  };

  template<typename Func>
  void Interpolate_Final(Func&& func, const real_t& Tstop, Matrix* M) const
  {
    *M = std::forward<Func>(func)(Pb_, Tstop);
  };

  // T 为Matrix或者std::vector<Matrix>
  template <typename T>
  void nodal_to_modal1D(const T& N, T* M) const {
    Matrix temp = v_u_diaginv_ * test_ref_ * wqua_diag_;
    *M = temp * N;
  };

  // T 为Matrix或者std::vector<Matrix>
  template <typename T>
  void modal_to_nodal1D(const T& M, T* N) const {
    Matrix temp = test_ref_.transpose();
    *N = temp * M;
  };

  template <typename Derived>
  void nodal_to_modal1D(const Eigen::MatrixBase<Derived>& N, Matrix* M) const {
    Matrix temp = v_u_diaginv_ * test_ref_ * wqua_diag_;
    *M = temp * N;
  };

  template <typename Derived>
  void modal_to_nodal1D(const Eigen::MatrixBase<Derived>& M, Matrix* N) const {
    *N = test_ref_.transpose() * M;
  }

  // T 为Matrix或者std::vector<Matrix>
  template <typename T>
  void Assemble_F(const T& function_f_quavalue, const int& test_dx_order, T* RHS) const {
    const Vector& Jx = mesh1D_->getJx();
    const Matrix* test_temp;
    switch (test_dx_order)
    {
    case 0: test_temp = &test_ref_;  break;
    case 1: test_temp = &test_dx_ref_;  break;
    // other: basis_->dx(0).Map(xqua, &test_ref_);
    default: QUEST_ERROR("unsupported test_dx_order ! "); break;
    }
    Matrix test_qua = std::pow(Jx(0),test_dx_order) * (*test_temp) * wqua_diag_;
    *RHS = test_qua * function_f_quavalue;
  };

  void Assemble_Flux(const std::vector<Matrix>& flux_int, const std::vector<Matrix>& flux_ext, 
                      std::vector<Matrix>* RHS) const;

  void Assemble_Flux(const Matrix& flux_int, const Matrix& flux_ext, 
                      Matrix* RHS) const;

  void Assemble_Flux_int(const Matrix& flux_int, Matrix* RHS) const;
  void Assemble_Flux_bc(const Matrix& flux_ext, Matrix* RHS) const;
  
  void DisplayResult(const Vector& numerical, const Vector& exact,  
                    const std::string& title, std::ostream& Outfile) const;

  void DisplayResult(const std::vector<Vector>& output, 
                    const std::vector<std::string>& outputname,
                    const std::string title, std::ostream& Outfile) const;

  virtual void computerrorL1(const Matrix& u_numerical_nodal, const Matrix& u_exact_nodal, 
                          real_t* error) const;
  virtual void computerrorL2(const Matrix& u_numerical_nodal, const Matrix& u_exact_nodal, 
                          real_t* error) const;
  virtual void computerrorLinf(const Matrix& u_numerical_nodal, const Matrix& u_exact_nodal, 
                          real_t* error) const;
   
private:
  
  const TensorMesh1D* mesh1D_;
  BasisFunction1D* basis_;
  const int& qua_order1D_;
  const QuadratureType& quatype_;
  const int num_equations_;
  const int dim_ = 1;
  Eigen::MatrixXi Tm_;
  int NTdofs_;

  std::vector<Matrix> flux_u_v_;
  Matrix dv_u_;
  Matrix v_u_;
  DiagnalMatrix v_u_diag_;
  DiagnalMatrix v_u_diaginv_; 
  std::vector<Vector> boundary_u_;
  Matrix test_ref_;
  Matrix test_ref_T_;
  Matrix test_dx_ref_;
  DiagnalMatrix wqua_diag_;
  Matrix Pb_;

  void generateTm();
  void generateDG();
  void generatePb();

};

class fespace2D 
{
public:
  // **** 标量有限元函数空间 ***** // 
  fespace2D(const TensorMesh2D* mesh2D, BasisFunction2D* basis, 
            const int& qua_order1D, const QuadratureType& quatype);

  // **** 向量有限元函数空间 ***** // 
  fespace2D(const TensorMesh2D* mesh2D, BasisFunction2D* basis, 
            const int& qua_order1D, const QuadratureType& quatype, 
            const int num_equations);
  virtual ~fespace2D() = default;

  int NTH_;

  virtual void setNTH(const int& NTH);
  
  virtual void init();
  const int& getdim() const;
  const TensorMesh2D* getmesh2D() const;
  const BasisFunction2D* getbasis() const;
  const int& getnum_equations() const;
  const int& getnumqua() const;
  const int& getnumqua1d() const;
  const IntMatrix& getTm() const;
  const int& getNTdofs() const;
  const std::vector<std::vector<Matrix>>& getflux_u_v() const;
  const Matrix& getdxv_u() const;
  const Matrix& getdyv_u() const;
  const Matrix& getv_u() const;
  const DiagnalMatrix& getv_u_diag() const;
  const DiagnalMatrix& getv_u_diaginv() const;
  const Matrix& getCoorRef() const;
  const std::vector<Matrix>& getCoorBdrRef() const;
  const std::vector<Matrix>& getboundary_u() const;
  const Matrix& gettest_ref() const;
  const Matrix& gettest_ref_T() const;
  const Matrix& gettest_dx_ref() const;
  const Matrix& gettest_dy_ref() const;
  const DiagnalMatrix& getwqua_diag() const;
  const DiagnalMatrix& getwqua_diag1d() const;
  const std::vector<Matrix>& getPb() const;
  const int& getNTH() const;
  
  template<typename Func>
  void Project_Initial(Func&& func, Matrix* M) const
  {
    // ****** 传入相关变量  ******* //
    const int& ncell = mesh2D_->getncell();
    const int& polydim = basis_->getpolydim();
    // ************************** //
    Matrix f_quavalue(numqua_, ncell);
    f_quavalue = std::forward<Func>(func)(Pb_[0], Pb_[1]);
    nodal_to_modal2D(f_quavalue, M);
  };

  template<typename Func>
  void Project_Final(Func&& func, const real_t& Tstop, Matrix* M) const
  {
    // ****** 传入相关变量  ******* //
    const int& ncell = mesh2D_->getncell();
    const int& polydim = basis_->getpolydim();
    // ************************** //
    Matrix f_quavalue(numqua_, ncell);
    f_quavalue = std::forward<Func>(func)(Pb_[0], Pb_[1], Tstop);
    nodal_to_modal2D(f_quavalue, M);
  };

  template<typename Func>
  void Interpolate_Initial(Func&& func, Matrix* M) const
  {
    *M = std::forward<Func>(func)(Pb_[0], Pb_[1]);
  };

  template<typename Func>
  void Interpolate_Initial_Bdr(Func&& func, Matrix* M) const
  {
    // ****** 传入相关变量  ******* //
    const real_t& hx = mesh2D_->gethx();
    const real_t& hy = mesh2D_->gethy();
    const int& extboundaryNum = mesh2D_->getextboundaryNum();
    const std::vector<int>& ExtBTypeIndex = 
      mesh2D_->getextboundarytypeindex();
    const std::vector<int>& extNei = 
      mesh2D_->getextboundaryneighbors();
    const std::vector<Eigen::Vector2d> extboundarycenter = 
      mesh2D_->getextboundarycenter();
    const Matrix& CellCenter = mesh2D_->getCellCenter();
    // *************************** //

    M->resize(qua_order1D_, extboundaryNum);
    M->setZero();
    Matrix temp;
    real_t x0, y0;
    int cell_index;
    for (int i = 0; i < extboundaryNum; i++)
    {
      temp = CoorBdrRef_[ExtBTypeIndex[i]];
      cell_index = extNei[i];
      for (int q = 0; q < qua_order1D_; q++)
      {
        x0 = CellCenter(0, cell_index) + temp(0, q) * hx;
        y0 = CellCenter(1, cell_index) + temp(1, q) * hy;
        (*M)(q, i) = std::forward<Func>(func)(x0, y0);
      }
    }
  };


  template<typename Func>
  void Interpolate_Final(Func&& func, const real_t& Tstop, Matrix* M) const
  {
    *M = std::forward<Func>(func)(Pb_[0], Pb_[1], Tstop);
  };

  // T 为Matrix或者std::vector<Matrix>
  template <typename T>
  void nodal_to_modal2D(const T& N, T* M) const {
    Matrix temp = v_u_diaginv_ * test_ref_ * wqua_diag_;
    *M = temp * N;
  };

  // T 为Matrix或者std::vector<Matrix>
  template <typename T>
  void modal_to_nodal2D(const T& M, T* N) const {
    Matrix temp = test_ref_.transpose();
    *N = temp * M;
  };

  template <typename Derived>
  void nodal_to_modal2D(const Eigen::MatrixBase<Derived>& N, Matrix* M) const {
    Matrix temp = v_u_diaginv_ * test_ref_ * wqua_diag_;
    *M = temp * N;
  };

  template <typename Derived>
  void modal_to_nodal2D(const Eigen::MatrixBase<Derived>& M, Matrix* N) const {
    *N = test_ref_.transpose() * M;
  }

  // T 为Matrix或者std::vector<Matrix>
  template <typename T>
  void Assemble_F(const T& function_f_quavalue, 
                  const int& test_dx_order, const int& test_dy_order,
                  T* RHS) const {
    const Vector& Jx = mesh2D_->getJx();
    const Vector& Jy = mesh2D_->getJy();
    const Matrix* test_temp;
    if (test_dx_order == 0 && test_dy_order == 0)
    {
      test_temp = &test_ref_;
    } else if (test_dx_order == 1 && test_dy_order == 0)
    {
      test_temp = &test_dx_ref_;
    } else if (test_dx_order == 0 && test_dy_order == 1)
    {
      test_temp = &test_dy_ref_;
    } else 
    {
      QUEST_ERROR("unsupported test_dx_order and test_dy_order ! ");
    }
    real_t sx = (test_dx_order == 0) ? 1 : Jx(0);
    real_t sy = (test_dy_order == 0) ? 1 : Jy(0);
    Matrix test_qua = sx * sy * (*test_temp) * wqua_diag_;
    *RHS = test_qua * function_f_quavalue;
  };

  void Assemble_Flux(const std::vector<Matrix>& flux_int, const std::vector<Matrix>& flux_ext, 
                      std::vector<Matrix>* RHS) const;

  void Assemble_Flux(const Matrix& flux_int, const Matrix& flux_ext, 
                      Matrix* RHS) const;

  void Assemble_Flux_int(const Matrix& flux_int, Matrix* RHS) const;
  void Assemble_Flux_bc(const Matrix& flux_ext, Matrix* RHS) const;
  
  void DisplayResult(const Vector& numerical, const Vector& exact,  
                    const std::string& title, std::ostream& Outfile) const;

  void DisplayResult(const std::vector<Vector>& output, 
                    const std::vector<std::string>& outputname,
                    const std::string title, std::ostream& Outfile) const;

  virtual void computerrorL1(const Matrix& u_numerical_nodal, const Matrix& u_exact_nodal, 
                          real_t* error) const;
  virtual void computerrorL2(const Matrix& u_numerical_nodal, const Matrix& u_exact_nodal, 
                          real_t* error) const;
  virtual void computerrorLinf(const Matrix& u_numerical_nodal, const Matrix& u_exact_nodal, 
                          real_t* error) const;
   
private:
  
  const TensorMesh2D* mesh2D_;
  BasisFunction2D* basis_;
  const int& qua_order1D_;
  const QuadratureType& quatype_;
  const int num_equations_;
  const int dim_ = 2;
  Eigen::MatrixXi Tm_;
  int NTdofs_;
  
  Matrix CoorRef_;
  Vector wqua2d_;
  std::vector<Matrix> CoorBdrRef_;
  std::vector<std::vector<Matrix>> flux_u_v_;
  Matrix dxv_u_;
  Matrix dyv_u_;
  Matrix v_u_;
  DiagnalMatrix v_u_diag_;
  DiagnalMatrix v_u_diaginv_; 
  std::vector<Matrix> boundary_u_;
  int numqua_;
  Matrix test_ref_;
  Matrix test_ref_T_;
  Matrix test_dx_ref_;
  Matrix test_dy_ref_;
  DiagnalMatrix wqua_diag_;
  DiagnalMatrix wqua_diag1d_;
  std::vector<Matrix> Pb_;

  void generateTm();
  void generateDG();
  void generatePb();

};

}; // namespace QUEST 

#endif // QUEST_FESPACE_HPP 

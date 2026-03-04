#ifndef QUEST_DRIFT_DIFFUSION1D_IMEX_ACCTEST_HPP
#define QUEST_DRIFT_DIFFUSION1D_IMEX_ACCTEST_HPP

#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "hyperbolicproblemsbase.hpp"
#include "config.hpp"
#include "error.hpp"
#include "poissonsolver1D.hpp"

namespace QUEST
{

class DriftDiffusion1DIMEX: public HyperbolicProblems1DBase
{
protected:

  static real_t cfl_;
  static real_t T0_;
  static real_t vbias_;

  void init() override;
  void setdt(real_t* dt) override;
  void time_stepping(real_t& dt) override;
  void update() override;

public:
  DriftDiffusion1DIMEX(const fespace1D* fe, const int& t_order,
                      PoissonSolver1D* poisol, ImplicitLinearDiffusionSolver* diff);
  ~DriftDiffusion1DIMEX() = default;

  static real_t rho_init(const real_t& x);
  static real_t phi_init(const real_t& x);
  static real_t rhod_init(const real_t& x);  // doping profile

  virtual void fluxint_compute();
  virtual void fluxext_compute();
  virtual void fk_compute();
};


// r_t  = D * ∇·q
// ∫_T r_t:v = D * ( - ∫_T q:∇v + ∫_e q⋅[[v]] )
// ∫_T q:w = D * ( - ∫_T r∇·w + ∫_e r⋅[[w]] )
// M1 * r_t = A * q
// M2 * q = B * r
// M1 * r_t = S * r = ( A * M2^{-1} * B ) * r
enum DiffusionType {
  LeftRight, // \q = q^-
  Center, // \q = 1/2 * (q^- + q^+), 
  RightLeft,  //  \q = q^+
}

inline real_t sign(real_t x);

// 这个是针对迪利克雷边界条件
class LinearDiffusionSolver1D
{
protected:
  
  int num_equations_;
  int intboundaryNum_;
  int extboundaryNum_;
  int dim_;
  int polydim_;
  int numqua_;
  int ncell_;
  int hx_;
  const fespace1D* fe_;
  const real_t& D_;
  const DiffusionType& diftype_;
  SparseMatrix M1_;
  SparseMatrix A_;
  SparseMatrix M2_;
  SparseMatrix M2inv_;
  SparseMatrix B_;
  SparseMatrix S_;

  void generateM1();
  virtual void generateA();
  void generateM2();
  void generateM2inv();
  virtual void generateB();
  void generateS();

  void init(const PoissonSolver1DType& poitype,
      const real_t& tol);
public:
  LinearDiffusionSolver1D(const fespace1D* fe, const real_t& D,
                        const DiffusionType& diftype);

  init();
  ~LinearDiffusionSolver1D() = default;
};

LinearDiffusionSolver1D::LinearDiffusionSolver1D(const fespace1D* fe, 
  const real_t& D, const DiffusionType& diftype)
  : fe_(fe), D_(D), diftype_(diftype) {
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
  generateM1();
  generateA();
  generateM2();
  generateM2inv();
  generateB();
  generateS();
}

real_t sign(real_t x) {
  if (x > 0) return 1.e0;
  if (x < 0) return -1.e0;
  return 0.e0;
};

void LinearDiffusionSolver1D::generateM1() {
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

  M1_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_M1;
  int estimatedNonZeros = NTdofs_;
  tripletList_M1.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u_(basis_index,basis_index) * JacobiDet_(i);
      int alpha = alpha_start + basis_index;
      tripletList_M1.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  M1_.setFromTriplets(tripletList_M1.begin(), tripletList_M1.end());
}

void LinearDiffusionSolver1D::generateM2() {
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

  M2_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_M2;
  int estimatedNonZeros = NTdofs_;
  tripletList_M2.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u_(basis_index,basis_index) * JacobiDet_(i);
      int alpha = alpha_start + basis_index;
      tripletList_M2.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  M2_.setFromTriplets(tripletList_M2.begin(), tripletList_M2.end());
}

void LinearDiffusionSolver1D::generateM2inv() {
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

  M2inv_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_M2inv;
  int estimatedNonZeros = NTdofs_;
  tripletList_M2inv.reserve(estimatedNonZeros);

  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int basis_index = 0; basis_index < polydim_; basis_index++) {
      real_t ini_value = v_u_(basis_index,basis_index) * JacobiDet_(i);
      ini_value = 1.e0 / ini_value;
      int alpha = alpha_start + basis_index;
      tripletList_M2inv.push_back(Eigen::Triplet<real_t>(alpha, alpha, ini_value));
    };
  };

  M2_.setFromTriplets(tripletList_M2.begin(), tripletList_M2.end());
}

void LinearDiffusionSolver1D::generateA() {
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
  const std::vector<Eigen::Vector2i>& extNei_period_ = 
    fe_->getmesh1D()->getextboundaryneighbors_period();
  // ******************************** //

  A_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_A;
  int estimatedNonZeros = NTdofs_;
  tripletList_A.reserve(estimatedNonZeros);

  Matrix Aref;
  Aref = dv_u_;
  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        real_t ini_value = Aref(test_basis_index, trial_basis_index) * JacobiDet_(i) * Jx_(i);
        ini_value = - ini_value;
        tripletList_A.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };

  real_t v;
  switch (diftype_) {
    case LeftRight: v = 1.e0;  break;
    case Center: v = 0.e0;  break;
    case RightLeft: v = -1.e0;  break;
    default: break;
  }

  for (int i = 0; i < intboundaryNum_; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = intNei_[i](test_cell);
      Vector trial_qua_value = boundary_u_[test_cell];
      real_t test_normal = intnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = intNei_[i](trial_cell);
        Vector test_qua_value = boundary_u_[trial_cell];
        real_t trial_normal = intnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value = (0.5e0 * trial_qua_value(trial_basis_index) + 
                0.5e0 * sign(v * trial_normal) * trial_qua_value(trial_basis_index))
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = D_ * ini_value;
            tripletList_A.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_; i++) {

    int test_cell_Index = extNei_[i];
    Vector trial_qua_value = boundary_u_[1-i];
    real_t test_normal = extnormals_[i];

    int trial_cell_Index = extNei_[i];
    Vector test_qua_value = boundary_u_[1-i];
    real_t trial_normal = extnormals_[i];

    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      int alpha = Tm_(test_basis_index, test_cell_Index);
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        int beta = Tm_(trial_basis_index, trial_cell_Index);
        real_t ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) + 
                0.5e0 * sign(v * trial_normal) * trial_qua_value(trial_basis_index))
                * test_qua_value(test_basis_index) * test_normal;
        ini_value = D_ * ini_value;
        tripletList_A.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };
  A_.setFromTriplets(tripletList_A.begin(), tripletList_A.end());
};

void LinearDiffusionSolver1D::generateB() {
  // ********* 传入相关变量 ********** //
  const int& NTdofs_ = fe_->getNTdofs();
  const DiagnalMatrix& wqua_diag_ = fe_->getwqua_diag();
  const Matrix& dv_u_ = fe_->getdv_u();
  const Matrix& v_u_ = fe_->getv_u();
  const Vector& JacobiDet_ = fe_->getmesh1D()->getJacobiDet();
  const Vector& Jx_ = fe_->getmesh1D()->getJx();
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const IntMatrix& Tm_ = fe_->getTm();
  const std::vector<real_t>& intnormals_ = fe_->getmesh1D()->getintboundarynormal();
  const std::vector<Eigen::Vector2i>& intNei_ = fe_->getmesh1D()->getintboundaryneighbors();
  const std::vector<real_t>& extnormals_ = fe_->getmesh1D()->getextboundarynormal();
  const std::vector<int>& extNei_ = fe_->getmesh1D()->getextboundaryneighbors();
  const std::vector<Eigen::Vector2i>& extNei_period_ = fe_->getmesh1D()->getextboundaryneighbors_period();
  // ******************************** //

  B_.resize(NTdofs_, NTdofs_);
  std::vector<Eigen::Triplet<real_t>> tripletList_B;
  int estimatedNonZeros = NTdofs_;
  tripletList_B.reserve(estimatedNonZeros);

  Matrix Bref;
  Bref = dv_u_;
  for (int i = 0; i < ncell_; i++) {
    int alpha_start = i * polydim_;
    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        real_t ini_value = Bref(test_basis_index, trial_basis_index) * JacobiDet_(i) * Jx_(i);
        ini_value = - ini_value;
        tripletList_B.push_back(Eigen::Triplet<real_t>(alpha_start + test_basis_index, 
                        alpha_start + trial_basis_index, ini_value));
      };
    };
  };

  real_t v;
  switch (diftype_) {
    case LeftRight: v = - 1.e0;  break;
    case Center: v = 0.e0;  break;
    case RightLeft: v = 1.e0;  break;
    default: break;
  }

  for (int i = 0; i < intboundaryNum_; i++) {
    for (int test_cell = 0; test_cell < 2; test_cell++) {
      int test_cell_Index = intNei_[i](test_cell);
      Vector trial_qua_value = boundary_u_[test_cell];
      real_t test_normal = intnormals_[i] * std::pow(-1.e0,test_cell);
      for (int trial_cell = 0; trial_cell < 2; trial_cell++) {
        int trial_cell_Index = intNei_[i](trial_cell);
        Vector test_qua_value = boundary_u_[trial_cell];
        real_t trial_normal = intnormals_[i] * std::pow(-1.e0,trial_cell);
        for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
          int alpha = Tm_(test_basis_index, test_cell_Index);
          for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
            int beta = Tm_(trial_basis_index, trial_cell_Index);
            real_t ini_value = (0.5e0 * trial_qua_value(trial_basis_index) + 
                0.5e0 * sign(v * trial_normal) * trial_qua_value(trial_basis_index))
                * test_qua_value(test_basis_index) * test_normal;
            ini_value = D_ * ini_value;
            tripletList_A.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
          };
        };
      };
    };
  };

  for (int i = 0; i < extboundaryNum_; i++) {

    int test_cell_Index = extNei_[i];
    Vector trial_qua_value = boundary_u_[1-i];
    real_t test_normal = extnormals_[i];

    int trial_cell_Index = extNei_[i];
    Vector test_qua_value = boundary_u_[1-i];
    real_t trial_normal = extnormals_[i];

    for (int test_basis_index = 0; test_basis_index < polydim_; test_basis_index++) {
      int alpha = Tm_(test_basis_index, test_cell_Index);
      for (int trial_basis_index = 0; trial_basis_index < polydim_; trial_basis_index++) {
        int beta = Tm_(trial_basis_index, trial_cell_Index);
        real_t ini_value =  (0.5e0 * trial_qua_value(trial_basis_index) + 
                0.5e0 * sign(v * trial_normal) * trial_qua_value(trial_basis_index))
                * test_qua_value(test_basis_index) * test_normal;
        ini_value = D_ * ini_value;
        tripletList_C.push_back(Eigen::Triplet<real_t>(alpha, beta, ini_value));
      };
    };
  };
};

} // namespace QUEST

#endif  // QUEST_DRIFT_DIFFUSION1D_IMEX_ACCTEST_HPP

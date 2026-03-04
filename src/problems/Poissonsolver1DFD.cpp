#include "Poissonsolver1DFD.hpp"

namespace QUEST
{

Poissonsolver1DFD_Base::Poissonsolver1DFD_Base(const FDmesh* mesh1D,
                    const int& x_order)
  : mesh1D_(mesh1D), x_order_(x_order) {};

const SparseMatrix& Poissonsolver1DFD_Base::getS() const {
  return S_;
};

void Poissonsolver1DFD_Base::generateS() {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***

  S_.resize(Nx, Nx);
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = 1.0 / (hx * hx);

  for (int i = 0; i < Nx; i++) {
    tripletList_S.emplace_back(i, i, 2.0 * coef);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, -1.0 * coef);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, -1.0 * coef);
    } 
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());
};

void Poissonsolver1DFD_Base::init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol) {
  generateS();
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::CG:
    std::cout << "  Preparation for Conjugate Gradient the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    cg_.compute(S_);
    cg_.setTolerance(tol);
    std::cout << "  The end of preparation for Conjugate Gradient method ! " << std::endl;
    break;

  case PoissonSolver1DTypeFD::PCG:
    QUEST_ERROR(" The PCG method is not been implemented ! ");
    break;
  
  case PoissonSolver1DTypeFD::LDLT:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LDLT the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LDLT the matrix by Pardiso ......\n";
#endif 
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    chol_.analyzePattern(S_);
    chol_.factorize(S_);
    std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  TOC;
};

void Poissonsolver1DFD_Base::Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***
  const real_t coef = 1.0 / (hx * hx);
  *bf = f;
  (*bf)(0) += coef * bc_L;
  (*bf)(Nx - 1) += coef * bc_R;
};

void Poissonsolver1DFD_Base::Solve(const Vector& bf,
            Vector* solution) {
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::CG:
    std::cout << "  Solve the Poisson by CG ......\n";
    *solution = cg_.solve(bf);
    std::cout << "  iterations = " << cg_.iterations() << std::endl;
    std::cout << "  eerror = " << cg_.error()      << std::endl;
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DTypeFD::LDLT:
    std::cout << "  Solve the Poisson by LDLT ......\n";
    *solution = chol_.solve(bf);
    std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;
}

void Poissonsolver1DFD_Base::SolveAll(const real_t& bc_L,
                                const real_t& bc_R,
                                const Vector& f,
                                Vector* solution,
                                Vector* solution_dx) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const Vector& xpoints = mesh1D_->getxpoints();
  // ***
  Vector bf;
  Assemble_bf(bc_L, bc_R, f, &bf);
  Solve(bf, solution);
  solution_dx->resize(Nx + 2);
  solution_dx->setZero();
  
  real_t chx = 2.e0 * hx;
  for (int i = 2; i <= Nx - 1; i++) {
    (*solution_dx)(i) = ((*solution)(i) - (*solution)(i - 2)) / chx;
  }
  (*solution_dx)(1) = ((*solution)(1) - bc_L) / chx;
  (*solution_dx)(Nx) = (bc_R - (*solution)(Nx - 2)) / chx;

  if (x_order_ == 1) {
    (*solution_dx)(0) = ((*solution)(0) - bc_L) / hx;
    (*solution_dx)(Nx + 1) = (bc_R - (*solution)(Nx - 1)) / hx;
  } else {
    QUEST_ERROR(" The 2nd-order or higher-order approximation for u_x has not been implemented ! ");
  }
};

Poissonsolver1DFD_Base_period::Poissonsolver1DFD_Base_period(const FDmesh_period* mesh1D, 
                    const int& x_order) :
  Poissonsolver1DFD_Base(mesh1D, x_order) {};

void Poissonsolver1DFD_Base_period::generateS() {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***

  S_.resize(Nx, Nx);
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 4 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = 1.0 / (hx * hx);

  for (int i = 0; i < Nx; i++) {
    tripletList_S.emplace_back(i, i, 2.0 * coef);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, -1.0 * coef);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, -1.0 * coef);
    } 
  }
  tripletList_S.emplace_back(0, Nx-1, -1.0 * coef);
  tripletList_S.emplace_back(Nx-1, 0, -1.0 * coef);

  int center = Nx / 2;
  real_t value = 1.e0 / hx;
  for (int i = 0; i < Nx; i++) {
    tripletList_S.emplace_back(center, i, value);
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());
};

void Poissonsolver1DFD_Base_period::init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol) {
  generateS();
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::BICG:
    std::cout << "  Preparation for BICG the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(S_);
    bicg_.setTolerance(tol);
    std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  
  case PoissonSolver1DTypeFD::LU:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LU the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
#endif 
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.factorize(S_);
    std::cout << "  The end of preparation for LU method ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  TOC;
};

void Poissonsolver1DFD_Base_period::Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***
  *bf = f;
};

void Poissonsolver1DFD_Base_period::Solve(const Vector& bf,
            Vector* solution) {
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::BICG:
    // std::cout << "  Solve the Poisson by CG ......\n";
    *solution = bicg_.solve(bf);
    // std::cout << "  iterations = " << bicg_.iterations() << std::endl;
    // std::cout << "  eerror = " << bicg_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DTypeFD::LU:
    // std::cout << "  Solve the Poisson by LDLT ......\n";
    *solution = lu_.solve(bf);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;
};

void Poissonsolver1DFD_Base_period::SolveAll(const real_t& bc_L,
                                const real_t& bc_R,
                                const Vector& f,
                                Vector* solution,
                                Vector* solution_dx) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const Vector& xpoints = mesh1D_->getxpoints();
  // ***
  Vector bf;
  Assemble_bf(bc_L, bc_R, f, &bf);
  Solve(bf, solution);
  solution_dx->resize(Nx);
  solution_dx->setZero();
  
  real_t chx = 2.e0 * hx;
  for (int i = 1; i < Nx-1; i++) {
    (*solution_dx)(i) = ((*solution)(i+1) - (*solution)(i-1)) / chx;
  }
  (*solution_dx)(0) = ((*solution)(1) - (*solution)(Nx-1)) / chx;
  (*solution_dx)(Nx-1) = ((*solution)(0) - (*solution)(Nx-2)) / chx;

};

Poissonsolver1DFD::Poissonsolver1DFD(const KineticFDmesh* mesh1D,
                    const int& x_order)
  : mesh1D_(mesh1D), x_order_(x_order) {};

const SparseMatrix& Poissonsolver1DFD::getS() const {
  return S_;
};

void Poissonsolver1DFD::generateS() {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***

  S_.resize(Nx, Nx);
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = 1.0 / (hx * hx);

  for (int i = 0; i < Nx; i++) {
    tripletList_S.emplace_back(i, i, 2.0 * coef);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, -1.0 * coef);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, -1.0 * coef);
    } 
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());
};

void Poissonsolver1DFD::init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol) {
  generateS();
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::CG:
    std::cout << "  Preparation for Conjugate Gradient the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    cg_.compute(S_);
    cg_.setTolerance(tol);
    std::cout << "  The end of preparation for Conjugate Gradient method ! " << std::endl;
    break;

  case PoissonSolver1DTypeFD::PCG:
    QUEST_ERROR(" The PCG method is not been implemented ! ");
    break;
  
  case PoissonSolver1DTypeFD::LDLT:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LDLT the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LDLT the matrix by Pardiso ......\n";
#endif 
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    chol_.analyzePattern(S_);
    chol_.factorize(S_);
    std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  TOC;
};

void Poissonsolver1DFD::Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***
  const real_t coef = 1.0 / (hx * hx);
  *bf = f;
  (*bf)(0) += coef * bc_L;
  (*bf)(Nx - 1) += coef * bc_R;
};

void Poissonsolver1DFD::Solve(const Vector& bf,
            Vector* solution) {
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::CG:
    // std::cout << "  Solve the Poisson by CG ......\n";
    *solution = cg_.solve(bf);
    // std::cout << "  iterations = " << cg_.iterations() << std::endl;
    // std::cout << "  eerror = " << cg_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DTypeFD::LDLT:
    // std::cout << "  Solve the Poisson by LDLT ......\n";
    *solution = chol_.solve(bf);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;
}

void Poissonsolver1DFD::SolveAll(const real_t& bc_L,
                                const real_t& bc_R,
                                const Vector& f,
                                Vector* solution,
                                Vector* solution_dx) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const Vector& xpoints = mesh1D_->getxpoints();
  // ***
  Vector bf;
  Assemble_bf(bc_L, bc_R, f, &bf);
  Solve(bf, solution);
  solution_dx->resize(Nx + 2);
  solution_dx->setZero();
  
  real_t chx = 2.e0 * hx;
  for (int i = 2; i <= Nx - 1; i++) {
    (*solution_dx)(i) = ((*solution)(i) - (*solution)(i - 2)) / chx;
  }
  (*solution_dx)(1) = ((*solution)(1) - bc_L) / chx;
  (*solution_dx)(Nx) = (bc_R - (*solution)(Nx - 2)) / chx;

  if (x_order_ == 1) {
    (*solution_dx)(0) = ((*solution)(0) - bc_L) / hx;
    (*solution_dx)(Nx + 1) = (bc_R - (*solution)(Nx - 1)) / hx;
  } else {
    QUEST_ERROR(" The 2nd-order or higher-order approximation for u_x has not been implemented ! ");
  }
};

Poissonsolver1DFD_period::Poissonsolver1DFD_period(const KineticFDmesh_period* mesh1D, 
                    const int& x_order) :
  Poissonsolver1DFD(mesh1D, x_order) {};

void Poissonsolver1DFD_period::generateS() {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***

  S_.resize(Nx, Nx);
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 4 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = 1.0 / (hx * hx);

  for (int i = 0; i < Nx; i++) {
    tripletList_S.emplace_back(i, i, 2.0 * coef);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, -1.0 * coef);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, -1.0 * coef);
    } 
  }
  tripletList_S.emplace_back(0, Nx-1, -1.0 * coef);
  tripletList_S.emplace_back(Nx-1, 0, -1.0 * coef);

  int center = Nx / 2;
  real_t value = 1.e0 / hx;
  for (int i = 0; i < Nx; i++) {
    tripletList_S.emplace_back(center, i, value);
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());
};

void Poissonsolver1DFD_period::init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol) {
  generateS();
  poitype_ = poitype;
  TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::BICG:
    // std::cout << "  Preparation for BICG the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(S_);
    bicg_.setTolerance(tol);
    // std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  
  case PoissonSolver1DTypeFD::LU:
// #ifndef QUEST_USE_MKL
//     std::cout << "  Preparation for LU the matrix by Eigen ......\n";
// #else 
//     std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
// #endif 
//     std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.factorize(S_);
    // std::cout << "  The end of preparation for LU method ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  TOC;
};

void Poissonsolver1DFD_period::Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  // ***
  *bf = f;
};

void Poissonsolver1DFD_period::Solve(const Vector& bf,
            Vector* solution) {
  // TIC;
  switch (poitype_)
  {
  case PoissonSolver1DTypeFD::BICG:
    // std::cout << "  Solve the Poisson by CG ......\n";
    *solution = bicg_.solve(bf);
    // std::cout << "  iterations = " << bicg_.iterations() << std::endl;
    // std::cout << "  eerror = " << bicg_.error()      << std::endl;
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  case PoissonSolver1DTypeFD::LU:
    // std::cout << "  Solve the Poisson by LDLT ......\n";
    *solution = lu_.solve(bf);
    // std::cout << "  The end of Solving ! " << std::endl;
    break;
  
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
  // TOC;
};

void Poissonsolver1DFD_period::SolveAll(const real_t& bc_L,
                                const real_t& bc_R,
                                const Vector& f,
                                Vector* solution,
                                Vector* solution_dx) {
  // *** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const Vector& xpoints = mesh1D_->getxpoints();
  // ***
  Vector bf;
  Assemble_bf(bc_L, bc_R, f, &bf);
  Solve(bf, solution);
  solution_dx->resize(Nx);
  solution_dx->setZero();
  
  real_t chx = 2.e0 * hx;
  for (int i = 1; i < Nx-1; i++) {
    (*solution_dx)(i) = ((*solution)(i+1) - (*solution)(i-1)) / chx;
  }
  (*solution_dx)(0) = ((*solution)(1) - (*solution)(Nx-1)) / chx;
  (*solution_dx)(Nx-1) = ((*solution)(0) - (*solution)(Nx-2)) / chx;

};

std::string name(const PoissonSolver1DTypeFD& poitype) {
  switch (poitype)
  {
  case PoissonSolver1DTypeFD::CG: return "CG"; break;
  case PoissonSolver1DTypeFD::PCG: return "PCG"; break;
  case PoissonSolver1DTypeFD::LDLT: return "LDLT"; break;
  case PoissonSolver1DTypeFD::LU: return "LU"; break;
  case PoissonSolver1DTypeFD::GMRES: return "GMRES"; break;
  case PoissonSolver1DTypeFD::BICG: return "BICG"; break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    return "Unknown !";
    break;
  } 
};

};
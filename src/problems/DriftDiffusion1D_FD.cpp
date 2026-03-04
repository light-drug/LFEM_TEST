#include "DriftDiffusion1D_FD.hpp"

namespace QUEST
{
DriftDiffusion1D_FD::DriftDiffusion1D_FD(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order) :
    mesh1D_(mesh1D), poisol_(poisol), x_order_(x_order), t_order_(t_order) {};

void DriftDiffusion1D_FD::setcfl(const real_t& cfl) {
  cfl_ = cfl;
};

void DriftDiffusion1D_FD::setpi(const real_t& pi) {
  pi_ = pi;
};

void DriftDiffusion1D_FD::setgamma(const real_t& gamma) {
  gamma_ = gamma;
};

void DriftDiffusion1D_FD::setsparsetol(const real_t& sparsetol) {
  sparsetol_ = sparsetol;
};
  
void DriftDiffusion1D_FD::setiterationtol(const real_t& iterationtol) {
  iterationtol_ = iterationtol;
};

const real_t& DriftDiffusion1D_FD::getcfl() const {
  return cfl_;
};

const real_t& DriftDiffusion1D_FD::getgamma() const {
  return gamma_;
};


real_t DriftDiffusion1D_FD::rho_init(const real_t& x) {
  return 1.e0;
};

Vector DriftDiffusion1D_FD::rho_init(const Vector& x) {
  Vector rho = x;
  rho.setConstant(1.e0);
  return rho;
};

real_t DriftDiffusion1D_FD::phi_bc(const real_t& x, const real_t & t) {
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

real_t DriftDiffusion1D_FD::rho_bc(const real_t& x, const real_t & t) {
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t rho;
  if (std::abs(x - x1) < 1.e-10) {
    return 1.e0;
  } else if (std::abs(x - x2) < 1.e-10) {
    return 1.e0;
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return rho;
};

real_t DriftDiffusion1D_FD::rho_d(const real_t& x) {
  return 1.e0 - (1.e0 - 0.001)/2 * (std::tanh((x - 0.3e0)/0.02e0) - std::tanh((x - 0.7e0)/0.02e0));
};

Vector DriftDiffusion1D_FD::rho_d(const Vector& x) {
  Vector x1 = (x.array() - 0.3) / 0.02;
  Vector x2 = (x.array() - 0.7) / 0.02;
  Vector rhod = 1.0 - (1.0 - 0.001) / 2.0 * (x1.array().tanh() - x2.array().tanh());
  return rhod;
};

real_t DriftDiffusion1D_FD::source(const real_t& x,  const real_t& t) {
  return 0.e0;
};

Vector DriftDiffusion1D_FD::source(const Vector& x,  const real_t& t) {
  Vector S = x;
  S.setZero();
  return S;
};

const Vector& DriftDiffusion1D_FD::getrho() const {
  return rho_;
};

const Vector& DriftDiffusion1D_FD::getE() const {
  return E_;
};

const Vector& DriftDiffusion1D_FD::getrhod() const {
  return rhod_;
};

const Vector& DriftDiffusion1D_FD::getphi() const {
  return phi_;
};

Vector DriftDiffusion1D_FD::getrhoinit() {
  const Vector& X = mesh1D_->getX();

  Vector rho = rho_init(X);
  return rho;
};

void DriftDiffusion1D_FD::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  rho_ = rho_init(X);
  rhod_ = rho_d(X);
  soltype_ = soltype;
  std::cout << " gamma_ = " << gamma_ << std::endl;
  std::cout << " cfl_ = " << cfl_ << std::endl;
};

void DriftDiffusion1D_FD::setdt(real_t* dt) {
  // ** 传入相关变量
  const real_t& hx = mesh1D_->gethx();
  // **
  *dt = cfl_ * hx;
};

void DriftDiffusion1D_FD::update_E(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  Vector bf = - (rho_ - rhod_) / gamma_;
  real_t bc_L = phi_bc(x1, Trun + dt);
  real_t bc_R = phi_bc(x2, Trun + dt);
  poisol_->SolveAll(bc_L, bc_R, bf, &phi_, &E_);
};

void DriftDiffusion1D_FD::update_rho(const real_t& Trun, 
                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  rho_L_ = rho_bc(x1, Trun + dt);
  rho_R_ = rho_bc(x2, Trun + dt);

  Source_ = source(X, Trun + dt) * dt;

  // 开始迭代
  real_t err = 1.e0;
  Vector bf;
  rhon_ = rho_;
  iter_ = 0;
  Vector prerhotemp = rho_;
  while (err > iterationtol_) {
    update_E(Trun, dt);
    generateS(Trun, dt);
    generatebf(Trun, dt, &bf);
    Vector rhotemp = rho_;
    SolveS(Trun, bf);  
    if (!rho_.allFinite() || !bf.allFinite()) {
      std::cout << "[ERROR] rho_ contains NaN or Inf after SolveS!" << std::endl;
      std::ofstream Sout("result/DD_acctest/S.txt");
      for (int k = 0; k < S_.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(S_, k); it; ++it) {
          Sout << it.row() << " " 
               << it.col() << " "
               << it.value() << "\n";
        }
      }
      std::cout << " rho_.allFinite() = " << rho_.allFinite() << std::endl;
      std::cout << " bf.allFinite() = " << bf.allFinite() << std::endl;
      Sout.close();
      std::abort();
    }
    iter_++;
    rhotemp = rhotemp - rho_;
    err = rhotemp.cwiseAbs().maxCoeff();
    std::cout <<" size = " << rho_.size() << ", err = " << std::setprecision(10) << err << ", iter = " << iter_ << std::endl;
    prerhotemp = rhotemp;
    // if (iter_ > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort." << std::endl;
    //   break;
    // }
  }
  std::cout << "err = " << err << std::endl;
};

void DriftDiffusion1D_FD::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + 2.e0 * coef - phi_(i) * coef;
    phi_temp_L = (i == 0) ? phi_L : phi_(i - 1);
    val = val + phi_temp_L / 2.e0 * coef;
    phi_temp_R = (i == Nx-1) ? phi_R : phi_(i + 1);
    val = val + phi_temp_R / 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - coef - (phi_(i) - phi_(i-1)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - coef + (phi_(i+1) - phi_(i)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    std::cout << "  Preparation for BICG the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(S_);
    bicg_.setTolerance(sparsetol_);
    std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  case Solver1DTypeFD::LU:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LU the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
#endif 
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.factorize(S_);
    std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  case Solver1DTypeFD::GMRES:
    gmres_.setTolerance(sparsetol_);
    gmres_.set_restart(30);
    gmres_.setMaxIterations(Nx);
    gmres_.compute(S_);
    break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
};

void DriftDiffusion1D_FD::SolveS(const real_t& Trun, const Vector& bf) {
  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    rho_ = bicg_.solve(bf);
    if (bicg_.info() != Eigen::Success) {
      std::cout << "[ERROR] BICG failed to converge!\n";
    }
    break;

  case Solver1DTypeFD::LU:
// #ifndef QUEST_USE_MKL
//     std::cout << "  Preparation for LU the matrix by Eigen ......\n";
// #else 
//     std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
// #endif 
    rho_ = lu_.solve(bf);
    break;

  case Solver1DTypeFD::GMRES:
    rho_ = gmres_.solve(bf);
    std::cout << "Iterations: " << gmres_.iterations() << std::endl;
    std::cout << "Estimated error: " << gmres_.error() << std::endl;
    if (gmres_.info() != Eigen::Success) {
      std::cout << "[ERROR] GMRES failed to converge!\n";
    }
    break;
  default:
    QUEST_ERROR(" The Dift Diffusion Solver is not been implemented ! ");
    break;
  }
};

// the right hand side of the drift-diffusion equation
void DriftDiffusion1D_FD::generatebf(const real_t& Trun, 
                                    const real_t& dt,
                                    Vector* bf) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t coef = dt / (hx * hx);
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);

  *bf = rhon_ + Source_;
  (*bf)(0) = (*bf)(0) + coef * rho_L_;
  (*bf)(0) = (*bf)(0) + (phi_(0) - phi_L) / 2.e0 * rho_L_ * coef;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) - (phi_R - phi_(Nx-1)) / 2.e0 * rho_R_ * coef;
};

void DriftDiffusion1D_FD::updateAll(const real_t& Trun, const real_t& dt) {
  update_rho(Trun, dt);
};

DriftDiffusion1D_FD_PNjunction::DriftDiffusion1D_FD_PNjunction(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order) :
  DriftDiffusion1D_FD(mesh1D, poisol, x_order, t_order) {};

real_t DriftDiffusion1D_FD_PNjunction::rho_init(const real_t& x) {
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector DriftDiffusion1D_FD_PNjunction::rho_init(const Vector& x)
{
  // ****** 传入相关变量 ****** //
  const int Nx = x.size();
  // ************************ //
  Vector rho(Nx);
  for (int i = 0; i < Nx; i++) 
  {
    rho(i) = (x(i) < 0.5e0) ? 0.9e0 : 0.1e0;
  }
  return rho;
};

real_t DriftDiffusion1D_FD_PNjunction::phi_bc(const real_t& x, const real_t & t)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t phi;
  if (std::abs(x - x1) < 1.e-10) {
    return 0.e0;
  } else if (std::abs(x - x2) < 1.e-10) {
    return 4.e0;
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return phi;
};

real_t DriftDiffusion1D_FD_PNjunction::rho_bc(const real_t& x, const real_t & t)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t rho;
  if (std::abs(x - x1) < 1.e-10) {
    return 0.9e0;
  } else if (std::abs(x - x2) < 1.e-10) {
    return 0.1e0;
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return rho;
};

real_t DriftDiffusion1D_FD_PNjunction::rho_d(const real_t& x)
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector DriftDiffusion1D_FD_PNjunction::rho_d(const Vector& x)
{
  // ****** 传入相关变量 ****** //
  const int Nx = x.size();
  // ************************ //
  Vector rho(Nx);
  for (int i = 0; i < Nx; i++) 
  {
    rho(i) = (x(i) < 0.5e0) ? 0.9e0 : 0.1e0;
  }
  return rho;
};

real_t DriftDiffusion1D_FD_PNjunction::source(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector DriftDiffusion1D_FD_PNjunction::source(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};



DriftDiffusion1D_FD_penaltyIter::DriftDiffusion1D_FD_penaltyIter(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order)
  : DriftDiffusion1D_FD(mesh1D, poisol, x_order, t_order) {};

void DriftDiffusion1D_FD_penaltyIter::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + 2.e0 * coef - phi_(i) * coef;
    phi_temp_L = (i == 0) ? phi_L : phi_(i - 1);
    val = val + phi_temp_L / 2.e0 * coef;
    phi_temp_R = (i == Nx-1) ? phi_R : phi_(i + 1);
    val = val + phi_temp_R / 2.e0 * coef;

    val = val + mu_ / gamma_;   // penalty term
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - coef - (phi_(i) - phi_(i-1)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - coef + (phi_(i+1) - phi_(i)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    std::cout << "  Preparation for BICG the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(S_);
    bicg_.setTolerance(sparsetol_);
    std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  case Solver1DTypeFD::LU:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LU the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
#endif 
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.factorize(S_);
    std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  case Solver1DTypeFD::GMRES:
    gmres_.setTolerance(sparsetol_);
    gmres_.set_restart(30);
    gmres_.setMaxIterations(Nx);
    gmres_.compute(S_);
    break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
};

void DriftDiffusion1D_FD_penaltyIter::generatebf(const real_t& Trun, 
                                    const real_t& dt,
                                    Vector* bf) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t coef = dt / (hx * hx);
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);

  *bf = rhon_ + Source_;
  (*bf)(0) = (*bf)(0) + coef * rho_L_;
  (*bf)(0) = (*bf)(0) + (phi_(0) - phi_L) / 2.e0 * rho_L_ * coef;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) - (phi_R - phi_(Nx-1)) / 2.e0 * rho_R_ * coef;

  // penalty term
  (*bf) = (*bf) + mu_ / gamma_ * rho_;
};

void DriftDiffusion1D_FD_penaltyIter::setmu(const real_t& mu) {
  mu_ = mu;
};

const real_t& DriftDiffusion1D_FD_penaltyIter::getmu() const {
  return mu_;
};

DriftDiffusion1D_FD_penaltyImplicit::DriftDiffusion1D_FD_penaltyImplicit(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order) :
  DriftDiffusion1D_FD(mesh1D, poisol, x_order, t_order) {};

void DriftDiffusion1D_FD_penaltyImplicit::init(const Solver1DTypeFD& soltype) 
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  rho_ = rho_init(X);
  rhod_ = rho_d(X);
  soltype_ = soltype;

  // 需要求解初始泊松方程
  Vector bf = - (rho_ - rhod_) / gamma_;
  real_t bc_L = phi_bc(x1, 0.e0);
  real_t bc_R = phi_bc(x2, 0.e0);
  poisol_->SolveAll(bc_L, bc_R, bf, &phi_, &E_);

  std::cout << " gamma_ = " << gamma_ << std::endl;
  std::cout << " cfl_ = " << cfl_ << std::endl;
}

void DriftDiffusion1D_FD_penaltyImplicit::update_rho(const real_t& Trun, 
                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const SparseMatrix& Poi = poisol_->getS();
  // *** 
  rho_L_ = rho_bc(x1, Trun + dt);
  rho_R_ = rho_bc(x2, Trun + dt);

  Source_ = source(X, Trun + dt) * dt;
  rhon_ = rho_;

  generateS(Trun, dt);
  generateB(Trun, dt);
  update_E(Trun, dt);
  
  real_t coef = 1.e0 / (hx * hx);
  rho_ = rhod_ - gamma_ * Poi * phi_;
  rho_(0) += gamma_ * phi_bc(x1, Trun + dt) * coef;
  rho_(Nx - 1) += gamma_ * phi_bc(x2, Trun + dt) * coef;

}

void DriftDiffusion1D_FD_penaltyImplicit::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - coef;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - coef;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());
};

void DriftDiffusion1D_FD_penaltyImplicit::generatebf(const real_t& Trun, 
                                    const real_t& dt,
                                    Vector* bf) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t coef = dt / (hx * hx);
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);

  *bf = rhon_ + Source_;
  real_t phi_temp_L, phi_temp_R;
  real_t rho_temp_L, rho_temp_R;
  real_t rhod_temp_L, rhod_temp_R;
  for (int i = 0; i < Nx; i++)
  { 
    phi_temp_L = (i == 0) ? phi_L : phi_(i - 1);
    rho_temp_L = (i == 0) ? rho_L_ : rho_(i - 1);
    rhod_temp_L = (i == 0) ? rho_d(x1) : rhod_(i - 1);

    phi_temp_R = (i == Nx-1) ? phi_R : phi_(i + 1);
    rho_temp_R = (i == Nx-1) ? rho_R_ : rho_(i + 1);
    rhod_temp_R = (i == Nx-1) ? rho_d(x2) : rhod_(i + 1);
    (*bf)(i) = (*bf)(i) - ((phi_temp_R - phi_(i)) * (rho_temp_R + rho_(i) - rhod_(i) - rhod_temp_R)
        - (phi_(i) - phi_temp_L) * (rho_(i) + rho_temp_L - rhod_(i) - rhod_temp_L)) / 2.e0 * coef;
  };
  
  (*bf)(0) = (*bf)(0) + coef * rho_L_;
  (*bf)(0) = (*bf)(0) - (rhod_(0) + rho_d(x1)) / 2.e0 * phi_L * coef;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) - (rhod_(Nx-1) + rho_d(x2)) / 2.e0 * phi_R * coef;

};

void DriftDiffusion1D_FD_penaltyImplicit::generateB(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  B_.resize(Nx, Nx);  B_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_B;
  int estimatedNonZeros = 3 * Nx;
  tripletList_B.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  // real_t phi_L = phi_bc(x1, Trun + dt);
  // real_t phi_R = phi_bc(x2, Trun + dt);
  real_t rhod_temp_L, rhod_temp_R;
  for (int i = 0; i < Nx; i++) {
    rhod_temp_L = (i == 0) ? rho_d(x1) : rhod_(i - 1);
    rhod_temp_R = (i == Nx - 1) ? rho_d(x2) : rhod_(i + 1);
    val = - (rhod_(i) + rhod_temp_L) / 2.e0 * coef - (rhod_temp_R + rhod_(i)) / 2.e0 * coef;
    tripletList_B.emplace_back(i, i, val);
    if (i > 0) {
      val = (rhod_(i) + rhod_(i-1)) / 2.e0 * coef;
      tripletList_B.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = (rhod_(i+1) + rhod_(i)) / 2.e0 * coef;
      tripletList_B.emplace_back(i, i + 1, val);
    } 
  }

  B_.setFromTriplets(tripletList_B.begin(), tripletList_B.end());
};

void DriftDiffusion1D_FD_penaltyImplicit::update_E(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const SparseMatrix& Poi = poisol_->getS();
  // *** 
  Vector bf;
  real_t bc_L = gamma_ * phi_bc(x1, Trun + dt);
  real_t bc_R = gamma_ * phi_bc(x2, Trun + dt);
  poisol_->Assemble_bf(bc_L, bc_R, rhod_, &bf);

  Vector g;
  generatebf(Trun, dt, &g);

  // Schur complement
  g = g - S_ * bf;
  g = - g;

  SparseMatrix Schur = - B_ + gamma_ * S_ * Poi;
  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    std::cout << "  Preparation for BICG the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(Schur);
    bicg_.setTolerance(sparsetol_);
    std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  case Solver1DTypeFD::LU:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LU the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
#endif 
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(Schur);
    lu_.factorize(Schur);
    std::cout << "  The end of preparation for LU method ! " << std::endl;
    break;
  case Solver1DTypeFD::GMRES:
    gmres_.setTolerance(sparsetol_);
    gmres_.set_restart(30);
    gmres_.setMaxIterations(Nx);
    gmres_.compute(Schur);
    break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    phi_ = bicg_.solve(g);
    if (bicg_.info() != Eigen::Success) {
      std::cout << "[ERROR] BICG failed to converge!\n";
    }
    break;

  case Solver1DTypeFD::LU:
// #ifndef QUEST_USE_MKL
//     std::cout << "  Preparation for LU the matrix by Eigen ......\n";
// #else 
//     std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
// #endif 
    phi_ = lu_.solve(g);
    break;

  case Solver1DTypeFD::GMRES:
    phi_ = gmres_.solve(g);
    std::cout << "Iterations: " << gmres_.iterations() << std::endl;
    std::cout << "Estimated error: " << gmres_.error() << std::endl;
    if (gmres_.info() != Eigen::Success) {
      std::cout << "[ERROR] GMRES failed to converge!\n";
    }
    break;
  default:
    QUEST_ERROR(" The Dift Diffusion Solver is not been implemented ! ");
    break;
  }

  if (!phi_.allFinite() || !g.allFinite()) {
    std::cout << "[ERROR] phi_ contains NaN or Inf after SolveS!" << std::endl;
    std::string OutputDir = "result/DD_FD_unipolar_penaltyImplicit";
    OutputDir = OutputDir + "_gamma" + std::to_string(gamma_);
    std::ofstream Sout(OutputDir + "/S.txt");
    for (int k = 0; k < Schur.outerSize(); ++k) {
      for (SparseMatrix::InnerIterator it(Schur, k); it; ++it) {
        Sout << it.row() << " " 
              << it.col() << " "
              << it.value() << "\n";
      }
    }
    std::cout << " phi.allFinite() = " << phi_.allFinite() << std::endl;
    std::cout << " g.allFinite() = " << g.allFinite() << std::endl;
    Sout.close();
    std::abort();
  }
};

DriftDiffusion1D_FD_PNjunction_penaltyImplicit::DriftDiffusion1D_FD_PNjunction_penaltyImplicit(const FDmesh* mesh1D,
                        Poissonsolver1DFD_Base* poisol,
                        const int& x_order,
                        const int& t_order) : 
  DriftDiffusion1D_FD(mesh1D, poisol, x_order, t_order),
  DriftDiffusion1D_FD_penaltyImplicit(mesh1D, poisol, x_order, t_order),
  DriftDiffusion1D_FD_PNjunction(mesh1D, poisol, x_order, t_order) {};

DriftDiffusion1D_FD_period::DriftDiffusion1D_FD_period(const FDmesh_period* mesh1D,
                        Poissonsolver1DFD_Base_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : DriftDiffusion1D_FD(mesh1D, poisol, x_order, t_order) {};

real_t DriftDiffusion1D_FD_period::rho_init(const real_t& x) 
{
  return rho_d(x) + gamma_ * std::sin(x);
};

Vector DriftDiffusion1D_FD_period::rho_init(const Vector& x) 
{
  // return gamma_ * x.array().sin();
  Vector rho = rho_d(x).array() + gamma_ * x.array().sin();
  return rho;
};

real_t DriftDiffusion1D_FD_period::rho_d(const real_t& x) 
{
  // return 1.e0;
  return 1.e0;
};

Vector DriftDiffusion1D_FD_period::rho_d(const Vector& x) 
{
  Vector rhod = x;
  rhod.setConstant(1.e0);
  return rhod;
};

real_t DriftDiffusion1D_FD_period::source(const real_t& x, const real_t& t) 
{
  real_t S = - std::exp(- 2.e0 * t) * std::cos(2.e0 * x);
  S = S * gamma_;
  S = S + std::exp(- t) * std::sin(x);
  return S;
};

Vector DriftDiffusion1D_FD_period::source(const Vector& x, const real_t& t) 
{
  Vector S = 2.e0 * x;
  S = - std::exp(- 2.e0 * t) * S.array().cos() * gamma_;
  S = S.array() + std::exp(- t) * x.array().sin();
  return S;
};

real_t DriftDiffusion1D_FD_period::rho_real(const real_t& x, const real_t& t)
{
  return rho_d(x) + gamma_ * std::exp(- t) * std::sin(x);
};

Vector DriftDiffusion1D_FD_period::rho_real(const Vector& x, const real_t& t)
{
  Vector rho = rho_d(x).array() + gamma_ * std::exp(- t) * x.array().sin();
  return rho;
};

Vector DriftDiffusion1D_FD_period::E_real(const Vector& x, const real_t& t) 
{
  Vector E = - std::exp(- t) * x.array().cos();
  return E;
};

void DriftDiffusion1D_FD_period::setrhofinal(const real_t& t) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  // *** 
  rho_final_ = rho_real(X, t);
};

const Vector& DriftDiffusion1D_FD_period::getrhofinal() const {
  return rho_final_;
};

Vector DriftDiffusion1D_FD_period::getEfinal(const real_t& t) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  // *** 
  Vector E = E_real(X, t);
  return E;
};

Vector DriftDiffusion1D_FD_period::getphifinal(const real_t& t) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  // *** 
  Vector phi = - std::exp(- t) * X.array().sin();
  return phi;
};

void DriftDiffusion1D_FD_period::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + 2.e0 * coef - phi_(i) * coef;
    phi_temp_L = (i == 0) ? phi_(Nx-1) : phi_(i - 1);
    val = val + phi_temp_L / 2.e0 * coef;
    phi_temp_R = (i == Nx-1) ? phi_(0) : phi_(i + 1);
    val = val + phi_temp_R / 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - coef - (phi_(i) - phi_(i-1)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - coef + (phi_(i+1) - phi_(i)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }
  // i = 0
  tripletList_S.emplace_back(0, Nx-1, - coef - (phi_(0) - phi_(Nx-1)) / 2.e0 * coef);
  // i = Nx-1
  tripletList_S.emplace_back(Nx-1, 0, - coef + (phi_(0) - phi_(Nx-1)) / 2.e0 * coef);

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    std::cout << "  Preparation for BICG the matrix ......\n";
    std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(S_);
    bicg_.setTolerance(sparsetol_);
    std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  case Solver1DTypeFD::LU:
// #ifndef QUEST_USE_MKL
//     std::cout << "  Preparation for LU the matrix by Eigen ......\n";
// #else 
//     std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
// #endif 
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    lu_.analyzePattern(S_);
    lu_.factorize(S_);
    // std::cout << "  The end of preparation for LU method ! " << std::endl;
    break;
  case Solver1DTypeFD::GMRES:
    gmres_.setTolerance(sparsetol_);
    gmres_.set_restart(30);
    gmres_.setMaxIterations(Nx);
    gmres_.compute(S_);
    break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
};

void DriftDiffusion1D_FD_period::generatebf(const real_t& Trun, 
                                    const real_t& dt,
                                    Vector* bf) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  *bf = rhon_ + Source_;
};

void DriftDiffusion1D_FD_period::update_E(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  Vector bf = - (rho_ - rhod_) / gamma_;
  poisol_->SolveAll(0.e0, 0.e0, bf, &phi_, &E_);
};

} // namespace QUEST

#include "KineticDriftDiffusion1D_SL_Newton.hpp"
namespace QUEST
{

KineticDriftDiffusion1D_SL_Newton::KineticDriftDiffusion1D_SL_Newton(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order):
  KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order) {};

void KineticDriftDiffusion1D_SL_Newton::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const real_t& hx = mesh1D_->gethx();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  iter_ = 0.e0;
  rho_ = rho_init(X);
  rhod_ = rho_d(X);
  f_.resize(Nv);
  xback_.resize(Nv);
  vback_.resize(Nv);
  finter_.resize(Nv);
  Fx_.resize(Nv);
  Source_f_.resize(Nv);
  SLflux_ = Vector::Zero(Nx);
  Source_v_ = Vector::Zero(Nx);
  soltype_ = soltype;
  for (int j = 0; j < Nv; j++) {
    f_[j] = f_init(X, V(j));
    xback_[j] = Vector::Zero(Nx);
    vback_[j] = Vector::Zero(Nx);
    Fx_[j] = Vector::Zero(Nx);
    finter_[j] = Vector::Zero(Nx);
  };

  De_.resize(Nx+2, Nx);
  int estimatedNonZeros = 2 * Nx;
  std::vector<Eigen::Triplet<real_t>> tripletList_De;
  tripletList_De.reserve(estimatedNonZeros);
  real_t val = 0.5e0 / hx;
  for (int i = 0; i < Nx; i++) {
    tripletList_De.emplace_back(i, i, val);
    tripletList_De.emplace_back(i+2, i, - val);
  }
  tripletList_De.emplace_back(0, 0, val);
  tripletList_De.emplace_back(Nx+1, Nx-1, - val);
  De_.setFromTriplets(tripletList_De.begin(), tripletList_De.end());
}

void KineticDriftDiffusion1D_SL_Newton::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  const SparseMatrix& Sp = poisol_->getS();
  // *** 

  real_t omega = std::exp( - dt / (eps_ * eps_));
  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (1.e0 - omega) * D_ * 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, - (1.e0 - omega) * D_  * coef);
      val = - (1.e0 - omega) * D_ / theta_ * E_(i) * dt / (2.e0 * hx);
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, - (1.e0 - omega) * D_  * coef);
      val = (1.e0 - omega) * D_ / theta_ * E_(i + 2) * dt / (2.e0 * hx);
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  A_.resize(Nx, Nx+2); A_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_A;
  estimatedNonZeros = 2 * Nx;
  tripletList_A.reserve(estimatedNonZeros);
  for (int i = 0; i < Nx; i++) {
    if (i > 0) {
      val = - (1.e0 - omega) * D_ / theta_ * rho_(i - 1) * dt / (2.e0 * hx);
      tripletList_A.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = (1.e0 - omega) * D_ / theta_ * rho_(i + 1) * dt / (2.e0 * hx);
      tripletList_A.emplace_back(i, i + 1, val);
    } 
  }
  A_.setFromTriplets(tripletList_A.begin(), tripletList_A.end());
  S_ = gamma_ * S_ * Sp - A_ * De_;

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

void KineticDriftDiffusion1D_SL_Newton::update_E(const real_t& Trun) {
  // ** 传入相关变量
  const real_t& hx = mesh1D_->gethx();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  real_t bc_L = phi_bc(x1, Trun);
  real_t bc_R = phi_bc(x2, Trun);
  real_t chx = 2.e0 * hx;
  for (int i = 2; i <= Nx - 1; i++) {
    E_(i) = (phi_(i) - phi_(i - 2)) / chx;
  }
  E_(1) = (phi_(1) - bc_L) / chx;
  E_(Nx) = (bc_R - phi_(Nx - 2)) / chx;

  if (x_order_ == 1) {
    E_(0) = (phi_(0) - bc_L) / hx;
    E_(Nx + 1) = (bc_R - phi_(Nx - 1)) / hx;
  } else {
    QUEST_ERROR(" The 2nd-order or higher-order approximation for u_x has not been implemented ! ");
  }
};

void KineticDriftDiffusion1D_SL_Newton::generatebf(const real_t& Trun, 
                          const real_t& dt,
                          Vector* bf) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  Source_v_.resize(Nx); Source_v_.setZero();
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);

  *bf = rhon_ - omega * dt / eps_ * SLflux_ + Source_ + Source_v_;
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ * coef * rho_L_;
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ / theta_ * E_(0) * dt / (2.e0 * hx) * rho_L_;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + (1.e0 - omega) * D_ * coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) - (1.e0 - omega) * D_ / theta_ * E_(Nx+1) * dt / (2.e0 * hx) * rho_R_;
  // for (int i = 1; i < Nx-1; i++) {
  //   (*bf)(i) = (*bf)(i) + (1.e0 - omega) * D_ * coef * (rho_(i-1) - 2.e0 * rho_(i) + rho_(i+1))
  //                       - (1.e0 - omega) * D_ / theta_ * dt / (2.e0 * hx) * 
  //                       (E_(i+2) * rho_(i+1) - E_(i) * rho_(i-1));
  // }
  *bf = *bf - S_ * rho_;
}


void KineticDriftDiffusion1D_SL_Newton::SolveS(const real_t& Trun, const Vector& bf) {
  // ** 传入相关变量
  const SparseMatrix& Sp = poisol_->getS();
  // *** 

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    deltaphi_ = bicg_.solve(bf);
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
    deltaphi_ = lu_.solve(bf);
    break;

  case Solver1DTypeFD::GMRES:
    deltaphi_ = gmres_.solve(bf);
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
  deltaphi_ = - deltaphi_;
  deltarho_ = - gamma_ * Sp * deltaphi_;
  deltaE_ = De_ * deltaphi_;

  phi_ = phi_ + deltaphi_;
  rho_ = rho_ + deltarho_;
  E_ = E_ + deltaE_;
}

void KineticDriftDiffusion1D_SL_Newton::update_rho(const real_t& Trun, 
                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  const real_t& hv = mesh1D_->gethv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  // *** 
  rho_L_ = rho_numericalbc(x1, Trun);
  rho_R_ = rho_numericalbc(x2, Trun);
  SLflux_.setZero();
  real_t area = hx * hv;
#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t vtemp, xtemp, ftemp;
  int is, js, ii;
  real_t xita, vita;
  real_t M;
  real_t dvdx;
  Vector fx(4), fv(4), rhotemp(4);
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    Fx_[j].setZero();
    M = Maxwell(vtemp);
    for (int i = 0; i < Nx; i++) {
      vtemp = vback_[j](i);
      xtemp = xback_[j](i);
      // xtemp\in [x_{is-1},x_{is}]; vtemp\in [v_{js-1},v_{js}];
      // 0 -> is-2, 1->is-1, 2->is, 3->is+1
      is = floor((xtemp - X(0)) / hx) + 1;
      js = floor((vtemp - V(0)) / hv) + 1;
      xita = (is * hx + X(0) - xtemp) / hx;
      vita = (js * hv + V(0) - vtemp) / hv;
      interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
      if (x_order_ == 1) {
        dvdx = - dt * (rho_(i) - rhod_(i)) / (eps_ * gamma_);
        if (V(j) >= 0) {
          Fx_[j](i) = (fx(1) - fx(0)) / hx + dvdx * (fv(1) - fv(0)) / hv
                      - (rhotemp(3) - rhotemp(2)) / hx * M
                      + dvdx * vtemp / theta_ * M;
          Fx_[j](i) = Fx_[j](i) / eps_;
        } else if (V(j) < 0) {
          Fx_[j](i) = (fx(3) - fx(2)) / hx + dvdx * (fv(3) - fv(2)) / hv
                      - (rhotemp(1) - rhotemp(0)) / hx * M
                      + dvdx * vtemp / theta_ * M;
          Fx_[j](i) = Fx_[j](i) / eps_;
        }
        if (js - 1 < 0 || js - 1 >= Nv) {
          ftemp= 0.e0;
        } else {
          ftemp = (is - 1 < 0) ? fL_bc(js-1, Trun) : 
                    (is - 1 >= Nx) ? fR_bc(js-1, Trun) : f_[js-1](is - 1);
        }
        finter_[j](i) = xita * vita * ftemp + (1.e0 - xita) * vita * fv(1) 
                        + xita * (1.e0 - vita) * fx(1) + (1.e0 - xita) * (1.e0 - vita) * fx(2);
      } else if (x_order_ == 2) {
        QUEST_ERROR(" The second SL method has not been implemented ! ");
      }
    }
    #pragma omp critical
    SLflux_ += Fx_[j] * Vweights(j) * V(j);
  }
}

  // 开始迭代
  real_t err = 1.e0;
  Vector bf;
  rhon_ = rho_;
  iter_ = 0;
  Vector prerhotemp = rho_;
  update_E(Trun);
  while (err > iterationtol_) {
    generateS(Trun, dt);
    generatebf(Trun, dt, &bf);
    Vector rhotemp = rho_;
    SolveS(Trun, bf);  // update rho_, phi_, E_
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
      
      Sout.close();
      std::abort();
    }
    std::cout << " rho_.allFinite() = " << rho_.allFinite() << std::endl;
    std::cout << " bf.allFinite() = " << bf.allFinite() << std::endl;
    iter_++;
    rhotemp = rhotemp - rho_;
    err = rhotemp.cwiseAbs().maxCoeff();
    std::cout <<" size = " << rho_.size() << ", err = " << std::setprecision(10) << err << ", iter = " << iter_ << std::endl;
    real_t cosang = rhotemp.dot(prerhotemp) / (rhotemp.norm() * rhotemp.norm() + 1e-300);
    std::cout << "cos(angle(delta_k, delta_{k-1})) = " << cosang << "\n"; 
    std::cout << " rhotemp.norm() / prerhotemp.norm() = " << rhotemp.norm() / prerhotemp.norm() << std::endl;
    prerhotemp = rhotemp;
    // if (iter_ > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort." << std::endl;
    //   break;
    // }
  }
  std::cout << "err = " << err << std::endl;
};

} // namespace QUEST

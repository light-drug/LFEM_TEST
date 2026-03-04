#include "APDriftDiffusion_SL.hpp"

namespace QUEST
{
  
APDriftDiffusion_SL::APDriftDiffusion_SL(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order) {};

void APDriftDiffusion_SL::setquasipoi_tol(const real_t& quasipoi_tol) 
{
  quasipoi_tol_ = quasipoi_tol;
}

void APDriftDiffusion_SL::setsoltype_poi(const Solver1DTypeFD& soltype_poi) 
{
  soltype_poi_ = soltype_poi;
}

void APDriftDiffusion_SL::update_rho(const real_t& Trun, 
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
  switch (x_order_)
  {
  case 1:
  {
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
      for (int i = 0; i < Nx; i++) {
        vtemp = vback_[j](i);
        xtemp = xback_[j](i);
        M = Maxwell(vtemp);
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
  }
    break;
  
  default:
    break;
  }

  Source_v_.resize(Nx); Source_v_.setZero();
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);

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
    SolveS(Trun, bf);  // update rho_
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
    prerhotemp = rhotemp;
    // if (iter_ > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort." << std::endl;
    //   break;
    // }
  }
  std::cout << "err = " << err << std::endl;
};


// 为了满足电中性极限的AP格式需要对泊松方程进行改写
// quasipoi_为修正后的泊松矩阵
void APDriftDiffusion_SL::update_E(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  const SparseMatrix& poimatrix = poisol_->getS();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  real_t omega = std::exp( - dt / (eps_ * eps_));
  quasipoi_.resize(Nx, Nx);  quasipoi_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_quasipoi_;
  int estimatedNonZeros = 5 * Nx;
  tripletList_quasipoi_.reserve(estimatedNonZeros);

  const real_t coef = 1.e0 / (hx * hx);
  real_t val, rhotemp;
  for (int i = 0; i < Nx; i++) {
    val = (2.e0 * gamma_ + (1 - omega) * dt * D_ / theta_ * rho_(i)) * coef;
    tripletList_quasipoi_.emplace_back(i, i, val);

    rhotemp = (i == 0) ? rho_L_ : rho_(i - 1);
    val = (1 - omega) * dt * D_ / theta_ * rhotemp / 2.e0 * coef;
    tripletList_quasipoi_.emplace_back(i, i, val);

    rhotemp = (i == Nx-1) ? rho_R_ : rho_(i + 1);
    val = (1 - omega) * dt * D_ / theta_ * rhotemp / 2.e0 * coef;
    tripletList_quasipoi_.emplace_back(i, i, val);

    if (i > 0) {
      val = (- gamma_ - (1 - omega) * dt * D_ / theta_ * (rho_(i - 1) + rho_(i)) / 2.e0) * coef;
      tripletList_quasipoi_.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = (- gamma_ - (1 - omega) * dt * D_ / theta_ * (rho_(i + 1) + rho_(i)) / 2.e0) * coef;
      tripletList_quasipoi_.emplace_back(i, i + 1, val);
    } 
  }

  quasipoi_.setFromTriplets(tripletList_quasipoi_.begin(), tripletList_quasipoi_.end());

  switch (soltype_poi_) {
    case Solver1DTypeFD::LU:
#ifndef QUEST_USE_MKL
      std::cout << "  Preparation for LU the matrix by Eigen ......\n";
#else 
      std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
#endif 
      lupoi_.analyzePattern(quasipoi_);
      lupoi_.factorize(quasipoi_);
      std::cout << "  The end of preparation for LU method ! " << std::endl;
      if (lupoi_.info() != Eigen::Success) {
        std::cout << "[ERROR] LU decomposition failed!\n";
        QUEST_ABORT(" The LU decomposition failed ! ");
      }
      break;
    case Solver1DTypeFD::GMRES:
      std::cout << "quasipoi_tol_ = " << quasipoi_tol_ << std::endl; 
      gmrespoi_.setTolerance(quasipoi_tol_);
      gmrespoi_.set_restart(Nx / 5);
      gmrespoi_.setMaxIterations(Nx);
      gmrespoi_.compute(quasipoi_);
      break;
    default:
      QUEST_ERROR(" The solver type for quasi-neutral Poisson equation has not been implemented ! ");
  }
  

  Vector bf = - (rho_ - rhod_ - omega * dt * SLflux_ + Source_ + Source_v_);
  bf = bf + (1 - omega) * dt * D_ * poimatrix * rho_;
  real_t bc_L = phi_bc(x1, Trun + dt);
  real_t bc_R = phi_bc(x2, Trun + dt);
  bf(0) = bf(0) - (1 - omega) * dt * D_ * rho_L_ * coef;
  bf(Nx - 1) = bf(Nx - 1) - (1 - omega) * dt * D_ * rho_R_ * coef;

  bf(0) = bf(0) + (gamma_ + (1 - omega) * dt * D_ / theta_ * (rho_L_ + rho_(0)) / 2.e0) * bc_L * coef;
  bf(Nx - 1) = bf(Nx - 1) + (gamma_ + (1 - omega) * dt * D_ / theta_ * (rho_R_ + rho_(Nx - 1)) / 2.e0) * bc_R * coef;

  // update phi_
  std::cout << " hsdhadsfuabl " << std::endl;
  switch (soltype_poi_)
  {
  case Solver1DTypeFD::LU:
    phi_ = lupoi_.solve(bf);
    break;
  case Solver1DTypeFD::GMRES:
    phi_ = gmrespoi_.solve(bf);
    std::cout << "Iterations: " << gmrespoi_.iterations() << std::endl;
    std::cout << "Estimated error: " << gmrespoi_.error() << std::endl;
    if (gmrespoi_.info() != Eigen::Success) {
      std::cout << "[ERROR] GMRES failed to converge!\n";
      QUEST_ABORT(" The GMRES failed to converge ! ");
    }
  default:
    break;
  }


  // update E_
  E_.resize(Nx + 2);
  E_.setZero();
  real_t hx2 = 2.e0 * hx;
  for (int i = 2; i <= Nx - 1; i++) {  
    E_(i) = (phi_(i) - phi_(i - 2)) / hx2;
  }
  E_(1) = (phi_(1) - bc_L) / hx2;
  E_(Nx) = (bc_R - phi_(Nx - 2)) / hx2;

  if (x_order_ == 1) {
    E_(0) = (phi_(0) - bc_L) / hx;
    E_(Nx + 1) = (bc_R - phi_(Nx - 1)) / hx;
  } else {
    QUEST_ERROR(" The 2nd-order or higher-order approximation for u_x has not been implemented ! ");
  }
}

void APDriftDiffusion_SL::generateS(const real_t& Trun, const real_t& dt) {
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

  real_t omega = std::exp( - dt / (eps_ * eps_));
  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (1.e0 - omega) * D_ * 2.e0 * coef 
          - (1.e0 - omega) * D_ / theta_ * phi_(i) * coef;
    phi_temp_L = (i == 0) ? phi_L : phi_(i - 1);
    val = val + (1.e0 - omega) * D_ / theta_ * phi_temp_L / 2.e0 * coef;
    phi_temp_R = (i == Nx-1) ? phi_R : phi_(i + 1);
    val = val + (1.e0 - omega) * D_ / theta_ * phi_temp_R / 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - (1.e0 - omega) * D_ * coef 
            - (1.e0 - omega) * D_ / theta_ * (phi_(i) - phi_(i-1)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - (1.e0 - omega) * D_ * coef
            + (1.e0 - omega) * D_ / theta_ * (phi_(i+1) - phi_(i)) / 2.e0 * coef;
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
    std::cout << "  The end of preparation for LU method ! " << std::endl;
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

// drift-diffusion方程的右端项
void APDriftDiffusion_SL::generatebf(const real_t& Trun, 
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
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  real_t phi_L = phi_bc(x1, Trun + dt);
  real_t phi_R = phi_bc(x2, Trun + dt);

  *bf = rhon_ - omega * dt * SLflux_ + Source_ + Source_v_;
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ * coef * rho_L_;
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ / theta_ * (phi_(0) - phi_L) / 2.e0 * rho_L_ * coef;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + (1.e0 - omega) * D_ * coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) - (1.e0 - omega) * D_ / theta_ 
                    * (phi_R - phi_(Nx-1)) / 2.e0 * rho_R_ * coef;
};


// 一维 PN 结问题
APDriftDiffusion1D_SL_PNjunction::APDriftDiffusion1D_SL_PNjunction(const KineticFDmesh* mesh1D,
                      Poissonsolver1DFD* poisol,
                      const int& x_order,
                      const int& t_order) : 
  APDriftDiffusion_SL(mesh1D, poisol, x_order, t_order) { };

real_t APDriftDiffusion1D_SL_PNjunction::rho_init(const real_t& x) 
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector APDriftDiffusion1D_SL_PNjunction::rho_init(const Vector& x)
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

real_t APDriftDiffusion1D_SL_PNjunction::rho_d(const real_t& x)
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector APDriftDiffusion1D_SL_PNjunction::rho_d(const Vector& x)
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

real_t APDriftDiffusion1D_SL_PNjunction::source(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector APDriftDiffusion1D_SL_PNjunction::source(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t APDriftDiffusion1D_SL_PNjunction::fsource(const real_t& x, const real_t& v, const real_t& t)
{
  return 0.e0;
};

Vector APDriftDiffusion1D_SL_PNjunction::fsource(const Vector& x, const real_t& v, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t APDriftDiffusion1D_SL_PNjunction::fsource_dx(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector APDriftDiffusion1D_SL_PNjunction::fsource_dx(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t APDriftDiffusion1D_SL_PNjunction::f_init(const real_t& x, const real_t& v)
{
  return rho_init(x) * Maxwell(v) + eps_ * g_init(x, v);
};

Vector APDriftDiffusion1D_SL_PNjunction::f_init(const Vector& x, const real_t& v)
{
  return rho_init(x) * Maxwell(v) + eps_ * g_init(x, v);
};

real_t APDriftDiffusion1D_SL_PNjunction::g_init(const real_t& x, const real_t& v)
{
  return 0.e0;
};

Vector APDriftDiffusion1D_SL_PNjunction::g_init(const Vector& x, const real_t& v)
{
  return Vector::Zero(x.size());
};

real_t APDriftDiffusion1D_SL_PNjunction::fL_bc(const real_t& v, const real_t& t) 
{
  QUEST_ERROR(" This fL_bc(real_t, real_t) is not valid ! ");
  return 0.e0;
};

real_t APDriftDiffusion1D_SL_PNjunction::fL_bc(const int& v_index, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(v_index) >= 0) {
    f = 0.9e0 * Maxwell(V(v_index));   // 入流边界条件
  } else {
    f = 2.e0 * f_[v_index](0) - f_[v_index](1);
  }
  return f;
};

real_t APDriftDiffusion1D_SL_PNjunction::fR_bc(const real_t& v, const real_t& t)
{
  QUEST_ERROR(" This fR_bc(real_t, real_t) is not valid ! ");
  return 0.e0;
};

real_t APDriftDiffusion1D_SL_PNjunction::fR_bc(const int& v_index, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  // *** 
  real_t f;
  if (V(v_index) <= 0) {
    f = 0.1e0 * Maxwell(V(v_index));
  } else {
    f = 2.e0 * f_[v_index](Nx-1) - f_[v_index](Nx-2);
  }
  return f;
};

real_t APDriftDiffusion1D_SL_PNjunction::phi_bc(const real_t& x, const real_t & t)
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




// 周期边界条件，多用于精度测试
APDriftDiffusion1D_SL_period::APDriftDiffusion1D_SL_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order) {};

void APDriftDiffusion1D_SL_period::setquasipoi_tol(const real_t& quasipoi_tol) 
{ 
  quasipoi_tol_ = quasipoi_tol;
}

void APDriftDiffusion1D_SL_period::update_rho(const real_t& Trun, const real_t& dt) 
{
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
    for (int i = 0; i < Nx; i++) {
      vtemp = vback_[j](i);
      xtemp = xback_[j](i);
      M = Maxwell(vtemp);
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
        ii = (((is - 1) % Nx) + Nx) % Nx;
        ftemp = (js - 1 >= 0 && js - 1 < Nv) ? f_[js - 1](ii) : 0.0;
        finter_[j](i) = xita * vita * ftemp + (1.e0 - xita) * vita * fv(1) 
                        + xita * (1.e0 - vita) * fx(1) + (1.e0 - xita) * (1.e0 - vita) * fx(2) ;
      } else if (x_order_ == 2) {
        QUEST_ERROR(" The second SL method has not been implemented ! ");
      }
    }
    #pragma omp critical
    SLflux_ += Fx_[j] * Vweights(j) * V(j);
  }
}
  Source_v_.resize(Nx); Source_v_.setZero();
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);

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
    SolveS(Trun, bf);  // update rho_
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

void APDriftDiffusion1D_SL_period::generateS(const real_t& Trun, const real_t& dt) 
{
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

  real_t omega = std::exp( - dt / (eps_ * eps_));
  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (1.e0 - omega) * D_ * 2.e0 * coef 
          - (1.e0 - omega) * D_ / theta_ * phi_(i) * coef;
    phi_temp_L = (i == 0) ? phi_(Nx - 1) : phi_(i - 1);
    val = val + (1.e0 - omega) * D_ / theta_ * phi_temp_L / 2.e0 * coef;
    phi_temp_R = (i == Nx - 1) ? phi_(0) : phi_(i + 1);
    val = val + (1.e0 - omega) * D_ / theta_ * phi_temp_R / 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - (1.e0 - omega) * D_ * coef 
            - (1.e0 - omega) * D_ / theta_ * (phi_(i) - phi_(i-1)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - (1.e0 - omega) * D_ * coef
            + (1.e0 - omega) * D_ / theta_ * (phi_(i+1) - phi_(i)) / 2.e0 * coef;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }

  // i = 0
  val = - (1.e0 - omega) * D_ * coef - (1.e0 - omega) * D_ / theta_ * (phi_(0) - phi_(Nx-1)) / 2.e0 * coef;
  tripletList_S.emplace_back(0, Nx-1, val);
  // i = Nx-1
  val = - (1.e0 - omega) * D_ * coef + (1.e0 - omega) * D_ / theta_ * (phi_(0) - phi_(Nx-1)) / 2.e0 * coef;
  tripletList_S.emplace_back(Nx-1, 0, val);

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    // std::cout << "  Preparation for BICG the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    bicg_.compute(S_);
    bicg_.setTolerance(sparsetol_);
    // std::cout << "  The end of preparation for BICG method ! " << std::endl;
    break;
  case Solver1DTypeFD::LU:
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

// rho更新矩阵右端项发生改变
void APDriftDiffusion1D_SL_period::generatebf(const real_t& Trun, 
                        const real_t& dt,
                        Vector* bf) 
{
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
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);

  *bf = rhon_ - omega * dt * SLflux_ + Source_ + Source_v_;
  // 周期边界下不需要处理边界的那些自由度
};

// 电场更新方程变为重写后的泊松方程
void APDriftDiffusion1D_SL_period::update_E(const real_t& Trun, const real_t& dt) 
{
   // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  const SparseMatrix& poimatrix = poisol_->getS();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 

  real_t omega = std::exp( - dt / (eps_ * eps_));
  quasipoi_.resize(Nx, Nx);  quasipoi_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_quasipoi_;
  int estimatedNonZeros = 6 * Nx;
  tripletList_quasipoi_.reserve(estimatedNonZeros);

  const real_t coef = 1.e0 / (hx * hx);
  real_t val, rhotemp;
  for (int i = 0; i < Nx; i++) {
    val = (2.e0 * gamma_ + (1 - omega) * dt * D_ / theta_ * rho_(i)) * coef;
    tripletList_quasipoi_.emplace_back(i, i, val);
    
    rhotemp = (i == 0) ? rho_(Nx-1) : rho_(i - 1);
    val = (1 - omega) * dt * D_ / theta_ * rhotemp / 2.e0 * coef;
    tripletList_quasipoi_.emplace_back(i, i, val);

    rhotemp = (i == Nx-1) ? rho_(0) : rho_(i + 1);
    val = (1 - omega) * dt * D_ / theta_ * rhotemp / 2.e0 * coef;
    tripletList_quasipoi_.emplace_back(i, i, val);

    if (i > 0) {
      val = (- gamma_ - (1 - omega) * dt * D_ / theta_ * (rho_(i - 1) + rho_(i)) / 2.e0) * coef;
      tripletList_quasipoi_.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = (- gamma_ - (1 - omega) * dt * D_ / theta_ * (rho_(i + 1) + rho_(i)) / 2.e0) * coef;
      tripletList_quasipoi_.emplace_back(i, i + 1, val);
    } 
  };
  val = (- gamma_ - (1 - omega) * dt * D_ / theta_ * (rho_(Nx - 1) + rho_(0)) / 2.e0) * coef;
  tripletList_quasipoi_.emplace_back(0, Nx-1, val);
  tripletList_quasipoi_.emplace_back(Nx-1, 0, val);

  int center = Nx / 2;
  val = 1.e0 / hx;
  for (int i = 0; i < Nx; i++) {
    tripletList_quasipoi_.emplace_back(center, i, val);
  }

  quasipoi_.setFromTriplets(tripletList_quasipoi_.begin(), tripletList_quasipoi_.end());

  std::cout << " APDriftDiffusion1D_SL_period quasipoi_tol_ = " << quasipoi_tol_ << std::endl;
  gmrespoi_.setTolerance(quasipoi_tol_);
  gmrespoi_.set_restart(Nx);
  gmrespoi_.setMaxIterations(5 * Nx);
  gmrespoi_.compute(quasipoi_);

  Vector bf = - (rho_ - rhod_ - omega * dt * SLflux_ + Source_ + Source_v_);
  real_t rho_Ltemp, rho_Rtemp;
  for (int i = 0; i < Nx; i++) {
    rho_Ltemp = (i == 0) ? rho_(Nx - 1) : rho_(i - 1);
    rho_Rtemp = (i == Nx - 1) ? rho_(0) : rho_(i + 1);
    bf(i) = bf(i) + (1 - omega) * dt * D_ * (- rho_Ltemp + 2.e0 * rho_(i) - rho_Rtemp) * coef;
  }
  

  // update phi_
  // std::cout << "quasipoi_ = " << quasipoi_tol_ << std::endl;
  phi_ = gmrespoi_.solve(bf);
  // std::cout << "Iterations: " << gmrespoi_.iterations() << std::endl;
  // std::cout << "Estimated error: " << gmrespoi_.error() << std::endl;
  if (gmrespoi_.info() != Eigen::Success) 
  {
    std::cout << "[ERROR] GMRES failed to converge!\n";
    QUEST_ABORT(" The GMRES failed to converge ! ");
  };

  // update E_
  E_.resize(Nx);
  E_.setZero();
  real_t hx2 = 2.e0 * hx;
  for (int i = 1; i < Nx - 1; i++) 
  {
    E_(i) = (phi_(i + 1) - phi_(i - 1)) / hx2;
  };
  E_(0) = (phi_(1) - phi_(Nx - 1)) / hx2;
  E_(Nx-1) = (phi_(0) - phi_(Nx - 2)) / hx2;

};

real_t APDriftDiffusion1D_SL_period::rho_init(const real_t& x) 
{
  return rho_d(x) + gamma_ * std::sin(x);
};

Vector APDriftDiffusion1D_SL_period::rho_init(const Vector& x) 
{
  Vector rho = rho_d(x).array() + gamma_ * x.array().sin();
  return rho;
};

real_t APDriftDiffusion1D_SL_period::rho_d(const real_t& x) 
{
  return 1.e0;
};

Vector APDriftDiffusion1D_SL_period::rho_d(const Vector& x) 
{
  Vector rhod = x;
  rhod.setConstant(1.e0);
  return rhod;
};

real_t APDriftDiffusion1D_SL_period::source(const real_t& x, const real_t& t) 
{
  real_t S = - std::exp(- 2.e0 * D_  * t) * std::cos(2.e0 * x);
  S = S * gamma_;
  S = S + std::exp(- D_ * t) * std::sin(x);
  return S;
};

Vector APDriftDiffusion1D_SL_period::source(const Vector& x, const real_t& t) 
{
  Vector S = 2.e0 * x;
  S = - std::exp(- 2.e0 * D_ * t) * S.array().cos() * gamma_;
  S = S.array() + std::exp(- D_ * t) * x.array().sin();
  return S;
};

real_t APDriftDiffusion1D_SL_period::fsource(const real_t& x, const real_t& v, const real_t& t) 
{
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  real_t S;
  real_t E = - std::exp(- D_ * t) * std::cos(x);
  S = std::exp(- D_ * t) * M * (eps_ * v * D_ * std::cos(x) - 
                                D_ * sin(x) + v2 * sin(x)) +
      std::exp(- 2 * D_ * t) * M * (eps_ * v * D_ * std::sin(2*x)/theta_ 
                                - v2 / theta_ * std::cos(2*x) + 
                                std::cos(x) * std::cos(x) * (1.e0 - v2 / theta_)) +
      std::exp(- 3 * D_ * t) * M * (std::cos(x) * std::cos(x)) * std::sin(x) * 
                            (1.e0 / theta_  - v2 / (theta_ * theta_));
  S = S * gamma_ + v2 * M / theta_ * std::exp(- D_ * t) * std::sin(x);
  S = S + std::pow(E, 2) * (M / theta_ - v2 * M / std::pow(theta_, 2)) + eps_ * v * M / theta_ * (- D_ * E);
  return S;
};

Vector APDriftDiffusion1D_SL_period::fsource(const Vector& x, const real_t& v, const real_t& t) 
{
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  Vector S;
  Vector sinx = x.array().sin();
  Vector cosx = x.array().cos();
  Vector sin2x = (2.e0 * x).array().sin();
  Vector cos2x = (2.e0 * x).array().cos();
  Vector E = - std::exp(- D_ * t) * cosx;
  S = std::exp(- D_ * t) * M * (eps_ * v * D_ * cosx.array()
                                - D_ * sinx.array() + v2 * sinx.array()) +
      std::exp(- 2 * D_ * t) * M * (eps_ * v * D_ * sin2x.array()/theta_ 
                                - 1.e0 / theta_ * v2 * cos2x.array() + 
                                cosx.array() * cosx.array() * (1.e0 - v2 / theta_)) +
      std::exp(- 3 * D_ * t) * M * (cosx.array() * cosx.array()) * sinx.array() * 
                            (1.e0 / theta_  - v2 / (theta_ * theta_));
  S = S * gamma_ + v2 * M / theta_ * std::exp(- D_ * t) * sinx;
  S = S.array() + E.array().square() * (M / theta_ - v2 * M / std::pow(theta_, 2)) + eps_ * v * M / theta_ * (- D_ * E.array());
  return S;
};

real_t APDriftDiffusion1D_SL_period::fsource_dx(const real_t& x, const real_t& t) 
{
  // <vS> = eps_ * D_ * (D_ * std::cos(x) * std::exp(- D_ * t)
  //     + (D_) / theta_ * std::exp(- 2.e0 * D_ * t) * std::sin(2*x))
  // <vS> = <vS> * gamma_ + eps_ * D_ * std::exp(- D_ * t) * std::cos(x)
  real_t S = eps_ * D_ * ( - D_ * std::sin(x) * std::exp(- D_ * t)
      + (2.e0 * D_) / theta_ * std::exp(- 2.e0 * D_ * t) * std::cos(2.e0 * x));
  S = S * gamma_ - eps_ * D_ * std::exp(- D_ * t) * std::sin(x);
  return S;
};

Vector APDriftDiffusion1D_SL_period::fsource_dx(const Vector& x, const real_t& t) 
{
  Vector S = eps_ * D_ * ( - D_ * x.array().sin() * std::exp(- D_ * t)
      + (2.e0 * D_) / theta_ * std::exp(- 2.e0 * D_ * t) * (2.e0 * x).array().cos());
  S = S.array() * gamma_ - eps_ * D_ * std::exp(- D_ * t) * x.array().sin();
  return S;
};

real_t APDriftDiffusion1D_SL_period::f_init(const real_t& x, const real_t& v) 
{
  real_t f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Vector APDriftDiffusion1D_SL_period::f_init(const Vector& x, const real_t& v) 
{
  Vector f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t APDriftDiffusion1D_SL_period::g_init(const real_t& x, const real_t& v) 
{
  real_t M = Maxwell(v);
  real_t g = - v * M * (std::cos(x) + 1.e0 / theta_ * std::cos(x) * std::sin(x));
  g = g * gamma_ - v * M * std::cos(x) / theta_;
  return g;
};

Vector APDriftDiffusion1D_SL_period::g_init(const Vector& x, const real_t& v) 
{
  real_t M = Maxwell(v);
  Vector g = - v * M * (x.array().cos() + 1.e0 / theta_ * x.array().cos() * x.array().sin());
  g = g.array() * gamma_ - v * M * x.array().cos() / theta_;
  return g;
};

real_t APDriftDiffusion1D_SL_period::rho_real(const real_t& x, const real_t& t)
{
  return rho_d(x) + gamma_ * std::exp(- D_ * t) * std::sin(x);
};

Vector APDriftDiffusion1D_SL_period::rho_real(const Vector& x, const real_t& t)
{
  Vector rho = rho_d(x).array() + gamma_ * std::exp(- D_ * t) * x.array().sin();
  return rho;
};

Vector APDriftDiffusion1D_SL_period::E_real(const Vector& x, const real_t& t) 
{
  Vector E = - std::exp(- D_ * t) * x.array().cos();
  return E;
};

real_t APDriftDiffusion1D_SL_period::g_real(const real_t& x, const real_t& v, const real_t& t) 
{
  real_t M = Maxwell(v);
  real_t g =  - v * M * (std::exp(- D_ * t) * std::cos(x) + 
            std::exp(- 2.e0 * D_ * t) / theta_ * std::cos(x) * std::sin(x));
  g = g * gamma_ - std::exp(- D_ * t) * v * M * std::cos(x) / theta_;
  return g;
};

Vector APDriftDiffusion1D_SL_period::g_real(const Vector& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  Vector g =  - v * M * (std::exp(- D_ * t) * x.array().cos() + 
            std::exp(- 2.e0 * D_ * t) / theta_ * x.array().cos() * x.array().sin());
  g = g.array() * gamma_ - std::exp(- D_ * t) * v * M * x.array().cos() / theta_;
  return g;
};

void APDriftDiffusion1D_SL_period::update_E_final(const real_t& Tstop)
{
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

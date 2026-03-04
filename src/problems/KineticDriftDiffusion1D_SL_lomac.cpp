#include "KineticDriftDiffusion1D_SL_lomac.hpp"

namespace QUEST
{

KineticDriftDiffusion1D_SL_lomac::KineticDriftDiffusion1D_SL_lomac(
                        const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    time_stepping_(TimeSteppingType::Rktype) {};

void KineticDriftDiffusion1D_SL_lomac::setTimeSteppingType(const TimeSteppingType& time_stepping)
{
  time_stepping_ = time_stepping;
};

void KineticDriftDiffusion1D_SL_lomac::setsigmas(const real_t& sigmas)
{
  sigmas_ = sigmas;
};

const Vector& KineticDriftDiffusion1D_SL_lomac::getphi()
{
  return phi_;
}

void KineticDriftDiffusion1D_SL_lomac::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const int& xDiv = mesh1D_->getxDiv();
  // *** 
  QUEST_VERIFY(x_order_ == t_order_, "x_order must be equal to t_order (space order == time order) !");
  QUEST_VERIFY(xDiv != Nx, " The mesh is not the periodical mesh ! ");

  pi_ = 3.14159265358979323846264338327;
  iter_ = 0;
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
    Fx_[j] = Vector::Zero(xDiv);
    finter_[j] = Vector::Zero(Nx);
  };
  std::cout << " eps_ = " << eps_ << std::endl;
  std::cout << " gamma_ = " << gamma_ << std::endl;
  std::cout << " D_ = " << D_ << std::endl;
  std::cout << " theta_ = " << theta_ << std::endl;
  std::cout << " cfl_ = " << cfl_ << std::endl;
};

void KineticDriftDiffusion1D_SL_lomac::update_E(const real_t& Trun) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  Vector bf = - (rho_ - rhod_) / gamma_;
  real_t bc_L = phi_bc(x1, Trun);
  real_t bc_R = phi_bc(x2, Trun);
  poisol_->SolveAll(bc_L, bc_R, bf, &phi_, &E_);
};

void KineticDriftDiffusion1D_SL_lomac::update_rho(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const int& xDiv = mesh1D_->getxDiv();
  const real_t& hx = mesh1D_->gethx();
  const real_t& hv = mesh1D_->gethv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const Matrix& Cellcenter = mesh1D_->getCellCenter();
  // *** 
  phi_L_ = phi_bc(x1, Trun);
  phi_R_ = phi_bc(x2, Trun);
  rho_L_ = rho_numericalbc(x1, Trun);
  rho_R_ = rho_numericalbc(x2, Trun);
  // std::cout << "rho_L_ = " << rho_L_ << std::endl;
  // std::cout << "rho_R_ = " << rho_R_ << std::endl;
  rho_L_ = 0.9;
  rho_R_ = 0.1;
  // PAUSE();
  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));
  real_t coef = dt / (hx * hx);

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
      for (int i = 0; i < xDiv; i++) {
        xtemp = Cellcenter(0, i) - 1.e0 / eps_ * V(j) * dt;
        if (i == 0)
        {
          vtemp = V(j) - 1.e0 / eps_ * (phi_(0) - phi_L_) / hx * dt;
        } else if (i == xDiv -1)
        {
          vtemp = V(j) - 1.e0 / eps_ * (phi_R_ - phi_(Nx-1)) / hx * dt;
        } else 
        {
          vtemp = V(j) - 1.e0 / eps_ * (phi_(i) - phi_(i-1)) / hx * dt;
        };
        M = Maxwell(vtemp);
        // xtemp\in [x_{is-1},x_{is}]; vtemp\in [v_{js-1},v_{js}];
        // 0 -> is-2, 1->is-1, 2->is, 3->is+1
        is = floor((xtemp - X(0)) / hx) + 1;
        js = floor((vtemp - V(0)) / hv) + 1;
        xita = (is * hx + X(0) - xtemp) / hx;
        vita = (js * hv + V(0) - vtemp) / hv;
        interpolateSL_x(is, js, Cellcenter(0, i), V(j), Trun, &fx, &fv, &rhotemp);
        if (js - 1 < 0 || js - 1 >= Nv) {
          ftemp = 0.e0;
        } else {
          ftemp = (is - 1 < 0) ? fL_bc(js-1, Trun) : 
                    (is - 1 >= Nx) ? fR_bc(js-1, Trun) : f_[js-1](is - 1);
        }
        switch (x_order_)
        {
        case 1:
          Fx_[j](i) = (xita * vita * ftemp 
                    + (1.e0 - xita) * vita * fv(1) 
                    + xita * (1.e0 - vita) * fx(1) 
                    + (1.e0 - xita) * (1.e0 - vita) * fx(2)
                    - xita * rhotemp(1) * M - (1.e0 - xita) * rhotemp(2) * M) / eps_;
          Fx_[j](i) = omega * V(j) * Fx_[j](i)
                    + eps_ / sigmas_ * (1.e0 - omega) *  
                      V(j) * fsource(Cellcenter(0, i), V(j), Trun + dt);
          break;
        default:
          QUEST_ERROR(" The higher order (>= 1) SL method has not been implemented ! ");
          break;
        }
      }
      // #pragma omp critical
      // SLflux_ += Fx_[j] * Vweights(j) * V(j);
    }
  };

  Source_v_.resize(Nx);
  Source_ = source(X, Trun + dt);
  SLflux_ = Vector::Zero(Nx);
  for (int j = 0; j < Nv; j++) {
    for (int i = 0; i < Nx; i++) {
      SLflux_(i) += (Fx_[j](i + 1) - Fx_[j](i)) * Vweights(j) / hx;
    }
  }

  // 开始迭代
  real_t err = 1.e0;
  Vector bf;
  rhon_ = rho_;
  iter_ = 0;
  Vector prerhotemp = rho_;
  while (err > iterationtol_) {
    update_E(Trun);
    generateS(Trun, dt);
    generatebf(Trun, dt, &bf);
    Vector rhotemp = rho_;
    SolveS(Trun, bf);  // update rho_
    if (!rho_.allFinite() || !bf.allFinite()) {
      std::cout << "[ERROR] rho_ contains NaN or Inf after SolveS!" << std::endl;
      std::cout << " rho_.allFinite() = " << rho_.allFinite() << std::endl;
      std::cout << " bf.allFinite() = " << bf.allFinite() << std::endl;
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

void KineticDriftDiffusion1D_SL_lomac::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));
  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  real_t phi_temp_L, phi_temp_R;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (2.e0 - phi_(i) / theta_) * (1.e0 - omega) * coef * D_ / sigmas_;
    phi_temp_L = (i == 0) ? phi_L_ : phi_(i - 1);
    val = val + phi_temp_L / 2.e0 * (1.e0 - omega) * coef * D_ / sigmas_;
    phi_temp_R = (i == Nx-1) ? phi_R_ : phi_(i + 1);
    val = val + phi_temp_R / 2.e0 * (1.e0 - omega) * coef * D_ / sigmas_;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = (- 1.e0 - (phi_(i) - phi_(i-1)) / 2.e0 / theta_)   
          * (1.e0 - omega) * coef * D_ / sigmas_;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = (- 1.e0 + (phi_(i + 1) - phi_(i)) / 2.e0 / theta_)   
          * (1.e0 - omega) * coef * D_ / sigmas_;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  };

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());
  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    // std::cout << "  Preparation for BICG the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
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

// the right hand side of the drift-diffusion equation
void KineticDriftDiffusion1D_SL_lomac::generatebf(const real_t& Trun, 
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
  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));
  real_t coef = dt / (hx * hx);

  *bf = rhon_ - dt * SLflux_ + dt * Source_;
  // (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ * coef * rho_L_;
  // (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ / theta_ 
  //          * (phi_(0) - phi_L_) / 2.e0 * coef * rho_L_;
  (*bf)(0) = (*bf)(0) + 
    (1.e0 + (phi_(0) - phi_L_) / 2.e0 / theta_)   
      * (1.e0 - omega) * coef * D_ / sigmas_ * rho_L_;

  // (*bf)(Nx - 1) = (*bf)(Nx - 1) + (1.e0 - omega) * D_ * coef * rho_R_;
  // (*bf)(Nx - 1) = (*bf)(Nx - 1) - (1.e0 - omega) * D_ / theta_ 
  //               * (phi_R_ - phi_(Nx-1)) / 2.e0 * coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) + 
    (1.e0 - (phi_R_ - phi_(Nx-1)) / 2.e0 / theta_)   
      * (1.e0 - omega) * coef * D_ / sigmas_ * rho_R_;
};

void KineticDriftDiffusion1D_SL_lomac::SolveS(const real_t& Trun, const Vector& bf)
{
  KineticDriftDiffusion1D_SL::SolveS(Trun, bf);
}

void KineticDriftDiffusion1D_SL_lomac::updateAll(const real_t& Trun, 
                                    const real_t& dt)
{
  // ** 传入相关变量
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  // ***

  update_E(Trun);
  // character_tracing_half(Trun, dt);
  update_rho(Trun, dt);
  character_tracing(Trun, dt);
  update_f(Trun, dt);
  Vector rho_new = Vector::Zero(Nx);
  for (int j = 0; j < Nv; j++) {
    rho_new = rho_new + Vweights(j) * f_[j];
  };

  real_t Mj;
  for (int j = 0; j < Nv; j++) {
    Mj = Maxwell(V(j));
    f_[j] = f_[j] - rho_new * Mj + rho_ * Mj;
  }; // 为了守恒而作的后处理
  
};

void KineticDriftDiffusion1D_SL_lomac::
  character_tracing(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // dconst int& xDiv = mesh1D_->getxDiv();
  const real_t& hx = mesh1D_->gethx();
  const real_t& hv = mesh1D_->gethv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const Matrix& Cellcenter = mesh1D_->getCellCenter();
  // ***
#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t vtemp, xtemp, ftemp;
  int is, js, ii;
  real_t xita, vita;
  Vector fx(4), fv(4), rhotemp(4);
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    finter_[j].setZero();
    for (int i = 0; i < Nx; i++) {
      xtemp = X(i) - 1.e0 / eps_ * V(j) * dt;
      vtemp = V(j) - 1.e0 / eps_ * E_(i+1) * dt;
      is = floor((xtemp - X(0)) / hx) + 1;
      js = floor((vtemp - V(0)) / hv) + 1;
      xita = (is * hx + X(0) - xtemp) / hx;
      vita = (js * hv + V(0) - vtemp) / hv;
      interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
      if (js - 1 < 0 || js - 1 >= Nv) {
        ftemp = 0.e0;
      } else {
        ftemp = (is - 1 < 0) ? fL_bc(js-1, Trun) : 
                  (is - 1 >= Nx) ? fR_bc(js-1, Trun) : f_[js-1](is - 1);
      }
      switch (x_order_) {
        case 1:
          finter_[j](i) = xita * vita * ftemp 
                      + (1.e0 - xita) * vita * fv(1) 
                      + xita * (1.e0 - vita) * fx(1) 
                      + (1.e0 - xita) * (1.e0 - vita) * fx(2);
          break;
        default:
          QUEST_ERROR(" The higher order (>= 1) SL method has not been implemented ! ");
          break;
      }
    }
  }
}
};

void KineticDriftDiffusion1D_SL_lomac::update_f(const real_t& Trun, const real_t& dt)
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  // std::vector<Vector> ftemp = f_;
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  real_t ee2 = 2 * eps2 + dt;
  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));

#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t M;
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    M = Maxwell(V(j));
    // Vector rhotemp = rho_real(X, Trun + dt);
    switch (time_stepping_)
    {
    case TimeSteppingType::Exponential:
      if (t_order_ == 1) {
        f_[j] = omega * finter_[j] + (1 - omega) * rho_ * M;
        Source_f_[j] = fsource(X, V(j), Trun + dt);
        f_[j] = f_[j] + (1 - omega) * eps2 / sigmas_ * Source_f_[j];
      } else 
      {
        QUEST_ERROR(" The high order (>= 1) SL method has not been implemented ! ");
      };
      break;

    case TimeSteppingType::Rktype:
      if (t_order_ == 1) {
        f_[j] = (eps2 / ee) * finter_[j] + (dt / ee) * rho_ * M;
        Source_f_[j] = fsource(X, V(j), Trun + dt);
        f_[j] = f_[j] + dt * (eps2 / ee) * Source_f_[j];
      } else 
      {
        QUEST_ERROR(" The high order (>= 1) SL method has not been implemented ! ");
      };
      break;
    
    default:
      QUEST_ERROR(" The time stepping type has not been implemented ! ");
      break;
    }
  };
}
};

void KineticDriftDiffusion1D_SL_lomac::interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp)
{
  KineticDriftDiffusion1D_SL::interpolateSL_x(is, js, x0, v0, Trun, fx, fv, rhotemp);
};

KineticDriftDiffusion1D_SL_lomac_period::KineticDriftDiffusion1D_SL_lomac_period(
                        const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_lomac(mesh1D, poisol, x_order, t_order) {};

real_t KineticDriftDiffusion1D_SL_lomac_period::rho_init(const real_t& x)
{
  return KineticDriftDiffusion1D_SL_period::rho_init(x);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::rho_init(const Vector& x)
{
  return KineticDriftDiffusion1D_SL_period::rho_init(x);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::rho_d(const real_t& x)
{
  return KineticDriftDiffusion1D_SL_period::rho_d(x);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::rho_d(const Vector& x)
{
  return KineticDriftDiffusion1D_SL_period::rho_d(x);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::source(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::source(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::source(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::source(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::fsource(const real_t& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::fsource(x, v, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::fsource(const Vector& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::fsource(x, v, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::fsource_dx(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::fsource_dx(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::fsource_dx(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::fsource_dx(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::f_init(const real_t& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period::f_init(x, v);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::f_init(const Vector& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period::f_init(x, v);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::g_init(const real_t& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period::g_init(x, v);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::g_init(const Vector& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period::g_init(x, v);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::
  rho_real(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::rho_real(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::
  rho_real(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::rho_real(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::
  E_real(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::E_real(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period::
  g_real(const real_t& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::g_real(x, v, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period::
  g_real(const Vector& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period::g_real(x, v, t);
};

void KineticDriftDiffusion1D_SL_lomac_period::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const int& xDiv = mesh1D_->getxDiv();
  // *** 
  QUEST_VERIFY(x_order_ == t_order_, "x_order must be equal to t_order (space order == time order) !");
  QUEST_VERIFY(xDiv == Nx, " The mesh must be the periodical mesh ! ");

  pi_ = 3.14159265358979323846264338327;
  iter_ = 0;
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
  std::cout << " eps_ = " << eps_ << std::endl;
  std::cout << " gamma_ = " << gamma_ << std::endl;
  std::cout << " D_ = " << D_ << std::endl;
  std::cout << " theta_ = " << theta_ << std::endl;
  std::cout << " cfl_ = " << cfl_ << std::endl;
  std::cout << " NTH_ = " << NTH_ << std::endl;
};

void KineticDriftDiffusion1D_SL_lomac_period::
  character_tracing(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // dconst int& xDiv = mesh1D_->getxDiv();
  const real_t& hx = mesh1D_->gethx();
  const real_t& hv = mesh1D_->gethv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const Matrix& Cellcenter = mesh1D_->getCellCenter();
  // ***
#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t vtemp, xtemp, ftemp;
  int is, js, ii;
  real_t xita, vita;
  Vector fx(4), fv(4), rhotemp(4);
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    finter_[j].setZero();
    for (int i = 0; i < Nx; i++) {
      xtemp = X(i) - 1.e0 / eps_ * V(j) * dt;
      vtemp = V(j) - 1.e0 / eps_ * E_(i) * dt;
      is = floor((xtemp - X(0)) / hx) + 1;
      js = floor((vtemp - V(0)) / hv) + 1;
      xita = (is * hx + X(0) - xtemp) / hx;
      vita = (js * hv + V(0) - vtemp) / hv;
      interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
      ii = (((is - 1) % Nx) + Nx) % Nx;
      ftemp = (js - 1 >= 0 && js - 1 < Nv) ? f_[js - 1](ii) : 0.0;
      switch (x_order_) {
        case 1:
          finter_[j](i) = xita * vita * ftemp 
                      + (1.e0 - xita) * vita * fv(1) 
                      + xita * (1.e0 - vita) * fx(1) 
                      + (1.e0 - xita) * (1.e0 - vita) * fx(2);
          break;
        default:
          QUEST_ERROR(" The higher order (>= 1) SL method has not been implemented ! ");
          break;
      }
    }
  }
}
};

void KineticDriftDiffusion1D_SL_lomac_period::update_rho(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // dconst int& xDiv = mesh1D_->getxDiv();
  const real_t& hx = mesh1D_->gethx();
  const real_t& hv = mesh1D_->gethv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const Matrix& Cellcenter = mesh1D_->getCellCenter();
  // ***
  // Vector temp = Cellcenter.row(0).transpose();
  // std::cout << " temp = \n" << temp.transpose() << std::endl;
  // std::cout << " X = \n" << X.transpose() << std::endl;
  // std::cout << " temp.rows() = " << temp.rows() << std::endl;
  // std::cout << " X.rows() = " << X.rows() << std::endl;
  // PAUSE();
  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));
  #pragma omp parallel num_threads(NTH_), default(shared)
  {
    real_t vtemp, xtemp, ftemp;
    int is, js, ii;
    real_t xita, vita;
    real_t M;
    Vector fx(4), fv(4), rhotemp(4);
  #pragma omp for schedule(static)
    for (int j = 0; j < Nv; j++) {
      Fx_[j].setZero();
      for (int i = 0; i < Nx; i++) {
        xtemp = Cellcenter(0, i) - 1.e0 / eps_ * V(j) * dt;
        if (i == Nx - 1)
        {
          vtemp = V(j) - 1.e0 / eps_ * (phi_(0) - phi_(Nx - 1)) / hx * dt;
        } else 
        {
          vtemp = V(j) - 1.e0 / eps_ * (phi_(i+1) - phi_(i)) / hx * dt;
        };
        M = Maxwell(vtemp);
        // xtemp\in [x_{is-1},x_{is}]; vtemp\in [v_{js-1},v_{js}];
        // 0 -> is-2, 1->is-1, 2->is, 3->is+1
        // fx(1) --(1-xita) (xita)-- fx(2),fv(2)
        //   |                         |
        //  vita                       |
        // 1-vita                      |
        //   |                         |
        // ftemp -------------------- fv(1)
        is = floor((xtemp - X(0)) / hx) + 1;
        js = floor((vtemp - V(0)) / hv) + 1;
        xita = (is * hx + X(0) - xtemp) / hx;
        vita = (js * hv + V(0) - vtemp) / hv;
        interpolateSL_x(is, js, Cellcenter(0, i), V(j), Trun, &fx, &fv, &rhotemp);
        ii = (((is - 1) % Nx) + Nx) % Nx;
        ftemp = (js - 1 >= 0 && js - 1 < Nv) ? f_[js - 1](ii) : 0.0;
        switch (x_order_) {
        case 1:
          Fx_[j](i) = (xita * vita * ftemp 
                    + (1.e0 - xita) * vita * fv(1) 
                    + xita * (1.e0 - vita) * fx(1) 
                    + (1.e0 - xita) * (1.e0 - vita) * fx(2)
                    - xita * rhotemp(1) * M - (1.e0 - xita) * rhotemp(2) * M) / eps_;
          Fx_[j](i) = omega * V(j) * Fx_[j](i)
                    + eps_ / sigmas_ * (1.e0 - omega) *  
                      V(j) * fsource(Cellcenter(0, i), V(j), Trun + dt);
          break;
        default:
          QUEST_ERROR(" The higher order (>= 1) SL method has not been implemented ! ");
          break;
        }
      }
      // #pragma omp critical
      // SLflux_ += Fx_[j] * Vweights(j) * V(j);
    }
  };

  Source_v_.resize(Nx);
  Source_ = source(X, Trun + dt);
  SLflux_ = Vector::Zero(Nx);
  for (int j = 0; j < Nv; j++) {
    for (int i = 0; i < Nx; i++) {
      if (i == 0)
      {
        SLflux_(i) += (Fx_[j](0) - Fx_[j](Nx - 1)) * Vweights(j) / hx;
      } else 
      {
        SLflux_(i) += (Fx_[j](i) - Fx_[j](i - 1)) * Vweights(j) / hx;
      };
    }
  }

  // 开始迭代
  real_t err = 1.e0;
  Vector bf;
  rhon_ = rho_;
  iter_ = 0;
  Vector prerhotemp = rho_;
  generatebf(Trun, dt, &bf);

  while (err > iterationtol_) {
    update_E(Trun);
    generateS(Trun, dt);
    Vector rhotemp = rho_;
    SolveS(Trun, bf);  // update rho_
    if (!rho_.allFinite() || !bf.allFinite()) 
    {
      std::cout << " rho_.allFinite() = " << rho_.allFinite() << std::endl;
      std::cout << " bf.allFinite() = " << bf.allFinite() << std::endl;
      // std::ofstream Sout("result/DD_acctest/S.txt");
      for (int k = 0; k < S_.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(S_, k); it; ++it) {
          std::cout << it.row() << " " 
               << it.col() << " "
               << it.value() << "\n";
        }
      }
      // Sout.close();
      QUEST_ERROR(" rho_ contains NaN or Inf after SolveS ! ");
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

void KineticDriftDiffusion1D_SL_lomac_period::generateS(const real_t& Trun, const real_t& dt)
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));
  S_.resize(Nx, Nx);  S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t phi_temp_L, phi_temp_R;
  real_t val;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (2.e0 - phi_(i) / theta_) * (1.e0 - omega) * coef * D_ / sigmas_;
    phi_temp_L = (i == 0) ? phi_(Nx-1) : phi_(i - 1);
    val = val + phi_temp_L / theta_ / 2.e0 * (1.e0 - omega) * coef * D_ / sigmas_;
    phi_temp_R = (i == Nx-1) ? phi_(0) : phi_(i + 1);
    val = val + phi_temp_R / theta_ / 2.e0 * (1.e0 - omega) * coef * D_ / sigmas_;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = (- 1.e0 - (phi_(i) - phi_(i-1)) / 2.e0 / theta_)   
          * (1.e0 - omega) * coef * D_ / sigmas_;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = (- 1.e0 + (phi_(i + 1) - phi_(i)) / 2.e0 / theta_)   
          * (1.e0 - omega) * coef * D_ / sigmas_;
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  };

  // i = 0
  val = (- 1.e0 - (phi_(0) - phi_(Nx-1)) / 2.e0 / theta_)   
      * (1.e0 - omega) * coef * D_ / sigmas_;
  tripletList_S.emplace_back(0, Nx-1, val);
  // i = Nx-1
  val = (- 1.e0 + (phi_(0) - phi_(Nx-1)) / 2.e0 / theta_)   
      * (1.e0 - omega) * coef * D_ / sigmas_;
  tripletList_S.emplace_back(Nx-1, 0, val);

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  switch (soltype_) {
  case Solver1DTypeFD::BICG:
    // std::cout << "  Preparation for BICG the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
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

void KineticDriftDiffusion1D_SL_lomac_period::generatebf(const real_t& Trun, 
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
  *bf = rhon_ - dt * SLflux_ + dt * Source_;
};

void KineticDriftDiffusion1D_SL_lomac_period::interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp)
{
  KineticDriftDiffusion1D_SL_period::interpolateSL_x(is, js, x0, v0, Trun, fx, fv, rhotemp);
};

void KineticDriftDiffusion1D_SL_lomac_period::SolveS(const real_t& Trun, const Vector& bf)
{
  KineticDriftDiffusion1D_SL_period::SolveS(Trun, bf);
};

void KineticDriftDiffusion1D_SL_lomac_period::updateAll(const real_t& Trun,
                                    const real_t& dt)
{
  // ** 传入相关变量
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  // ***

  update_E(Trun);
  // character_tracing_half(Trun, dt);
  update_rho(Trun, dt);
  character_tracing(Trun, dt);
  update_f(Trun, dt);
  Vector rho_new = Vector::Zero(Nx);
  for (int j = 0; j < Nv; j++) {
    rho_new = rho_new + Vweights(j) * f_[j];
  };
  // rho_ = rho_new;
  real_t Mj;
  for (int j = 0; j < Nv; j++) {
    Mj = Maxwell(V(j));
    f_[j] = f_[j] - rho_new * Mj + rho_ * Mj;
  }; // 为了守恒而作的后处理
};

void KineticDriftDiffusion1D_SL_lomac_period::update_E(const real_t& Trun)
{
  Vector bf = - (rho_ - rhod_) / gamma_;
  poisol_->SolveAll(0.e0, 0.e0, bf, &phi_, &E_);
};

void KineticDriftDiffusion1D_SL_lomac_period::update_f(const real_t& Trun, 
                                                const real_t& dt)
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  // std::vector<Vector> ftemp = f_;
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  real_t ee2 = 2 * eps2 + dt;
  real_t omega = std::exp( - dt * sigmas_ / (eps_ * eps_));

#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t M;
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    M = Maxwell(V(j));
    // Vector rhotemp = rho_real(X, Trun + dt);
    switch (time_stepping_)
    {
    case TimeSteppingType::Exponential:
      if (t_order_ == 1) {
        f_[j] = omega * finter_[j] + (1 - omega) * rho_ * M;
        Source_f_[j] = fsource(X, V(j), Trun + dt);
        f_[j] = f_[j] + (1 - omega) * eps2 / sigmas_ * Source_f_[j];
      } else 
      {
        QUEST_ERROR(" The high order (>= 1) SL method has not been implemented ! ");
      };
      break;

    case TimeSteppingType::Rktype:
      if (t_order_ == 1) {
        f_[j] = (eps2 / ee) * finter_[j] + (dt / ee) * rho_ * M;
        Source_f_[j] = fsource(X, V(j), Trun + dt);
        f_[j] = f_[j] + dt * (eps2 / ee) * Source_f_[j];
      } else 
      {
        QUEST_ERROR(" The high order (>= 1) SL method has not been implemented ! ");
      };
      break;
    
    default:
      QUEST_ERROR(" The time stepping type has not been implemented ! ");
      break;
    }
  };
};
};

KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma(
                        const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_lomac(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_lomac_period(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period_Different_Gamma(mesh1D, poisol, x_order, t_order) {};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  init(const Solver1DTypeFD& soltype)
{
  KineticDriftDiffusion1D_SL_lomac_period::init(soltype);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  character_tracing(const real_t& Trun, const real_t& dt)
{
  KineticDriftDiffusion1D_SL_lomac_period::character_tracing(Trun, dt);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  update_rho(const real_t& Trun, const real_t& dt) 
{
  KineticDriftDiffusion1D_SL_lomac_period::
    update_rho(Trun, dt);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  updateAll(const real_t& Trun, const real_t& dt)
{
  KineticDriftDiffusion1D_SL_lomac_period::
    updateAll(Trun, dt);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  generateS(const real_t& Trun, const real_t& dt)
{
  KineticDriftDiffusion1D_SL_lomac_period::
    generateS(Trun, dt);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  generatebf(const real_t& Trun, 
                const real_t& dt,
                Vector* bf) 
{
  KineticDriftDiffusion1D_SL_lomac_period::
    generatebf(Trun, dt, bf);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  SolveS(const real_t& Trun, const Vector& bf)
{
  KineticDriftDiffusion1D_SL_lomac_period::
    SolveS(Trun, bf);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  update_E(const real_t& Trun)
{
  KineticDriftDiffusion1D_SL_lomac_period::
    update_E(Trun);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  update_f(const real_t& Trun, const real_t& dt) 
{
  KineticDriftDiffusion1D_SL_lomac_period::
    update_f(Trun, dt);
};

void KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  interpolateSL_x(const int& is, const int& js,
                const real_t& x0, const real_t& v0,
                const real_t& Trun,
                Vector* fx, Vector* fv, Vector* rhotemp)
{
  KineticDriftDiffusion1D_SL_lomac_period::
    interpolateSL_x(is, js, x0, v0, Trun, fx, fv, rhotemp);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  rho_init(const real_t& x) 
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_init(x);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  rho_init(const Vector& x)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_init(x);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  rho_d(const real_t& x)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_d(x);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  rho_d(const Vector& x)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_d(x);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  source(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::source(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  source(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::source(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  fsource(const real_t& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource(x, v, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  fsource(const Vector& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource(x, v, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  fsource_dx(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource_dx(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  fsource_dx(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource_dx(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  f_init(const real_t& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::f_init(x, v);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  f_init(const Vector& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::f_init(x, v);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  g_init(const real_t& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::g_init(x, v);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  g_init(const Vector& x, const real_t& v)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::g_init(x, v);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  rho_real(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_real(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  rho_real(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_real(x, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  E_real(const Vector& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::E_real(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  g_real(const real_t& x, const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::g_real(x, v, t);
};

Vector KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma::
  g_real(const Vector& x, const real_t& v, const real_t& t) 
{
  return KineticDriftDiffusion1D_SL_period_Different_Gamma::g_real(x, v, t);
};


// KineticDriftDiffusion1D_SL_lomac_PNjunction
KineticDriftDiffusion1D_SL_lomac_PNjunction::
  KineticDriftDiffusion1D_SL_lomac_PNjunction(const KineticFDmesh* mesh1D,
                          Poissonsolver1DFD* poisol,
                          const int& x_order,
                          const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_lomac(mesh1D, poisol, x_order, t_order) {};;

void KineticDriftDiffusion1D_SL_lomac_PNjunction::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const int& xDiv = mesh1D_->getxDiv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  QUEST_VERIFY(x_order_ == t_order_, "x_order must be equal to t_order (space order == time order) !");
  QUEST_VERIFY(xDiv != Nx, " The mesh is not the periodical mesh ! ");
  QUEST_VERIFY(x1 == 0.e0, " The left boundary is equal to one ! ");
  QUEST_VERIFY(x2 == 1.e0, " The right boundary is equal to one ! ");

  pi_ = 3.14159265358979323846264338327;
  iter_ = 0;
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
    Fx_[j] = Vector::Zero(xDiv);
    finter_[j] = Vector::Zero(Nx);
  };
  std::cout << " eps_ = " << eps_ << std::endl;
  std::cout << " gamma_ = " << gamma_ << std::endl;
  std::cout << " D_ = " << D_ << std::endl;
  std::cout << " theta_ = " << theta_ << std::endl;
  std::cout << " cfl_ = " << cfl_ << std::endl;
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::rho_init(const real_t& x)
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::rho_init(const Vector& x)
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

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  rho_numericalbc(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::rho_numericalbc(x, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  phi_bc(const real_t& x, const real_t & t)
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

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  rho_d(const real_t& x)
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::
  rho_d(const Vector& x)
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

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  source(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::
  source(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fsource(const real_t& x, const real_t& v, const real_t& t)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fsource(const Vector& x, const real_t& v, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fsource_dx(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fsource_dx(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  f_init(const real_t& x, const real_t& v)
{
  return rho_init(x) * Maxwell(v) + eps_ * g_init(x, v);
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::
  f_init(const Vector& x, const real_t& v)
{
  return rho_init(x) * Maxwell(v) + eps_ * g_init(x, v);
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  g_init(const real_t& x, const real_t& v)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_lomac_PNjunction::
  g_init(const Vector& x, const real_t& v)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fL_bc(const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::fL_bc(v, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fL_bc(const int& v_index, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(v_index) >= 0) {
    f = 0.9e0 * Maxwell(V(v_index));
  } else {
    f = 2.e0 * f_[v_index](0) - f_[v_index](1);
  }
  return f;
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fR_bc(const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::fR_bc(v, t);
};

real_t KineticDriftDiffusion1D_SL_lomac_PNjunction::
  fR_bc(const int& v_index, const real_t& t)
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

} // namespace QUEST

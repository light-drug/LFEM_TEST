#include "KineticLinearDiffusion_SL.hpp"

namespace QUEST
{

KineticLinearDiffusion1D_SL::KineticLinearDiffusion1D_SL(const KineticFDmesh* mesh1D,
                      Poissonsolver1DFD* poisol,
                      const int& x_order,
                      const int& t_order)
: KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order) {};

void KineticLinearDiffusion1D_SL::setsigmas(const real_t& sigmas)
{
  sigmas_ = sigmas; 
};

real_t KineticLinearDiffusion1D_SL::rho_init(const real_t& x)
{
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL::rho_init(const Vector& x)
{
  Vector rho = x;
  rho.setZero();
  return rho;
};

real_t KineticLinearDiffusion1D_SL::rho_numericalbc(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::rho_numericalbc(x, t);
};  

real_t KineticLinearDiffusion1D_SL::source(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL::source(const Vector& x, const real_t& t)
{
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticLinearDiffusion1D_SL::fsource(const real_t& x, const real_t& v, const real_t& t)
{
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL::fsource(const Vector& x, const real_t& v, const real_t& t)
{
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticLinearDiffusion1D_SL::fsource_dx(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL::fsource_dx(const Vector& x, const real_t& t)
{
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticLinearDiffusion1D_SL::f_init(const real_t& x, const real_t& v)
{
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL::f_init(const Vector& x, const real_t& v)
{
  Vector f = x;
  f.setZero();
  return f;
};

real_t KineticLinearDiffusion1D_SL::g_init(const real_t& x, const real_t& v)
{
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL::g_init(const Vector& x, const real_t& v)
{
  Vector g = x;
  g.setZero();
  return g;
};

real_t KineticLinearDiffusion1D_SL::fL_bc(const real_t& v, const real_t& t)
{
  QUEST_ERROR(" This fL_bc(real_t, real_t) is not valid ! ");
  return 0.e0;
};

real_t KineticLinearDiffusion1D_SL::fL_bc(const int& v_index, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(v_index) >= 0) {
    f = 1.e0;
  } else {
    f = 2.e0 * f_[v_index](0) - f_[v_index](1);
  }
  return f;
};

real_t KineticLinearDiffusion1D_SL::fR_bc(const real_t& v, const real_t& t)
{
  QUEST_ERROR(" This fR_bc(real_t, real_t) is not valid ! ");
  return 0.e0;
};

real_t KineticLinearDiffusion1D_SL::fR_bc(const int& v_index, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  // *** 
  real_t f;
  if (V(v_index) <= 0) {
    f = 0.e0;
  } else {
    f = 2.e0 * f_[v_index](Nx-1) - f_[v_index](Nx-2);
  }
  return f;
};

void KineticLinearDiffusion1D_SL::init(const Solver1DTypeFD& soltype)
{
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  // *** 
  QUEST_VERIFY(std::abs(x1 - 0.e0) < 1.e-13, " x1 has to be equal to 0.e0 !");
  QUEST_VERIFY(std::abs(x2 - 1.e0) < 1.e-13, " x2 has to be equal to 1.e0 !");
  QUEST_VERIFY(std::abs(v1 + 1.e0) < 1.e-13, " v1 has to be equal to -1.e0 !");
  QUEST_VERIFY(std::abs(v2 - 1.e0) < 1.e-13, " v2 has to be equal to 1.e0 !");
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
};

void KineticLinearDiffusion1D_SL::character_tracing(const real_t& Trun, 
                                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
  for (int j = 0; j < Nv; j++) {
    for (int i = 0; i < Nx; i++) {
      xback_[j](i) = X(i) - 1.e0 / eps_ * V(j) * dt;
      vback_[j](i) = V(j);
    }
  }
};

void KineticLinearDiffusion1D_SL::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  real_t omega = std::exp( - dt / (eps_ * eps_));
  // real_t omega = 0.e0;
  S_.resize(Nx, Nx); S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + 2.e0 * (1.e0 - omega) * coef / sigmas_;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      val = - (1.e0 - omega) * coef / sigmas_;
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      val = - (1.e0 - omega) * coef / sigmas_;
      tripletList_S.emplace_back(i, i + 1, val);
    }
  };
  
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

void KineticLinearDiffusion1D_SL::generatebf(const real_t& Trun, 
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

  *bf = rhon_ - omega * dt * SLflux_ + Source_ + Source_v_;
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * coef / sigmas_ * rho_L_;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + (1.e0 - omega) * coef / sigmas_ * rho_R_;
};

void KineticLinearDiffusion1D_SL::SolveS(const real_t& Trun, const Vector& bf) {
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

void KineticLinearDiffusion1D_SL::update_rho(const real_t& Trun, 
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
  SLflux_.setZero();
  rho_L_ = rho_numericalbc(x1, Trun);
  rho_R_ = rho_numericalbc(x2, Trun);

  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;

#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t vtemp, xtemp, ftemp;
  int is, js, ii;
  real_t xita;
  Vector fx(4), fv(4), rhotemp(4);
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    Fx_[j].setZero();
    for (int i = 0; i < Nx; i++) {
      vtemp = vback_[j](i);
      xtemp = xback_[j](i);
      // xtemp\in [x_{is-1},x_{is}];
      // 0 -> is-2, 1->is-1, 2->is, 3->is+1
      if (x_order_ == 1) 
      {
        if (V(j) >= 0) 
        {
          is = floor((xtemp - X(0)) / hx) + 1;
          js = j;
          // xtemp\in [x_{is-1},x_{is});
          interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
          xita = (is * hx + X(0) - xtemp) / hx;
          Fx_[j](i) = (fx(1) - fx(0)) / hx - (rhotemp(3) - rhotemp(2)) / hx;
          Fx_[j](i) = Fx_[j](i) / eps_;
        } else if (V(j) < 0) 
        {
          is = ceil((xtemp - X(0)) / hx);
          js = j;
          // xtemp\in (x_{is-1},x_{is}];
          xita = (is * hx + X(0) - xtemp) / hx;
          interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
          Fx_[j](i) = (fx(3) - fx(2)) / hx - (rhotemp(1) - rhotemp(0)) / hx;
          Fx_[j](i) = Fx_[j](i) / eps_;
        }
        finter_[j](i) = xita * fx(1) + (1.e0 - xita) * fx(2);
      } else if (x_order_ == 2) 
      {
        QUEST_ERROR(" The second SL method has not been implemented ! ");
      }
    }
    #pragma omp critical
    SLflux_ = SLflux_ + Fx_[j] * Vweights(j) * V(j);
  };
}

  Source_v_.resize(Nx); Source_v_.setZero();
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);

  Vector bf;
  rhon_ = rho_;
  // PAUSE();
  generateS(Trun, dt);
  generatebf(Trun, dt, &bf);
  SolveS(Trun, bf);  // update rho_
  if (!rho_.allFinite() || !bf.allFinite()) {
    std::cout << "[ERROR] rho_ contains NaN or Inf after SolveS!" << std::endl;
    std::ofstream Sout("result/LD_acctest/S.txt");
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
  };
};

void KineticLinearDiffusion1D_SL::updateAll(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  // *** 

  character_tracing(Trun, dt);
  update_rho(Trun, dt);
  update_f(Trun, dt);
  rho_.setZero();
  for (int j = 0; j < Nv; j++) {
    rho_ = rho_ + Vweights(j) * f_[j];
  }
};

void KineticLinearDiffusion1D_SL::update_f(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  std::vector<Vector> ftemp = f_;
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  real_t ee2 = 2 * eps2 + dt;
#pragma omp parallel num_threads(NTH_), default(shared)
{
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    // Vector rhotemp = rho_real(X, Trun + dt);
    if (t_order_ == 1) {
      ftemp[j] = (eps2 / ee) * finter_[j] + (dt / ee) * rho_;
      Source_f_[j] = fsource(X, V(j), Trun + dt);
      f_[j] = ftemp[j] + dt * (eps2 / ee) * Source_f_[j];
    } 
    else {
      QUEST_ERROR(" The high order (>= 2) SL method has not been implemented ! ");
    }
    
  };
};
};

void KineticLinearDiffusion1D_SL::interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) {
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
  real_t fL, fR;
  fv->setZero();
  if (x_order_ == 1) {
    if (js >= 0 && js < Nv) { // v方向拉出边界
      fL = fL_bc(js, Trun); // x方向拉出
      fR = fR_bc(js, Trun);
      (*fx)(0) = (is - 2 < 0) ? fL : 
                  (is - 2 >= Nx) ? fR : f_[js](is - 2);
      (*fx)(1) = (is - 1 < 0) ? fL : 
                  (is - 1 >= Nx) ? fR : f_[js](is - 1);
      (*fx)(2) = (is < 0) ? fL : 
                  (is >= Nx) ? fR : f_[js](is);
      (*fx)(3) = (is + 1 < 0) ? fL : 
                  (is + 1 >= Nx) ? fR : f_[js](is + 1);
    } else {
      fx->setZero();
    };

    (*rhotemp)(0) = (is - 2 < 0) ? rho_L_ : 
                (is - 2 >= Nx) ? rho_R_ : rho_(is - 2);
    (*rhotemp)(1) = (is - 1 < 0) ? rho_L_ : 
                (is - 1 >= Nx) ? rho_R_ : rho_(is - 1);
    (*rhotemp)(2) = (is < 0) ? rho_L_ : 
                (is >= Nx) ? rho_R_ : rho_(is);
    (*rhotemp)(3) = (is + 1 < 0) ? rho_L_ : 
                (is + 1 >= Nx) ? rho_R_ : rho_(is + 1);
  } else if (x_order_ == 2) {
    QUEST_ERROR(" The second SL method has not been implemented ! ");
  }
}

KineticLinearDiffusion1D_SL_period::KineticLinearDiffusion1D_SL_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order) {};


real_t KineticLinearDiffusion1D_SL_period::rho_init(const real_t& x) {
  return std::sin(x);
};

Vector KineticLinearDiffusion1D_SL_period::rho_init(const Vector& x) {
  Vector rho = x.array().sin();
  return rho;
};

real_t KineticLinearDiffusion1D_SL_period::rho_real(const real_t& x, const real_t& t) {
  return std::exp(- D_ * t) * std::sin(x);
};

Vector KineticLinearDiffusion1D_SL_period::rho_real(const Vector& x, const real_t& t) {
  Vector rho = std::exp(- D_ * t) * x.array().sin();
  return rho;
};

real_t KineticLinearDiffusion1D_SL_period::source(const real_t& x,  const real_t& t) {
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL_period::source(const Vector& x,  const real_t& t) {
  Vector rhod = 0.e0 * x;
  return rhod;
};

real_t KineticLinearDiffusion1D_SL_period::fsource(const real_t& x, 
              const real_t& v, 
              const real_t& t) {
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  real_t S;
  S = std::exp(- D_ * t) * M * (eps_ * v * D_ * std::cos(x) - D_ * sin(x) + v2 * sin(x));
  return S;
};

Vector KineticLinearDiffusion1D_SL_period::fsource(const Vector& x, 
              const real_t& v, 
              const real_t& t) {
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  Vector S;
  Vector sinx = x.array().sin();
  Vector cosx = x.array().cos();
  S = std::exp(- D_ * t) * M * (eps_ * v * D_ * cosx - D_ * sinx + v2 * sinx);
  return S;
};

real_t KineticLinearDiffusion1D_SL_period::fsource_dx(const real_t& x, 
              const real_t& t) {
  // <vS> = eps_ * D_ * D_ * std::cos(x) * std::exp(- D_ * t)
  real_t S = - eps_ * D_ * D_ * std::sin(x) * std::exp(- D_ * t);
  return S;
};

Vector KineticLinearDiffusion1D_SL_period::fsource_dx(const Vector& x, 
              const real_t& t) {
  Vector S = - eps_ * D_ * D_ * x.array().sin() * std::exp(- D_ * t);
  return S;
};


real_t KineticLinearDiffusion1D_SL_period::f_init(const real_t& x, const real_t& v) {
  real_t f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Vector KineticLinearDiffusion1D_SL_period::f_init(const Vector& x, const real_t& v) {
  Vector f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t KineticLinearDiffusion1D_SL_period::g_init(const real_t& x, const real_t& v) {
  real_t M = Maxwell(v);
  real_t g =  - v * M * std::cos(x);
  return g;
};

Vector KineticLinearDiffusion1D_SL_period::g_init(const Vector& x, const real_t& v) {
  real_t M = Maxwell(v);
  Vector g =  - v * M * x.array().cos();
  return g;
};

real_t KineticLinearDiffusion1D_SL_period::f_real(const real_t& x, const real_t& v, const real_t& t) {
  real_t f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

Vector KineticLinearDiffusion1D_SL_period::f_real(const Vector& x, const real_t& v, const real_t& t) {
  Vector f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

real_t KineticLinearDiffusion1D_SL_period::g_real(const real_t& x, const real_t& v, const real_t& t) {
  real_t M = Maxwell(v);
  real_t g =  - v * M * (std::exp(- D_ * t) * std::cos(x));
  return g;
};

Vector KineticLinearDiffusion1D_SL_period::g_real(const Vector& x, const real_t& v, const real_t& t) {
  real_t M = Maxwell(v);
  Vector g =  - v * M * (std::exp(- D_ * t) * x.array().cos());
  return g;
};

void KineticLinearDiffusion1D_SL_period::character_tracing(const real_t& Trun, 
                                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
  for (int j = 0; j < Nv; j++) {
    for (int i = 0; i < Nx; i++) {
      xback_[j](i) = X(i) - 1.e0 / eps_ * V(j) * dt;
      vback_[j](i) = V(j);
    }
  }
};

void KineticLinearDiffusion1D_SL_period::generateS(const real_t& Trun, 
                                    const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const real_t& hx = mesh1D_->gethx();
  // *** 

  real_t omega = std::exp( - dt / (eps_ * eps_));
  // real_t omega = 0.e0;
  S_.resize(Nx, Nx); S_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_S;
  int estimatedNonZeros = 3 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (1.e0 - omega) * D_ * 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, - (1.e0 - omega) * D_  * coef);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, - (1.e0 - omega) * D_  * coef);
    } 
  }
  // i = 0
  tripletList_S.emplace_back(0, Nx-1, - (1.e0 - omega) * D_  * coef);

  // i = Nx-1
  tripletList_S.emplace_back(Nx-1, 0, - (1.e0 - omega) * D_  * coef);

  S_.setFromTriplets(tripletList_S.begin(), tripletList_S.end());

  switch (soltype_) {
  case Solver1DTypeFD::CG:
    // std::cout << "  Preparation for CG the matrix ......\n";
    // std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    cg_.compute(S_);
    cg_.setTolerance(sparsetol_);
    // std::cout << "  The end of preparation for CG method ! " << std::endl;
    break;
  case Solver1DTypeFD::LDLT:
// #ifndef QUEST_USE_MKL
//     std::cout << "  Preparation for LDLT the matrix by Eigen ......\n";
// #else 
//     std::cout << "  Preparation for LU the matrix by Pardiso ......\n";
// #endif 
//     std::cout << "  The size of Matrix = " << S_.rows() << std::endl;
    ldlt_.analyzePattern(S_);
    ldlt_.factorize(S_);
    // std::cout << "  The end of preparation for LDLT method ! " << std::endl;
    break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
};

void KineticLinearDiffusion1D_SL_period::generatebf(const real_t& Trun, 
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
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);
  // *bf = rhon_ + Source_;
  *bf = rhon_ - omega * dt * SLflux_ + Source_v_ + Source_;
  // std::cout << " ===== bf ===== \n" << *bf <<  std::endl;
  // PAUSE();
};

void KineticLinearDiffusion1D_SL_period::SolveS(const real_t& Trun, const Vector& bf) {
  switch (soltype_) {
  case Solver1DTypeFD::CG:
    rho_ = cg_.solve(bf);
    if (cg_.info() != Eigen::Success) {
      std::cout << "[ERROR] CG failed to converge!\n";
    }
    break;
  case Solver1DTypeFD::LDLT:
#ifndef QUEST_USE_MKL
    std::cout << "  Preparation for LDLT the matrix by Eigen ......\n";
#else 
    std::cout << "  Preparation for LDLT the matrix by Pardiso ......\n";
#endif 
    rho_ = ldlt_.solve(bf);
    break;
  default:
    QUEST_ERROR(" The Dift Diffusion Solver is not been implemented ! ");
    break;
  }
};

void KineticLinearDiffusion1D_SL_period::update_rho(const real_t& Trun, 
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
  SLflux_.setZero();
  Vector fx(4), fv(4), rhotemp(4);
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  real_t xita;
  for (int j = 0; j < Nv; j++) {
    Fx_[j].setZero();
    // std::cout << f_[j] << std::endl;
    // PAUSE();
    for (int i = 0; i < Nx; i++) {
      real_t vtemp = vback_[j](i);
      real_t xtemp = xback_[j](i);
      // xtemp\in [x_{is-1},x_{is}];
      // 0 -> is-2, 1->is-1, 2->is, 3->is+1
      // std::cout << " ================== " << std::endl;
      // std::cout << "X(i) = " << X(i) << ", xtemp = " << xtemp << ",vtemp = " << vtemp << std::endl;
      
      if (x_order_ == 1) {
        real_t M = Maxwell(vtemp);
        if (V(j) >= 0) {
          int is = floor((xtemp - X(0)) / hx) + 1;
          int js = j;
          // xtemp\in [x_{is-1},x_{is});
          interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
          xita = is * hx + X(0) - xtemp;
          Fx_[j](i) = (fx(1) - fx(0)) / hx - (rhotemp(3) - rhotemp(2)) / hx * M;
          Fx_[j](i) = Fx_[j](i) / eps_;
        } else if (V(j) < 0) {
          int is = ceil((xtemp - X(0)) / hx);
          int js = j;
          // xtemp\in (x_{is-1},x_{is}];
          xita = is * hx + X(0) - xtemp;
          interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
          Fx_[j](i) = (fx(3) - fx(2)) / hx - (rhotemp(1) - rhotemp(0)) / hx * M;
          Fx_[j](i) = Fx_[j](i) / eps_;
        }
        // std::cout << " V(j) = " << V(j) << ", fx = " << fx << std::endl;
        finter_[j](i) = xita / hx * fx(1) + (1.e0 - xita / hx) * fx(2);
      } else if (x_order_ == 2) {
        QUEST_ERROR(" The second SL method has not been implemented ! ");
      }
    }
    SLflux_ = SLflux_ + Fx_[j] * Vweights(j) * V(j);
  }
  Vector bf;
  rhon_ = rho_;
  // PAUSE();
  generateS(Trun, dt);
  generatebf(Trun, dt, &bf);
  SolveS(Trun, bf);  // update rho_
  if (!rho_.allFinite() || !bf.allFinite()) {
    std::cout << "[ERROR] rho_ contains NaN or Inf after SolveS!" << std::endl;
    std::ofstream Sout("result/LD_acctest/S.txt");
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
};

void KineticLinearDiffusion1D_SL_period::update_f(const real_t& Trun, const real_t& dt) {
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
  std::vector<Vector> ftemp = f_;
  int is, js;
  int ii1, ii2;
  real_t fv_temp;
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  for (int j = 0; j < Nv; j++) {
    real_t M = Maxwell(V(j));
    if (x_order_ == 1) {
      ftemp[j] = (eps2 / ee) * (finter_[j]) + (dt / ee) * rho_  * M;
    } else if (x_order_ == 2) {
      QUEST_ERROR(" The second SL method has not been implemented ! ");
    }
    Source_f_[j] = fsource(X, V(j), Trun + dt);
    f_[j] = ftemp[j] + dt * (eps2 / ee) * Source_f_[j];
  };
};

void KineticLinearDiffusion1D_SL_period::interpolateSL_x(const int& is, const int& js,
                  const real_t& x0, const real_t& v0,
                  const real_t& Trun,
                  Vector* fx, Vector* fv, Vector* rhotemp) {
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
  // *** 
  int ii;
  fv->setZero();
  if (x_order_ == 1) {
    ii = (((is - 2) % Nx) + Nx) % Nx;
    (*fx)(0) = f_[js](ii);
    (*rhotemp)(0) = rho_(ii);
    ii = (((is - 1) % Nx) + Nx) % Nx;
    (*fx)(1) = f_[js](ii);
    (*rhotemp)(1) = rho_(ii);
    ii = ((is % Nx) + Nx) % Nx;
    (*fx)(2) = f_[js](ii);
    (*rhotemp)(2) = rho_(ii);
    ii = (((is + 1) % Nx) + Nx) % Nx;
    (*fx)(3) = f_[js](ii);
    (*rhotemp)(3) = rho_(ii);
  } else if (x_order_ == 2) {
    QUEST_ERROR(" The second SL method has not been implemented ! ");
  }
};

void KineticLinearDiffusion1D_SL_period::updateAll(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  // *** 

  character_tracing(Trun, dt);
  update_rho(Trun, dt);
  update_f(Trun, dt);
  rho_.setZero();
  for (int j = 0; j < Nv; j++) {
    rho_ = rho_ + Vweights(j) * f_[j];
  }
};


KineticLinearDiffusion1D_SL_period_twovel::KineticLinearDiffusion1D_SL_period_twovel(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticLinearDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order) {};


real_t KineticLinearDiffusion1D_SL_period_twovel::rho_init(const real_t& x) {
  return std::sin(x) / gamma_;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::rho_init(const Vector& x) {
  Vector rho = x.array().sin() / gamma_;
  return rho;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::rho_real(const real_t& x, const real_t& t) {
  return std::exp(gamma_ * t) * std::sin(x) / gamma_;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::rho_real(const Vector& x, const real_t& t) {
  Vector rho = std::exp(gamma_ * t) * x.array().sin() / gamma_;
  return rho;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::source(const real_t& x,  const real_t& t) {
  return 0.e0;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::source(const Vector& x,  const real_t& t) {
  Vector rhod = 0.e0 * x;
  return rhod;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::fsource(const real_t& x, 
              const real_t& v, 
              const real_t& t) {
  real_t S = 0.e0;
  return S;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::fsource(const Vector& x, 
              const real_t& v, 
              const real_t& t) {
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::fsource_dx(const real_t& x, 
              const real_t& t) {
  real_t S = 0.e0;
  return S;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::fsource_dx(const Vector& x, 
              const real_t& t) {
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::f_init(const real_t& x, const real_t& v) {
  real_t f = rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::f_init(const Vector& x, const real_t& v) {
  Vector f = rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::g_init(const real_t& x, const real_t& v) {
  real_t g =  v * std::cos(x);
  return g;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::g_init(const Vector& x, const real_t& v) {
  Vector g =  v * x.array().cos();
  return g;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::f_real(const real_t& x, const real_t& v, const real_t& t) {
  real_t f = rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::f_real(const Vector& x, const real_t& v, const real_t& t) {
  Vector f = rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

real_t KineticLinearDiffusion1D_SL_period_twovel::g_real(const real_t& x, const real_t& v, const real_t& t) {
  real_t g = v * (std::exp(gamma_ * t) * std::cos(x));
  return g;
};

Vector KineticLinearDiffusion1D_SL_period_twovel::g_real(const Vector& x, const real_t& v, const real_t& t) {
  Vector g = v * (std::exp(gamma_ * t) * x.array().cos());
  return g;
};

void KineticLinearDiffusion1D_SL_period_twovel::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  gamma_ = - 2.e0 / (1.e0 + std::sqrt(1 - 4 * std::pow(eps_, 2.e0)));
  iter_ = 0;
  rho_ = rho_init(X);
  rhod_ = rho_d(X);
  f_.resize(Nv);
  xback_.resize(Nv);
  vback_.resize(Nv);
  finter_.resize(Nv);
  rho_inter_.resize(Nv);
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
    rho_inter_[j] = Vector::Zero(Nx);
  };
  std::cout << " eps_ = " << eps_ << std::endl;
  std::cout << " gamma_ = " << gamma_ << std::endl;
  std::cout << " D_ = " << D_ << std::endl;
  std::cout << " theta_ = " << theta_ << std::endl;
  std::cout << " cfl_ = " << cfl_ << std::endl;
};

void KineticLinearDiffusion1D_SL_period_twovel::setdt(real_t* dt) {
  // ** 传入相关变量
  const real_t& hx = mesh1D_->gethx();
  // **
  *dt = cfl_ * hx;
};

void KineticLinearDiffusion1D_SL_period_twovel::update_rho(const real_t& Trun, 
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
  SLflux_.setZero();
  Vector fx(4), fv(4), rhotemp(4);
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  real_t xita;
  for (int j = 0; j < Nv; j++) {
    Fx_[j].setZero();
    for (int i = 0; i < Nx; i++) {
      real_t vtemp = vback_[j](i);
      real_t xtemp = xback_[j](i);
      // xtemp\in [x_{is-1},x_{is}];
      // 0 -> is-2, 1->is-1, 2->is, 3->is+1
      // std::cout << " ================== " << std::endl;
      // std::cout << "X(i) = " << X(i) << ", xtemp = " << xtemp << ",vtemp = " << vtemp << std::endl;
      if (x_order_ == 1) {
        real_t M = Maxwell(vtemp);
        if (V(j) >= 0) {
          int is = floor((xtemp - X(0)) / hx) + 1;
          int js = j;
          // xtemp\in [x_{is-1},x_{is});
          // std::cout << "  v = " << V(j) << std::endl;
          // std::cout << " (is-1) * hx = " << (is-1) * hx <<  ", xtemp =  " << xtemp << ", is * hx = " << is*hx << std::endl;
          xita = is * hx + X(0) - xtemp;
          interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
          Fx_[j](i) = (fx(1) - fx(0)) / hx - (rhotemp(3) - rhotemp(2)) / hx;
          Fx_[j](i) = Fx_[j](i) / eps_;
        } else if (V(j) < 0) {
          int is = ceil((xtemp - X(0)) / hx);
          int js = j;
          // xtemp\in (x_{is-1},x_{is}];
          // std::cout << "  v = " << V(j) << std::endl;
          // std::cout << " (is-1) * hx = " << (is-1) * hx <<  ", xtemp =  " << xtemp << ", is*hx = " << is*hx << std::endl;
          xita = is * hx + X(0) - xtemp;
          interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);
          Fx_[j](i) = (fx(3) - fx(2)) / hx - (rhotemp(1) - rhotemp(0)) / hx;
          Fx_[j](i) = Fx_[j](i) / eps_;
        }
        finter_[j](i) = xita / hx * fx(1) + (1.e0 - xita / hx) * fx(2);
        rho_inter_[j](i) = xita / hx * rhotemp(1) + (1.e0 - xita / hx) * rhotemp(2);
      } else if (x_order_ == 2) {
        QUEST_ERROR(" The second SL method has not been implemented ! ");
      }
      // std::cout << "is = " << is << std::endl;
      // std::cout << "xita = " << xita << std::endl;
      // std::cout << " is * hx = " << is * hx << ", (is-1) * hx = " << (is-1) * hx << std::endl;
    }
    SLflux_ = SLflux_ + Fx_[j] * Vweights(j) * V(j);
    // std::cout << " V(j) = " << V(j) << ", Vweights(j) = " << Vweights(j) << std::endl;
  }
  Vector bf;
  rhon_ = rho_;
  
  generateS(Trun, dt);
  generatebf(Trun, dt, &bf);
  SolveS(Trun, bf);  // update rho_
  if (!rho_.allFinite() || !bf.allFinite()) {
    std::cout << "[ERROR] rho_ contains NaN or Inf after SolveS!" << std::endl;
    std::ofstream Sout("result/LD_acctest/S.txt");
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
};

void KineticLinearDiffusion1D_SL_period_twovel::update_f(const real_t& Trun, const real_t& dt) {
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
  std::vector<Vector> ftemp = f_;
  int is, js;
  int ii1, ii2;
  real_t fv_temp;
  real_t eps2 = eps_ * eps_;
  real_t eps22 = 2 * eps2; 
  real_t ee = eps2 + dt;
  real_t ee2 = dt + eps22;
  // Vector rho_exact = rho_real(X, Trun + dt);
  for (int j = 0; j < Nv; j++) {
    // if (t_order_ == 1) {
    //   ftemp[j] = (eps2 / ee) * (finter_[j]) + (dt / ee) * rho_;
    // } else if (t_order_ == 2) {
      ftemp[j] = (eps22 / ee2) * (finter_[j]) + (dt / ee2) * rho_
               - (dt / ee2) * (finter_[j] - rho_inter_[j]);
    // }
    Source_f_[j] = fsource(X, V(j), Trun + dt);
    f_[j] = ftemp[j] + dt * (eps2 / ee) * Source_f_[j];
  };
};

KineticLinearDiffusion1D_SL_period_twovel_v2::KineticLinearDiffusion1D_SL_period_twovel_v2(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticLinearDiffusion1D_SL_period_twovel(mesh1D, poisol, x_order, t_order) {};

void KineticLinearDiffusion1D_SL_period_twovel_v2::update_f(const real_t& Trun, const real_t& dt) {
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
  std::vector<Vector> ftemp = f_;
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;
  std::cout << " The implicit finitue difference method for f update " << std::endl;
  Sv1_.resize(Nx, Nx); Sv1_.setZero();
  bfv1_.resize(Nx); bfv1_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_Sv1; // v = -1
  int estimatedNonZeros = 2 * Nx;
  tripletList_Sv1.reserve(estimatedNonZeros);

  Sv2_.resize(Nx, Nx); Sv2_.setZero();
  bfv2_.resize(Nx); bfv2_.setZero();
  std::vector<Eigen::Triplet<real_t>> tripletList_Sv2; // v = 1
  estimatedNonZeros = 2 * Nx;
  tripletList_Sv2.reserve(estimatedNonZeros);

  real_t val = 1.e0 + eps_ * dt / hx / ee;
  real_t val1 = - eps_ * dt / hx / ee;
  for (int i = 0; i < Nx; i++) {
    tripletList_Sv1.emplace_back(i, i, val);
    tripletList_Sv2.emplace_back(i, i, val);
    if (i < Nx - 1) {
      tripletList_Sv1.emplace_back(i, i + 1, val1);
    }
    if (i > 0) {
      tripletList_Sv2.emplace_back(i, i - 1, val1);
    }
  }
  // i = 0
  tripletList_Sv1.emplace_back(Nx-1, 0, val1);

  // i = Nx-1
  tripletList_Sv2.emplace_back(0, Nx-1, val1);

  Sv1_.setFromTriplets(tripletList_Sv1.begin(), tripletList_Sv1.end());
  Sv2_.setFromTriplets(tripletList_Sv2.begin(), tripletList_Sv2.end());

  Source_f_[0] = fsource(X, V(0), Trun + dt);
  bfv1_ = dt / ee * rho_ + eps2 / ee * f_[0] + eps2 * dt / ee * Source_f_[0];
  Source_f_[1] = fsource(X, V(1), Trun + dt);
  bfv2_ = dt / ee * rho_ + eps2 / ee * f_[1] + eps2 * dt / ee * Source_f_[1];

  gmresv1_.setTolerance(sparsetol_);
  gmresv1_.set_restart(30);
  gmresv1_.setMaxIterations(Nx);
  gmresv1_.compute(Sv1_);

  gmresv2_.setTolerance(sparsetol_);
  gmresv2_.set_restart(30);
  gmresv2_.setMaxIterations(Nx);
  gmresv2_.compute(Sv2_);

  f_[0] = gmresv1_.solve(bfv1_);
  if (gmresv1_.info() != Eigen::Success) {
    std::cout << "[ERROR] GMRESv1 failed to converge!\n";
  }
  // std::cout << "bfv1_ = " << bfv1_.transpose() << std::endl;
  // std::cout << f_[0].transpose() << std::endl;
  std::cout << "Iterations: " << gmresv1_.iterations() << std::endl;
  std::cout << "Estimated error: " << gmresv1_.error() << std::endl;
  f_[1] = gmresv2_.solve(bfv2_);
  if (gmresv2_.info() != Eigen::Success) {
    std::cout << "[ERROR] GMRESv2 failed to converge!\n";
  }
  std::cout << "Iterations: " << gmresv2_.iterations() << std::endl;
  std::cout << "Estimated error: " << gmresv2_.error() << std::endl;
};

} // namespace QUEST

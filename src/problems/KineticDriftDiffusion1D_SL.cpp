#include "KineticDriftDiffusion1D_SL.hpp"

namespace QUEST
{

KineticDriftDiffusion1D_SL::KineticDriftDiffusion1D_SL(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order)
  : mesh1D_(mesh1D), x_order_(x_order), poisol_(poisol), t_order_(t_order) {};

void KineticDriftDiffusion1D_SL::settheta(const real_t& theta) {
  theta_ = theta;
};

void KineticDriftDiffusion1D_SL::setcfl(const real_t& cfl) {
  cfl_ = cfl;
};

void KineticDriftDiffusion1D_SL::setpi(const real_t& pi) {
  pi_ = pi;
};

void KineticDriftDiffusion1D_SL::setgamma(const real_t& gamma) {
  gamma_ = gamma;
};

void KineticDriftDiffusion1D_SL::seteps(const real_t& eps) {
  eps_ = eps;
};

void KineticDriftDiffusion1D_SL::setsparsetol(const real_t& sparsetol) {
  sparsetol_ = sparsetol;
};
  
void KineticDriftDiffusion1D_SL::setiterationtol(const real_t& iterationtol) {
  iterationtol_ = iterationtol;
};

void KineticDriftDiffusion1D_SL::setD(const real_t& D) {
  D_ = D;
};

void KineticDriftDiffusion1D_SL::setNTH(const int& NTH) {
  NTH_ = NTH;
};

const real_t& KineticDriftDiffusion1D_SL::getcfl() const {
  return cfl_;
};

const int& KineticDriftDiffusion1D_SL::getNTH() const {
  return NTH_;
};

real_t KineticDriftDiffusion1D_SL::rho_init(const real_t& x) {
  return 1.e0;
};

Vector KineticDriftDiffusion1D_SL::rho_init(const Vector& x) {
  Vector rho = x;
  rho.setConstant(1.e0);
  return rho;
};

real_t KineticDriftDiffusion1D_SL::rho_numericalbc(const real_t& x, const real_t& t) {
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  // *** 
  real_t rho = 0.e0;
  if (std::abs(x - x1) < 1.e-10) {
    for (int j = 0; j < Nv; j++) {
      rho += Vweights(j) * fL_bc(j, t);
    }
  } else if (std::abs(x - x2) < 1.e-10) {
    for (int j = 0; j < Nv; j++) {
      rho += Vweights(j) * fR_bc(j, t);
    }
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return rho;
};

real_t KineticDriftDiffusion1D_SL::phi_bc(const real_t& x, const real_t& t) {
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

real_t KineticDriftDiffusion1D_SL::rho_d(const real_t& x) {
  return 1.e0 - (1.e0 - 0.001)/2 * (std::tanh((x - 0.3e0)/0.02e0) - std::tanh((x - 0.7e0)/0.02e0));
};

Vector KineticDriftDiffusion1D_SL::rho_d(const Vector& x) {
  Vector x1 = (x.array() - 0.3) / 0.02;
  Vector x2 = (x.array() - 0.7) / 0.02;
  Vector rhod = 1.0 - (1.0 - 0.001) / 2.0 * (x1.array().tanh() - x2.array().tanh());
  return rhod;
};

real_t KineticDriftDiffusion1D_SL::f_init(const real_t& x, const real_t& v) {
  return Maxwell(v);
};

Vector KineticDriftDiffusion1D_SL::f_init(const Vector& x, const real_t& v) {
  Vector f = x;
  f.setConstant(1.e0);
  f = f * Maxwell(v);
  return f;
};

real_t KineticDriftDiffusion1D_SL::g_init(const real_t& x, const real_t& v) {
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL::g_init(const Vector& x, const real_t& v) {
  Vector g = x;
  g.setZero();
  return g;
};

real_t KineticDriftDiffusion1D_SL::fL_bc(const real_t& v, const real_t& t) {
  QUEST_ERROR(" This fL_bc(real_t, real_t) is not valid ! ");
  return 0.e0;
};

real_t KineticDriftDiffusion1D_SL::fL_bc(const int& v_index, const real_t& t) {
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t f;
  if (V(v_index) >= 0) {
    f = Maxwell(V(v_index));
  } else {
    f = 2.e0 * f_[v_index](0) - f_[v_index](1);
  }
  return f;
};

real_t KineticDriftDiffusion1D_SL::fR_bc(const real_t& v, const real_t& t) {
  QUEST_ERROR(" This fR_bc(real_t, real_t) is not valid ! ");
  return 0.e0;
};

real_t KineticDriftDiffusion1D_SL::fR_bc(const int& v_index, const real_t& t) {
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  // *** 
  real_t f;
  if (V(v_index) <= 0) {
    f = Maxwell(V(v_index));
  } else {
    f = 2.e0 * f_[v_index](Nx-1) - f_[v_index](Nx-2);
  }
  return f;
};

real_t KineticDriftDiffusion1D_SL::source(const real_t& x,  const real_t& t) {
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL::source(const Vector& x,  const real_t& t) {
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticDriftDiffusion1D_SL::fsource(const real_t& x, 
              const real_t& v, 
              const real_t& t) {
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL::fsource(const Vector& x, 
              const real_t& v, 
              const real_t& t) {
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticDriftDiffusion1D_SL::fsource_dx(const real_t& x, 
              const real_t& t) {
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL::fsource_dx(const Vector& x, 
              const real_t& t) {
  Vector S = x;
  S.setZero();
  return S;
};

real_t KineticDriftDiffusion1D_SL::Maxwell(const real_t& v) {
  real_t M = 1.e0 / std::sqrt(2.e0 * pi_ * theta_) * std::exp( - v * v / ( 2 * theta_));
  return M;
};

const Vector& KineticDriftDiffusion1D_SL::getrho() const {
  return rho_;
};

const Vector& KineticDriftDiffusion1D_SL::getE() const {
  return E_;
};

const Vector KineticDriftDiffusion1D_SL::getrhod() const {
  return rhod_;
};

const std::vector<Vector>& KineticDriftDiffusion1D_SL::getf() const {
  return f_;
};

const int& KineticDriftDiffusion1D_SL::getiter() const {
  return iter_;
};

Vector KineticDriftDiffusion1D_SL::getrhoinit() {
  const Vector& X = mesh1D_->getX();

  Vector rho = rho_init(X);
  return rho;
}

std::vector<Vector> KineticDriftDiffusion1D_SL::getfinit() {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // ***
  std::vector<Vector> f(Nv);
  for (int j = 0; j < Nv; j++) {
    f[j] = f_init(X, V(j));
  }
  return f;
}

void KineticDriftDiffusion1D_SL::init(const Solver1DTypeFD& soltype) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
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

void KineticDriftDiffusion1D_SL::setdt(real_t* dt) {
  // ** 传入相关变量
  const real_t& hx = mesh1D_->gethx();
  // **
  *dt = cfl_ * hx;
}

void KineticDriftDiffusion1D_SL::update_E(const real_t& Trun) {
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
}

void KineticDriftDiffusion1D_SL::character_tracing(const real_t& Trun, 
                                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
#pragma omp parallel num_threads(NTH_), default(shared)
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    for (int i = 0; i < Nx; i++) {
      xback_[j](i) = X(i) - 1.e0 / eps_ * V(j) * dt;
      vback_[j](i) = V(j) - 1.e0 / eps_ * E_(i+1) * dt;
    }
  }
};

void KineticDriftDiffusion1D_SL::update_rho(const real_t& Trun, 
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
        
      }
      #pragma omp critical
      SLflux_ += Fx_[j] * Vweights(j) * V(j);
    }
  }
  }
  break;
  
  // case 2:
  // {
  // #pragma omp parallel num_threads(NTH_), default(shared)
  // {
  //   real_t vtemp, xtemp, ftemp;
  //   int is, js, ii;
  //   real_t xita, vita;
  //   real_t M;
  //   real_t dvdx;
  //   Vector fx(4), fv(4), rhotemp(4);
  // #pragma omp for schedule(static)
  //   for (int j = 0; j < Nv; j++) {
  //     Fx_[j].setZero();
  //     for (int i = 0; i < Nx; i++) {
  //       vtemp = vback_[j](i);
  //       xtemp = xback_[j](i);
  //       M = Maxwell(vtemp);
  //       // xtemp\in [x_{is-1},x_{is}]; vtemp\in [v_{js-1},v_{js}];
  //       // 0 -> is-2, 1->is-1, 2->is, 3->is+1
  //       is = floor((xtemp - X(0)) / hx) + 1;
  //       js = floor((vtemp - V(0)) / hv) + 1;
  //       xita = (is * hx + X(0) - xtemp) / hx;
  //       vita = (js * hv + V(0) - vtemp) / hv;
  //       interpolateSL_x(is, js, X(i), V(j), Trun, &fx, &fv, &rhotemp);

  //       dvdx = - dt * (rho_(i) - rhod_(i)) / (eps_ * gamma_);
  //       if (V(j) >= 0) {
  //         Fx_[j](i) = (fx(1) - fx(0)) / hx + dvdx * (fv(1) - fv(0)) / hv
  //                     - (rhotemp(3) - rhotemp(2)) / hx * M
  //                     + dvdx * vtemp / theta_ * M;
  //         Fx_[j](i) = Fx_[j](i) / eps_;
  //       } else if (V(j) < 0) {
  //         Fx_[j](i) = (fx(3) - fx(2)) / hx + dvdx * (fv(3) - fv(2)) / hv
  //                     - (rhotemp(1) - rhotemp(0)) / hx * M
  //                     + dvdx * vtemp / theta_ * M;
  //         Fx_[j](i) = Fx_[j](i) / eps_;
  //       }
  //       if (js - 1 < 0 || js - 1 >= Nv) {
  //         ftemp= 0.e0;
  //       } else {
  //         ftemp = (is - 1 < 0) ? fL_bc(js-1, Trun) : 
  //                   (is - 1 >= Nx) ? fR_bc(js-1, Trun) : f_[js-1](is - 1);
  //       }
  //       finter_[j](i) = xita * vita * ftemp + (1.e0 - xita) * vita * fv(1) 
  //                       + xita * (1.e0 - vita) * fx(1) + (1.e0 - xita) * (1.e0 - vita) * fx(2);
        
  //     }
  //     #pragma omp critical
  //     SLflux_ += Fx_[j] * Vweights(j) * V(j);
  //   }
  // }
  // }
  // break;

  default:
  {
    QUEST_ERROR(" The second SL method has not been implemented ! ");
  }
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
    update_E(Trun);
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
  
}

void KineticDriftDiffusion1D_SL::interpolateSL_x(const int& is, const int& js,
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
  if (x_order_ == 1) {
    if (js >= 0 && js < Nv) { // v方向拉出边界
      real_t fL = fL_bc(js, Trun); // x方向拉出
      real_t fR = fR_bc(js, Trun);
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
    if (js - 2 < 0 || js - 2 >= Nv) {
      (*fv)(0) = 0.e0;
    } else {
      (*fv)(0) = (is < 0) ? fL_bc(js - 2, Trun) : 
                (is >= Nx) ? fR_bc(js - 2, Trun) : f_[js-2](is);
    }
    if (js - 1 < 0 || js - 1 >= Nv) {
      (*fv)(1) = 0.e0;
    } else {
      (*fv)(1) = (is < 0) ? fL_bc(js - 1, Trun) : 
                (is >= Nx) ? fR_bc(js - 1, Trun) : f_[js-1](is);
    }
    if (js < 0 || js >= Nv) {
      (*fv)(2) = 0.e0;
    } else {
      (*fv)(2) = (is < 0) ? fL_bc(js, Trun) : 
                (is >= Nx) ? fR_bc(js, Trun) : f_[js](is);
    }
    if (js + 1 < 0 || js + 1 >= Nv) {
      (*fv)(3) = 0.e0;
    } else {
      (*fv)(3) = (is < 0) ? fL_bc(js + 1, Trun) : 
                (is >= Nx) ? fR_bc(js + 1, Trun) : f_[js+1](is);
    }
  } else if (x_order_ == 2) {
    QUEST_ERROR(" The second SL method has not been implemented ! ");
  }
}

void KineticDriftDiffusion1D_SL::update_f(const real_t& Trun, 
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
  real_t M;
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    M = Maxwell(V(j));
    // Vector rhotemp = rho_real(X, Trun + dt);
    if (t_order_ == 1) {
      ftemp[j] = (eps2 / ee) * finter_[j] + (dt / ee) * rho_ * M;
      Source_f_[j] = fsource(X, V(j), Trun + dt);
      f_[j] = ftemp[j] + dt * (eps2 / ee) * Source_f_[j];
    } 
    // else if (t_order_ == 2) {
    //   ftemp[j] = (2 * eps2 / ee2) * finter_[j] + (dt / ee2) * rho_ * M;
    //   for (int i = 0; i < Nx; i++)
    //   { 
    //     Source_f_[j] = 1.e0/ 2.e0 * fsource(X, V(j), Trun + dt) 
    //                   + 1.e0/ 2.e0 * fsource(xback_[j](i), vback_[j](i), Trun);
    //     Source_f_[j] = 1.e0/ 2.e0 * fsource(X, V(j), Trun + dt) 
    //                   + 1.e0/ 2.e0 * fsource(xback_[j](i), vback_[j](i), Trun);
        
    //   }
    // } 
    else {
      QUEST_ERROR(" The high order (>= 2) SL method has not been implemented ! ");
    }
    
  };
};
};

void KineticDriftDiffusion1D_SL::generateS(const real_t& Trun, 
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

void KineticDriftDiffusion1D_SL::SolveS(const real_t& Trun, const Vector& bf) {
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
void KineticDriftDiffusion1D_SL::generatebf(const real_t& Trun, 
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
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ * coef * rho_L_;
  (*bf)(0) = (*bf)(0) + (1.e0 - omega) * D_ / theta_ * E_(0) * dt / (2.e0 * hx) * rho_L_;

  (*bf)(Nx - 1) = (*bf)(Nx - 1) + (1.e0 - omega) * D_ * coef * rho_R_;
  (*bf)(Nx - 1) = (*bf)(Nx - 1) - (1.e0 - omega) * D_ / theta_ * E_(Nx+1) * dt / (2.e0 * hx) * rho_R_;
};

void KineticDriftDiffusion1D_SL::updateAll(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nv = mesh1D_->getNv();
  // *** 

  update_E(Trun);
  character_tracing(Trun, dt);
  update_rho(Trun, dt);
  update_f(Trun, dt);
  rho_.setZero();
  for (int j = 0; j < Nv; j++) {
    rho_ = rho_ + Vweights(j) * f_[j];
  }
}

// ** 周期边界系列 ** // 
KineticDriftDiffusion1D_SL_period::KineticDriftDiffusion1D_SL_period(const KineticFDmesh_period* mesh1D,
                        Poissonsolver1DFD_period* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order) {};

real_t KineticDriftDiffusion1D_SL_period::rho_init(const real_t& x) {
  return std::sin(x);
};

Vector KineticDriftDiffusion1D_SL_period::rho_init(const Vector& x) {
  Vector rho = x.array().sin();
  return rho;
};

real_t KineticDriftDiffusion1D_SL_period::rho_real(const real_t& x, const real_t& t) {
  return std::exp(- D_ * t) * std::sin(x);
};

Vector KineticDriftDiffusion1D_SL_period::rho_real(const Vector& x, const real_t& t) {
  Vector rho = std::exp(- D_ * t) * x.array().sin();
  return rho;
};

Vector KineticDriftDiffusion1D_SL_period::E_real(const Vector& x, const real_t& t) {
  Vector E = - std::exp(- D_ * t) * x.array().cos();
  return E;
};

Vector KineticDriftDiffusion1D_SL_period::getEfinal(const real_t& t) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  // *** 
  Vector E = E_real(X, t);
  return E;
};

std::vector<Vector> KineticDriftDiffusion1D_SL_period::getffinal(const real_t& t) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nv = mesh1D_->getNv();
  const int& Nx = mesh1D_->getNx();
  // *** 
  std::vector<Vector> f(Nv);
  for (int j = 0; j < Nv; j++) {
    f[j] = f_real(X, V(j), t);
  }
  return f;
};

void KineticDriftDiffusion1D_SL_period::setrhofinal(const real_t& t) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  // *** 
  rho_final_ = rho_real(X, t);
};

const Vector& KineticDriftDiffusion1D_SL_period::getrhofinal() const {
  return rho_final_;
};

real_t KineticDriftDiffusion1D_SL_period::rho_d(const real_t& x) {
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_period::rho_d(const Vector& x) {
  Vector rhod = 0.e0 * x;
  return rhod;
};

real_t KineticDriftDiffusion1D_SL_period::source(const real_t& x,  const real_t& t) {
  return - std::exp(- 2.e0 * D_  * t) * std::cos(2.e0 * x);
};

Vector KineticDriftDiffusion1D_SL_period::source(const Vector& x,  const real_t& t) {
  Vector rhod = 2.e0 * x;
  rhod = - std::exp(- 2.e0 * D_ * t) * rhod.array().cos();
  return rhod;
};

real_t KineticDriftDiffusion1D_SL_period::fsource(const real_t& x, 
              const real_t& v, 
              const real_t& t) {
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  real_t S;
  S = std::exp(- D_ * t) * M * (eps_ * v * D_ * std::cos(x) - 
                                D_ * sin(x) + v2 * sin(x)) +
      std::exp(- 2 * D_ * t) * M * (eps_ * v * D_ * std::sin(2*x)/theta_ 
                                - v2 / theta_ * std::cos(2*x) + 
                                std::cos(x) * std::cos(x) * (1.e0 - v2 / theta_)) +
      std::exp(- 3 * D_ * t) * M * (std::cos(x) * std::cos(x)) * std::sin(x) * 
                            (1.e0 / theta_  - v2 / (theta_ * theta_));
  return S;
};

Vector KineticDriftDiffusion1D_SL_period::fsource(const Vector& x, 
              const real_t& v, 
              const real_t& t) {
  real_t M = Maxwell(v);
  real_t v2 = v * v;
  Vector S;
  Vector sinx = x.array().sin();
  Vector cosx = x.array().cos();
  Vector sin2x = (2.e0 * x).array().sin();
  Vector cos2x = (2.e0 * x).array().cos();
  S = std::exp(- D_ * t) * M * (eps_ * v * D_ * cosx.array()
                                - D_ * sinx.array() + v2 * sinx.array()) +
      std::exp(- 2 * D_ * t) * M * (eps_ * v * D_ * sin2x.array()/theta_ 
                                - 1.e0 / theta_ * v2 * cos2x.array() + 
                                cosx.array() * cosx.array() * (1.e0 - v2 / theta_)) +
      std::exp(- 3 * D_ * t) * M * (cosx.array() * cosx.array()) * sinx.array() * 
                            (1.e0 / theta_  - v2 / (theta_ * theta_));
  return S;
};

real_t KineticDriftDiffusion1D_SL_period::fsource_dx(const real_t& x, 
              const real_t& t) {
  // <vS> = eps_ * D_ * (D_ * std::cos(x) * std::exp(- D_ * t)
  //     + (D_) / theta_ * std::exp(- 2.e0 * D_ * t) * std::sin(2*x))
  real_t S = eps_ * D_ * ( - D_ * std::sin(x) * std::exp(- D_ * t)
      + (2.e0 * D_) / theta_ * std::exp(- 2.e0 * D_ * t) * std::cos(2.e0 * x));
  return S;
};

Vector KineticDriftDiffusion1D_SL_period::fsource_dx(const Vector& x, 
              const real_t& t) {
  Vector S = eps_ * D_ * ( - D_ * x.array().sin() * std::exp(- D_ * t)
      + (2.e0 * D_) / theta_ * std::exp(- 2.e0 * D_ * t) * (2.e0 * x).array().cos());
  return S;
};


real_t KineticDriftDiffusion1D_SL_period::f_init(const real_t& x, const real_t& v) {
  real_t f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Vector KineticDriftDiffusion1D_SL_period::f_init(const Vector& x, const real_t& v) {
  Vector f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t KineticDriftDiffusion1D_SL_period::g_init(const real_t& x, const real_t& v) {
  real_t M = Maxwell(v);
  real_t g = - v * M * (std::cos(x) + 1.e0 / theta_ * std::cos(x) * std::sin(x));
  return g;
};

Vector KineticDriftDiffusion1D_SL_period::g_init(const Vector& x, const real_t& v) {
  real_t M = Maxwell(v);
  Vector g =  - v * M * (x.array().cos() + 1.e0 / theta_ * x.array().cos() * x.array().sin());
  return g;
};

real_t KineticDriftDiffusion1D_SL_period::f_real(const real_t& x, const real_t& v, const real_t& t) {
  real_t f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

Vector KineticDriftDiffusion1D_SL_period::f_real(const Vector& x, const real_t& v, const real_t& t) {
  Vector f = Maxwell(v) * rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

real_t KineticDriftDiffusion1D_SL_period::g_real(const real_t& x, const real_t& v, const real_t& t) {
  real_t M = Maxwell(v);
  real_t g =  - v * M * (std::exp(- D_ * t) * std::cos(x) + 
            std::exp(- 2.e0 * D_ * t) / theta_ * std::cos(x) * std::sin(x));
  return g;
};

Vector KineticDriftDiffusion1D_SL_period::g_real(const Vector& x, const real_t& v, const real_t& t) {
  real_t M = Maxwell(v);
  Vector g =  - v * M * (std::exp(- D_ * t) * x.array().cos() + 
            std::exp(- 2.e0 * D_ * t) / theta_ * x.array().cos() * x.array().sin());
  return g;
};

void KineticDriftDiffusion1D_SL_period::character_tracing(const real_t& Trun, 
                                                const real_t& dt) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  const int& Nv = mesh1D_->getNv();
  // *** 
  // E_ = E_real(X, Trun + dt);
#pragma omp parallel num_threads(NTH_), default(shared)
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    for (int i = 0; i < Nx; i++) {
      xback_[j](i) = X(i) - 1.e0 / eps_ * V(j) * dt;
      vback_[j](i) = V(j) - 1.e0 / eps_ * E_(i) * dt;
    }
  }
};

void KineticDriftDiffusion1D_SL_period::generateS(const real_t& Trun, 
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
  int estimatedNonZeros = 5 * Nx;
  tripletList_S.reserve(estimatedNonZeros);

  const real_t coef = dt / (hx * hx);
  real_t val;
  for (int i = 0; i < Nx; i++) {
    val = 1.e0 + (1.e0 - omega) * D_ * 2.e0 * coef;
    tripletList_S.emplace_back(i, i, val);
    if (i > 0) {
      tripletList_S.emplace_back(i, i - 1, - (1.e0 - omega) * D_  * coef);
      val = - (1.e0 - omega) * D_ / theta_ * E_(i - 1) * dt / (2.e0 * hx);
      tripletList_S.emplace_back(i, i - 1, val);
    }
    if (i < Nx - 1) {
      tripletList_S.emplace_back(i, i + 1, - (1.e0 - omega) * D_  * coef);
      val = (1.e0 - omega) * D_ / theta_ * E_(i + 1) * dt / (2.e0 * hx);
      tripletList_S.emplace_back(i, i + 1, val);
    } 
  }
  // i = 0
  tripletList_S.emplace_back(0, Nx-1, - (1.e0 - omega) * D_  * coef);
  val = - (1.e0 - omega) * D_ / theta_ * E_(Nx-1) * dt / (2.e0 * hx);
  tripletList_S.emplace_back(0, Nx-1, val);
  // i = Nx-1
  tripletList_S.emplace_back(Nx-1, 0, - (1.e0 - omega) * D_  * coef);
  val = (1.e0 - omega) * D_ / theta_ * E_(0) * dt / (2.e0 * hx);
  tripletList_S.emplace_back(Nx-1, 0, val);

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
    gmres_.set_restart(100);
    gmres_.setMaxIterations(Nx);
    gmres_.compute(S_);
    break;
  default:
    QUEST_ERROR(" The Poisson Solver is not been implemented ! ");
    break;
  }
};

void KineticDriftDiffusion1D_SL_period::generatebf(const real_t& Trun, 
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
  // *bf = rhon_ + Source_;
  *bf = rhon_ - omega * dt * SLflux_ + Source_ + Source_v_;
};

void KineticDriftDiffusion1D_SL_period::SolveS(const real_t& Trun, const Vector& bf) {
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

void KineticDriftDiffusion1D_SL_period::update_E(const real_t& Trun) {
  // ** 传入相关变量
  const Vector& X = mesh1D_->getX();
  const int& Nx = mesh1D_->getNx();
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  Vector bf = - (rho_ - rhod_) / gamma_;
  poisol_->SolveAll(0.e0, 0.e0, bf, &phi_, &E_);
};

void KineticDriftDiffusion1D_SL_period::update_rho(const real_t& Trun, 
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
  real_t area = hx * hv;

  Source_v_.resize(Nx); Source_v_.setZero();
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);

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
      interpolateSL_x(is, js, xtemp, vtemp, Trun, &fx, &fv, &rhotemp);
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
                        + xita * (1.e0 - vita) * fx(1) + (1.e0 - xita) * (1.e0 - vita) * fx(2);
      } else if (x_order_ == 2) {
        QUEST_ERROR(" The second SL method has not been implemented ! ");
      }
    }
    #pragma omp critical
    SLflux_ += Fx_[j] * Vweights(j) * V(j);
  }
}

  real_t err = 1.e0;
  Vector bf;
  rhon_ = rho_;
  iter_ = 0;
  while (err > iterationtol_) {
    update_E(Trun);
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
    std::cout <<" size = " << rho_.size() << ", err = " << err << ", iter = " << iter_ << std::endl;
    if (iter_ > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort." << std::endl;
      break;
    }
  }
  std::cout << "err = " << err << std::endl;
};

void KineticDriftDiffusion1D_SL_period::update_f(const real_t& Trun, const real_t& dt) {
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
#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t M;
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    M = Maxwell(V(j));
    // Vector rhotemp = rho_real(X, Trun + dt);
    if (x_order_ == 1) {
      ftemp[j] = (eps2 / ee) * finter_[j] + (dt / ee) * rho_ * M;
    } else if (x_order_ == 2) {
      QUEST_ERROR(" The second SL method has not been implemented ! ");
    }
    Source_f_[j] = fsource(X, V(j), Trun + dt);
    f_[j] = ftemp[j] + dt * (eps2 / ee) * Source_f_[j];
  };
}
  
};

void KineticDriftDiffusion1D_SL_period::interpolateSL_x(const int& is, const int& js,
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
  int ii;
  if (x_order_ == 1) {
    if (js <= 0 || js >= Nv) { // v方向拉出边界
      fx->setZero(); fv->setZero();
      rhotemp->setZero();
    } else {
      ii = (((is - 2) % Nx) + Nx) % Nx;
      (*fx)(0) = f_[js](ii);
      (*rhotemp)(0) = rho_(ii);
      ii = (((is - 1) % Nx) + Nx) % Nx;
      (*fx)(1) = f_[js](ii);
      (*rhotemp)(1) = rho_(ii);
      ii = ((is % Nx) + Nx) % Nx;
      (*fx)(2) = f_[js](ii);
      (*rhotemp)(2) = rho_(ii);
      (*fv)(0) = (js - 2 >= 0 && js - 2 < Nv) ? f_[js - 2](ii) : 0.0;
      (*fv)(1) = (js - 1 >= 0 && js - 1 < Nv) ? f_[js - 1](ii) : 0.0;
      (*fv)(2) = (js >= 0 && js < Nv) ? f_[js](ii) : 0.0;
      (*fv)(3) = (js + 1 >= 0 && js + 1 < Nv) ? f_[js + 1](ii) : 0.0;
      ii = (((is + 1) % Nx) + Nx) % Nx;
      (*fx)(3) = f_[js](ii);
      (*rhotemp)(3) = rho_(ii);
    }
  } else if (x_order_ == 2) {
    QUEST_ERROR(" The second SL method has not been implemented ! ");
  };
};

void KineticDriftDiffusion1D_SL_period::updateAll(const real_t& Trun, const real_t& dt) {
  // ** 传入相关变量
  const Vector& Vweights = mesh1D_->getvweights();
  const int& Nv = mesh1D_->getNv();
  // *** 

  update_E(Trun);
  character_tracing(Trun, dt);
  update_rho(Trun, dt);
  update_f(Trun, dt);

  rho_.setZero();
  for (int j = 0; j < Nv; j++) {
    rho_ = rho_ + Vweights(j) * f_[j];
  }
};


// 对不同gamma算例的精度测试
KineticDriftDiffusion1D_SL_period_Different_Gamma::KineticDriftDiffusion1D_SL_period_Different_Gamma(const KineticFDmesh_period* mesh1D,
                          Poissonsolver1DFD_period* poisol,
                          const int& x_order,
                          const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order) {};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_init(const real_t& x) 
{ 
  // return gamma_ * std::sin(x);
  return rho_d(x) + gamma_ * std::sin(x);
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_init(const Vector& x) 
{
  // return gamma_ * x.array().sin();
  Vector rho = rho_d(x).array() + gamma_ * x.array().sin();
  return rho;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_d(const real_t& x) 
{
  // return 1.e0;
  return 1.e0;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_d(const Vector& x) 
{
  Vector rhod = x;
  rhod.setConstant(1.e0);
  return rhod;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::source(const real_t& x, const real_t& t) 
{
  real_t S = - std::exp(- 2.e0 * D_  * t) * std::cos(2.e0 * x);
  S = S * gamma_;
  S = S + std::exp(- D_ * t) * std::sin(x);
  return S;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::source(const Vector& x, const real_t& t) 
{
  Vector S = 2.e0 * x;
  S = - std::exp(- 2.e0 * D_ * t) * S.array().cos() * gamma_;
  S = S.array() + std::exp(- D_ * t) * x.array().sin();
  return S;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource(const real_t& x, const real_t& v, const real_t& t) 
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

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource(const Vector& x, const real_t& v, const real_t& t) 
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

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource_dx(const real_t& x, const real_t& t) 
{
  // <vS> = eps_ * D_ * (D_ * std::cos(x) * std::exp(- D_ * t)
  //     + (D_) / theta_ * std::exp(- 2.e0 * D_ * t) * std::sin(2*x))
  // <vS> = <vS> * gamma_ + eps_ * D_ * std::exp(- D_ * t) * std::cos(x)
  real_t S = eps_ * D_ * ( - D_ * std::sin(x) * std::exp(- D_ * t)
      + (2.e0 * D_) / theta_ * std::exp(- 2.e0 * D_ * t) * std::cos(2.e0 * x));
  S = S * gamma_ - eps_ * D_ * std::exp(- D_ * t) * std::sin(x);
  return S;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::fsource_dx(const Vector& x, const real_t& t) 
{
  Vector S = eps_ * D_ * ( - D_ * x.array().sin() * std::exp(- D_ * t)
      + (2.e0 * D_) / theta_ * std::exp(- 2.e0 * D_ * t) * (2.e0 * x).array().cos());
  S = S.array() * gamma_ - eps_ * D_ * std::exp(- D_ * t) * x.array().sin();
  return S;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::f_init(const real_t& x, const real_t& v) 
{
  real_t f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::f_init(const Vector& x, const real_t& v) 
{
  Vector f = Maxwell(v) * rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::g_init(const real_t& x, const real_t& v) 
{
  real_t M = Maxwell(v);
  real_t g = - v * M * (std::cos(x) + 1.e0 / theta_ * std::cos(x) * std::sin(x));
  g = g * gamma_ - v * M * std::cos(x) / theta_;
  return g;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::g_init(const Vector& x, const real_t& v) 
{
  real_t M = Maxwell(v);
  Vector g = - v * M * (x.array().cos() + 1.e0 / theta_ * x.array().cos() * x.array().sin());
  g = g.array() * gamma_ - v * M * x.array().cos() / theta_;
  return g;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_real(const real_t& x, const real_t& t)
{
  return rho_d(x) + gamma_ * std::exp(- D_ * t) * std::sin(x);
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::rho_real(const Vector& x, const real_t& t)
{
  Vector rho = rho_d(x).array() + gamma_ * std::exp(- D_ * t) * x.array().sin();
  return rho;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::E_real(const Vector& x, const real_t& t) 
{
  Vector E = - std::exp(- D_ * t) * x.array().cos();
  return E;
};

real_t KineticDriftDiffusion1D_SL_period_Different_Gamma::g_real(const real_t& x, const real_t& v, const real_t& t) 
{
  real_t M = Maxwell(v);
  real_t g =  - v * M * (std::exp(- D_ * t) * std::cos(x) + 
            std::exp(- 2.e0 * D_ * t) / theta_ * std::cos(x) * std::sin(x));
  g = g * gamma_ - std::exp(- D_ * t) * v * M * std::cos(x) / theta_;
  return g;
};

Vector KineticDriftDiffusion1D_SL_period_Different_Gamma::g_real(const Vector& x, const real_t& v, const real_t& t)
{
  real_t M = Maxwell(v);
  Vector g =  - v * M * (std::exp(- D_ * t) * x.array().cos() + 
            std::exp(- 2.e0 * D_ * t) / theta_ * x.array().cos() * x.array().sin());
  g = g.array() * gamma_ - std::exp(- D_ * t) * v * M * x.array().cos() / theta_;
  return g;
};

// KineticDriftDiffusion1D_FVSL::KineticDriftDiffusion1D_FVSL(const KineticFVmesh* mesh1D,
//                         Poissonsolver1DFD* poisol,
//                         const int& x_order,
//                         const int& t_order)
//   : mesh1D_(mesh1D), poisol_(poisol), x_order_(x_order), t_order_(t_order) {};

// void KineticDriftDiffusion1D_FVSL::init(const Solver1DTypeFD& soltype) 
// {
//   // ** 传入相关变量
//   const Vector& X = mesh1D_->getX();
//   const Vector& V = mesh1D_->getV();
//   const int& Nx = mesh1D_->getNx();
//   const int& Nv = mesh1D_->getNv();
//   const real_t& hx = mesh1D->gethx();
//   const real_t& hv = mesh1D->gethv();
//   const int& qua_num = mesh1D->getquanum();
//   // *** 
//   // *** 
//   pi_ = 3.14159265358979323846264338327;
//   iter_ = 0;
//   rho_ = rho_init(X);
//   rhod_ = rho_d(X);
//   f_.resize(Nv);
//   xback_.resize(Nv);
//   vback_.resize(Nv);
//   finter_.resize(Nv);
//   Fx_.resize(Nv);
//   Source_f_.resize(Nv);
//   SLflux_ = Vector::Zero(Nx);
//   Source_v_ = Vector::Zero(Nx);
//   soltype_ = soltype;
//   for (int j = 0; j < Nv; j++) {
//     f_[j] = f_init(X, V(j));
//     xback_[j] = Vector::Zero(Nx);
//     vback_[j] = Vector::Zero(Nx);
//     Fx_[j] = Vector::Zero(Nx);
//     finter_[j] = Vector::Zero(Nx);
//   };
//   std::cout << " eps_ = " << eps_ << std::endl;
//   std::cout << " gamma_ = " << gamma_ << std::endl;
//   std::cout << " D_ = " << D_ << std::endl;
//   std::cout << " theta_ = " << theta_ << std::endl;
//   std::cout << " cfl_ = " << cfl_ << std::endl;
// };


KineticDriftDiffusion1D_SL_PNjunction::KineticDriftDiffusion1D_SL_PNjunction(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D_, poisol, x_order, t_order) {};

real_t KineticDriftDiffusion1D_SL_PNjunction::rho_init(const real_t& x)
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector KineticDriftDiffusion1D_SL_PNjunction::rho_init(const Vector& x)
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

real_t KineticDriftDiffusion1D_SL_PNjunction::rho_numericalbc(const real_t& x, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::rho_numericalbc(x, t);
};

real_t KineticDriftDiffusion1D_SL_PNjunction::phi_bc(const real_t& x, const real_t & t)
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

real_t KineticDriftDiffusion1D_SL_PNjunction::
  rho_d(const real_t& x)
{
  return (x < 0.5e0) ? 0.9e0 : 0.1e0;
};

Vector KineticDriftDiffusion1D_SL_PNjunction::
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

real_t KineticDriftDiffusion1D_SL_PNjunction::
  source(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_PNjunction::
  source(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_PNjunction::
  fsource(const real_t& x, const real_t& v, const real_t& t)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_PNjunction::
  fsource(const Vector& x, const real_t& v, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_PNjunction::
  fsource_dx(const real_t& x, const real_t& t)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_PNjunction::
  fsource_dx(const Vector& x, const real_t& t)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_PNjunction::
  f_init(const real_t& x, const real_t& v)
{
  return rho_init(x) * Maxwell(v) + eps_ * g_init(x, v);
};

Vector KineticDriftDiffusion1D_SL_PNjunction::f_init(const Vector& x, const real_t& v)
{
  return rho_init(x) * Maxwell(v) + eps_ * g_init(x, v);
};

real_t KineticDriftDiffusion1D_SL_PNjunction::g_init(const real_t& x, const real_t& v)
{
  return 0.e0;
};

Vector KineticDriftDiffusion1D_SL_PNjunction::g_init(const Vector& x, const real_t& v)
{
  return Vector::Zero(x.size());
};

real_t KineticDriftDiffusion1D_SL_PNjunction::fL_bc(const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::fL_bc(v, t);
};

real_t KineticDriftDiffusion1D_SL_PNjunction::fL_bc(const int& v_index, const real_t& t)
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

real_t KineticDriftDiffusion1D_SL_PNjunction::fR_bc(const real_t& v, const real_t& t)
{
  return KineticDriftDiffusion1D_SL::fR_bc(v, t);
};

real_t KineticDriftDiffusion1D_SL_PNjunction::fR_bc(const int& v_index, const real_t& t)
{
  // ** 传入相关变量
  const Vector& V = mesh1D_->getV();
  const int& Nx = mesh1D_->getNx();
  // *** 
  real_t f;
  if (V(v_index) <= 0) 
  {
    f = 0.1e0 * Maxwell(V(v_index));
  } else 
  {
    f = 2.e0 * f_[v_index](Nx-1) - f_[v_index](Nx-2);
  };
  return f;
};

void KineticDriftDiffusion1D_SL_PNjunction::init(const Solver1DTypeFD& soltype) {
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
  QUEST_VERIFY(std::abs(x1 - 0.e0) < 1.e-13, " x1 has to be equal to 0.e0 !");
  QUEST_VERIFY(std::abs(x2 - 1.e0) < 1.e-13, " x2 has to be equal to 1.e0 !");

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

// 10-27所做的一些小尝试
KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative::
  KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative(
                          const KineticFDmesh_period* mesh1D,
                          Poissonsolver1DFD_period* poisol,
                          const int& x_order,
                          const int& t_order)
  : KineticDriftDiffusion1D_SL(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period(mesh1D, poisol, x_order, t_order),
    KineticDriftDiffusion1D_SL_period_Different_Gamma(mesh1D, poisol, x_order, t_order) {};

void KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative::
                      update_rho(const real_t& Trun, 
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

  Source_v_.resize(Nx); Source_v_.setZero();
  real_t omega = std::exp( - dt / (eps_ * eps_));
  real_t coef = dt / (hx * hx);
  Source_ = source(X, Trun + dt) * dt;
  Source_v_ = (1.e0 - omega) * eps_ * (- 1.e0) * dt * fsource_dx(X, Trun + dt);

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
      // std::cout << "M = " << M << std::endl;
      // xtemp\in [x_{is-1},x_{is}]; vtemp\in [v_{js-1},v_{js}];
      // fx:: 0 -> is-2, 1->is-1, 2->is, 3->is+1, js
      // fv:: 0 -> js-2, ... , is
      is = floor((xtemp - X(0)) / hx) + 1;
      js = floor((vtemp - V(0)) / hv) + 1;
      xita = (is * hx + X(0) - xtemp) / hx;
      vita = (js * hv + V(0) - vtemp) / hv;
      
      interpolateSL_x(is, js, xtemp, vtemp, Trun, &fx, &fv, &rhotemp);
      if (x_order_ == 1) {
        // dvdx = - dt * (rho_(i) - rhod_(i)) / (eps_ * gamma_);
        dvdx = 0.e0;
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
        // std::cout << " fx(2) - fv(2) = " << fx(2) - fv(2) << std::endl;
        // if (j == 10)
        // {
        //   std::cout << " i = " << i << " j = " << j 
        //             << " is = " << is << " js = " << js 
        //             << " xita = " << xita
        //             << " vita = " << vita << std::endl;
        // };
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
  real_t temp1 = 0.e0, temp2 = 0.e0;
  temp1 = SLflux_.sum() * hx;
  temp2 = Source_.sum() * hx;
  std::cout << " temp1 = " << temp1 << std::endl;
  std::cout << " temp2 = " << temp2 << std::endl;

  real_t err = 1.e0;
  Vector bf;
  rhon_ = rho_;
  iter_ = 0;
  while (err > iterationtol_) {
    update_E(Trun);
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
    std::cout <<" size = " << rho_.size() << ", err = " << err << ", iter = " << iter_ << std::endl;
    if (iter_ > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort." << std::endl;
      break;
    }
  }
  std::cout << "err = " << err << std::endl;
}

void KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative::
  update_f(const real_t& Trun, const real_t& dt)
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
  std::vector<Vector> ftemp = f_;
  real_t eps2 = eps_ * eps_;
  real_t ee = eps2 + dt;

  real_t area = hx * hv;
  real_t temp3i = 0.e0, temp3 = 0.e0, temp4 = 0.e0;
  for(int j = 0; j < Nv; j++)
  { 
    temp3i = temp3i + f_[j].sum() * area;
    temp3 = temp3 + finter_[j].sum() * area;
    temp4 = temp4 + Source_f_[j].sum() * area;
  };
  std::cout << " temp3i = " << temp3i << std::endl;
  std::cout << " temp3 = " << temp3 << std::endl;
  std::cout << " temp4 = " << temp4 << std::endl;

#pragma omp parallel num_threads(NTH_), default(shared)
{
  real_t M;
#pragma omp for schedule(static)
  for (int j = 0; j < Nv; j++) {
    M = Maxwell(V(j));
    // Vector rhotemp = rho_real(X, Trun + dt);
    if (x_order_ == 1) {
      ftemp[j] = (eps2 / ee) * finter_[j] + (dt / ee) * rho_ * M;
    } else if (x_order_ == 2) {
      QUEST_ERROR(" The second SL method has not been implemented ! ");
    }
    Source_f_[j] = fsource(X, V(j), Trun + dt);
    f_[j] = ftemp[j] + dt * (eps2 / ee) * Source_f_[j];
  };
}
  
};

} // namespace QUEST

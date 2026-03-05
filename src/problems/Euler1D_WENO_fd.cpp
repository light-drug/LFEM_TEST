#include "Euler1D_WENO_fd.hpp"

#include <cmath>

#include "error.hpp"

namespace QUEST
{
namespace
{
constexpr real_t kWenoEps = 1.e-6;
}

Euler1D_WENO_FD::Euler1D_WENO_FD(const FDmesh* mesh1D,
                                 const EX_TVDRK* rk_table,
                                 const int x_order)
  : mesh1D_(mesh1D), rk_table_(rk_table), gamma_(1.4e0), x_order_(x_order),
    pi_(3.141592653589793238462643383279), CFL_(0.4e0)
{
  if (x_order_ != 3 && x_order_ != 5) {
    QUEST_ERROR("Euler1D_WENO_FD only supports WENO3 or WENO5 in space.");
  }
}

void Euler1D_WENO_FD::setgamma(const real_t& gamma)
{
  QUEST_VERIFY(gamma > 1.e0, "gamma should be > 1 for Euler equation.");
  gamma_ = gamma;
}

void Euler1D_WENO_FD::setCFL(const real_t& cfl)
{
  CFL_ = cfl;
}

const std::vector<Vector>& Euler1D_WENO_FD::getu() const { return u_; }
const Vector& Euler1D_WENO_FD::getrho() const { return u_[0]; }
const Vector& Euler1D_WENO_FD::getrhovx() const { return u_[1]; }
const Vector& Euler1D_WENO_FD::getE() const { return u_[2]; }

Vector Euler1D_WENO_FD::getpre(const std::vector<Vector>& u) const
{
  const Vector vx = u[1].array() / u[0].array();
  return ((u[2].array() - 0.5e0 * u[0].array() * vx.array().square()) * (gamma_ - 1.e0)).matrix();
}

void Euler1D_WENO_FD::init()
{
  const Vector& xc = mesh1D_->getCellCenter_vec();
  u_.resize(3);
  u_[0] = rho_init(xc);
  const Vector vx = vx_init(xc);
  const Vector pre = pre_init(xc);
  u_[1] = (u_[0].array() * vx.array()).matrix();
  u_[2] = (pre.array() / (gamma_ - 1.e0) + 0.5e0 * u_[0].array() * vx.array().square()).matrix();

  const int stages = rk_table_->getstages();
  u_stages_.assign(stages, std::vector<Vector>(3, Vector::Zero(u_[0].size())));
  Lu_stages_.assign(stages, std::vector<Vector>(3, Vector::Zero(u_[0].size())));
}

void Euler1D_WENO_FD::setdt(real_t* dt) const
{
  const Vector vx = u_[1].array() / u_[0].array();
  const Vector pre = getpre(u_);
  const Vector cs = (gamma_ * pre.array() / u_[0].array()).sqrt();
  const real_t max_speed = (vx.array().abs() + cs.array()).maxCoeff();
  *dt = CFL_ * mesh1D_->gethx() / max_speed;
}

void Euler1D_WENO_FD::extend_left_ghost(const int ghost_id, const real_t Trun,
                                        std::vector<Vector>* ue) const
{
  const real_t xg = mesh1D_->getx1() - (real_t(ghost_id) + 0.5e0) * mesh1D_->gethx();
  const int idx = 2 - ghost_id;
  const real_t rho = rho_bc(xg, Trun);
  const real_t vx = vx_bc(xg, Trun);
  const real_t pre = pre_bc(xg, Trun);
  ue->at(0)(idx) = rho;
  ue->at(1)(idx) = rho * vx;
  ue->at(2)(idx) = pre / (gamma_ - 1.e0) + 0.5e0 * rho * vx * vx;
}

void Euler1D_WENO_FD::extend_right_ghost(const int ghost_id, const real_t,
                                         std::vector<Vector>* ue) const
{
  const int N = static_cast<int>(u_[0].size());
  const int idx = 3 + N + ghost_id;
  ue->at(0)(idx) = u_[0](N - 1 - ghost_id);
  ue->at(1)(idx) = u_[1](N - 1 - ghost_id);
  ue->at(2)(idx) = u_[2](N - 1 - ghost_id);
}

std::vector<Vector> Euler1D_WENO_FD::extend_with_ghost(const std::vector<Vector>& u, const real_t& Trun) const
{
  const int ng = 3;
  const int N = static_cast<int>(u[0].size());
  std::vector<Vector> ue(3, Vector::Zero(N + 2 * ng));
  for (int k = 0; k < 3; ++k) {
    ue[k].segment(ng, N) = u[k];
  }
  for (int g = 0; g < ng; ++g) {
    extend_left_ghost(g, Trun, &ue);
    const int idx = ng + N + g;
    for (int k = 0; k < 3; ++k) {
      ue[k](idx) = u[k](N - 1 - g);
    }
  }
  return ue;
}

real_t Euler1D_WENO_FD::weno3_left_biased(const Vector& ue, const int iface) const
{
  const real_t um1 = ue(iface - 1);
  const real_t u0 = ue(iface);
  const real_t u1 = ue(iface + 1);

  const real_t p0 = 0.5e0 * (-um1 + 3.e0 * u0);
  const real_t p1 = 0.5e0 * (u0 + u1);

  const real_t beta0 = (u0 - um1) * (u0 - um1);
  const real_t beta1 = (u1 - u0) * (u1 - u0);

  const real_t a0 = (1.e0 / 3.e0) / ((kWenoEps + beta0) * (kWenoEps + beta0));
  const real_t a1 = (2.e0 / 3.e0) / ((kWenoEps + beta1) * (kWenoEps + beta1));
  return (a0 * p0 + a1 * p1) / (a0 + a1);
}

real_t Euler1D_WENO_FD::weno5_left_biased(const Vector& ue, const int iface) const
{
  const real_t um2 = ue(iface - 2);
  const real_t um1 = ue(iface - 1);
  const real_t u0 = ue(iface);
  const real_t u1 = ue(iface + 1);
  const real_t u2 = ue(iface + 2);

  const real_t p0 = (2.e0 * um2 - 7.e0 * um1 + 11.e0 * u0) / 6.e0;
  const real_t p1 = (-um1 + 5.e0 * u0 + 2.e0 * u1) / 6.e0;
  const real_t p2 = (2.e0 * u0 + 5.e0 * u1 - u2) / 6.e0;

  const real_t beta0 = (13.e0 / 12.e0) * std::pow(um2 - 2.e0 * um1 + u0, 2) +
                       0.25e0 * std::pow(um2 - 4.e0 * um1 + 3.e0 * u0, 2);
  const real_t beta1 = (13.e0 / 12.e0) * std::pow(um1 - 2.e0 * u0 + u1, 2) +
                       0.25e0 * std::pow(um1 - u1, 2);
  const real_t beta2 = (13.e0 / 12.e0) * std::pow(u0 - 2.e0 * u1 + u2, 2) +
                       0.25e0 * std::pow(3.e0 * u0 - 4.e0 * u1 + u2, 2);

  const real_t a0 = 0.1e0 / ((kWenoEps + beta0) * (kWenoEps + beta0));
  const real_t a1 = 0.6e0 / ((kWenoEps + beta1) * (kWenoEps + beta1));
  const real_t a2 = 0.3e0 / ((kWenoEps + beta2) * (kWenoEps + beta2));
  return (a0 * p0 + a1 * p1 + a2 * p2) / (a0 + a1 + a2);
}

real_t Euler1D_WENO_FD::weno3_right_biased(const Vector& ue, const int iface) const
{
  const real_t u0 = ue(iface);
  const real_t u1 = ue(iface + 1);
  const real_t u2 = ue(iface + 2);

  const real_t p0 = 0.5e0 * (-u2 + 3.e0 * u1);
  const real_t p1 = 0.5e0 * (u0 + u1);

  const real_t beta0 = (u2 - u1) * (u2 - u1);
  const real_t beta1 = (u1 - u0) * (u1 - u0);

  const real_t a0 = (1.e0 / 3.e0) / ((kWenoEps + beta0) * (kWenoEps + beta0));
  const real_t a1 = (2.e0 / 3.e0) / ((kWenoEps + beta1) * (kWenoEps + beta1));
  return (a0 * p0 + a1 * p1) / (a0 + a1);
}

real_t Euler1D_WENO_FD::weno5_right_biased(const Vector& ue, const int iface) const
{
  const real_t um1 = ue(iface - 1);
  const real_t u0 = ue(iface);
  const real_t u1 = ue(iface + 1);
  const real_t u2 = ue(iface + 2);
  const real_t u3 = ue(iface + 3);

  const real_t p0 = (2.e0 * u3 - 7.e0 * u2 + 11.e0 * u1) / 6.e0;
  const real_t p1 = (-u2 + 5.e0 * u1 + 2.e0 * u0) / 6.e0;
  const real_t p2 = (2.e0 * u1 + 5.e0 * u0 - um1) / 6.e0;

  const real_t beta0 = (13.e0 / 12.e0) * std::pow(u3 - 2.e0 * u2 + u1, 2) +
                       0.25e0 * std::pow(u3 - 4.e0 * u2 + 3.e0 * u1, 2);
  const real_t beta1 = (13.e0 / 12.e0) * std::pow(u2 - 2.e0 * u1 + u0, 2) +
                       0.25e0 * std::pow(u2 - u0, 2);
  const real_t beta2 = (13.e0 / 12.e0) * std::pow(u1 - 2.e0 * u0 + um1, 2) +
                       0.25e0 * std::pow(3.e0 * u1 - 4.e0 * u0 + um1, 2);

  const real_t a0 = 0.1e0 / ((kWenoEps + beta0) * (kWenoEps + beta0));
  const real_t a1 = 0.6e0 / ((kWenoEps + beta1) * (kWenoEps + beta1));
  const real_t a2 = 0.3e0 / ((kWenoEps + beta2) * (kWenoEps + beta2));
  return (a0 * p0 + a1 * p1 + a2 * p2) / (a0 + a1 + a2);
}

Vector Euler1D_WENO_FD::flux_eval(const Vector& u) const
{
  const real_t rho = u(0);
  const real_t vx = u(1) / rho;
  const real_t pre = (u(2) - 0.5e0 * rho * vx * vx) * (gamma_ - 1.e0);
  Vector f(3);
  f << u(1), rho * vx * vx + pre, (u(2) + pre) * vx;
  return f;
}

void Euler1D_WENO_FD::eigensystem_from_state(const Vector& u,
                                             Matrix* R,
                                             Matrix* L,
                                             real_t* max_lambda) const
{
  const real_t rho = u(0);
  const real_t v = u(1) / rho;
  const real_t p = (u(2) - 0.5e0 * rho * v * v) * (gamma_ - 1.e0);
  const real_t c = std::sqrt(gamma_ * p / rho);
  const real_t H = (u(2) + p) / rho;

  R->resize(3, 3);
  (*R)(0, 0) = 1.e0;      (*R)(0, 1) = 1.e0;     (*R)(0, 2) = 1.e0;
  (*R)(1, 0) = v - c;     (*R)(1, 1) = v;        (*R)(1, 2) = v + c;
  (*R)(2, 0) = H - v * c; (*R)(2, 1) = 0.5e0 * v * v; (*R)(2, 2) = H + v * c;

  *L = R->inverse();
  *max_lambda = std::max(std::abs(v - c), std::max(std::abs(v), std::abs(v + c)));
}

void Euler1D_WENO_FD::Lu_compute(const std::vector<Vector>& u, const real_t& Trun,
                                 std::vector<Vector>* Lu)
{
  const int N = static_cast<int>(u[0].size());
  const real_t hx = mesh1D_->gethx();
  const std::vector<Vector> ue = extend_with_ghost(u, Trun);

  std::vector<Vector> flux(3, Vector::Zero(N + 1));
  for (int i = 0; i <= N; ++i) {
    const int iface = i + 2;

    Vector u_avg(3);
    for (int k = 0; k < 3; ++k) {
      u_avg(k) = 0.5e0 * (ue[k](iface) + ue[k](iface + 1));
    }

    Matrix R(3, 3), L(3, 3);
    real_t alpha = 0.e0;
    eigensystem_from_state(u_avg, &R, &L, &alpha);

    std::vector<Vector> gplus(3, Vector::Zero(ue[0].size()));
    std::vector<Vector> gminus(3, Vector::Zero(ue[0].size()));

    for (int s = iface - 3; s <= iface + 3; ++s) {
      Vector uq(3);
      uq << ue[0](s), ue[1](s), ue[2](s);
      const Vector fq = flux_eval(uq);

      const Vector w = L * uq;
      const Vector g = L * fq;
      const Vector gp = 0.5e0 * (g.array() + alpha * w.array()).matrix();
      const Vector gm = 0.5e0 * (g.array() - alpha * w.array()).matrix();

      for (int m = 0; m < 3; ++m) {
        gplus[m](s) = gp(m);
        gminus[m](s) = gm(m);
      }
    }

    Vector ghat(3);
    for (int m = 0; m < 3; ++m) {
      const real_t gp_hat = (x_order_ == 3) ? weno3_left_biased(gplus[m], iface)
                                             : weno5_left_biased(gplus[m], iface);
      const real_t gm_hat = (x_order_ == 3) ? weno3_right_biased(gminus[m], iface)
                                             : weno5_right_biased(gminus[m], iface);
      ghat(m) = gp_hat + gm_hat;
    }

    const Vector fhat = R * ghat;
    for (int k = 0; k < 3; ++k) {
      flux[k](i) = fhat(k);
    }
  }

  Lu->assign(3, Vector::Zero(N));
  for (int k = 0; k < 3; ++k) {
    (*Lu)[k].array() = -(flux[k].tail(N).array() - flux[k].head(N).array()) / hx;
  }
}

void Euler1D_WENO_FD::updateAll(const real_t& Trun, const real_t& dt)
{
  const Matrix& A = rk_table_->getA();
  const Matrix& Be = rk_table_->getBe();
  const Vector& c = rk_table_->getc();
  const int stages = rk_table_->getstages();

  for (int s = 0; s < stages; ++s) {
    for (int k = 0; k < 3; ++k) {
      u_stages_[s][k].setZero();
      for (int j = 0; j < s; ++j) {
        u_stages_[s][k] += A(s, j) * u_stages_[j][k] + dt * Be(s, j) * Lu_stages_[j][k];
      }
      if (s == 0) {
        u_stages_[s][k] = u_[k];
      }
    }
    if (s == stages - 1) {
      break;
    }
    Lu_compute(u_stages_[s], Trun + c(s) * dt, &Lu_stages_[s]);
  }
  u_ = u_stages_[stages - 1];
}

real_t Euler1D_WENO_FD::rho_init(const real_t& x) const { return 1.e0 + 0.2e0 * std::sin(2.e0 * pi_ * x); }
Vector Euler1D_WENO_FD::rho_init(const Vector& x) const { return (1.e0 + 0.2e0 * (2.e0 * pi_ * x.array()).sin()).matrix(); }
real_t Euler1D_WENO_FD::rho_bc(const real_t& x, const real_t&) const { return rho_init(x); }
Vector Euler1D_WENO_FD::rho_bc(const Vector& x, const real_t&) const { return rho_init(x); }

real_t Euler1D_WENO_FD::vx_init(const real_t&) const { return 1.e0; }
Vector Euler1D_WENO_FD::vx_init(const Vector& x) const { return Vector::Ones(x.size()); }
real_t Euler1D_WENO_FD::vx_bc(const real_t& x, const real_t&) const { return vx_init(x); }
Vector Euler1D_WENO_FD::vx_bc(const Vector& x, const real_t&) const { return vx_init(x); }

real_t Euler1D_WENO_FD::pre_init(const real_t&) const { return 1.e0; }
Vector Euler1D_WENO_FD::pre_init(const Vector& x) const { return Vector::Ones(x.size()); }
real_t Euler1D_WENO_FD::pre_bc(const real_t& x, const real_t&) const { return pre_init(x); }
Vector Euler1D_WENO_FD::pre_bc(const Vector& x, const real_t&) const { return pre_init(x); }

Euler1D_WENO_FD_period::Euler1D_WENO_FD_period(const FDmesh_period* mesh1D,
                                               const EX_TVDRK* rk_table,
                                               const int x_order)
  : Euler1D_WENO_FD(mesh1D, rk_table, x_order), mesh1D_period_(mesh1D)
{
}

std::vector<Vector> Euler1D_WENO_FD_period::extend_with_ghost(const std::vector<Vector>& u,
                                                              const real_t&) const
{
  const int ng = 3;
  const int N = static_cast<int>(u[0].size());
  std::vector<Vector> ue(3, Vector::Zero(N + 2 * ng));
  for (int k = 0; k < 3; ++k) {
    ue[k].segment(ng, N) = u[k];
    ue[k].head(ng) = u[k].tail(ng);
    ue[k].tail(ng) = u[k].head(ng);
  }
  return ue;
}

void Euler1D_WENO_FD_period::extend_left_ghost(const int, const real_t, std::vector<Vector>*) const {}
void Euler1D_WENO_FD_period::extend_right_ghost(const int, const real_t, std::vector<Vector>*) const {}

real_t Euler1D_WENO_FD_period::rho_init(const real_t& x) const { return 1.e0 + 0.2e0 * std::sin(2.e0 * pi_ * x); }
Vector Euler1D_WENO_FD_period::rho_init(const Vector& x) const { return (1.e0 + 0.2e0 * (2.e0 * pi_ * x.array()).sin()).matrix(); }
real_t Euler1D_WENO_FD_period::rho_bc(const real_t& x, const real_t&) const { return rho_init(x); }
Vector Euler1D_WENO_FD_period::rho_bc(const Vector& x, const real_t&) const { return rho_init(x); }

real_t Euler1D_WENO_FD_period::vx_init(const real_t&) const { return 1.e0; }
Vector Euler1D_WENO_FD_period::vx_init(const Vector& x) const { return Vector::Ones(x.size()); }
real_t Euler1D_WENO_FD_period::vx_bc(const real_t& x, const real_t&) const { return vx_init(x); }
Vector Euler1D_WENO_FD_period::vx_bc(const Vector& x, const real_t&) const { return vx_init(x); }

real_t Euler1D_WENO_FD_period::pre_init(const real_t&) const { return 1.e0; }
Vector Euler1D_WENO_FD_period::pre_init(const Vector& x) const { return Vector::Ones(x.size()); }
real_t Euler1D_WENO_FD_period::pre_bc(const real_t& x, const real_t&) const { return pre_init(x); }
Vector Euler1D_WENO_FD_period::pre_bc(const Vector& x, const real_t&) const { return pre_init(x); }

} // namespace QUEST

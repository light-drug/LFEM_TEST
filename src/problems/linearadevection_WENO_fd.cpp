#include "linearadevection_WENO_fd.hpp"

#include <cmath>

namespace QUEST
{
namespace
{
constexpr real_t kWenoEps = 1.e-6;
}

LinearAdvection_WENO_FD::LinearAdvection_WENO_FD(const FDmesh* mesh1D,
                                                 const EX_TVDRK* rk_table,
                                                 const real_t a,
                                                 const int x_order)
  : mesh1D_(mesh1D), rk_table_(rk_table), a_(a), x_order_(x_order),
    pi_(3.141592653589793238462643383279), CFL_(0.4e0)
{
  if (x_order_ != 3 && x_order_ != 5) {
    QUEST_ERROR("LinearAdvection_WENO_FD only supports WENO3 or WENO5 in space.");
  }
}

void LinearAdvection_WENO_FD::setCFL(const real_t& cfl)
{
  CFL_ = cfl;
}

const Vector& LinearAdvection_WENO_FD::getu() const
{
  return u_;
}

void LinearAdvection_WENO_FD::init()
{
  const Vector& xc = mesh1D_->getCellCenter_vec();
  u_ = u_init(xc);
  u_stages_.assign(rk_table_->getstages(), Vector::Zero(u_.size()));
  Lu_stages_.assign(rk_table_->getstages(), Vector::Zero(u_.size()));
}

void LinearAdvection_WENO_FD::setdt(real_t* dt) const
{
  const int& t_order = rk_table_->gett_order();
  *dt = CFL_ * mesh1D_->gethx() / std::abs(a_);
  *dt = std::pow(*dt, real_t(x_order_) / real_t(t_order));
}

real_t LinearAdvection_WENO_FD::left_ghost_value(const int ghost_id, const real_t& Trun) const
{
  const real_t xg = mesh1D_->getx1() - (real_t(ghost_id) + 0.5e0) * mesh1D_->gethx();
  return u_bc(xg, Trun);
}

real_t LinearAdvection_WENO_FD::right_ghost_value(const Vector& u, const int ghost_id, const real_t&) const
{
  const int N = static_cast<int>(u.size());
  return u(N - 1 - ghost_id);
}

Vector LinearAdvection_WENO_FD::extend_with_ghost(const Vector& u, const real_t& Trun) const
{
  const int ng = 3;
  const int N = static_cast<int>(u.size());
  Vector ue(N + 2 * ng);
  ue.segment(ng, N) = u;
  for (int g = 0; g < ng; ++g) {
    ue(ng - 1 - g) = left_ghost_value(g, Trun);
    ue(ng + N + g) = right_ghost_value(u, g, Trun);
  }
  return ue;
}

real_t LinearAdvection_WENO_FD::weno3_left_biased(const Vector& ue, const int iface) const
{
  const real_t um1 = ue(iface - 1);
  const real_t u0 = ue(iface);
  const real_t u1 = ue(iface + 1);

  const real_t p0 = 0.5e0 * (- um1 + 3.e0 * u0);
  const real_t p1 = 0.5e0 * (u0 + u1);

  const real_t beta0 = (u0 - um1) * (u0 - um1);
  const real_t beta1 = (u1 - u0) * (u1 - u0);

  const real_t d0 = 1.e0 / 3.e0;
  const real_t d1 = 2.e0 / 3.e0;
  return d0 * p0 + d1 * p1;
  // const real_t a0 = d0 / ((kWenoEps + beta0) * (kWenoEps + beta0));
  // const real_t a1 = d1 / ((kWenoEps + beta1) * (kWenoEps + beta1));
  // const real_t asum = a0 + a1;

  // return (a0 * p0 + a1 * p1) / asum;
}

real_t LinearAdvection_WENO_FD::weno5_left_biased(const Vector& ue, const int iface) const
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

  const real_t d0 = 0.1e0;
  const real_t d1 = 0.6e0;
  const real_t d2 = 0.3e0;

  const real_t a0 = d0 / ((kWenoEps + beta0) * (kWenoEps + beta0));
  const real_t a1 = d1 / ((kWenoEps + beta1) * (kWenoEps + beta1));
  const real_t a2 = d2 / ((kWenoEps + beta2) * (kWenoEps + beta2));
  const real_t asum = a0 + a1 + a2;

  return (a0 * p0 + a1 * p1 + a2 * p2) / asum;
}

void LinearAdvection_WENO_FD::Lu_compute(const Vector& u, const real_t& Trun, Vector* Lu)
{
  const int N = static_cast<int>(u.size());
  const real_t hx = mesh1D_->gethx();

  Vector ue = extend_with_ghost(u, Trun);
  Vector flux(N + 1);
  for (int i = 0; i <= N; ++i) {
    const int iface = i + 2;
    const real_t uL = (x_order_ == 3) ? weno3_left_biased(ue, iface) : weno5_left_biased(ue, iface);
    flux(i) = a_ * uL;
  }

  Lu->resize(N);
  Lu->array() = - (flux.tail(N).array() - flux.head(N).array()) / hx;
}

void LinearAdvection_WENO_FD::updateAll(const real_t& Trun, const real_t& dt)
{
  const Matrix& A = rk_table_->getA();
  const Matrix& Be = rk_table_->getBe();
  const Vector& c = rk_table_->getc();
  const int& stages = rk_table_->getstages();
  const int& Nx = mesh1D_->getxDiv();

  Vector zeros = Vector::Zero(Nx);
  for (int s = 0; s < stages; ++s) {
    u_stages_[s] = zeros;
    for (int j = 0; j < s; ++j) {
      u_stages_[s] += A(s, j) * u_stages_[j] + dt * Be(s, j) * Lu_stages_[j];
    }
    if (s == 0)
    {
      u_stages_[s] = u_;
    };
    if (s == stages - 1) break;
    Lu_compute(u_stages_[s], Trun + c(s) * dt, &Lu_stages_[s]);
  }
  u_ = u_stages_[stages - 1];
}

real_t LinearAdvection_WENO_FD::u_init(const real_t& x) const
{
  return std::sin(2.e0 * pi_ * x);
}

Vector LinearAdvection_WENO_FD::u_init(const Vector& x) const
{
  return (2.e0 * pi_ * x.array()).sin().matrix();
}

real_t LinearAdvection_WENO_FD::u_bc(const real_t& x, const real_t& t) const
{
  return u_real(x, t);
}

Vector LinearAdvection_WENO_FD::u_bc(const Vector& x, const real_t& t) const
{
  Vector y = x.array() - a_ * t;
  return u_init(y);
}

real_t LinearAdvection_WENO_FD::u_real(const real_t& x, const real_t& t) const
{
  return u_init(x - a_ * t);
}

Vector LinearAdvection_WENO_FD::u_real(const Vector& x, const real_t& t) const
{
  Vector y = x.array() - a_ * t;
  return u_init(y);
}

LinearAdvection_WENO_FD_period::LinearAdvection_WENO_FD_period(const FDmesh_period* mesh1D,
                                                               const EX_TVDRK* rk_table,
                                                               const real_t a,
                                                               const int x_order)
  : LinearAdvection_WENO_FD(mesh1D, rk_table, a, x_order), mesh1D_period_(mesh1D)
{
}

Vector LinearAdvection_WENO_FD_period::extend_with_ghost(const Vector& u, const real_t&) const
{
  const int ng = 3;
  const int N = static_cast<int>(u.size());
  Vector ue(N + 2 * ng);
  ue.segment(ng, N) = u;
  ue.head(ng) = u.tail(ng);
  ue.tail(ng) = u.head(ng);
  return ue;
}

real_t LinearAdvection_WENO_FD_period::left_ghost_value(const int, const real_t&) const
{
  return 0.e0;
}

real_t LinearAdvection_WENO_FD_period::right_ghost_value(const Vector&, const int, const real_t&) const
{
  return 0.e0;
}

real_t LinearAdvection_WENO_FD_period::u_init(const real_t& x) const
{
  return std::sin(2.e0 * pi_ * x);
}

Vector LinearAdvection_WENO_FD_period::u_init(const Vector& x) const
{
  return (2.e0 * pi_ * x.array()).sin().matrix();
}

real_t LinearAdvection_WENO_FD_period::u_real(const real_t& x, const real_t& t) const
{
  const real_t x1 = mesh1D_period_->getx1();
  const real_t x2 = mesh1D_period_->getx2();
  const real_t L = x2 - x1;
  real_t y = x - a_ * t;
  y = y - std::floor((y - x1) / L) * L;
  return u_init(y);
}

Vector LinearAdvection_WENO_FD_period::u_real(const Vector& x, const real_t& t) const
{
  const real_t x1 = mesh1D_period_->getx1();
  const real_t x2 = mesh1D_period_->getx2();
  const real_t L = x2 - x1;
  Vector y = x.array() - a_ * t;
  y = (y.array() - x1 - ((y.array() - x1) / L).floor() * L + x1).matrix();
  return u_init(y);
}

} // namespace QUEST

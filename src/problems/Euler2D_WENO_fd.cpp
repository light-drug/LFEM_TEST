#include "Euler2D_WENO_fd.hpp"

#include <cmath>

#include "error.hpp"

namespace QUEST
{
namespace
{
constexpr real_t kWenoEps = 1.e-6;
}

Euler2D_WENO_FD::Euler2D_WENO_FD(const FDmesh2D* mesh2D,
                                 const EX_TVDRK* rk_table,
                                 const int x_order)
  : mesh2D_(mesh2D), rk_table_(rk_table), x_order_(x_order),
    gamma_(1.4e0), CFL_(0.4e0), num_equations_(4)
{
  if (x_order_ != 3 && x_order_ != 5) {
    QUEST_ERROR("Euler2D_WENO_FD only supports WENO3/WENO5.");
  }
}

void Euler2D_WENO_FD::setgamma(const real_t& gamma) { gamma_ = gamma; }
void Euler2D_WENO_FD::setCFL(const real_t& cfl) { CFL_ = cfl; }
const std::vector<Matrix>& Euler2D_WENO_FD::getu() const { return u_; }
const Matrix& Euler2D_WENO_FD::getrho() const { return u_[0]; }

void Euler2D_WENO_FD::init()
{
  const int Nx = mesh2D_->getxDiv();
  const int Ny = mesh2D_->getyDiv();
  const Vector& xc = mesh2D_->getCellCenter_xvec();
  const Vector& yc = mesh2D_->getCellCenter_yvec();

  u_.assign(num_equations_, Matrix::Zero(Nx, Ny));

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      const real_t x = xc(i);
      const real_t y = yc(j);
      const real_t rho = rho_init(x, y);
      const real_t vx = vx_init(x, y);
      const real_t vy = vy_init(x, y);
      const real_t pre = pre_init(x, y);

      u_[0](i, j) = rho;
      u_[1](i, j) = rho * vx;
      u_[2](i, j) = rho * vy;
      u_[3](i, j) = pre / (gamma_ - 1.e0) + 0.5e0 * rho * (vx * vx + vy * vy);
    }
  }

  const int stages = rk_table_->getstages();
  u_stages_.assign(stages, std::vector<Matrix>(num_equations_, Matrix::Zero(Nx, Ny)));
  Lu_stages_.assign(stages, std::vector<Matrix>(num_equations_, Matrix::Zero(Nx, Ny)));
}

real_t Euler2D_WENO_FD::max_wave_speed(const std::vector<Matrix>& u) const
{
  Matrix rho = u[0];
  Matrix vx = u[1].array() / rho.array();
  Matrix vy = u[2].array() / rho.array();
  Matrix p = (gamma_ - 1.e0)
           * (u[3].array() - 0.5e0 * rho.array() * (vx.array().square() + vy.array().square()));
  Matrix c = (gamma_ * p.array() / rho.array()).sqrt();
  Matrix s = (vx.array().square() + vy.array().square()).sqrt() + c.array();
  return s.maxCoeff();
}

void Euler2D_WENO_FD::setdt(real_t* dt) const
{
  const int& t_order = rk_table_->gett_order();
  const real_t hmin = std::min(mesh2D_->gethx(), mesh2D_->gethy());
  *dt = CFL_ * hmin / max_wave_speed(u_);
  *dt = std::pow(*dt, real_t(x_order_) / real_t(t_order));
}

Vector Euler2D_WENO_FD::weno3_left_biased(const std::vector<Vector>& unei) const
{
  Vector temp(num_equations_);
  for (int m = 0; m < num_equations_; ++m) {
    const real_t um1 = unei[0](m), u0 = unei[1](m), u1 = unei[2](m);
    const real_t p0 = 0.5e0 * (-um1 + 3.e0 * u0);
    const real_t p1 = 0.5e0 * (u0 + u1);
    const real_t b0 = (u0 - um1) * (u0 - um1);
    const real_t b1 = (u1 - u0) * (u1 - u0);
    const real_t a0 = (1.e0 / 3.e0) / ((kWenoEps + b0) * (kWenoEps + b0));
    const real_t a1 = (2.e0 / 3.e0) / ((kWenoEps + b1) * (kWenoEps + b1));
    temp(m) = (a0 * p0 + a1 * p1) / (a0 + a1);
  }
  return temp;
}

Vector Euler2D_WENO_FD::weno5_left_biased(const std::vector<Vector>& unei) const
{
  Vector temp(num_equations_);
  for (int m = 0; m < num_equations_; ++m) {
    const real_t um2 = unei[0](m), um1 = unei[1](m), u0 = unei[2](m), u1 = unei[3](m), u2 = unei[4](m);
    const real_t p0 = (2.e0 * um2 - 7.e0 * um1 + 11.e0 * u0) / 6.e0;
    const real_t p1 = (-um1 + 5.e0 * u0 + 2.e0 * u1) / 6.e0;
    const real_t p2 = (2.e0 * u0 + 5.e0 * u1 - u2) / 6.e0;

    const real_t b0 = (13.e0 / 12.e0) * std::pow(um2 - 2.e0 * um1 + u0, 2)
                    + 0.25e0 * std::pow(um2 - 4.e0 * um1 + 3.e0 * u0, 2);
    const real_t b1 = (13.e0 / 12.e0) * std::pow(um1 - 2.e0 * u0 + u1, 2)
                    + 0.25e0 * std::pow(um1 - u1, 2);
    const real_t b2 = (13.e0 / 12.e0) * std::pow(u0 - 2.e0 * u1 + u2, 2)
                    + 0.25e0 * std::pow(3.e0 * u0 - 4.e0 * u1 + u2, 2);

    const real_t a0 = 0.1e0 / ((kWenoEps + b0) * (kWenoEps + b0));
    const real_t a1 = 0.6e0 / ((kWenoEps + b1) * (kWenoEps + b1));
    const real_t a2 = 0.3e0 / ((kWenoEps + b2) * (kWenoEps + b2));

    temp(m) = (a0 * p0 + a1 * p1 + a2 * p2) / (a0 + a1 + a2);
  }
  return temp;
}

Vector Euler2D_WENO_FD::weno3_right_biased(const std::vector<Vector>&unei) const
{
  return weno3_left_biased({unei[2],unei[1],unei[0]});
}

Vector Euler2D_WENO_FD::weno5_right_biased(const std::vector<Vector>&unei) const
{
  return weno5_left_biased({unei[4],unei[3],unei[2],unei[1],unei[0]});
}

Vector Euler2D_WENO_FD::flux_x_eval(const Vector& u) const
{
  const real_t rho = u(0);
  const real_t vx = u(1) / rho;
  const real_t vy = u(2) / rho;
  const real_t p = (gamma_ - 1.e0) * (u(3) - 0.5e0 * rho * (vx * vx + vy * vy));
  Vector f(4);
  f << u(1), rho * vx * vx + p, rho * vx * vy, (u(3) + p) * vx;
  return f;
}

Vector Euler2D_WENO_FD::flux_y_eval(const Vector& u) const
{
  const real_t rho = u(0);
  const real_t vx = u(1) / rho;
  const real_t vy = u(2) / rho;
  const real_t p = (gamma_ - 1.e0) * (u(3) - 0.5e0 * rho * (vx * vx + vy * vy));
  Vector f(4);
  f << u(2), rho * vx * vy, rho * vy * vy + p, (u(3) + p) * vy;
  return f;
}

Vector Euler2D_WENO_FD::get_state(const std::vector<Matrix>& u,
                                  const int i,
                                  const int j,
                                  const real_t& Trun,
                                  const bool periodic) const
{
  const int Nx = mesh2D_->getxDiv();
  const int Ny = mesh2D_->getyDiv();
  const real_t hx = mesh2D_->gethx();
  const real_t hy = mesh2D_->gethy();

  int ii = i;
  int jj = j;
  if (periodic) {
    ii = (ii % Nx + Nx) % Nx;
    jj = (jj % Ny + Ny) % Ny;
  }

  Vector U(4);
  if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) {
    for (int m = 0; m < 4; ++m) {
      U(m) = u[m](ii, jj);
    }
    return U;
  }

  real_t x = mesh2D_->getx1() + (real_t(ii) + 0.5e0) * hx;
  real_t y = mesh2D_->gety1() + (real_t(jj) + 0.5e0) * hy;
  const real_t rho = rho_bc(x, y, Trun);
  const real_t vx = vx_bc(x, y, Trun);
  const real_t vy = vy_bc(x, y, Trun);
  const real_t pre = pre_bc(x, y, Trun);

  U(0) = rho;
  U(1) = rho * vx;
  U(2) = rho * vy;
  U(3) = pre / (gamma_ - 1.e0) + 0.5e0 * rho * (vx * vx + vy * vy);
  return U;
}

void Euler2D_WENO_FD::Lu_compute(const std::vector<Matrix>& u,
                                 const real_t& Trun,
                                 std::vector<Matrix>* Lu)
{
  const int Nx = mesh2D_->getxDiv();
  const int Ny = mesh2D_->getyDiv();
  const real_t hx = mesh2D_->gethx();
  const real_t hy = mesh2D_->gethy();
  const int k = (x_order_ - 1) / 2;

  const bool periodic = false;
  const real_t alpha = max_wave_speed(u);

  std::vector<Matrix> fluxx(num_equations_, Matrix::Zero(Nx + 1, Ny));
  std::vector<Matrix> fluxy(num_equations_, Matrix::Zero(Nx, Ny + 1));

  for (int i = 0; i <= Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      std::vector<Vector> gp(x_order_), gm(x_order_);
      int idx = 0;
      for (int s = i - k; s <= i + k; ++s) {
        Vector U = get_state(u, s, j, Trun, periodic);
        gp[idx] = 0.5e0 * (flux_x_eval(U) + alpha * U);
        ++idx;
      }
      idx = 0;
      for (int s = i - k + 1; s <= i + k + 1; ++s) {
        Vector U = get_state(u, s, j, Trun, periodic);
        gm[idx] = 0.5e0 * (flux_x_eval(U) - alpha * U);
        ++idx;
      }

      Vector fhat = (x_order_ == 3) ? (weno3_left_biased(gp) + weno3_right_biased(gm))
                                    : (weno5_left_biased(gp) + weno5_right_biased(gm));
      for (int m = 0; m < num_equations_; ++m) {
        fluxx[m](i, j) = fhat(m);
      }
    }
  }

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j <= Ny; ++j) {
      std::vector<Vector> gp(x_order_), gm(x_order_);
      int idx = 0;
      for (int s = j - k; s <= j + k; ++s) {
        Vector U = get_state(u, i, s, Trun, periodic);
        gp[idx] = 0.5e0 * (flux_y_eval(U) + alpha * U);
        ++idx;
      }
      idx = 0;
      for (int s = j - k + 1; s <= j + k + 1; ++s) {
        Vector U = get_state(u, i, s, Trun, periodic);
        gm[idx] = 0.5e0 * (flux_y_eval(U) - alpha * U);
        ++idx;
      }

      Vector fhat = (x_order_ == 3) ? (weno3_left_biased(gp) + weno3_right_biased(gm))
                                    : (weno5_left_biased(gp) + weno5_right_biased(gm));
      for (int m = 0; m < num_equations_; ++m) {
        fluxy[m](i, j) = fhat(m);
      }
    }
  }

  Lu->assign(num_equations_, Matrix::Zero(Nx, Ny));
  for (int m = 0; m < num_equations_; ++m) {
    (*Lu)[m] = -((fluxx[m].bottomRows(Nx) - fluxx[m].topRows(Nx)) / hx
               + (fluxy[m].rightCols(Ny) - fluxy[m].leftCols(Ny)) / hy);
  }
}

void Euler2D_WENO_FD::updateAll(const real_t& Trun, const real_t& dt)
{
  const Matrix& A = rk_table_->getA();
  const Matrix& Be = rk_table_->getBe();
  const Vector& c = rk_table_->getc();
  const int stages = rk_table_->getstages();

  for (int s = 0; s < stages; ++s) {
    for (int k = 0; k < num_equations_; ++k) {
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

real_t Euler2D_WENO_FD::rho_init(const real_t& x, const real_t& y) const
{
  return 1.e0 + 0.2e0 * std::sin(x + y);
}

real_t Euler2D_WENO_FD::vx_init(const real_t&, const real_t&) const { return 1.e0; }
real_t Euler2D_WENO_FD::vy_init(const real_t&, const real_t&) const { return 1.e0; }
real_t Euler2D_WENO_FD::pre_init(const real_t&, const real_t&) const { return 1.e0; }

real_t Euler2D_WENO_FD::rho_bc(const real_t& x, const real_t& y, const real_t& t) const
{
  return 1.e0 + 0.2e0 * std::sin(x + y - 2.e0 * t);
}
real_t Euler2D_WENO_FD::vx_bc(const real_t&, const real_t&, const real_t&) const { return 1.e0; }
real_t Euler2D_WENO_FD::vy_bc(const real_t&, const real_t&, const real_t&) const { return 1.e0; }
real_t Euler2D_WENO_FD::pre_bc(const real_t&, const real_t&, const real_t&) const { return 1.e0; }

Matrix Euler2D_WENO_FD::rho_real(const Matrix& x, const Matrix& y, const real_t& t) const
{
  return (1.e0 + 0.2e0 * (x.array() + y.array() - 2.e0 * t).sin()).matrix();
}

Euler2D_WENO_FD_period::Euler2D_WENO_FD_period(const FDmesh2D* mesh2D,
                                               const EX_TVDRK* rk_table,
                                               const int x_order)
  : Euler2D_WENO_FD(mesh2D, rk_table, x_order)
{
}


void Euler2D_WENO_FD_period::Lu_compute(const std::vector<Matrix>& u,
                                        const real_t& Trun,
                                        std::vector<Matrix>* Lu)
{
  const int Nx = mesh2D_->getxDiv();
  const int Ny = mesh2D_->getyDiv();
  const real_t hx = mesh2D_->gethx();
  const real_t hy = mesh2D_->gethy();
  const int k = (x_order_ - 1) / 2;

  const bool periodic = true;
  const real_t alpha = max_wave_speed(u);

  std::vector<Matrix> fluxx(num_equations_, Matrix::Zero(Nx + 1, Ny));
  std::vector<Matrix> fluxy(num_equations_, Matrix::Zero(Nx, Ny + 1));

  for (int i = 0; i <= Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      std::vector<Vector> gp(x_order_), gm(x_order_);
      int idx = 0;
      for (int s = i - k; s <= i + k; ++s) {
        Vector U = get_state(u, s, j, Trun, periodic);
        gp[idx] = 0.5e0 * (flux_x_eval(U) + alpha * U);
        ++idx;
      }
      idx = 0;
      for (int s = i - k + 1; s <= i + k + 1; ++s) {
        Vector U = get_state(u, s, j, Trun, periodic);
        gm[idx] = 0.5e0 * (flux_x_eval(U) - alpha * U);
        ++idx;
      }
      Vector fhat = (x_order_ == 3) ? (weno3_left_biased(gp) + weno3_right_biased(gm))
                                    : (weno5_left_biased(gp) + weno5_right_biased(gm));
      for (int m = 0; m < num_equations_; ++m) {
        fluxx[m](i, j) = fhat(m);
      }
    }
  }

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j <= Ny; ++j) {
      std::vector<Vector> gp(x_order_), gm(x_order_);
      int idx = 0;
      for (int s = j - k; s <= j + k; ++s) {
        Vector U = get_state(u, i, s, Trun, periodic);
        gp[idx] = 0.5e0 * (flux_y_eval(U) + alpha * U);
        ++idx;
      }
      idx = 0;
      for (int s = j - k + 1; s <= j + k + 1; ++s) {
        Vector U = get_state(u, i, s, Trun, periodic);
        gm[idx] = 0.5e0 * (flux_y_eval(U) - alpha * U);
        ++idx;
      }
      Vector fhat = (x_order_ == 3) ? (weno3_left_biased(gp) + weno3_right_biased(gm))
                                    : (weno5_left_biased(gp) + weno5_right_biased(gm));
      for (int m = 0; m < num_equations_; ++m) {
        fluxy[m](i, j) = fhat(m);
      }
    }
  }

  Lu->assign(num_equations_, Matrix::Zero(Nx, Ny));
  for (int m = 0; m < num_equations_; ++m) {
    (*Lu)[m] = -((fluxx[m].bottomRows(Nx) - fluxx[m].topRows(Nx)) / hx
               + (fluxy[m].rightCols(Ny) - fluxy[m].leftCols(Ny)) / hy);
  }
}

real_t Euler2D_WENO_FD_period::rho_init(const real_t& x, const real_t& y) const
{
  return 1.e0 + 0.2e0 * std::sin(x + y);
}

real_t Euler2D_WENO_FD_period::rho_bc(const real_t& x, const real_t& y, const real_t& t) const
{
  QUEST_ERROR("rho_bc should not be called for periodic case.");
  return rho_init(x, y);
}

Matrix Euler2D_WENO_FD_period::rho_real(const Matrix& x, const Matrix& y, const real_t& t) const
{
  return (1.e0 + 0.2e0 * (x.array() + y.array() - 2.e0 * t).sin()).matrix();
}

} // namespace QUEST

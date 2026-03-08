#include "Euler1D_acctest.hpp"

#include <algorithm>
#include <cmath>

namespace QUEST
{

Euler1D_DG_TVDRK::Euler1D_DG_TVDRK(const TensorMesh1D* mesh1D,
                                  const fespace1D* fe,
                                  const EX_TVDRK* rk_table)
  : HyperbolicProblems1DBase(mesh1D, fe, rk_table)
{
  gamma_ = 1.4e0;
  cfl_ = 0.2e0;
  pi_ = std::acos(-1.e0);
};

Matrix Euler1D_DG_TVDRK::rho_init(const Matrix& x) {
  Matrix result;
  result = pi_ * x;
  result = 1 + 0.2e0 * result.array().sin();
  return result;
}

real_t Euler1D_DG_TVDRK::rho_init(const real_t& x) {
  real_t result;
  result = pi_ * x;
  result = 1 + 0.2e0 * std::sin(result);
  return result;
};

real_t Euler1D_DG_TVDRK::rho_bc(const real_t& x, const real_t& t) {
  real_t result;
  result = pi_ * (x - t);
  result = 1 + 0.2e0 * std::sin(result);
  return result;
};

Matrix Euler1D_DG_TVDRK::vx_init(const Matrix& x) {
  Matrix result;
  result = Matrix::Constant(x.rows(), x.cols(), 1.e0);
  return result;
}

real_t Euler1D_DG_TVDRK::vx_init(const real_t& x) {
  return 1.e0;
};

real_t Euler1D_DG_TVDRK::vx_bc(const real_t& x, const real_t& t) {
  return 1.e0;
};

Matrix Euler1D_DG_TVDRK::pre_init(const Matrix& x) {
  Matrix result;
  result = Matrix::Constant(x.rows(), x.cols(), 1.e0);
  return result;
};

real_t Euler1D_DG_TVDRK::pre_init(const real_t& x) {
  return 1.e0;
};

real_t Euler1D_DG_TVDRK::pre_bc(const real_t& x, const real_t& t) {
  return 1.e0;
};

void Euler1D_DG_TVDRK::init()
{
  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &u_modal_[0]);
  Matrix rho_nodal;
  fe_->Interpolate_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &rho_nodal);
  Matrix vx_nodal;
  fe_->Interpolate_Initial(
      [this](const Matrix& x) { return this->vx_init(x); }, &vx_nodal);
  Matrix rhovx_nodal = rho_nodal.array() * vx_nodal.array();
  fe_->nodal_to_modal1D(rhovx_nodal, &(u_modal_[1]));
  Matrix pre_nodal;
  fe_->Interpolate_Initial(
      [this](const Matrix& x) { return this->pre_init(x); }, &pre_nodal);
  Matrix erg_nodal =
    pre_nodal.array() / (gamma_ - 1.e0) +
    0.5e0 * vx_nodal.array().square() * rho_nodal.array();
  fe_->nodal_to_modal1D(erg_nodal, &(u_modal_[2]));
};

void Euler1D_DG_TVDRK::setdt(real_t* dt) {
  const int& k1D = fe_->getbasis()->getk1D();
  const real_t& hx = mesh1D_->gethx();

  real_t max_speed = max_speed_compute(u_modal_);
  *dt = cfl_ * hx / max_speed / real_t(2 * k1D - 1);
};

real_t Euler1D_DG_TVDRK::max_speed_compute(const std::vector<Matrix>& u_modal)
{
  std::vector<Matrix> u_nodal;
  fe_->modal_to_nodal1D(u_modal, &u_nodal);
  Matrix pre_nodal = (u_nodal[2].array() - 0.5e0 *
    (u_nodal[1].array().square() / u_nodal[0].array())) * (gamma_ - 1.e0);
  Matrix sound_nodal = (gamma_ * (pre_nodal.array() / u_nodal[0].array())).sqrt();
  Matrix vx_nodal = u_nodal[1].array() / u_nodal[0].array();
  real_t max_speed = (vx_nodal.cwiseAbs().array() + sound_nodal.array()).maxCoeff();
  return max_speed;
}

void Euler1D_DG_TVDRK::updateAll(const real_t& Trun, const real_t& dt)
{
  const Matrix& A = rk_table_->getA();
  const Matrix& Be = rk_table_->getBe();
  const Vector& c = rk_table_->getc();
  const int stages = rk_table_->getstages();

  for (int s = 0; s < stages; ++s) {
    for (int k = 0; k < num_equations_; ++k) {
      u_modal_stages_[s][k].setZero();
      for (int j = 0; j < s; ++j) {
        u_modal_stages_[s][k] += A(s, j) * u_modal_stages_[j][k]
                              + dt * Be(s, j) * Lu_stages_[j][k];
      }
      if (s == 0) {
        u_modal_stages_[s][k] = u_modal_[k];
      }
    }

    if (use_limiter_) {
      tvd_limiter_characteristic(&u_modal_stages_[s]);
    }

    if (s == stages - 1) {
      break;
    }
    Lu_compute(u_modal_stages_[s], Trun + c(s) * dt, dt, &Lu_stages_[s]);
  }
  u_modal_ = u_modal_stages_[stages - 1];
};

void Euler1D_DG_TVDRK::fu_compute(const std::vector<Matrix>& u_modal,
                                  std::vector<Matrix>* fu_nodal)
{
  std::vector<Matrix> u_nodal;
  fe_->modal_to_nodal1D(u_modal, &u_nodal);
  Matrix kie_nodal = u_nodal[1].array().square() / u_nodal[0].array();
  Matrix pre_nodal = (u_nodal[2] - 0.5e0 * kie_nodal) * (gamma_ - 1);
  Matrix vx_nodal = u_nodal[1].array() / u_nodal[0].array();

  fu_nodal->resize(num_equations_);
  fu_nodal->at(0) = u_nodal[1];
  fu_nodal->at(1) = kie_nodal + pre_nodal;
  fu_nodal->at(2) = (u_nodal[2].array() + pre_nodal.array()) * vx_nodal.array();
};

void Euler1D_DG_TVDRK::Lu_compute(const std::vector<Matrix>& u_modal,
                                  const real_t& Trun,
                                  const real_t& dt,
                                  std::vector<Matrix>* Lu)
{
  std::vector<Matrix> fu_nodal;
  fu_compute(u_modal, &fu_nodal);

  std::vector<Matrix> rhs(num_equations_);
  for (int k = 0; k < num_equations_; ++k) {
    fe_->Assemble_F(fu_nodal[k], 1, &rhs[k]);
  }

  std::vector<Matrix> dirichlet(num_equations_, Matrix::Zero(1, 2));
  const real_t xL = mesh1D_->getx1();
  const real_t xR = mesh1D_->getx2();
  dirichlet[0](0, 0) = rho_bc(xL, Trun);
  dirichlet[0](0, 1) = rho_bc(xR, Trun);
  dirichlet[1](0, 0) = dirichlet[0](0, 0) * vx_bc(xL, Trun);
  dirichlet[1](0, 1) = dirichlet[0](0, 1) * vx_bc(xR, Trun);
  dirichlet[2](0, 0) = pre_bc(xL, Trun) / (gamma_ - 1.e0)
                     + 0.5e0 * dirichlet[0](0, 0) * std::pow(vx_bc(xL, Trun), 2);
  dirichlet[2](0, 1) = pre_bc(xR, Trun) / (gamma_ - 1.e0)
                     + 0.5e0 * dirichlet[0](0, 1) * std::pow(vx_bc(xR, Trun), 2);

  fluxint_compute(u_modal, &flux_int_);
  fluxext_compute(u_modal, dirichlet, &flux_ext_);

  std::vector<Matrix> flux_rhs;
  fe_->Assemble_Flux(flux_int_, flux_ext_, &flux_rhs);

  const DiagnalMatrix& Minv = fe_->getv_u_diaginv();
  Lu->resize(num_equations_);
  for (int k = 0; k < num_equations_; ++k) {
    // NOTE:
    // Assemble_F(test_dx_order=1) and Assemble_Flux already contain the
    // geometric 1/hx scaling from physical mapping, but they do NOT apply
    // the modal mass inverse. Therefore Lu must multiply M^{-1} explicitly.
    (*Lu)[k] = Minv * (-rhs[k] + flux_rhs[k]);
  }
};

void Euler1D_DG_TVDRK::numerical_flux(const Vector& u_L,
                                      const Vector& u_R,
                                      const real_t& normal,
                                      const real_t& max_speed,
                                      Vector* flux)
{
  real_t rho_L = u_L(0), rho_R = u_R(0);
  real_t rhovx_L = u_L(1), rhovx_R = u_R(1);
  real_t vx_L = rhovx_L / rho_L, vx_R = rhovx_R / rho_R;
  real_t E_L = u_L(2), E_R = u_R(2);

  real_t kie_L = rho_L * vx_L * vx_L;
  real_t kie_R = rho_R * vx_R * vx_R;
  real_t pre_L = (E_L - 0.5e0 * kie_L) * (gamma_ - 1.e0);
  real_t pre_R = (E_R - 0.5e0 * kie_R) * (gamma_ - 1.e0);

  real_t phi_L = (E_L + pre_L) * vx_L;
  real_t phi_R = (E_R + pre_R) * vx_R;

  flux->resize(num_equations_);
  (*flux)(0) = 0.5e0 * (rhovx_L + rhovx_R) * normal - max_speed * 0.5e0 * (rho_R - rho_L);
  (*flux)(1) = 0.5e0 * (kie_L + pre_L + kie_R + pre_R)  * normal - max_speed * 0.5e0 * (rhovx_R - rhovx_L);
  (*flux)(2) = 0.5e0 * (phi_L + phi_R) * normal - max_speed * 0.5e0 * (E_R - E_L);
};

void Euler1D_DG_TVDRK::fluxint_compute(const std::vector<Matrix>& u_modal,
                                        std::vector<Matrix>* flux_int)
{
  const int intboundarynum = mesh1D_->getintboundaryNum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();

  const real_t max_speed = max_speed_compute(u_modal);
  flux_int->assign(num_equations_, Matrix::Zero(1, intboundarynum));

  for (int i = 0; i < intboundarynum; ++i)
  {
    const int cL = IntBNei[i](0);
    const int cR = IntBNei[i](1);
    Vector uL(num_equations_), uR(num_equations_), flux_face;
    for (int k = 0; k < num_equations_; ++k) {
      uL(k) = u_modal[k].col(cL).dot(boundary_u[1]);
      uR(k) = u_modal[k].col(cR).dot(boundary_u[0]);
    }
    numerical_flux(uL, uR, IntBNormal[i], max_speed, &flux_face);
    for (int k = 0; k < num_equations_; ++k) {
      (*flux_int)[k](0, i) = flux_face(k);
    }
  }
}

void Euler1D_DG_TVDRK::fluxext_compute(const std::vector<Matrix>& u_modal,
                                      const std::vector<Matrix>& Dirichlet,
                                      std::vector<Matrix>* flux_ext)
{
  const int extboundarynum = mesh1D_->getextboundaryNum();
  const std::vector<int>& ExtBNei = mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();

  const real_t max_speed = max_speed_compute(u_modal);
  flux_ext->assign(num_equations_, Matrix::Zero(1, extboundarynum));

  for (int i = 0; i < extboundarynum; ++i)
  {
    const int cell = ExtBNei[i];
    const int side = 1 - i;
    Vector u_in(num_equations_), u_bc(num_equations_), flux_face;
    for (int k = 0; k < num_equations_; ++k) {
      u_in(k) = u_modal[k].col(cell).dot(boundary_u[side]);
      u_bc(k) = Dirichlet[k](0, i);
    }
    numerical_flux(u_in, u_bc, ExtBNormal[i], max_speed, &flux_face);
    for (int k = 0; k < num_equations_; ++k) {
      (*flux_ext)[k](0, i) = flux_face(k);
    }
  }
}

void Euler1D_DG_TVDRK::setgamma(const real_t& gamma)
{
  gamma_ = gamma;
};

void Euler1D_DG_TVDRK::setlimiter(const bool enable_limiter, const real_t theta)
{
  use_limiter_ = enable_limiter;
  limiter_theta_ = theta;
}

real_t Euler1D_DG_TVDRK::minmod3(const real_t a, const real_t b, const real_t c) const
{
  if (a > 0.e0 && b > 0.e0 && c > 0.e0) {
    return std::min(a, std::min(b, c));
  }
  if (a < 0.e0 && b < 0.e0 && c < 0.e0) {
    return std::max(a, std::max(b, c));
  }
  return 0.e0;
}

void Euler1D_DG_TVDRK::tvd_limiter_characteristic(std::vector<Matrix>* u_modal)
{
  const int k1D = fe_->getbasis()->getk1D();
  if (k1D < 2) {
    return;
  }

  const int ncell = mesh1D_->getncell();
  const real_t eps = 1.e-12;

  for (int i = 0; i < ncell; ++i)
  {
    const int iL = (i == 0) ? ncell - 1 : i - 1;
    const int iR = (i == ncell - 1) ? 0 : i + 1;

    Vector U0(3), UL0(3), UR0(3);
    for (int k = 0; k < 3; ++k) {
      U0(k) = (*u_modal)[k](0, i);
      UL0(k) = (*u_modal)[k](0, iL);
      UR0(k) = (*u_modal)[k](0, iR);
    }

    const real_t rho = std::max(U0(0), eps);
    const real_t vx = U0(1) / rho;
    const real_t E = U0(2);
    const real_t pre = std::max((gamma_ - 1.e0) * (E - 0.5e0 * U0(1) * vx), eps);
    const real_t c = std::sqrt(gamma_ * pre / rho);
    const real_t H = (E + pre) / rho;

    Matrix R(3, 3), L(3, 3);
    R << 1.e0, 1.e0, 1.e0,
         vx - c, vx, vx + c,
         H - vx * c, 0.5e0 * vx * vx, H + vx * c;

    const real_t gm1 = gamma_ - 1.e0;
    const real_t cc = c * c;
    L << 0.5e0 * ((0.5e0 * gm1 * vx * vx + vx * c) / cc),
         -0.5e0 * ((gm1 * vx + c) / cc),
         0.5e0 * gm1 / cc,
         1.e0 - (gm1 * vx * vx) / cc,
         gm1 * vx / cc,
         -gm1 / cc,
         0.5e0 * ((0.5e0 * gm1 * vx * vx - vx * c) / cc),
         -0.5e0 * ((gm1 * vx - c) / cc),
         0.5e0 * gm1 / cc;

    Vector dU1(3), dUm(3), dUp(3);
    for (int k = 0; k < 3; ++k) {
      dU1(k) = (*u_modal)[k](1, i);
      dUm(k) = U0(k) - UL0(k);
      dUp(k) = UR0(k) - U0(k);
    }

    Vector a1 = L * dU1;
    Vector am = L * dUm;
    Vector ap = L * dUp;

    Vector a_lim(3);
    for (int m = 0; m < 3; ++m) {
      const real_t c1 = limiter_theta_ * am(m);
      const real_t c2 = a1(m);
      const real_t c3 = limiter_theta_ * ap(m);
      a_lim(m) = minmod3(c1, c2, c3);
    }

    Vector dU1_lim = R * a_lim;
    for (int k = 0; k < 3; ++k) {
      (*u_modal)[k](1, i) = dU1_lim(k);
      for (int j = 2; j < k1D; ++j) {
        (*u_modal)[k](j, i) = 0.e0;
      }
    }
  }
}

Euler1D_DG_TVDRK_period::Euler1D_DG_TVDRK_period(const TensorMesh1D* mesh1D,
                                                const fespace1D* fe,
                                                const EX_TVDRK* rk_table)
  : Euler1D_DG_TVDRK(mesh1D, fe, rk_table) {};

void Euler1D_DG_TVDRK_period::fluxext_compute(const std::vector<Matrix>& u_modal,
                                              const std::vector<Matrix>& Dirichlet,
                                              std::vector<Matrix>* flux_ext)
{
  QUEST_VERIFY(mesh1D_->IsPeriodBoundary(), "The mesh is not periodical !");
  const std::vector<Eigen::Vector2i>& PNei = mesh1D_->getextboundaryneighbors_period();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();

  const real_t max_speed = max_speed_compute(u_modal);
  flux_ext->assign(num_equations_, Matrix::Zero(1, 2));

  for (int i = 0; i < 2; ++i)
  {
    const int cIn = PNei[i](0);
    const int cOut = PNei[i](1);
    const int inSide = 1 - i;
    const int outSide = i;
    Vector u_in(num_equations_), u_out(num_equations_), flux_face;
    for (int k = 0; k < num_equations_; ++k) {
      u_in(k) = u_modal[k].col(cIn).dot(boundary_u[inSide]);
      u_out(k) = u_modal[k].col(cOut).dot(boundary_u[outSide]);
    }
    numerical_flux(u_in, u_out, ExtBNormal[i], max_speed, &flux_face);
    for (int k = 0; k < num_equations_; ++k) {
      (*flux_ext)[k](0, i) = flux_face(k);
    }
  }
}

real_t Euler1D_DG_TVDRK_period::rho_bc(const real_t& x, const real_t& t) {
  return rho_real(x, t);
}

real_t Euler1D_DG_TVDRK_period::vx_bc(const real_t& x, const real_t& t) {
  return vx_real(x, t);
}

real_t Euler1D_DG_TVDRK_period::pre_bc(const real_t& x, const real_t& t) {
  return pre_real(x, t);
}

Matrix Euler1D_DG_TVDRK_period::rho_real(const Matrix& x, const real_t& t) {
  Matrix result = pi_ * (x.array() - t).matrix();
  result = 1.e0 + 0.2e0 * result.array().sin();
  return result;
}

real_t Euler1D_DG_TVDRK_period::rho_real(const real_t& x, const real_t& t) {
  return 1.e0 + 0.2e0 * std::sin(pi_ * (x - t));
}

Matrix Euler1D_DG_TVDRK_period::vx_real(const Matrix& x, const real_t& t) {
  return Matrix::Ones(x.rows(), x.cols());
}

real_t Euler1D_DG_TVDRK_period::vx_real(const real_t& x, const real_t& t) {
  return 1.e0;
}

Matrix Euler1D_DG_TVDRK_period::pre_real(const Matrix& x, const real_t& t) {
  return Matrix::Ones(x.rows(), x.cols());
}

real_t Euler1D_DG_TVDRK_period::pre_real(const real_t& x, const real_t& t) {
  return 1.e0;
}

} // namespace QUEST

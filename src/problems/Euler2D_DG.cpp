#include "Euler2D_DG.hpp"

#include <cmath>

#include "error.hpp"

namespace QUEST
{

Euler2D_DG_TVDRK::Euler2D_DG_TVDRK(const TensorMesh2D* mesh2D,
                                   const fespace2D* fe,
                                   const EX_TVDRK* rk_table)
  : HyperbolicProblems2DBase(mesh2D, fe, rk_table), gamma_(1.4e0), cfl_(0.2e0),
    pi_(std::acos(-1.e0))
{
}

Matrix Euler2D_DG_TVDRK::rho_init(const Matrix& x, const Matrix& y)
{
  return (1.e0 + 0.2e0 * (x.array() + y.array()).sin()).matrix();
}

Matrix Euler2D_DG_TVDRK::vx_init(const Matrix& x, const Matrix& y)
{
  return Matrix::Ones(x.rows(), x.cols());
}

Matrix Euler2D_DG_TVDRK::vy_init(const Matrix& x, const Matrix& y)
{
  return Matrix::Ones(x.rows(), x.cols());
}

Matrix Euler2D_DG_TVDRK::pre_init(const Matrix& x, const Matrix& y)
{
  return Matrix::Ones(x.rows(), x.cols());
}

real_t Euler2D_DG_TVDRK::rho_bc(const real_t& x, const real_t& y, const real_t& t)
{
  return 1.e0 + 0.2e0 * std::sin(x + y - 2.e0 * t);
}

real_t Euler2D_DG_TVDRK::vx_bc(const real_t&, const real_t&, const real_t&)
{
  return 1.e0;
}

real_t Euler2D_DG_TVDRK::vy_bc(const real_t&, const real_t&, const real_t&)
{
  return 1.e0;
}

real_t Euler2D_DG_TVDRK::pre_bc(const real_t&, const real_t&, const real_t&)
{
  return 1.e0;
}

void Euler2D_DG_TVDRK::init()
{
  fe_->Project_Initial([this](const Matrix& x, const Matrix& y) { return this->rho_init(x, y); },
                       &u_modal_[0]);
  Matrix rho_nodal, vx_nodal, vy_nodal, pre_nodal;
  fe_->Interpolate_Initial([this](const Matrix& x, const Matrix& y) { return this->rho_init(x, y); },
                           &rho_nodal);
  fe_->Interpolate_Initial([this](const Matrix& x, const Matrix& y) { return this->vx_init(x, y); },
                           &vx_nodal);
  fe_->Interpolate_Initial([this](const Matrix& x, const Matrix& y) { return this->vy_init(x, y); },
                           &vy_nodal);
  fe_->Interpolate_Initial([this](const Matrix& x, const Matrix& y) { return this->pre_init(x, y); },
                           &pre_nodal);

  fe_->nodal_to_modal2D((rho_nodal.array() * vx_nodal.array()).matrix(), &u_modal_[1]);
  fe_->nodal_to_modal2D((rho_nodal.array() * vy_nodal.array()).matrix(), &u_modal_[2]);
  Matrix erg_nodal = pre_nodal.array() / (gamma_ - 1.e0)
                   + 0.5e0 * rho_nodal.array() * (vx_nodal.array().square() + vy_nodal.array().square());
  fe_->nodal_to_modal2D(erg_nodal, &u_modal_[3]);
}

void Euler2D_DG_TVDRK::setdt(real_t* dt)
{
  const int& k1D = fe_->getbasis()->getk1D();
  const real_t hmin = std::min(mesh2D_->gethx(), mesh2D_->gethy());
  *dt = cfl_ * hmin / max_speed_compute(u_modal_) / real_t(2 * k1D - 1);
}

real_t Euler2D_DG_TVDRK::max_speed_compute(const std::vector<Matrix>& u_modal)
{
  std::vector<Matrix> u_nodal;
  fe_->modal_to_nodal2D(u_modal, &u_nodal);
  Matrix rho = u_nodal[0];
  Matrix vx = u_nodal[1].array() / rho.array();
  Matrix vy = u_nodal[2].array() / rho.array();
  Matrix pre = (gamma_ - 1.e0)
             * (u_nodal[3].array() - 0.5e0 * (u_nodal[1].array().square() + u_nodal[2].array().square()) / rho.array());
  Matrix cs = (gamma_ * pre.array() / rho.array()).sqrt();
  Matrix speed = (vx.array().square() + vy.array().square()).sqrt() + cs.array();
  return speed.maxCoeff();
}

void Euler2D_DG_TVDRK::updateAll(const real_t& Trun, const real_t& dt)
{
  const Matrix& A = rk_table_->getA();
  const Matrix& Be = rk_table_->getBe();
  const Vector& c = rk_table_->getc();
  const int stages = rk_table_->getstages();

  for (int s = 0; s < stages; ++s) {
    for (int k = 0; k < num_equations_; ++k) {
      u_modal_stages_[s][k].setZero();
      for (int j = 0; j < s; ++j) {
        u_modal_stages_[s][k] += A(s, j) * u_modal_stages_[j][k] + dt * Be(s, j) * Lu_stages_[j][k];
      }
      if (s == 0) {
        u_modal_stages_[s][k] = u_modal_[k];
      }
    }
    if (s == stages - 1) {
      break;
    }
    Lu_compute(u_modal_stages_[s], Trun + c(s) * dt, dt, &Lu_stages_[s]);
  }
  u_modal_ = u_modal_stages_[stages - 1];
}

void Euler2D_DG_TVDRK::flux_compute(const std::vector<Matrix>& u_modal,
                                    std::vector<Matrix>* fu_nodal,
                                    std::vector<Matrix>* fv_nodal)
{
  std::vector<Matrix> u_nodal;
  fe_->modal_to_nodal2D(u_modal, &u_nodal);

  Matrix rho = u_nodal[0];
  Matrix rhou = u_nodal[1];
  Matrix rhov = u_nodal[2];
  Matrix E = u_nodal[3];
  Matrix vx = rhou.array() / rho.array();
  Matrix vy = rhov.array() / rho.array();
  Matrix pre = (gamma_ - 1.e0)
             * (E.array() - 0.5e0 * rho.array() * (vx.array().square() + vy.array().square()));

  fu_nodal->resize(num_equations_);
  fv_nodal->resize(num_equations_);
  (*fu_nodal)[0] = rhou;
  (*fu_nodal)[1] = rhou.array() * vx.array() + pre.array();
  (*fu_nodal)[2] = rhov.array() * vx.array();
  (*fu_nodal)[3] = (E.array() + pre.array()) * vx.array();

  (*fv_nodal)[0] = rhov;
  (*fv_nodal)[1] = rhou.array() * vy.array();
  (*fv_nodal)[2] = rhov.array() * vy.array() + pre.array();
  (*fv_nodal)[3] = (E.array() + pre.array()) * vy.array();
}

void Euler2D_DG_TVDRK::Lu_compute(const std::vector<Matrix>& u_modal,
                                  const real_t& Trun,
                                  const real_t& dt,
                                  std::vector<Matrix>* Lu)
{
  std::vector<Matrix> fu_nodal, fv_nodal;
  flux_compute(u_modal, &fu_nodal, &fv_nodal);

  std::vector<Matrix> rhsx(num_equations_), rhsy(num_equations_);
  fe_->Assemble_F(fu_nodal, 1, 0, &rhsx);
  fe_->Assemble_F(fv_nodal, 0, 1, &rhsy);

  const int extboundarynum = mesh2D_->getextboundaryNum();
  std::vector<Matrix> dirichlet(num_equations_, Matrix::Zero(fe_->getnumqua1d(), extboundarynum));
  const std::vector<Eigen::Vector2d>& bc_center = mesh2D_->getextboundarycenter();
  for (int i = 0; i < extboundarynum; ++i) {
    const real_t x = bc_center[i](0);
    const real_t y = bc_center[i](1);
    const real_t rho = rho_bc(x, y, Trun);
    const real_t vx = vx_bc(x, y, Trun);
    const real_t vy = vy_bc(x, y, Trun);
    const real_t pre = pre_bc(x, y, Trun);
    for (int q = 0; q < fe_->getnumqua1d(); ++q) {
      dirichlet[0](q, i) = rho;
      dirichlet[1](q, i) = rho * vx;
      dirichlet[2](q, i) = rho * vy;
      dirichlet[3](q, i) = pre / (gamma_ - 1.e0) + 0.5e0 * rho * (vx * vx + vy * vy);
    }
  }

  fluxint_compute(u_modal, &flux_int_);
  fluxext_compute(u_modal, dirichlet, &flux_ext_);

  std::vector<Matrix> flux_rhs;
  fe_->Assemble_Flux(flux_int_, flux_ext_, &flux_rhs);

  const DiagnalMatrix& Minv = fe_->getv_u_diaginv();
  *Lu = Minv * (rhsx + rhsy - flux_rhs);
}

void Euler2D_DG_TVDRK::numerical_flux(const Vector& u_L, const Vector& u_R,
                                      const Eigen::Vector2d& normal,
                                      const real_t& max_speed,
                                      Vector* flux)
{
  const real_t nx = normal(0);
  const real_t ny = normal(1);

  auto flux_dot_n = [&](const Vector& u) {
    const real_t rho = u(0);
    const real_t rhou = u(1);
    const real_t rhov = u(2);
    const real_t E = u(3);
    const real_t vx = rhou / rho;
    const real_t vy = rhov / rho;
    const real_t pre = (gamma_ - 1.e0) * (E - 0.5e0 * rho * (vx * vx + vy * vy));
    Vector fn(4);
    fn(0) = rhou * nx + rhov * ny;
    fn(1) = (rhou * vx + pre) * nx + (rhou * vy) * ny;
    fn(2) = (rhov * vx) * nx + (rhov * vy + pre) * ny;
    fn(3) = (E + pre) * (vx * nx + vy * ny);
    return fn;
  };

  *flux = 0.5e0 * (flux_dot_n(u_L) + flux_dot_n(u_R)) - 0.5e0 * max_speed * (u_R - u_L);
}

void Euler2D_DG_TVDRK::fluxint_compute(const std::vector<Matrix>& u_modal,
                                       std::vector<Matrix>* flux_int)
{
  const int intboundarynum = mesh2D_->getintboundaryNum();
  const int nq1d = fe_->getnumqua1d();
  const auto& IntBNei = mesh2D_->getintboundaryneighbors();
  const auto& IntBNormal = mesh2D_->getintboundarynormal();
  const auto& IntBType = mesh2D_->getintboundarytypeindex();
  const auto& boundary_u = fe_->getboundary_u();

  const real_t max_speed = max_speed_compute(u_modal);
  flux_int->assign(num_equations_, Matrix::Zero(nq1d, intboundarynum));

  for (int i = 0; i < intboundarynum; ++i) {
    const int cL = IntBNei[i](0);
    const int cR = IntBNei[i](1);
    const int sL = IntBType[i](0);
    const int sR = IntBType[i](1);
    Eigen::Vector2d normal = IntBNormal[i];

    for (int q = 0; q < nq1d; ++q) {
      Vector uL(4), uR(4), flux_face;
      for (int k = 0; k < num_equations_; ++k) {
        uL(k) = u_modal[k].col(cL).dot(boundary_u[sL].col(q));
        uR(k) = u_modal[k].col(cR).dot(boundary_u[sR].col(q));
      }
      numerical_flux(uL, uR, normal, max_speed, &flux_face);
      for (int k = 0; k < num_equations_; ++k) {
        (*flux_int)[k](q, i) = flux_face(k);
      }
    }
  }
}

void Euler2D_DG_TVDRK::fluxext_compute(const std::vector<Matrix>& u_modal,
                                       const std::vector<Matrix>& Dirichlet,
                                       std::vector<Matrix>* flux_ext)
{
  const int extboundarynum = mesh2D_->getextboundaryNum();
  const int nq1d = fe_->getnumqua1d();
  const auto& ExtBNei = mesh2D_->getextboundaryneighbors();
  const auto& ExtBNormal = mesh2D_->getextboundarynormal();
  const auto& ExtBType = mesh2D_->getextboundarytypeindex();
  const auto& boundary_u = fe_->getboundary_u();

  const real_t max_speed = max_speed_compute(u_modal);
  flux_ext->assign(num_equations_, Matrix::Zero(nq1d, extboundarynum));

  for (int i = 0; i < extboundarynum; ++i) {
    const int cell = ExtBNei[i];
    const int side = ExtBType[i];
    Eigen::Vector2d normal = ExtBNormal[i];

    for (int q = 0; q < nq1d; ++q) {
      Vector u_in(4), u_bc(4), flux_face;
      for (int k = 0; k < num_equations_; ++k) {
        u_in(k) = u_modal[k].col(cell).dot(boundary_u[side].col(q));
        u_bc(k) = Dirichlet[k](q, i);
      }
      numerical_flux(u_in, u_bc, normal, max_speed, &flux_face);
      for (int k = 0; k < num_equations_; ++k) {
        (*flux_ext)[k](q, i) = flux_face(k);
      }
    }
  }
}

void Euler2D_DG_TVDRK::setgamma(const real_t& gamma)
{
  gamma_ = gamma;
}

Euler2D_DG_TVDRK_period::Euler2D_DG_TVDRK_period(const TensorMesh2D* mesh2D,
                                                 const fespace2D* fe,
                                                 const EX_TVDRK* rk_table)
  : Euler2D_DG_TVDRK(mesh2D, fe, rk_table)
{
}

void Euler2D_DG_TVDRK_period::Lu_compute(const std::vector<Matrix>& u_modal,
                                         const real_t& Trun,
                                         const real_t& dt,
                                         std::vector<Matrix>* Lu)
{
  std::vector<Matrix> fu_nodal, fv_nodal;
  flux_compute(u_modal, &fu_nodal, &fv_nodal);
  std::vector<Matrix> rhsx(num_equations_), rhsy(num_equations_);
  fe_->Assemble_F(fu_nodal, 1, 0, &rhsx);
  fe_->Assemble_F(fv_nodal, 0, 1, &rhsy);

  std::vector<Matrix> dirichlet(num_equations_, Matrix::Zero(fe_->getnumqua1d(), mesh2D_->getextboundaryNum()));
  fluxint_compute(u_modal, &flux_int_);
  fluxext_compute(u_modal, dirichlet, &flux_ext_);

  std::vector<Matrix> flux_rhs;
  fe_->Assemble_Flux(flux_int_, flux_ext_, &flux_rhs);

  const DiagnalMatrix& Minv = fe_->getv_u_diaginv();
  *Lu = Minv * (rhsx + rhsy - flux_rhs);
}

void Euler2D_DG_TVDRK_period::fluxext_compute(const std::vector<Matrix>& u_modal,
                                              const std::vector<Matrix>& Dirichlet,
                                              std::vector<Matrix>* flux_ext)
{
  QUEST_VERIFY(mesh2D_->IsPeriodBoundary(), "The mesh is not periodical !");
  const int extboundarynum = mesh2D_->getextboundaryNum();
  const int nq1d = fe_->getnumqua1d();
  const auto& PNei = mesh2D_->getextboundaryneighbors_period();
  const auto& PType = mesh2D_->getextboundarytypeindex_period();
  const auto& ExtBNormal = mesh2D_->getextboundarynormal();
  const auto& boundary_u = fe_->getboundary_u();

  const real_t max_speed = max_speed_compute(u_modal);
  flux_ext->assign(num_equations_, Matrix::Zero(nq1d, extboundarynum));

  for (int i = 0; i < extboundarynum; ++i) {
    const int cIn = PNei[i](0);
    const int cOut = PNei[i](1);
    const int sIn = PType[i](0);
    const int sOut = PType[i](1);
    const Eigen::Vector2d normal = ExtBNormal[i];
    for (int q = 0; q < nq1d; ++q) {
      Vector u_in(4), u_out(4), flux_face;
      for (int k = 0; k < num_equations_; ++k) {
        u_in(k) = u_modal[k].col(cIn).dot(boundary_u[sIn].col(q));
        u_out(k) = u_modal[k].col(cOut).dot(boundary_u[sOut].col(q));
      }
      numerical_flux(u_in, u_out, normal, max_speed, &flux_face);
      for (int k = 0; k < num_equations_; ++k) {
        (*flux_ext)[k](q, i) = flux_face(k);
      }
    }
  }
}

real_t Euler2D_DG_TVDRK_period::rho_bc(const real_t& x, const real_t& y, const real_t& t)
{
  return rho_real(x, y, t);
}

real_t Euler2D_DG_TVDRK_period::vx_bc(const real_t& x, const real_t& y, const real_t& t)
{
  return vx_real(x, y, t);
}

real_t Euler2D_DG_TVDRK_period::vy_bc(const real_t& x, const real_t& y, const real_t& t)
{
  return vy_real(x, y, t);
}

real_t Euler2D_DG_TVDRK_period::pre_bc(const real_t& x, const real_t& y, const real_t& t)
{
  return pre_real(x, y, t);
}

Matrix Euler2D_DG_TVDRK_period::rho_real(const Matrix& x, const Matrix& y, const real_t& t)
{
  return (1.e0 + 0.2e0 * (x.array() + y.array() - 2.e0 * t).sin()).matrix();
}

real_t Euler2D_DG_TVDRK_period::rho_real(const real_t& x, const real_t& y, const real_t& t)
{
  return 1.e0 + 0.2e0 * std::sin(x + y - 2.e0 * t);
}

Matrix Euler2D_DG_TVDRK_period::vx_real(const Matrix& x, const Matrix& y, const real_t& t)
{
  return Matrix::Ones(x.rows(), x.cols());
}

real_t Euler2D_DG_TVDRK_period::vx_real(const real_t& x, const real_t& y, const real_t& t)
{
  return 1.e0;
}

Matrix Euler2D_DG_TVDRK_period::vy_real(const Matrix& x, const Matrix& y, const real_t& t)
{
  return Matrix::Ones(x.rows(), x.cols());
}

real_t Euler2D_DG_TVDRK_period::vy_real(const real_t& x, const real_t& y, const real_t& t)
{
  return 1.e0;
}

Matrix Euler2D_DG_TVDRK_period::pre_real(const Matrix& x, const Matrix& y, const real_t& t)
{
  return Matrix::Ones(x.rows(), x.cols());
}

real_t Euler2D_DG_TVDRK_period::pre_real(const real_t& x, const real_t& y, const real_t& t)
{
  return 1.e0;
}

} // namespace QUEST

#include "Euler1D_acctest.hpp"

namespace QUEST
{
real_t EulerPerodicalBase::cfl_ = 0.8e0;

real_t EulerPerodicalBase::gamma_ = 1.4e0;

real_t EulerPerodicalBase::pi_ = 3.141592653589793238462643383279;

EulerPerodicalBase::EulerPerodicalBase(const fespace1D* fe, const int& t_order)
  : HyperbolicProblems1DBase(fe, t_order) {
};

Matrix EulerPerodicalBase::rho_init(const Matrix& x) {
  Matrix result;
  result = pi_ * x;
  result = 1 + 0.2e0 * result.array().sin();
  return result;
}

Matrix EulerPerodicalBase::vx_init(const Matrix& x) {
  Matrix result;
  result = Matrix::Constant(x.rows(), x.cols(), 1.e0);
  return result;
}

Matrix EulerPerodicalBase::rhovx_init(const Matrix& x) {
  Matrix result = rho_init(x);
  Matrix temp = vx_init(x);
  result = result.array() * temp.array();
  return result;
}

Matrix EulerPerodicalBase::pre_init(const Matrix& x) {
  Matrix result;
  result = Matrix::Constant(x.rows(), x.cols(), 1.e0);
  return result;
};

Matrix EulerPerodicalBase::erg_init(const Matrix& x) {
  Matrix pre = pre_init(x);
  Matrix rho = rho_init(x);
  Matrix vx = vx_init(x);
  Matrix result = pre.array() / (gamma_ - 1.e0) + 0.5e0 * vx.array().square() * rho.array();
  return result;
}

void EulerPerodicalBase::init() {
  fe_->Project_Initial(&rho_init, &u_modal_[0]);
  fe_->Project_Initial(&rhovx_init, &u_modal_[1]);
  fe_->Project_Initial(&erg_init, &u_modal_[2]);

  fe_->modal_to_nodal1D(u_modal_, &u_nodal_);
  vx_nodal_ = u_nodal_[1].array() / u_nodal_[0].array();
  pre_nodal_ = (u_nodal_[2].array() - 0.5e0 * 
    vx_nodal_.array().square() * u_nodal_[0].array()) * (gamma_ - 1.e0);
  std::cout << " Euler 1D has been initialized !" << std::endl; 
};

void EulerPerodicalBase::setdt(real_t* dt) {
  Matrix sound_nodal_ = gamma_ * (pre_nodal_.array() / u_nodal_[0].array());
  sound_nodal_ = sound_nodal_.cwiseSqrt();

  Matrix vxtemp = vx_nodal_.cwiseAbs() + sound_nodal_;
  real_t alphax = vxtemp.maxCoeff();
  max_speed_ = alphax;
  *dt = cfl_ * hx_ / max_speed_ / real_t(2*polydim_ - 1);
};

void EulerPerodicalBase::time_stepping(real_t& dt) {
  std::vector<Matrix> u0 = u_modal_;
  switch (t_order_)
  {
  case 1:
    fk_compute();
    u_modal_ = u_modal_ + dt * fk_;
    update();
    break;
  
  case 2:
    fk_compute();
    u_modal_ = u0 + dt * fk_;
    update();

    fk_compute();
    u_modal_ = 0.5e0 * u0 + 0.5e0 *(u_modal_ + dt * fk_);
    update();
    break;

  case 3:
    fk_compute();
    u_modal_ = u_modal_ + dt * fk_;
    update();

    fk_compute();
    u_modal_ = 0.75e0 * u0 + 0.25e0 * (u_modal_ + dt * fk_);
    update();

    fk_compute();
    u_modal_ = 1.e0/3.e0 * u0 + 2.e0/3.e0 * (u_modal_ + dt * fk_);
    update();
    break;
  default:
    throw std::invalid_argument("Wrong TVD Time stepping order!");
    break;
  }
};

void EulerPerodicalBase::update() {
  fe_->modal_to_nodal1D(u_modal_, &u_nodal_);
  vx_nodal_ = u_nodal_[1].array() / u_nodal_[0].array();
  pre_nodal_ = (u_nodal_[2].array() - 0.5e0 * 
    vx_nodal_.array().square() * u_nodal_[0].array()) * (gamma_ - 1.e0);

  fu_nodal_[0][0] = u_nodal_[1];
  fu_nodal_[0][1] = u_nodal_[1].array().square() / u_nodal_[0].array() + pre_nodal_.array();
  fu_nodal_[0][2] = (u_nodal_[2].array() + pre_nodal_.array()) * vx_nodal_.array();
};

void EulerPerodicalBase::flux_compute (
    const real_t& rho_L, const real_t& vx_L, const real_t& pre_L, const real_t& E_L,
    const real_t& rho_R, const real_t& vx_R, const real_t& pre_R, const real_t& E_R,
    const real_t& normal, 
    Vector* flux) {
  flux->resize(3);
  real_t rhovx_L = rho_L * vx_L;
  real_t rhovx_R = rho_R * vx_R;

  real_t kie_L = rho_L * vx_L * vx_L;
  real_t kie_R = rho_R * vx_R * vx_R;

  real_t phi_L = (E_L + pre_L) * vx_L;
  real_t phi_R = (E_R + pre_R) * vx_R;

  (*flux)(0) = 0.5e0 * (rhovx_L + rhovx_R) * normal - max_speed_ * 0.5e0 * (rho_R - rho_L);
  (*flux)(1) = 0.5e0 * (kie_L + pre_L + kie_R + pre_R)  * normal - max_speed_ * 0.5e0 * (rhovx_R - rhovx_L);
  (*flux)(2) = 0.5e0 * (phi_L + phi_R) * normal - max_speed_ * 0.5e0 * (E_R - E_L);
};

void EulerPerodicalBase::fluxint_compute() {
  // ********************* 传入相关变量 ************************** // 
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const Vector From = boundary_u_[0];
  const Vector To = boundary_u_[1];
  const std::vector<Eigen::Vector2i>& NeiCell_ = fe_->getmesh1D()->getintboundaryneighbors();
  const std::vector<real_t>& NeiNormal_ = fe_->getmesh1D()->getintboundarynormal();
  // ************************************************************ // 

// #ifndef QUEST_USE_OPENMP
//   QUEST_ERROR(" flux int computation procedure must use OPENMP ! ");
// #else
#pragma omp parallel num_threads(NTH_), default(shared)
    {
      int Cellindex0, Cellindex1;
      real_t rho_L, vx_L, rhovx_L, pre_L, E_L;
      real_t rho_R, vx_R, rhovx_R, pre_R, E_R;
      Vector flux; 
      real_t normal;
      int i;
#pragma omp for schedule(static)
      for (i = 0; i < intboundaryNum_; i++) {
        Cellindex0 = NeiCell_[i](0);
        Cellindex1 = NeiCell_[i](1);
        normal = NeiNormal_[i];

        rho_L = u_modal_[0].col(Cellindex0).dot(From);
        rhovx_L = u_modal_[1].col(Cellindex0).dot(From);
        E_L = u_modal_[2].col(Cellindex0).dot(From);
        vx_L = rhovx_L / rho_L;
        pre_L = (E_L - 0.5e0 * (rho_L * vx_L * vx_L)) * (gamma_ - 1.e0);

        rho_R = u_modal_[0].col(Cellindex1).dot(To);
        rhovx_R = u_modal_[1].col(Cellindex1).dot(To);
        E_R = u_modal_[2].col(Cellindex1).dot(To);
        vx_R = rhovx_R / rho_R;
        pre_R = (E_R - 0.5e0 * (rho_R * vx_R * vx_R)) * (gamma_ - 1.e0);
        flux_compute(rho_L, vx_L, pre_L, E_L,
                      rho_R, vx_R, pre_R, E_R,
                      normal, &flux);
        flux_int_[0](0, i) = flux(0);
        flux_int_[1](0, i) = flux(1);
        flux_int_[2](0, i) = flux(2);
      }
    }
// #endif 
}

void EulerPerodicalBase::fluxext_compute() {
  // ********************* 传入相关变量 ************************** // 
  const std::vector<Vector>& boundary_u_ = fe_->getboundary_u();
  const Vector From = boundary_u_[0];
  const Vector To = boundary_u_[1];
  const std::vector<Eigen::Vector2i>& NeiCell_ = fe_->getmesh1D()->getextboundaryneighbors_period();
  const std::vector<real_t>& NeiNormal_ = fe_->getmesh1D()->getextboundarynormal();
  // ************************************************************ // 
// #ifndef QUEST_USE_OPENMP
//   QUEST_ERROR(" flux ext computation procedure must use OPENMP ! ");
// #else
#pragma omp parallel num_threads(NTH_), default(shared)
    {
      int Cellindex0, Cellindex1;
      real_t rho_L, vx_L, rhovx_L, pre_L, E_L;
      real_t rho_R, vx_R, rhovx_R, pre_R, E_R;
      Vector flux; 
      real_t normal;
      int i;
#pragma omp for schedule(static)
      for (i = 0; i < intboundaryNum_; i++) {
        Cellindex0 = NeiCell_[i](0);
        Cellindex1 = NeiCell_[i](1);
        normal = NeiNormal_[i];

        rho_L = u_modal_[0].col(Cellindex0).dot(From);
        rhovx_L = u_modal_[1].col(Cellindex0).dot(From);
        E_L = u_modal_[2].col(Cellindex0).dot(From);
        vx_L = rhovx_L / rho_L;
        pre_L = (E_L - 0.5e0 * (rho_L * vx_L * vx_L)) * (gamma_ - 1.e0);

        rho_R = u_modal_[0].col(Cellindex1).dot(To);
        rhovx_R = u_modal_[1].col(Cellindex1).dot(To);
        E_R = u_modal_[2].col(Cellindex1).dot(To);
        vx_R = rhovx_R / rho_R;
        pre_R = (E_R - 0.5e0 * (rho_R * vx_R * vx_R)) * (gamma_ - 1.e0);
        flux_compute(rho_L, vx_L, pre_L, E_L,
                      rho_R, vx_R, pre_R, E_R,
                      normal, &flux);
        flux_ext_[0](0, i) = flux(0);
        flux_ext_[1](0, i) = flux(1);
        flux_ext_[2](0, i) = flux(2);
      }
    }
// #endif 
}

void EulerPerodicalBase::fk_compute() {
  // ********************* 传入相关变量 ************************** // 
  const DiagnalMatrix& Minv = fe_->getv_u_diaginv();
  // ************************************************************ //
  fluxint_compute();
  fluxext_compute();

  std::vector<Matrix> fx_test_;
  fe_->Assemble_F(fu_nodal_[0], 1, &fx_test_);

  std::vector<Matrix> flux_test_;
  fe_->Assemble_Flux(flux_int_, flux_ext_, &flux_test_);

  fk_ = Minv * (fx_test_ - flux_test_);
};

} // namespace QUEST

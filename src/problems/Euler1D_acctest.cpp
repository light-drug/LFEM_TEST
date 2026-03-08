#include "Euler1D_acctest.hpp"

namespace QUEST
{

Euler1D_DG_TVDRK::Euler1D_DG_TVDRK(const TensorMesh1D* mesh1D,
                                  const fespace1D* fe, 
                                  const EX_TVDRK* rk_table)
  : HyperbolicProblems1DBase(mesh1D, fe, rk_table) {};

Matrix Euler1D_DG_TVDRK::rho_init(const Matrix& x) {
  Matrix result;
  result = pi_ * x;
  result = 1 + 0.2e0 * result.array().sin();
  return result;
}

Matrix Euler1D_DG_TVDRK::vx_init(const Matrix& x) {
  Matrix result;
  result = Matrix::Constant(x.rows(), x.cols(), 1.e0);
  return result;
}

Matrix Euler1D_DG_TVDRK::rhovx_init(const Matrix& x) {
  Matrix result = rho_init(x);
  Matrix temp = vx_init(x);
  result = result.array() * temp.array();
  return result;
}

Matrix Euler1D_DG_TVDRK::pre_init(const Matrix& x) {
  Matrix result;
  result = Matrix::Constant(x.rows(), x.cols(), 1.e0);
  return result;
};

Matrix Euler1D_DG_TVDRK::erg_init(const Matrix& x) {
  Matrix pre = pre_init(x);
  Matrix rho = rho_init(x);
  Matrix vx = vx_init(x);
  Matrix result = pre.array() / (gamma_ - 1.e0) + 0.5e0 * vx.array().square() * rho.array();
  return result;
}

void Euler1D_DG_TVDRK::init() 
{
  fe_->Project_Initial(&rho_init, &u_modal_[0]);
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
  fe_->Project_Initial(erg_nodal, &(u_modal_[2]));
  
};

void Euler1D_DG_TVDRK::setdt(real_t* dt) {
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& k1D = fe_->getbasis()->getk1D();
  real_t max_speed = max_speed_compute(u_modal_);
  *dt = cfl_ * hx_ / max_speed / real_t(2 * k1D - 1);
};

real_t Euler1D_DG_TVDRK::max_speed_compute(const std::vector<Matrix>& u_modal)
{
  std::vector<Matrix> u_nodal;
  fe_->modal_to_nodal1D(u_modal, &u_nodal);
  Matrix pre_nodal = u_nodal[2].array() - 0.5e0 * 
    (u_nodal[1].array().square() / u_nodal[0].array());
  Matrix sound_nodal = gamma_ * (pre_nodal.array() / u_nodal[0].array());
  Matrix vx_nodal = u_nodal[1].array() / u_nodal[0].array();
  real_t max_speed = vx_nodal.cwiseAbs() + sound_nodal;
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
};

void Euler1D_DG_TVDRK::fu_compute(const std::vector<Matrix>& u_modal,
                                  std::vector<Matrix>* fu_nodal) 
{
  std::vector<Matrix> u_nodal;
  fe_->modal_to_nodal1D(u_modal, &u_nodal);
  Matrix kie_nodal = u_nodal[1].array().square() / u_nodal[0].array();
  Matrix pre_nodal = (u_nodal[2].array() - 0.5e0 * kie_nodal) * (gamma_ - 1);
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
  // TODO
};

void Euler1D_DG_TVDRK::numerical_flux(const Vector& u_L,
                                      const Vector& u_R,
                                      const real_t& normal,
                                      const real_t& max_speed,
                                      Vector* flux) 
{
  // 经典Lax-Friedshes
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
  // TODO
}

void Euler1D_DG_TVDRK::fluxext_compute(const std::vector<Matrix>& u_modal,
                                      const std::vector<Matrix>& Dirichlet,
                                      std::vector<Matrix>* flux_ext) 
{
  // TODO
}

void Euler1D_DG_TVDRK::setgamma(const real_t& gamma)
{
  gamma_ = gamma;
};

} // namespace QUEST

#include "KineticLD_DG_IMEX_EX.hpp"

namespace QUEST
{

Kinetic1D_LD_DG_IMEX_EX::Kinetic1D_LD_DG_IMEX_EX(const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table)
    : mesh1D_(mesh1D), fe_(fe), rk_table_(rk_table) {};

void Kinetic1D_LD_DG_IMEX_EX::seteps(const real_t& knu)
{
  eps_ = knu;
  eps2_ = knu * knu;
};

void Kinetic1D_LD_DG_IMEX_EX::setChy(const real_t& Chy)
{
  Chy_ = Chy;
};

void Kinetic1D_LD_DG_IMEX_EX::setCdif(const real_t& Cdif)
{
  Cdif_ = Cdif;
};

void Kinetic1D_LD_DG_IMEX_EX::setCR(const real_t& CR)
{
  CR_ = CR;
};

void Kinetic1D_LD_DG_IMEX_EX::setbeta1(const real_t& beta1)
{
  beta1_ = beta1;
};

void Kinetic1D_LD_DG_IMEX_EX::setsigmas(const real_t& sigmas)
{
  sigmas_ = sigmas;
};

void Kinetic1D_LD_DG_IMEX_EX::setNTH(const int& NTH)
{
  NTH_ = NTH;
};

real_t Kinetic1D_LD_DG_IMEX_EX::rho_numericalbc(const real_t& x, const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  // *** 
  real_t rho = 0.e0;
  if (std::abs(x - x1) < 1.e-11) {
    for (int j = 0; j < Nv; j++) {
      rho += Vweights(j) * fL_bc(j, t, modal);
    }
  } else if (std::abs(x - x2) < 1.e-11) {
    for (int j = 0; j < Nv; j++) {
      rho += Vweights(j) * fR_bc(j, t, modal);
    }
  } else {
    QUEST_ERROR(" This is not the boundary x value ! ");
  }
  return rho;
};

real_t Kinetic1D_LD_DG_IMEX_EX::fL_bc(const int& j,  const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  real_t f;
  if (V(j) >= 0) {
    f = 2.e0;
  } else {
    Vector temp = modal.rho.col(0) + eps_ * modal.g[j].col(0);
    f = temp.dot(boundary_u[1]);
  }
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX::fR_bc(const int& j,  const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  real_t f;
  if (V(j) <= 0) 
  {
    f = 1.e0;
  } else {
    Vector temp = modal.rho.col(ncell - 1) + eps_ * modal.g[j].col(ncell - 1);
    f = temp.dot(boundary_u[0]);
  }
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX::gL_bc(const int& j,  const real_t& t,
      const model_data_& modal, const real_t& rho_L) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t g;
  g = (fL_bc(j, t, modal) - rho_L) / eps_;
  return g;
};

real_t Kinetic1D_LD_DG_IMEX_EX::gR_bc(const int& j,  const real_t& t,
      const model_data_& modal, const real_t& rho_R) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t g;
  g = (fR_bc(j, t, modal) - rho_R) / eps_;
  return g;
};

Matrix Kinetic1D_LD_DG_IMEX_EX::rho_init(const Matrix& x) 
{
  Matrix rho = x;
  rho.setConstant(1.e0);
  rho.array() = (x.array() < 0.0).select(2.0, rho.array());
  return rho;
};

real_t Kinetic1D_LD_DG_IMEX_EX::rho_init(const real_t& x) 
{
  real_t rho;
  if (x < 0.e0) {
    rho = 2.e0;
  } else {
    rho = 1.e0;
  }
  return rho;
};

Matrix Kinetic1D_LD_DG_IMEX_EX::f_init(const Matrix& x, const real_t& v) 
{
  Matrix f = x;
  f.setConstant(1.e0);
  f.array() = (x.array() < 0.0).select(2.0, f.array());
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX::f_init(const real_t& x, const real_t& v) 
{
  real_t f;
  if (x < 0.e0) {
    f = 2.e0;
  } else {
    f = 1.e0;
  }
  return f;
};

Matrix Kinetic1D_LD_DG_IMEX_EX::g_init(const Matrix& x, const real_t& v) 
{
  Matrix g = Matrix::Zero(x.rows(), x.cols());
  return g;
};

real_t Kinetic1D_LD_DG_IMEX_EX::g_init(const real_t& x, const real_t& v) 
{
  real_t g = 0.e0;
  return g;
};

void Kinetic1D_LD_DG_IMEX_EX::setdt(real_t* dt) 
{
  // ** 传入相关变量
  const real_t& hx = mesh1D_->gethx();
  // *** 
  *dt = Chy_ * eps_ * hx + Cdif_ * hx * hx;
};

void Kinetic1D_LD_DG_IMEX_EX::init() 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  // *** 
  QUEST_VERIFY(Nv == 2, " This is the two velocity Riemann problem test !");
  
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();

  ah_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_);
  // rho_stages_.resize(stages_);
  // g_stages_.resize(stages_);  

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  }
  kinetic_modal_stages_.resize(stages_);
};

const Matrix& Kinetic1D_LD_DG_IMEX_EX::getrho_modal() const 
{
  return kinetic_modal_.rho;
};

const std::vector<Matrix>& Kinetic1D_LD_DG_IMEX_EX::getg_modal() const 
{
  return kinetic_modal_.g;
};

void Kinetic1D_LD_DG_IMEX_EX::fluxint_compute(const Matrix& modal, 
                          const real_t& beta,
                          Matrix* flux_int)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& intboundarynum = fe_->getmesh1D()->getintboundaryNum();
  // const int& extboundarynum = fe_->getmesh1D()->getextboundarynum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  // const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  // const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  
  *flux_int = Matrix::Zero(1, intboundarynum);
  int Cellindex0, Cellindex1;
  Vector temp0, temp1;
  real_t val0, val1;
  real_t normal;
  for (int i = 0; i < intboundarynum; i++)
  {
    Cellindex0 = IntBNei[i](0);
    Cellindex1 = IntBNei[i](1);
    temp0 = modal.col(Cellindex0).transpose();
    temp1 = modal.col(Cellindex1).transpose();
    val0 = temp0.dot(boundary_u[0]);
    val1 = temp1.dot(boundary_u[1]);
    normal = IntBNormal[i];
    (*flux_int)(0, i) = (val0 + val1) / 2.e0 + beta * (val0 - val1) * normal;
  };
  
};

void Kinetic1D_LD_DG_IMEX_EX::fluxext_compute(const Matrix& modal, 
                          const real_t& beta,
                          const Matrix& Dirichlet,
                          Matrix* flux_ext)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  // const int& intboundarynum = fe_->getmesh1D()->getintboundarynum();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  // const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  // const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  *flux_ext = Matrix::Zero(1, extboundarynum);
  int Cellindex0;
  Vector temp0, temp1;
  real_t val0, val1;
  real_t normal;
  for (int i = 0; i < extboundarynum; i++)
  {
    Cellindex0 = ExtBNei[i];
    temp0 = modal.col(Cellindex0).transpose();
    val0 = temp0.dot(boundary_u[1 - i]);
    val1 = Dirichlet(0, i);
    normal = ExtBNormal[i];
    (*flux_ext)(0, i) = (val0 + val1) / 2.e0 + beta * (val0 - val1) * normal;
  };
};

void Kinetic1D_LD_DG_IMEX_EX::fluxint_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          Matrix* flux_int)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& intboundarynum = fe_->getmesh1D()->getintboundaryNum();
  // const int& extboundarynum = fe_->getmesh1D()->getextboundarynum();
  const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  // const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  // const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  
  *flux_int = Matrix::Zero(1, intboundarynum);
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal;
  for (int i = 0; i < intboundarynum; i++)
  {
    Cellindex0 = IntBNei[i](0);
    Cellindex1 = IntBNei[i](1);
    val0 = a * modal.col(Cellindex0).dot(boundary_u[0]);
    val1 = a * modal.col(Cellindex1).dot(boundary_u[1]);
    normal = IntBNormal[i];
    (*flux_int)(0, i) = (a * normal > 0) ? val0 : val1;
  };
};

void Kinetic1D_LD_DG_IMEX_EX::fluxext_upwind_compute(const Matrix& modal, 
                          const real_t& a,
                          const Matrix& Dirichlet,
                          Matrix* flux_ext)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  // const int& intboundarynum = fe_->getmesh1D()->getintboundarynum();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  // const std::vector<Eigen::Vector2i>& IntBNei = mesh1D_->getintboundaryneighbors();
  // const std::vector<real_t>& IntBNormal = mesh1D_->getintboundarynormal();
  const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  
  *flux_ext = Matrix::Zero(1, extboundarynum);
  int Cellindex0, Cellindex1;
  Vector temp0, temp1;
  real_t val0, val1;
  real_t normal;
  for (int i = 0; i < extboundarynum; i++)
  {
    Cellindex0 = ExtBNei[i];
    temp0 = modal.col(Cellindex0).transpose();
    val0 = a * temp0.dot(boundary_u[1 - i]);
    val1 = a * Dirichlet(0, i);
    normal = ExtBNormal[i];
    (*flux_ext)(0, i) = (a * normal > 0) ? val0 : val1;
  };
};

void Kinetic1D_LD_DG_IMEX_EX::ah_compute(const model_data_& modal, 
        const real_t& Trun, 
        const Matrix& boundary_flux,
        Matrix* ah)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  Matrix vg_modal = Matrix::Zero(polydim, ncell);
  for (int j = 0; j < Nv; j++) {
    vg_modal += vweights(j) * V(j) * modal.g[j];
  };
  Matrix vg_nodal;
  fe_->modal_to_nodal1D(vg_modal, &vg_nodal);
  Matrix rhs;
  fe_->Assemble_F(vg_nodal, 1, &rhs);
  Matrix flux_int;
  fluxint_compute(vg_modal, beta1_, &flux_int);
  // fluxext_compute(vg_modal, beta1_, boundary_flux, &flux_ext);  // 交错通量 beta
  Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  fe_->Assemble_Flux(flux_int, boundary_flux, &flux_rhs);
  *ah = - rhs + flux_rhs;
  // std::cout << "boundary_flux = " << boundary_flux << std::endl;
  // std::cout << "flux_int = " << flux_int << std::endl;
  // std::cout << "flux_rhs = " << flux_rhs << std::endl;
};

void Kinetic1D_LD_DG_IMEX_EX::ah_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux)
{ 
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const real_t& hx = mesh1D_->gethx();
  // *** 
  Matrix vg_modal = Matrix::Zero(polydim, ncell);
  for (int j = 0; j < Nv; j++) {
    vg_modal += vweights(j) * V(j) * modal.g[j];
  };
  boundary_flux->resize(1, 2);
  boundary_flux->setZero();
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal;

  // Peng Schur中提到的边界取法
  Cellindex0 = ExtBNei[0];
  normal = ExtBNormal[0];
  val0 = vg_modal.col(Cellindex0).dot(boundary_u[1]);
  real_t rho_L_in = modal.rho.col(Cellindex0).dot(boundary_u[1]);
  real_t rho_L = rho_numericalbc(x1, Trun, modal);
  (*boundary_flux)(0, 0) = val0 + CR_ * (rho_L_in - rho_L) * normal;
  // 
  Cellindex0 = ExtBNei[1];
  normal = ExtBNormal[1];
  val0 = vg_modal.col(Cellindex0).dot(boundary_u[0]);
  real_t rho_R_in = modal.rho.col(Cellindex0).dot(boundary_u[0]);
  real_t rho_R = rho_numericalbc(x2, Trun, modal);
  (*boundary_flux)(0, 1) = val0 + CR_ * (rho_R_in - rho_R) * normal;
  // 
};


void Kinetic1D_LD_DG_IMEX_EX::bh_compute(const model_data_& modal, 
                          const real_t& Trun,
                          const std::vector<Matrix>& boundary_flux,
                          std::vector<Matrix>* bh) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  std::vector<Matrix> g_nodal;
  fe_->modal_to_nodal1D(modal.g, &g_nodal);
  Matrix rhs;
  Matrix flux_int, flux_ext;
  Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  Matrix dh = Matrix::Zero(polydim, ncell);
  Matrix temp;
  bh->resize(Nv);
  for (int j = 0; j < Nv; j++) {
    temp = V(j) * g_nodal[j];
    fe_->Assemble_F(temp, 1, &rhs);
    fluxint_upwind_compute(modal.g[j], V(j), &flux_int);
    // fluxext_upwind_compute(modal.g[j], V(j), boundary_flux[j], &flux_ext);
    fe_->Assemble_Flux(flux_int, boundary_flux[j], &flux_rhs);
    (*bh)[j] = (- rhs + flux_rhs);
    dh = dh + vweights(j) * (*bh)[j];  // \Pi(vg_x)
  };
  // (I - \Pi)(vg_x)
  for (int j = 0; j < Nv; j++)
  {
    (*bh)[j] = (*bh)[j] - dh;
  };
  // std::cout << "vweights = " << vweights.transpose() << std::endl; 
  // PAUSE();
};

void Kinetic1D_LD_DG_IMEX_EX::bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  const std::vector<int>& ExtBNei =  mesh1D_->getextboundaryneighbors();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  boundary_flux->resize(Nv);
  real_t rho_L = rho_numericalbc(x1, Trun, modal);
  real_t rho_R = rho_numericalbc(x2, Trun, modal);

  int Cellindex0;
  real_t temp0, temp1;
  for (int j = 0; j < Nv; j++) 
  { 
    boundary_flux->at(j).resize(1, 2);
    boundary_flux->at(j).setZero();
    if (V(j) >= 0)
    {
      Cellindex0 = ExtBNei[1];
      temp0 = V(j) * modal.g[j].col(Cellindex0).dot(boundary_u[0]);
      boundary_flux->at(j)(0, 0) = V(j) * gL_bc(j, Trun, modal, rho_L);
      boundary_flux->at(j)(0, 1) = temp0; 
    } else if (V(j) < 0)
    {
      Cellindex0 = ExtBNei[0];
      temp0 = V(j) * modal.g[j].col(Cellindex0).dot(boundary_u[1]);
      boundary_flux->at(j)(0, 0) = temp0;
      boundary_flux->at(j)(0, 1) = V(j) * gR_bc(j, Trun, modal, rho_R);
    };
  };
};

void Kinetic1D_LD_DG_IMEX_EX::dh_compute(const model_data_& modal,  
                          const real_t& Trun, 
                          const Matrix& boundary_flux,
                          Matrix* dh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  
  Matrix rho_nodal;
  fe_->modal_to_nodal1D(modal.rho, &rho_nodal);
  Matrix rhs;
  fe_->Assemble_F(rho_nodal, 1, &rhs);
  Matrix flux_int;
  fluxint_compute(modal.rho, - beta1_, &flux_int);
  // fluxext_compute(modal.rho, - beta1_, boundary_flux, &flux_ext); // 交错通量 - beta_
  Matrix flux_rhs = Matrix::Zero(polydim, ncell);
  fe_->Assemble_Flux(flux_int, boundary_flux, &flux_rhs);
  
  *dh = rhs - flux_rhs;
};

void Kinetic1D_LD_DG_IMEX_EX::dh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux)
{ 
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  // *** 
  boundary_flux->resize(1, 2);
  boundary_flux->setZero();

  // Peng Schur中提到的边界取法
  real_t rho_L = rho_numericalbc(x1, Trun, modal);
  real_t rho_R = rho_numericalbc(x2, Trun, modal);

  (*boundary_flux)(0, 0) = rho_L;
  (*boundary_flux)(0, 1) = rho_R;
};

void Kinetic1D_LD_DG_IMEX_EX::sh_compute(const model_data_& modal,  
                          const real_t& Trun, 
                          std::vector<Matrix>* sh)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  // *** 
  
  std::vector<Matrix> g_nodal;
  fe_->modal_to_nodal1D(modal.g, &g_nodal);
  std::vector<Matrix> rhs;
  fe_->Assemble_F(g_nodal, 0, &rhs);
  *sh = rhs;
};

void Kinetic1D_LD_DG_IMEX_EX::updateAll(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Matrix& vweights = mesh1D_->getvweights();
  const Matrix& Ae = rk_table_->getA();
  const Matrix& Ai = rk_table_->getAi();
  const Vector& be = rk_table_->getb();
  const Vector& bi = rk_table_->getbi();
  const Vector& ce = rk_table_->getc();
  const Vector& ci = rk_table_->getci();
  const DiagnalMatrix& M = fe_->getv_u_diag();
  const DiagnalMatrix& Minv = fe_->getv_u_diaginv();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  // *** 
  Matrix rhoboundary_flux, vgboundary_flux;
  std::vector<Matrix> gboundary_flux;
  real_t temp0, temp1;
  for (int s = 0; s < stages_; s++)
  { 
    kinetic_modal_stages_[s].rho = M * kinetic_modal_.rho;
    kinetic_modal_stages_[s].g = eps2_ * M * kinetic_modal_.g;
    temp0 = Ai(s, s) * dt;
    temp1 = eps2_ + sigmas_ * Ai(s, s) * dt;
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - dt * Ae(s, i) * ah_[i];
      for (int j = 0; j < Nv; j++)
      {
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j] 
                                      - (dt * Ae(s, i) * eps_)  * bh_[i][j]
                                      + (dt * Ai(s, i) * V(j)) * dh_[i]
                                      - (dt * sigmas_ * Ai(s, i)) * sh_[i][j];
      };
    };
    kinetic_modal_stages_[s].rho = Minv * kinetic_modal_stages_[s].rho;
    if (s == 0)
    {
      dh_extflux_compute(kinetic_modal_, Trun, &rhoboundary_flux);
    } else 
    {
      dh_extflux_compute(kinetic_modal_stages_[s-1], Trun + ci(s-1) * dt, &rhoboundary_flux);
    }
    dh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, rhoboundary_flux, &(dh_[s]));

    for (int j = 0; j < Nv; j++)
    {
      kinetic_modal_stages_[s].g[j] = (1.e0 / temp1) * kinetic_modal_stages_[s].g[j] + (temp0 / temp1) * V(j) * dh_[s];
    };
    kinetic_modal_stages_[s].g = Minv * kinetic_modal_stages_[s].g;
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));

    ah_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &vgboundary_flux);
    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));

    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
  };
  // PAUSE();
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
};

Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::Kinetic1D_LD_DG_IMEX_EX_IsentropicBC(
                      const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table) {};

void Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::init() 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const real_t& v1 = mesh1D_->getv1();
  const real_t& v2 = mesh1D_->getv2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  // ***   
  QUEST_VERIFY(std::abs(v1 + 1.e0) < 1.e-12, " v1 must be equal to - 1.e0 !");
  QUEST_VERIFY(std::abs(v2 - 1.e0) < 1.e-12, " v2 must be equal to 1.e0 !");
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();

  ah_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_);
  // rho_stages_.resize(stages_);
  // g_stages_.resize(stages_);  

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  }
  kinetic_modal_stages_.resize(stages_);
};


Matrix Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::rho_init(const Matrix& x) 
{
  Matrix rho = x;
  rho.setZero();
  return rho;
};

real_t Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::rho_init(const real_t& x) 
{
  return 0.e0;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::f_init(const Matrix& x, const real_t& v) 
{
  Matrix f = x;
  f.setZero();
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::f_init(const real_t& x, const real_t& v) 
{
  return 0.e0;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::g_init(const Matrix& x, const real_t& v) 
{
  Matrix g = Matrix::Zero(x.rows(), x.cols());
  return g;
};

real_t Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::g_init(const real_t& x, const real_t& v) 
{
  real_t g = 0.e0;
  return g;
};

real_t Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::fL_bc(const int& j,  const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  real_t f;
  if (V(j) >= 0) {
    f = 1.e0;
  } else {
    Vector temp = modal.rho.col(0) + eps_ * modal.g[j].col(0);
    f = temp.dot(boundary_u[1]);
  }
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX_IsentropicBC::fR_bc(const int& j,  const real_t& t,
      const model_data_& modal) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& Vweights = mesh1D_->getvweights();
  const Vector& V = mesh1D_->getV();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  real_t f;
  if (V(j) <= 0) 
  {
    f = 0.e0;
  } else {
    Vector temp = modal.rho.col(ncell - 1) + eps_ * modal.g[j].col(ncell - 1);
    f = temp.dot(boundary_u[0]);
  }
  return f;
};

Kinetic1D_LD_DG_IMEX_EX_twovel_period::Kinetic1D_LD_DG_IMEX_EX_twovel_period(
                      const KineticTensorMesh1D* mesh1D,
                      const fespace1D* fe,
                      const IMEX_RK* rk_table)
  : Kinetic1D_LD_DG_IMEX_EX(mesh1D, fe, rk_table) {};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::init() 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const BoundaryType boundary_type = mesh1D_->getboundaryType();
  // *** 
  pi_ = 3.14159265358979323846264338327;
  stages_ = rk_table_->getstages();
  QUEST_VERIFY(eps_ <= 0.5e0, " eps (Knudsen number has to be less than 0.5 !)");
  QUEST_VERIFY(boundary_type == BoundaryType::PeriodBoundary, " must be periodical boundary condition !");
  r_ = 2.e0 / (1.e0 + std::sqrt(1.e0 - 4.e0 * eps2_));

  ah_.resize(stages_);
  bh_.resize(stages_);
  dh_.resize(stages_);
  sh_.resize(stages_); 

  fe_->Project_Initial(
      [this](const Matrix& x) { return this->rho_init(x); }, &(kinetic_modal_.rho));
  kinetic_modal_.g.resize(Nv);
  real_t vj;
  for (int j = 0; j < Nv; j++) {
    vj = V(j);
    fe_->Project_Initial(
      [this, vj](const Matrix& x) { return this->g_init(x, vj); }, &(kinetic_modal_.g[j]));
  };
  kinetic_modal_stages_.resize(stages_);
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::updateAll(const real_t& Trun, const real_t& dt) 
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Matrix& vweights = mesh1D_->getvweights();
  const Matrix& Ae = rk_table_->getA();
  const Matrix& Ai = rk_table_->getAi();
  const Vector& be = rk_table_->getb();
  const Vector& bi = rk_table_->getbi();
  const Vector& ce = rk_table_->getc();
  const Vector& ci = rk_table_->getci();
  const DiagnalMatrix& M = fe_->getv_u_diag();
  const DiagnalMatrix& Minv = fe_->getv_u_diaginv();
  const Vector& JacobiDet = fe_->getmesh1D()->getJacobiDet();
  // *** 
  Matrix rhoboundary_flux, vgboundary_flux;
  std::vector<Matrix> gboundary_flux;
  real_t temp0, temp1;
  for (int s = 0; s < stages_; s++)
  { 
    kinetic_modal_stages_[s].rho = M * kinetic_modal_.rho;
    kinetic_modal_stages_[s].g = eps2_ * M * kinetic_modal_.g;
    temp0 = Ai(s, s) * dt;
    temp1 = eps2_ + sigmas_ * Ai(s, s) * dt;
    for (int i = 0; i < s; i++)
    {
      kinetic_modal_stages_[s].rho = kinetic_modal_stages_[s].rho 
                                    - dt * Ae(s, i) * ah_[i];
      for (int j = 0; j < Nv; j++)
      {
        kinetic_modal_stages_[s].g[j] = kinetic_modal_stages_[s].g[j] 
                                      - (dt * Ae(s, i) * eps_)  * bh_[i][j]
                                      + (dt * Ai(s, i) * V(j)) * dh_[i]
                                      - (dt * sigmas_ * Ai(s, i)) * sh_[i][j];
      };
    };
    kinetic_modal_stages_[s].rho = Minv * kinetic_modal_stages_[s].rho;
    
    dh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &rhoboundary_flux);
    dh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, rhoboundary_flux, &(dh_[s]));
    // std::cout << "dh_[" << s << "] = \n" << dh_[s] << std::endl;
    for (int j = 0; j < Nv; j++)
    {
      kinetic_modal_stages_[s].g[j] = (1.e0 / temp1) * kinetic_modal_stages_[s].g[j] + (temp0 / temp1) * V(j) * dh_[s];
    };
    kinetic_modal_stages_[s].g = Minv * kinetic_modal_stages_[s].g;
    
    if(s == stages_ - 1) break;
    sh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &(sh_[s]));

    ah_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &vgboundary_flux);
    ah_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, vgboundary_flux, &(ah_[s]));
    // std::cout << "ah_[" << s << "] = \n" << ah_[s] << std::endl;
    bh_extflux_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, &gboundary_flux);
    bh_compute(kinetic_modal_stages_[s], Trun + ci(s) * dt, gboundary_flux, &(bh_[s]));
  };
  // PAUSE();
  kinetic_modal_.rho = kinetic_modal_stages_[stages_ - 1].rho;
  kinetic_modal_.g = kinetic_modal_stages_[stages_ - 1].g;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_init(const Matrix& x)
{
  Matrix rho = - 1.e0 / (4.e0 * r_) * x.array().sin() + 0.5e0;
  return rho;
};

real_t Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_init(const real_t& x) 
{
  real_t rho = - 1.e0 / (4.e0 * r_) * std::sin(x) + 0.5e0;
  return rho;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_init(const Matrix& x, const real_t& v) 
{
  Matrix f = rho_init(x) + eps_ * g_init(x, v);
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_init(const real_t& x, const real_t& v) 
{
  real_t f = rho_init(x) + eps_ * g_init(x, v);
  return f;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_init(const Matrix& x, const real_t& v) 
{
  Matrix g = 1.e0 / 4.e0 * v * x.array().cos();
  return g;
};

real_t Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_init(const real_t& x, const real_t& v) 
{
  real_t g = 1.e0 / 4.e0 * v * std::cos(x);
  return g;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_real(const Matrix& x, const real_t& t)
{
  Matrix rho = - 1.e0 / (4.e0 * r_) * x.array().sin() * std::exp(- r_ * t) + 0.5e0;
  return rho;
};

real_t Kinetic1D_LD_DG_IMEX_EX_twovel_period::rho_real(const real_t& x, const real_t& t)
{
  real_t rho = - 1.e0 / (4.e0 * r_) * std::sin(x) * std::exp(- r_ * t) + 0.5e0;
  return rho;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_real(const Matrix& x, const real_t& v, const real_t& t)
{
  Matrix f = rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

real_t Kinetic1D_LD_DG_IMEX_EX_twovel_period::f_real(const real_t& x, const real_t& v, const real_t& t)
{
  real_t f = rho_real(x, t) + eps_ * g_real(x, v, t);
  return f;
};

Matrix Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_real(const Matrix& x, const real_t& v, const real_t& t)
{
  Matrix g = 1.e0 / 4.e0 * v * x.array().cos() * std::exp(- r_ * t);
  return g;
};

real_t Kinetic1D_LD_DG_IMEX_EX_twovel_period::g_real(const real_t& x, const real_t& v, const real_t& t)
{
  real_t g = 1.e0 / 4.e0 * v * std::cos(x) * std::exp(- r_ * t);
  return g;
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::ah_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux)
{ 
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  const std::vector<Eigen::Vector2i>& ExtBNei_period =  mesh1D_->getextboundaryneighbors_period();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  Matrix vg_modal = Matrix::Zero(polydim, ncell);
  for (int j = 0; j < Nv; j++) {
    vg_modal += vweights(j) * V(j) * modal.g[j];
  };
  boundary_flux->resize(1, 2);
  boundary_flux->setZero();
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal;

  // 周期边界
  Cellindex0 = ExtBNei_period[0](0);
  Cellindex1 = ExtBNei_period[0](1);
  val0 = vg_modal.col(Cellindex0).dot(boundary_u[1]);
  val1 = vg_modal.col(Cellindex1).dot(boundary_u[0]);
  normal = ExtBNormal[0];
  // beta = 0.5 flux = val1"-";
  // std::cout << "val0 = " << val0 << ", Cellindex0 = " << Cellindex0 << std::endl;
  // std::cout << "val1 = " << val1 << ", Cellindex1 = " << Cellindex1 << std::endl;
  (*boundary_flux)(0, 0) = (val0 + val1) / 2.e0 + beta1_ * (val0 - val1) * normal;
  // std::cout << " boundary_flux = " <<(*boundary_flux)(0, 0) << std::endl;
  // PAUSE();
  // Cellindex0 = ExtBNei_period[1](0);
  // Cellindex1 = ExtBNei_period[1](1);
  // val0 = vg_modal.col(Cellindex0).transpose().dot(boundary_u[0]);
  // val1 = vg_modal.col(Cellindex1).transpose().dot(boundary_u[1]);
  // normal = ExtBNormal[1];
  // (*boundary_flux)(0, 1) = (val0 + val1) / 2.e0 + beta1_ * (val0 - val1) * normal;
  (*boundary_flux)(0, 1) = (*boundary_flux)(0, 0);
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::bh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun,
                          std::vector<Matrix>* boundary_flux)
{
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const int& ncell = mesh1D_->getncell();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  const Vector& vweights = mesh1D_->getvweights();
  const int& extboundarynum = fe_->getmesh1D()->getextboundaryNum();
  const std::vector<Eigen::Vector2i>& ExtBNei_period =  mesh1D_->getextboundaryneighbors_period();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  // *** 
  boundary_flux->resize(Nv);
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal = ExtBNormal[0];
  real_t v;
  for (int j = 0; j < Nv; j++) 
  { 
    v = V(j);
    boundary_flux->at(j).setZero();
    boundary_flux->at(j).resize(1, extboundarynum);
    Cellindex0 = ExtBNei_period[0](0);
    Cellindex1 = ExtBNei_period[0](1);
    val0 = v * modal.g[j].col(Cellindex0).dot(boundary_u[1]);
    val1 = v * modal.g[j].col(Cellindex1).dot(boundary_u[0]);
    // std::cout << "g[" << j << "] = " << modal.g[j] << std::endl;
    // std::cout << " val0 = " << val0 << std::endl;
    // std::cout << " val1 = " << val1 << std::endl;
    boundary_flux->at(j)(0, 0) = (v * normal > 0) ? val0 : val1;
    // std::cout << " boundary_flux->at(j)(0, 0) = " << boundary_flux->at(j)(0, 0) << std::endl;
    boundary_flux->at(j)(0, 1) = boundary_flux->at(j)(0, 0);
  };
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::dh_extflux_compute(const model_data_& modal, 
                          const real_t& Trun, 
                          Matrix* boundary_flux)
{ 
  // ** 传入相关变量
  const real_t& x1 = mesh1D_->getx1();
  const real_t& x2 = mesh1D_->getx2();
  const Matrix& rhotemp = modal.rho;
  const std::vector<Vector>& boundary_u = fe_->getboundary_u();
  const std::vector<Eigen::Vector2i>& ExtBNei_period =  mesh1D_->getextboundaryneighbors_period();
  const std::vector<real_t>& ExtBNormal = mesh1D_->getextboundarynormal();
  // *** 
  boundary_flux->resize(1, 2);
  boundary_flux->setZero();
  int Cellindex0, Cellindex1;
  real_t val0, val1;
  real_t normal;

  Cellindex0 = ExtBNei_period[0](0);
  Cellindex1 = ExtBNei_period[0](1);
  val0 = rhotemp.col(Cellindex0).dot(boundary_u[1]);
  val1 = rhotemp.col(Cellindex1).dot(boundary_u[0]);
  normal = ExtBNormal[0];
  (*boundary_flux)(0, 0) = (val0 + val1) / 2.e0 - beta1_ * (val0 - val1) * normal;

  (*boundary_flux)(0, 1) = (*boundary_flux)(0, 0);
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::getrho_real_modal(const real_t& Tstop,
                            Matrix* rho_real_modal)
{
  fe_->Project_Final(
      [this](const Matrix& x, const real_t& t) { return this->rho_real(x, t); }, 
      Tstop, 
      rho_real_modal);
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::getrho_real_nodal(const real_t& Tstop,
                            Matrix* rho_real_nodal)
{
  fe_->Interpolate_Final(
      [this](const Matrix& x, const real_t& t) { return this->rho_real(x, t); },
      Tstop, 
      rho_real_nodal);
};

void Kinetic1D_LD_DG_IMEX_EX_twovel_period::getg_real_nodal(const real_t& Tstop, 
                      std::vector<Matrix>* g_real_nodal)
{
  // ** 传入相关变量
  const int& Nv = mesh1D_->getNv();
  const Vector& V = mesh1D_->getV();
  // *** 
  real_t vj;
  g_real_nodal->resize(Nv);
  for (int j = 0; j < Nv; j++)
  {
    vj = V(j);
    fe_->Interpolate_Final(
      [this, vj](const Matrix& x, const real_t& t) { return this->g_real(x, vj, t); },
      Tstop,
      &(g_real_nodal->at(j)));
  };
};

} // namespace QUEST

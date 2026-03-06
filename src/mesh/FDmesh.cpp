#include "FDmesh.hpp"

namespace QUEST
{

FDmesh::FDmesh(const real_t& x1, const real_t& x2,
        const int& xDiv)
  : x1_(x1), x2_(x2), xDiv_(xDiv) {};

void FDmesh::init() {
  generatehx();
  generatexpoints();
  generateCellCenter(); 
};

const real_t& FDmesh::getx1() const 
{
  return x1_;
};

const real_t& FDmesh::getx2() const 
{
  return x2_;
};

const int& FDmesh::getxDiv() const 
{
  return xDiv_;
};

const real_t& FDmesh::gethx() const 
{
  return hx_;
};

const int& FDmesh::getNx() const 
{
  return Nx_;
};

const Vector& FDmesh::getxpoints() const 
{
  return xpoints_;
};

const Vector& FDmesh::getX() const 
{
  return X_;
};

const Vector& FDmesh::getCellCenter_vec() const
{
  return CellCenter_vec_;
};

const Matrix& FDmesh::getCellCenter() const 
{
  return CellCenter_;
};

void FDmesh::generatehx() {
  hx_ = (x2_ - x1_) / real_t(xDiv_);
  Nx_ = xDiv_;
};

void FDmesh::generatexpoints() {
  Nx_ = xDiv_ - 1;
  xpoints_.resize(xDiv_ + 1);
  for (int i = 0; i <= xDiv_; ++i) {
    xpoints_(i) = i * hx_ + x1_;
  }
  X_.resize(Nx_);
  for (int i = 0; i < Nx_; ++i) {
    X_(i) = (i + 1) * hx_ + x1_;
  }
};

void FDmesh::generateCellCenter() 
{
  CellCenter_.resize(1, xDiv_);
  for (int i = 0; i < xDiv_; ++i) {
    CellCenter_(0, i) = (real_t(i) + 0.5e0) * hx_ + x1_;
  }
  CellCenter_vec_ = CellCenter_.row(0).transpose();
};

void FDmesh::DisplayResult_withX(const Vector& xtemp,
    const Vector& numerical, const Vector& exact,  
    const std::string& title, std::ostream& Outfile) const 
{
  QUEST_VERIFY(numerical.size() == exact.size(), 
              " numerical and exact must have the same mesh size !");
  QUEST_VERIFY(numerical.size() == xtemp.size(), 
              " x and numerical must have the same mesh size !");
  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = \"X\", \"Numerical Solution\", \"Exact Solution\"" << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << Nx_ << " F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < xtemp.size(); i++) {
    real_t x0 = xtemp(i);
    Outfile << x0 << " " 
            << numerical(i) << " "
            << exact(i) << " " << std::endl;
  };
};


void FDmesh::DisplayResult(const Vector& numerical, const Vector& exact,  
    const std::string& title, std::ostream& Outfile) const 
{
  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = \"X\", \"Numerical Solution\", \"Exact Solution\"" << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << Nx_ << " F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < Nx_; i++) {
    real_t x0 = X_(i);
    Outfile << x0 << " " 
            << numerical(i) << " "
            << exact(i) << " " << std::endl;
  };
};

void FDmesh::DisplayResult_ex(const Vector& numerical,   
    const std::string& title, std::ostream& Outfile) const {

  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = \"X\", \"Numerical Solution\", \"Exact Solution\"" << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << Nx_ << " F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < xDiv_ + 1; i++) {
    real_t x0 = xpoints_(i);
    Outfile << x0 << " " 
            << numerical(i) << " " << std::endl;
  };
};

void FDmesh::computerrorL1(const Vector& numerical, const Vector& exact, 
                          real_t* error) {
  Vector err = numerical - exact;
  err = err.cwiseAbs();
  *error = err.sum() * hx_;
};

void FDmesh::computerrorL2(const Vector& numerical, const Vector& exact, 
                          real_t* error) {
  Vector err = numerical - exact;
  err = err.cwiseAbs2();
  *error = err.sum() * hx_;
  *error = std::sqrt(*error);
};

void FDmesh::computerrorLinf(const Vector& numerical, const Vector& exact, 
                          real_t* error) {
  Vector err = numerical - exact;
  err = err.cwiseAbs();
  *error = err.maxCoeff();
};

FDmesh_period::FDmesh_period(const real_t& x1, const real_t& x2,
        const int& xDiv) : 
    FDmesh(x1, x2, xDiv) {};


void FDmesh_period::generatexpoints() 
{
  Nx_ = xDiv_;
  xpoints_.resize(xDiv_ + 1);
  for (int i = 0; i <= xDiv_; ++i) {
    xpoints_(i) = i * hx_ + x1_;
  }
  X_.resize(Nx_);
  for (int i = 0; i < Nx_; ++i) {
    X_(i) = (i) * hx_ + x1_;
  }
}


KineticFDmesh::KineticFDmesh(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv)
  : FDmesh(x1, x2, xDiv), v1_(v1), v2_(v2), vDiv_(vDiv) {};

void KineticFDmesh::init() 
{
  generatehx();
  generatehv();
  generatexpoints();
  generatevpoints();
  generatevweights();
  generateCellCenter(); 
};

const real_t& KineticFDmesh::getv1() const {
  return v1_;
};

const real_t& KineticFDmesh::getv2() const {
  return v2_;
};

const real_t& KineticFDmesh::gethv() const {
  return hv_;
};

const int& KineticFDmesh::getvDiv() const {
  return vDiv_;
};

const int& KineticFDmesh::getNv() const {
  return Nv_;
};

const Vector& KineticFDmesh::getV() const {
  return V_;
};

const Vector& KineticFDmesh::getvpoints() const {
  return vpoints_;
};

const Vector& KineticFDmesh::getvweights() const {
  return vweights_;
};

void KineticFDmesh::DisplayResult(const std::vector<Vector>& numerical, 
    const std::vector<Vector>& exact,
    const std::string& title, std::ostream& Outfile) const {

  Outfile <<  "TITLE = \" "<< title << " \" " << std::endl;
  Outfile <<  "VARIABLES = \"X\", \"V\", \"Numerical Solution\", \"Exact Solution\"" << std::endl;
  Outfile <<  "ZONE T = \"Solution Zone\", I = " << Nx_ << ", J = " << Nv_ << ", F = POINT" << std::endl;
  Outfile << std::fixed << std::setprecision(16);
  for (int j = 0; j < Nv_; j++) {
    for (int i = 0; i < Nx_; i++) {
      real_t x0 = X_(i);
      real_t v0 = V_(j);
      Outfile << x0 << " " 
              << v0 << " " 
              << numerical[j](i) << " "
              << exact[j](i) << " " << std::endl;
    };
  }
};


void KineticFDmesh::computerrorL1(const std::vector<Vector>& numerical, 
                          const std::vector<Vector>& exact, 
                          real_t* error) {
  std::vector<Vector> err = numerical - exact;
  *error = 0.e0;
  for (int j = 0; j < Nv_ ; j++) {
    err[j] = err[j].cwiseAbs();
    *error = *error + vweights_(j) * err[j].sum() * hx_;
  }
};

void KineticFDmesh::computerrorL2(const std::vector<Vector>& numerical, 
                          const std::vector<Vector>& exact, 
                          real_t* error) {
  std::vector<Vector> err = numerical - exact;
  *error = 0.e0;
  for (int j = 0; j < Nv_ ; j++) {
    err[j] = err[j].cwiseAbs2();
    *error = *error + vweights_(j) * err[j].sum() * hx_;
  }
  *error = std::sqrt(*error);
};

void KineticFDmesh::computerrorLinf(const std::vector<Vector>& numerical, 
                          const std::vector<Vector>& exact, 
                          real_t* error) {
  std::vector<Vector> err = numerical - exact;
  *error = 0.e0;
  for (int j = 0; j < Nv_ ; j++) {
    *error = std::max(*error, err[j].cwiseAbs().maxCoeff());
  }
};

void KineticFDmesh::generatehv() {
  hv_ = (v2_ - v1_) / real_t(vDiv_);
};

void KineticFDmesh::generatevpoints() {
  Nv_ = vDiv_ + 1;
  vpoints_.resize(vDiv_ + 1);
  for (int i = 0; i <= vDiv_; ++i) {
    vpoints_(i) = i * hv_ + v1_;
  }
  V_ = vpoints_;
};

void KineticFDmesh::generatevweights() {
  vweights_.resize(vDiv_ + 1);
  for (int i = 0; i <= vDiv_; ++i) {
    vweights_(i) = hv_;
  }
};

// KineticFVmesh::KineticFVmesh(const real_t& x1, const real_t& x2,
//         const real_t& v1, const real_t& v2,
//         const int& xDiv, const int& vDiv,
//         const int& qua_num, const QuadratureType& qua_type)
//   : KineticFDmesh(x1, x2, v1, v2, xDiv, vDiv),
//     qua_num_(qua_num), qua_type_(qua_type) {};

// void KineticFVmesh::init()
// {
//   generatehx();
//   generatehv();
//   generatexpoints();
//   generatevpoints();
//   generatevweights();
//   generateCellCenter(); 
//   generatexpoints();
//   generatevpoints();
//   generatevweights();
//   generatexPb();
//   generatevPb();
// };

// const Matrix& KineticFVmesh::getxPb() const
// {
//   return xPb_;
// };

// const Matrix& KineticFVmesh::getvPb() const
// {
//   return vPb_;
// };

// const Vector& KineticFVmesh::getwqua() const
// {
//   return wqua_;
// };

// const Vector& KineticFVmesh::getxqua() const
// {
//   return xqua_;
// };

// void KineticFVmesh::generatexpoints() 
// {
//   xpoints_.resize(xDiv_ + 1);
//   for (int i = 0; i <= xDiv_; ++i) {
//     xpoints_(i) = real_t(i) * hx_ + x1_;
//   };
//   Nx_ = xDiv_;
//   X_.resize(Nx_);
//   for (int i = 0; i < Nx_; ++i) 
//   {
//     X_(i) = (real_t(i) + 0.5e0) * hx_ + x1_;
//   };
// };

// void KineticFVmesh::generatevpoints()
// {
//   Nv_ = vDiv_;
//   vpoints_.resize(vDiv_);
//   for (int i = 0; i < vDiv_; ++i) 
//   {
//     vpoints_(i) = (real_t(i) + 0.5e0) * hv_ + v1_;
//   };
//   V_ = vpoints_;
// };

// void KineticFVmesh::generatevweights() 
// {
//   vweights_.resize(vDiv_);
//   for (int i = 0; i < vDiv_; ++i) {
//     vweights_(i) = hv_;
//   }
// };

// void KineticFVmesh::generatexPb()
// { 
//   QUEST::initializeQuadradure(qua_num_, qua_type_, &xqua_, &wqua_);
//   xPb_.resize(qua_num_, Nx_);xPb_.setZero();
//   for (int i = 0; i < Nx_; i++) {
//     for (int qu = 0; qu < qua_num_; qu++) {
//       xPb_(qu, i) = X_(i) + xqua_(qu) * hx_;
//     }
//   };
// };

// void KineticFVmesh::generatevPb()
// {
//   vPb_.resize(qua_num_, Nv_);vPb_.setZero();
//   for (int i = 0; i < Nv_; i++) {
//     for (int qu = 0; qu < qua_num_; qu++) {
//       vPb_(qu, i) = V_(i) + xqua_(qu) * hv_;
//     }
//   };
// };

KineticFDmesh_Gauss::KineticFDmesh_Gauss(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) :
    KineticFDmesh(x1, x2, v1, v2, xDiv, vDiv) {};

void KineticFDmesh_Gauss::generatevpoints()
{
  Nv_ = vDiv_;
  QUEST::Internal::initializeGaussLegendre(vDiv_, &vpoints_, &vweights_);
  vpoints_ = vpoints_ * 2;
  V_ = vpoints_;
};

void KineticFDmesh_Gauss::generatevweights() {};

KineticFDmesh_period::KineticFDmesh_period(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) : 
    KineticFDmesh(x1, x2, v1, v2, xDiv, vDiv) {};

void KineticFDmesh_period::generatexpoints() {
  Nx_ = xDiv_;
  xpoints_.resize(xDiv_ + 1);
  for (int i = 0; i <= xDiv_; ++i) {
    xpoints_(i) = i * hx_ + x1_;
  }
  X_.resize(Nx_);
  for (int i = 0; i < Nx_; ++i) {
    X_(i) = (i) * hx_ + x1_;
  }
};

KineticFDmesh_period_Gauss::KineticFDmesh_period_Gauss(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) :
    KineticFDmesh_period(x1, x2, v1, v2, xDiv, vDiv) {};

void KineticFDmesh_period_Gauss::generatevpoints()
{
  Nv_ = vDiv_;
  QUEST::Internal::initializeGaussLegendre(vDiv_, &vpoints_, &vweights_);
  vpoints_ = vpoints_ * 2;
  V_ = vpoints_;
};

void KineticFDmesh_period_Gauss::generatevweights() {};

KineticFDmesh_twovel::KineticFDmesh_twovel
        ( const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) : 
    KineticFDmesh(x1, x2, v1, v2, xDiv, vDiv) {};

void KineticFDmesh_twovel::generatevpoints() 
{
  QUEST_VERIFY(vDiv_ == 2, " vDiv_ must equal to 2 ! ");
  Nv_ = 2;
  vpoints_.resize(2);
  vpoints_(0) = -1.e0;
  vpoints_(1) = 1.e0;
  V_ = vpoints_;
};

void KineticFDmesh_twovel::generatevweights() {
  vweights_.resize(2);
  vweights_(0) = 1.e0 / 2.e0;
  vweights_(1) = 1.e0 / 2.e0;
};

KineticFDmesh_period_twovel::KineticFDmesh_period_twovel(
        const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) : 
    KineticFDmesh_period(x1, x2, v1, v2, xDiv, vDiv) {};

void KineticFDmesh_period_twovel::generatevpoints() {
  Nv_ = 2;
  vpoints_.resize(2);
  vpoints_(0) = -1.e0;
  vpoints_(1) = 1.e0;
  V_ = vpoints_;
};

void KineticFDmesh_period_twovel::generatevweights() {
  vweights_.resize(2);
  vweights_(0) = 1.e0 / 2.e0;
  vweights_(1) = 1.e0 / 2.e0;
};
  


};

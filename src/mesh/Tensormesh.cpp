/*
 * Copyright (C) 2025 Xiamen University
 *
 * @Author: Liang Pan
 * @Date:   2025-04-12
 * @Last Modified by: Liang Pan
 * @Last Modified time: 2025-04-12
 */

#include "Tensormesh.hpp"

namespace QUEST
{

TensorMesh1D::TensorMesh1D(const real_t& x1,
    const real_t& x2,
    const int& xDiv): x1_(x1), x2_(x2), xDiv_(xDiv) {};

void TensorMesh1D::init() {
  generatencell();
  generatehx();
  generatePT();
  generateCellCenter();
  generateCellVol();
  generateJx();
  generateJacobiDet();
  generateintNeighbors();
  generateextNeighbors();
};

const real_t& TensorMesh1D::getx1() const {
  return x1_;
};

const real_t& TensorMesh1D::getx2() const {
  return x2_;
};

const int& TensorMesh1D::getxDiv() const {
  return xDiv_;
};

const int& TensorMesh1D::getncell() const {
  return ncell_;
};

const real_t& TensorMesh1D::getRefVol1D() const {
  return RefVol1D_;
};

const real_t& TensorMesh1D::gethx() const {
  return hx_;
};

const Vector& TensorMesh1D::getCellVol() const {
  return CellVol_;
};

const Vector& TensorMesh1D::getJx() const {
  return Jx_;
};

const Vector& TensorMesh1D::getJacobiDet() const {
  return JacobiDet_;
};

const DiagnalMatrix& TensorMesh1D::getCellVol_diag() const {
  return CellVol_diag_;
};

const DiagnalMatrix& TensorMesh1D::getJx_diag() const {
  return Jx_diag_;
};

const DiagnalMatrix& TensorMesh1D::getJacobiDet_diag() const {
  return JacobiDet_diag_;
};

const Matrix& TensorMesh1D::getP() const {
  return P_;
};

const Matrix& TensorMesh1D::getCellCenter() const {
  return CellCenter_;
};

const Eigen::MatrixXi& TensorMesh1D::getT() const {
  return T_;
};

const BoundaryType& TensorMesh1D::getboundaryType() const {
  return boundary_type_;
};

bool TensorMesh1D::IsPeriodBoundary() const {
  if (boundary_type_ == BoundaryType::PeriodBoundary)
  {
    return true;
  }
  return false;
};

const int& TensorMesh1D::getintboundaryNum() const {
  return intboundaryNum_;
};

const std::vector<Eigen::Vector2i>& TensorMesh1D::getintboundaryneighbors() const {
  return intboundaryneighbors_;
};

const std::vector<real_t>& TensorMesh1D::getintboundarynormal() const {
  return intboundarynormal_;
};

const int& TensorMesh1D::getextboundaryNum() const {
  return extboundaryNum_;
};

const std::vector<int>& TensorMesh1D::getextboundaryneighbors() const {
  return extboundaryneighbors_;
};

const std::vector<real_t>& TensorMesh1D::getextboundarynormal() const {
  return extboundarynormal_;
};

const std::vector<BoundaryType>& TensorMesh1D::getextboundarytype() const {
  return extboundarytype_;
};

const std::vector<Eigen::Vector2i>& TensorMesh1D::getextboundaryneighbors_period() const {
  return extboundaryneighbors_period_;
};

void TensorMesh1D::generatencell() {
  ncell_ = xDiv_;
};

void TensorMesh1D::generatehx() {
  hx_ = (x2_ - x1_) / real_t(xDiv_);
};

void TensorMesh1D::generatePT() {
  P_.resize(1, (xDiv_ + 1)); // Resize the matrix to fit all nodes
  for (int i = 0; i <= xDiv_; ++i) {
    P_(0, i) = i * hx_ + x1_;
  }

  // 1 -- 2
  T_.resize(2, ncell_);
  for (int i = 0; i < xDiv_; ++i) {
    int Left = i;
    int Right = i + 1;
    T_(0, i) = Left;
    T_(1, i) = Right;
  }
};

void TensorMesh1D::generateCellCenter() {
  CellCenter_.resize(1, ncell_);
  for (int i = 0; i < xDiv_; ++i) {
    CellCenter_(0, i) = (real_t(i) + 0.5e0) * hx_ + x1_;
  }
};

void TensorMesh1D::generateCellVol() {
  CellVol_.resize(ncell_);
  CellVol_.setConstant(hx_);

  CellVol_diag_ = DiagnalMatrix(CellVol_);
};

void TensorMesh1D::generateJx() {
  Jx_.resize(ncell_);
  Jx_.setConstant(RefVol1D_/hx_);

  Jx_diag_ = DiagnalMatrix(Jx_);

};

void TensorMesh1D::generateJacobiDet() {
  JacobiDet_.resize(ncell_);
  JacobiDet_.setConstant(hx_ / RefVol1D_);

  JacobiDet_diag_ = DiagnalMatrix(JacobiDet_);
};

void TensorMesh1D::generateintNeighbors() {
  intboundaryNum_ = xDiv_ - 1;
  intboundaryneighbors_.resize(intboundaryNum_);
  intboundarynormal_.resize(intboundaryNum_);
  for (int i = 0; i < intboundaryNum_; i++) {
    int CellIndex = i;
    Eigen::Vector2i neighborIndices;
    real_t normals;

    neighborIndices(0) = CellIndex;
    neighborIndices(1) = CellIndex + 1;
    normals = 1.e0;

    intboundaryneighbors_[i] = neighborIndices;
    intboundarynormal_[i] = normals;
  }
};

void TensorMesh1D::generateextNeighbors() {
  extboundaryNum_ = 2;
  extboundaryneighbors_.resize(extboundaryNum_);
  extboundarynormal_.resize(extboundaryNum_);


  extboundaryneighbors_[0] = 0;
  extboundarynormal_[0] = - 1.e0;

  extboundaryneighbors_[1] = xDiv_ - 1;
  extboundarynormal_[1] = 1.e0;

};

void TensorMesh1D::generateextNeighbors_period(const BoundaryType& boundarytype) {
  extboundarytype_.resize(2);
  extboundarytype_[0] = boundarytype;
  extboundarytype_[1] = boundarytype;
  boundary_type_ = boundarytype;

  if (boundary_type_ == BoundaryType::PeriodBoundary) {
    extboundaryneighbors_period_.resize(2);

    extboundaryneighbors_period_[0](0) = 0;
    extboundaryneighbors_period_[0](1) = xDiv_ - 1;

    extboundaryneighbors_period_[1](0) = xDiv_ - 1;
    extboundaryneighbors_period_[1](1) = 0;

  } 
}

KineticTensorMesh1D::KineticTensorMesh1D(
        const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) : 
    TensorMesh1D(x1, x2, xDiv), v1_(v1), v2_(v2), vDiv_(vDiv) {};

void KineticTensorMesh1D::init() {
  generatencell();
  generatehx();
  generatePT();
  generateCellCenter();
  generateCellVol();
  generateJx();
  generateJacobiDet();
  generateintNeighbors();
  generateextNeighbors();

  generatehx();
  generatehv();
  generatevpoints();
  generatevweights();
}

const real_t& KineticTensorMesh1D::getv1() const {
  return v1_;
};

const real_t& KineticTensorMesh1D::getv2() const {
  return v2_;
};

const real_t& KineticTensorMesh1D::gethv() const {
  return hv_;
};

const int& KineticTensorMesh1D::getvDiv() const {
  return vDiv_;
};

const int& KineticTensorMesh1D::getNv() const {
  return Nv_;
};

const Vector& KineticTensorMesh1D::getvpoints() const {
  return vpoints_;
};

const Vector& KineticTensorMesh1D::getV() const {
  return V_;
};

const Vector& KineticTensorMesh1D::getvweights() const {
  return vweights_;
};

void KineticTensorMesh1D::generatehv() {
  hv_ = (v2_ - v1_) / real_t(vDiv_);
};

void KineticTensorMesh1D::generatevpoints() {
  Nv_ = vDiv_ + 1;
  vpoints_.resize(vDiv_ + 1);
  for (int i = 0; i <= vDiv_; ++i) {
    vpoints_(i) = i * hv_ + v1_;
  }
  V_ = vpoints_;
};

void KineticTensorMesh1D::generatevweights() {
  vweights_.resize(vDiv_ + 1);
  for (int i = 0; i <= vDiv_; ++i) {
    vweights_(i) = hv_;
  }
};

KineticTensorMesh1D_twovel::KineticTensorMesh1D_twovel(
        const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) : 
    KineticTensorMesh1D(x1, x2, v1, v2, xDiv, vDiv) {};

// hv_ is not constant in Gauss points
void KineticTensorMesh1D_twovel::generatehv() {};

void KineticTensorMesh1D_twovel::generatevpoints() {
  QUEST_VERIFY(v1_ == -1.e0, " In twovel mesh, v1 must be -1.0 ");
  QUEST_VERIFY(v2_ == 1.e0, " In twovel mesh, v2 must be 1.0 ");
  QUEST_VERIFY(vDiv_ == 2, " In twovel mesh, vDiv must be 1.0 ");
  Nv_ = vDiv_;
  vpoints_.resize(vDiv_);
  vweights_.resize(vDiv_);

  vpoints_ << -1.e0, 1.e0;
  vweights_ << 0.5e0, 0.5e0;

  V_ = vpoints_;
};  

void KineticTensorMesh1D_twovel::generatevweights() {};

const real_t& KineticTensorMesh1D_twovel::gethv() const {
  QUEST_ERROR(" In Gauss Mesh on the velocity space, hv is not constant.");
  return hv_;
};


KineticTensorMesh1D_Gauss::KineticTensorMesh1D_Gauss(
        const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv) : 
    KineticTensorMesh1D(x1, x2, v1, v2, xDiv, vDiv) {};

// hv_ is not constant in Gauss points
void KineticTensorMesh1D_Gauss::generatehv() {};

void KineticTensorMesh1D_Gauss::generatevpoints() {
  Nv_ = vDiv_;
  QUEST::Internal::initializeGaussLegendre(vDiv_, &vpoints_, &vweights_);
  vpoints_ = vpoints_ * 2;
  V_ = vpoints_;
};  

void KineticTensorMesh1D_Gauss::generatevweights() {};

const real_t& KineticTensorMesh1D_Gauss::gethv() const {
  QUEST_ERROR(" In Gauss Mesh on the velocity space, hv is not constant.");
  return hv_;
};

TensorMesh2D::TensorMesh2D(const real_t& x1,
          const real_t& x2,
          const real_t& y1,
          const real_t& y2,
          const int& xDiv,
          const int& yDiv): x1_(x1), x2_(x2), y1_(y1), y2_(y2), xDiv_(xDiv), yDiv_(yDiv) {};

void TensorMesh2D::init() {
  generatencell();
  generatehx();
  generatePT();
  generateCellCenter();
  generateCellVol();
  generateJxJy();
  generateJacobiDet();
  generateintNeighbors();
  generateextNeighbors();
};


const real_t& TensorMesh2D::getx1() const {
  return x1_;
};

const real_t& TensorMesh2D::getx2() const {
  return x2_;
};

const real_t& TensorMesh2D::gety1() const {
  return y1_;
};

const real_t& TensorMesh2D::gety2() const {
  return y2_;
};

const int& TensorMesh2D::getxDiv() const {
  return xDiv_;
};

const int& TensorMesh2D::getyDiv() const {
  return yDiv_;
};

const int& TensorMesh2D::getncell() const {
  return ncell_;
};

const real_t& TensorMesh2D::getRefVol1D() const {
  return RefVol1D_;
};

const real_t& TensorMesh2D::gethx() const {
  return hx_;
};

const real_t& TensorMesh2D::gethy() const {
  return hy_;
};

const Vector& TensorMesh2D::getCellVol() const {
  return CellVol_;
};

const Vector& TensorMesh2D::getJx() const {
  return Jx_;
};

const Vector& TensorMesh2D::getJy() const {
  return Jy_;
};

const Vector& TensorMesh2D::getJacobiDet() const {
  return JacobiDet_;
};

const DiagnalMatrix& TensorMesh2D::getCellVol_diag() const {
  return CellVol_diag_;
};

const DiagnalMatrix& TensorMesh2D::getJx_diag() const {
  return Jx_diag_;
};

const DiagnalMatrix& TensorMesh2D::getJy_diag() const {
  return Jy_diag_;
};

const DiagnalMatrix& TensorMesh2D::getJacobiDet_diag() const {
  return JacobiDet_diag_;
};

const Matrix& TensorMesh2D::getP() const {
  return P_;
};

const Matrix& TensorMesh2D::getCellCenter() const {
  return CellCenter_;
};

const Eigen::MatrixXi& TensorMesh2D::getT() const {
  return T_;
};

const BoundaryType& TensorMesh2D::getboundaryType() const {
  return boundary_type_;
};  

bool TensorMesh2D::IsPeriodBoundary() const {
  if (boundary_type_ == BoundaryType::PeriodBoundary)
  {
    return true;
  }
  return false;
};

int TensorMesh2D::CellIndex(int i, int j) const
{
  return j * xDiv_ + i;
};

int TensorMesh2D::xCellIndex(int i) const
{
  QUEST_VERIFY(i >= 0 && i < ncell_, "xCellIndex out of range");
  return i % xDiv_;
};

int TensorMesh2D::yCellIndex(int j) const
{
  QUEST_VERIFY(j >= 0 && j < ncell_, "yCellIndex out of range");
  return j / xDiv_;
};

const int& TensorMesh2D::getintboundaryNum() const {
  return intboundaryNum_;
};

const std::vector<real_t>& TensorMesh2D::getintboundarylength() const {
  return intboundarylength_;
};

const std::vector<Eigen::Vector2i>& TensorMesh2D::getintboundaryneighbors() const {
  return intboundaryneighbors_;
};

const std::vector<Eigen::Vector2i>& TensorMesh2D::getintboundarytypeindex() const {
  return intboundarytypeindex_;
};

const std::vector<Eigen::Vector2d>& TensorMesh2D::getintboundarynormal() const {
  return intboundarynormal_;
};

const int& TensorMesh2D::getextboundaryNum() const {
  return extboundaryNum_;
};

const std::vector<real_t>& TensorMesh2D::getextboundarylength() const {
  return extboundarylength_;
};

const std::vector<int>& TensorMesh2D::getextboundaryneighbors() const {
  return extboundaryneighbors_;
};

const std::vector<int>& TensorMesh2D::getextboundarytypeindex() const {
  return extboundarytypeindex_;
};

const std::vector<Eigen::Vector2d>& TensorMesh2D::getextboundarynormal() const {
  return extboundarynormal_;
};

const std::vector<Eigen::Vector2d>& TensorMesh2D::getextboundarycenter() const {
  return extboundarycenter_;
};

const std::vector<BoundaryType>& TensorMesh2D::getextboundarytype() const {
  return extboundarytype_;
};

const std::vector<Eigen::Vector2i>& TensorMesh2D::getextboundaryneighbors_period() const {
  return extboundaryneighbors_period_;
};

const std::vector<Eigen::Vector2i>& TensorMesh2D::getextboundarytypeindex_period() const {
  return extboundarytypeindex_period_;
};

void TensorMesh2D::generatencell()
{
  ncell_ = xDiv_ * yDiv_;
};

void TensorMesh2D::generatehx()
{
  hx_ = (x2_ - x1_) / real_t(xDiv_);
  hy_ = (y2_ - y1_) / real_t(yDiv_);
};

void TensorMesh2D::generatePT()
{
  P_.resize(2, (xDiv_ + 1) * (yDiv_ + 1)); // Resize the matrix to fit all nodes
  int nodeIndex = 0;
  for (int j = 0; j <= yDiv_; ++j) 
  {
    for (int i = 0; i <= xDiv_; ++i) 
    {
      P_(0, nodeIndex) = i * hx_ + x1_;
      P_(1, nodeIndex) = j * hy_ + y1_;
      ++nodeIndex;
    }
  }

  T_.resize(4, ncell_);
  int cellIndex = 0;
  for (int j = 0; j < yDiv_; ++j) 
  {
    for (int i = 0; i < xDiv_; ++i) 
    {
      int bottomLeft = j * (xDiv_ + 1) + i;
      int bottomRight = bottomLeft + 1;
      int topLeft = bottomLeft + (xDiv_ + 1);
      int topRight = topLeft + 1;
      T_(0, cellIndex) = bottomLeft; //左下角点
      T_(1, cellIndex) = bottomRight; // 右下角点 
      T_(2, cellIndex) = topRight;  // 右上角点
      T_(3, cellIndex) = topLeft; // 左上角点
      ++cellIndex;
    }
  }
};

void TensorMesh2D::generateCellCenter()
{
  CellCenter_.resize(2, ncell_);
  int cellIndex = 0;
  for (int j = 0; j < yDiv_; ++j) 
  {
    for (int i = 0; i < xDiv_; ++i) 
    {
      CellCenter_(0, cellIndex) = (real_t(i) + 0.5e0) * hx_ + x1_;
      CellCenter_(1, cellIndex) = (real_t(j) + 0.5e0) * hy_ + y1_;
      ++cellIndex;
    }
  }
};

void TensorMesh2D::generateCellVol()
{
  CellVol_.resize(ncell_);
  for (int i = 0; i < ncell_; i++)
  {
    CellVol_(i) = hx_ * hy_;
  }
};

void TensorMesh2D::generateJxJy()
{
  Jx_.resize(ncell_);
  Jy_.resize(ncell_);
  for (int i = 0; i < ncell_; i++)
  {
    Jx_(i) = RefVol2D_ / hx_; // \xi对x求偏导
    Jy_(i) = RefVol2D_ / hy_; // \eta对y求偏导
  }
};

void TensorMesh2D::generateJacobiDet()
{
  JacobiDet_.resize(ncell_);
  JacobiDet_.setConstant(hx_ * hy_ / RefVol2D_);

  JacobiDet_diag_ = DiagnalMatrix(JacobiDet_);
};

void TensorMesh2D::generateintNeighbors()
{
  intboundaryNum_ = xDiv_ * (yDiv_ - 1) + (xDiv_ - 1) * yDiv_;
  intboundarytypeindex_.resize(intboundaryNum_);
  intboundaryneighbors_.resize(intboundaryNum_);
  intboundarylength_.resize(intboundaryNum_);
  intboundarynormal_.resize(intboundaryNum_);
  
  for (int j = 0; j < yDiv_; j++)
  {
    for (int i = 0; i < (xDiv_ - 1); i++)
    {
      int Index = j * (xDiv_ - 1) + i;
      int CellIndex = j * xDiv_ + i;
      Eigen::Vector2i neighborIndices;
      real_t length;
      Eigen::Vector2d normals;

      neighborIndices(0) = CellIndex;
      neighborIndices(1) = CellIndex + 1;
      length = hy_;
      normals << 1.e0, 0.e0;

      intboundaryneighbors_[Index] = neighborIndices;
      intboundarylength_[Index] = length;
      intboundarynormal_[Index] = normals;
      intboundarytypeindex_[Index](0) = 0;
      intboundarytypeindex_[Index](1) = 1;
    }
  }

  for (int j = 0; j < (yDiv_ - 1); j++)
  {
    for (int i = 0; i < xDiv_; i++)
    {
      int Index = j * xDiv_ + i;
      int CellIndex = j * xDiv_ + i;
      Eigen::Vector2i neighborIndices;
      real_t length;
      Eigen::Vector2d normals;

      neighborIndices(0) = CellIndex;
      neighborIndices(1) = CellIndex + xDiv_;
      length = hx_;
      normals << 0.e0, 1.e0;

      intboundaryneighbors_[Index + (xDiv_ - 1) * yDiv_] = neighborIndices;
      intboundarylength_[Index + (xDiv_ - 1) * yDiv_] = length;
      intboundarynormal_[Index + (xDiv_ - 1) * yDiv_] = normals;
      intboundarytypeindex_[Index + (xDiv_ - 1) * yDiv_](0) = 2;
      intboundarytypeindex_[Index + (xDiv_ - 1) * yDiv_](1) = 3;
    }
  }
};

void TensorMesh2D::generateextNeighbors()
{
  extboundaryNum_ = (xDiv_ + yDiv_) * 2;
  extboundaryneighbors_.resize(extboundaryNum_);
  extboundarylength_.resize(extboundaryNum_);
  extboundarynormal_.resize(extboundaryNum_);
  extboundarycenter_.resize(extboundaryNum_);
  extboundarytypeindex_.resize(extboundaryNum_);

  int Index = 0;
  for (int i = 0; i < xDiv_; i++)
  { //下边界
    int CellIndex = i;
    real_t length = hx_;
    Eigen::Vector2d normals;

    normals << 0.e0 , - 1.e0;
    
    extboundaryneighbors_[Index] = CellIndex;
    extboundarylength_[Index] = length;
    extboundarynormal_[Index] = normals;
    extboundarycenter_[Index](0) = (real_t(i) + 0.5e0) * hx_ + x1_;
    extboundarycenter_[Index](1) = y1_;
    extboundarytypeindex_[Index] = 3;
    Index++;
  };

  for (int j = 0; j < yDiv_; j++)
  { //右边界
    int CellIndex = j * xDiv_ + (xDiv_ - 1);
    real_t length = hy_;
    Eigen::Vector2d normals;
    
    normals << 1.e0 , 0.e0;
    
    extboundaryneighbors_[Index] = CellIndex;
    extboundarylength_[Index] = length;
    extboundarynormal_[Index] = normals;
    extboundarycenter_[Index](0) = x2_;
    extboundarycenter_[Index](1) = (real_t(j) + 0.5e0) * hy_ + y1_;
    extboundarytypeindex_[Index] = 0;
    Index++;
  };

  for (int i = 0; i < xDiv_; i++)
  { //上边界
    int CellIndex = (yDiv_ - 1) * xDiv_ + i;
    real_t length = hx_;
    Eigen ::Vector2d normals;
    
    normals << 0.e0 , 1.e0;
    
    extboundaryneighbors_[Index] = CellIndex;
    extboundarylength_[Index] = length;
    extboundarynormal_[Index] = normals;
    extboundarycenter_[Index](0) = (real_t(i) + 0.5e0) * hx_ + x1_;
    extboundarycenter_[Index](1) = y2_;
    extboundarytypeindex_[Index] = 2;
    Index++;
  };
  
  for (int j = 0; j < yDiv_; j++)
  { //左边界
    int CellIndex = j * xDiv_ ;
    real_t length = hy_;
    Eigen ::Vector2d normals;

    normals << - 1.e0 , 0.e0;
    
    extboundaryneighbors_[Index] = CellIndex;
    extboundarylength_[Index] = length;
    extboundarynormal_[Index] = normals;
    extboundarycenter_[Index](0) = x1_;
    extboundarycenter_[Index](1) = (real_t(j) + 0.5e0) * hy_ + y1_;
    extboundarytypeindex_[Index] = 1;
    Index++;
  };
};

void TensorMesh2D::generateextNeighbors_period(const BoundaryType& boundarytype)
{
  extboundaryneighbors_period_.resize(extboundaryNum_);
  extboundarytypeindex_period_.resize(extboundaryNum_);
  extboundarytype_.resize(extboundaryNum_);
  boundary_type_ = boundarytype;
  if (boundary_type_ == BoundaryType::PeriodBoundary) 
  {
    int Index = 0;
    for (int i = 0; i < xDiv_; i++) 
    { //下边界
      int CellIndex0 = i;
      int CellIndex1 = (yDiv_ - 1) * xDiv_ + i;
      Eigen::Vector2i neighborIndices;

      neighborIndices(0) = CellIndex0;
      neighborIndices(1) = CellIndex1;

      extboundaryneighbors_period_[Index] = neighborIndices;
      extboundarytype_[Index] = boundarytype;
      extboundarytypeindex_period_[Index](0) = 3;
      extboundarytypeindex_period_[Index](1) = 2;
      Index++;
    };

    for (int j = 0; j < yDiv_; j++)
    { //右边界
      int CellIndex0 = j * xDiv_ + (xDiv_ - 1);
      int CellIndex1 = j * xDiv_;
      Eigen::Vector2i neighborIndices;
      
      neighborIndices(0) = CellIndex0;
      neighborIndices(1) = CellIndex1;

      extboundaryneighbors_period_[Index] = neighborIndices;
      extboundarytype_[Index] = boundarytype;
      extboundarytypeindex_period_[Index](0) = 0;
      extboundarytypeindex_period_[Index](1) = 1;
      Index++;
    }

    for (int i = 0; i < xDiv_; i++)
    { //上边界
      int CellIndex0 = (yDiv_ - 1) * xDiv_ + i;
      int CellIndex1 = i;
      Eigen::Vector2i neighborIndices;
      
      neighborIndices(0) = CellIndex0;
      neighborIndices(1) = CellIndex1;
      
      extboundaryneighbors_period_[Index] = neighborIndices;
      extboundarytype_[Index] = boundarytype;
      extboundarytypeindex_period_[Index](0) = 2;
      extboundarytypeindex_period_[Index](1) = 3;
      Index++;
    }
    
    for (int j = 0; j < yDiv_; j++)
    { //左边界
      int CellIndex0 = j * xDiv_;
      int CellIndex1 = j * xDiv_ + (xDiv_ - 1); 
      Eigen::Vector2i neighborIndices;
      
      neighborIndices(0) = CellIndex0;
      neighborIndices(1) = CellIndex1;
      
      extboundaryneighbors_period_[Index] = neighborIndices;
      extboundarytype_[Index] = boundarytype;
      extboundarytypeindex_period_[Index](0) = 1;
      extboundarytypeindex_period_[Index](1) = 0;
      Index++;
    }
  }
};


} // namespace QUEST



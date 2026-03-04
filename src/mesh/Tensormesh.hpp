/*
 * Copyright (C) 2025 Xiamen University
 *
 * @Author: Liang Pan
 * @Date:   2025-04-12
 * @Last Modified by: Liang Pan
 * @Last Modified time: 2025-04-12
 */

#ifndef QUEST_TENSORMESH_HPP
#define QUEST_TENSORMESH_HPP

#include <vector>

#include "config.hpp"
#include "integralrule.hpp"

namespace QUEST {

enum class BoundaryType : int {
  DirichletBoundary = -1,
  PeriodBoundary = -2,
  InflowBoundary = -3,
  OutflowBoundary = -4,
  RefelectiveBoundary = -5,
};

class TensorMesh1D {
public:
  // 一维均匀网格构造
  TensorMesh1D(const real_t& x1,
          const real_t& x2,
          const int& xDiv);
  virtual ~TensorMesh1D() = default;
  virtual void init();
  
  virtual const real_t& getx1() const;
  virtual const real_t& getx2() const;
  virtual const int& getxDiv() const;
  virtual const int& getncell() const;
  virtual const real_t& getRefVol1D() const;
  virtual const real_t& gethx() const;
  virtual const Vector& getCellVol() const;
  virtual const Vector& getJx() const;
  virtual const Vector& getJacobiDet() const;
  virtual const DiagnalMatrix& getCellVol_diag() const;
  virtual const DiagnalMatrix& getJx_diag() const;
  virtual const DiagnalMatrix& getJacobiDet_diag() const;
  virtual const Matrix& getP() const;
  virtual const Matrix& getCellCenter() const;
  virtual const Eigen::MatrixXi& getT() const;
  virtual const BoundaryType& getboundaryType() const;
  virtual bool IsPeriodBoundary() const;

  virtual const int& getintboundaryNum() const;
  virtual const std::vector<Eigen::Vector2i>& getintboundaryneighbors() const;
  virtual const std::vector<real_t>& getintboundarynormal() const;

  virtual const int& getextboundaryNum() const;
  virtual const std::vector<int>& getextboundaryneighbors() const;
  virtual const std::vector<real_t>& getextboundarynormal() const;
  virtual const std::vector<QUEST::BoundaryType>& getextboundarytype() const;
  virtual const std::vector<Eigen::Vector2i>& getextboundaryneighbors_period() const;

  virtual void generateextNeighbors_period(const BoundaryType& boundarytype);

protected:
  
  const real_t& x1_;
  const real_t& x2_;
  const int& xDiv_;

  int ncell_;
  real_t hx_;
  real_t RefVol1D_ = 1.e0;
  Vector CellVol_;
  Vector Jx_;
  Vector JacobiDet_;
  DiagnalMatrix CellVol_diag_;
  DiagnalMatrix Jx_diag_;
  DiagnalMatrix JacobiDet_diag_;
  Matrix CellCenter_;
  Matrix P_;
  Eigen::MatrixXi T_;
  BoundaryType boundary_type_;

  int intboundaryNum_;
  std::vector<Eigen::Vector2i> intboundaryneighbors_;
  std::vector<real_t> intboundarynormal_;

  int extboundaryNum_;
  std::vector<int> extboundaryneighbors_;
  std::vector<real_t> extboundarynormal_;
  std::vector<QUEST::BoundaryType> extboundarytype_;

  std::vector<Eigen::Vector2i> extboundaryneighbors_period_;

  virtual void generatencell();
  virtual void generatehx();
  virtual void generatePT();
  virtual void generateCellCenter();
  virtual void generateCellVol();
  virtual void generateJx();
  virtual void generateJacobiDet();
  virtual void generateintNeighbors();
  virtual void generateextNeighbors();

};

class KineticTensorMesh1D : public TensorMesh1D 
{
public:
  KineticTensorMesh1D(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);
  ~KineticTensorMesh1D() override = default;
  void init() override;

  virtual const real_t& getv1() const;
  virtual const real_t& getv2() const;
  virtual const real_t& gethv() const;
  virtual const int& getvDiv() const;
  virtual const int& getNv() const;
  virtual const Vector& getvpoints() const;
  virtual const Vector& getV() const;
  virtual const Vector& getvweights() const;

protected:
  const real_t& v1_;
  const real_t& v2_;
  const int& vDiv_;
  
  real_t hv_;
  int Nv_;
  Vector vpoints_;
  Vector V_;
  Vector vweights_;

  virtual void generatehv();
  virtual void generatevpoints();
  virtual void generatevweights();
};

class KineticTensorMesh1D_twovel : public KineticTensorMesh1D
{
public:
  KineticTensorMesh1D_twovel(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);
  ~KineticTensorMesh1D_twovel() override = default;

protected:
  void generatehv() override;
  void generatevpoints() override;
  void generatevweights() override;

  const real_t& gethv() const override;
};

class KineticTensorMesh1D_Gauss : public KineticTensorMesh1D
{
public:
  KineticTensorMesh1D_Gauss(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);
  ~KineticTensorMesh1D_Gauss() override = default;

protected:
  void generatehv() override;
  void generatevpoints() override;
  void generatevweights() override;
  
  const real_t& gethv() const override;
};


class TensorMesh2D {
public:
  // 二维均匀网格构造
  TensorMesh2D(const real_t& x1,
          const real_t& x2,
          const real_t& y1,
          const real_t& y2,
          const int& xDiv,
          const int& yDiv);
  virtual ~TensorMesh2D() = default;
  virtual void init();
  
  virtual const real_t& getx1() const;
  virtual const real_t& getx2() const;
  virtual const real_t& gety1() const;
  virtual const real_t& gety2() const;
  virtual const int& getxDiv() const;
  virtual const int& getyDiv() const;
  virtual const int& getncell() const;
  virtual const real_t& getRefVol1D() const;
  virtual const real_t& gethx() const;
  virtual const real_t& gethy() const;
  virtual const Vector& getCellVol() const;
  virtual const Vector& getJx() const;
  virtual const Vector& getJy() const;
  virtual const Vector& getJacobiDet() const;
  virtual const DiagnalMatrix& getCellVol_diag() const;
  virtual const DiagnalMatrix& getJx_diag() const;
  virtual const DiagnalMatrix& getJy_diag() const;
  virtual const DiagnalMatrix& getJacobiDet_diag() const;
  virtual const Matrix& getP() const;
  virtual const Matrix& getCellCenter() const;
  virtual const Eigen::MatrixXi& getT() const;
  virtual const BoundaryType& getboundaryType() const;
  virtual bool IsPeriodBoundary() const;
  virtual int CellIndex(int i, int j) const;
  virtual int xCellIndex(int i) const;
  virtual int yCellIndex(int j) const;

  virtual const int& getintboundaryNum() const;
  virtual const std::vector<real_t>& getintboundarylength() const;
  virtual const std::vector<Eigen::Vector2i>& getintboundaryneighbors() const;
  virtual const std::vector<Eigen::Vector2i>& getintboundarytypeindex() const;
  virtual const std::vector<Eigen::Vector2d>& getintboundarynormal() const;

  virtual const int& getextboundaryNum() const;
  virtual const std::vector<real_t>& getextboundarylength() const;
  virtual const std::vector<int>& getextboundaryneighbors() const;
  virtual const std::vector<int>& getextboundarytypeindex() const;
  virtual const std::vector<Eigen::Vector2d>& getextboundarynormal() const;
  virtual const std::vector<Eigen::Vector2d>& getextboundarycenter() const;
  virtual const std::vector<QUEST::BoundaryType>& getextboundarytype() const;

  virtual const std::vector<Eigen::Vector2i>& getextboundaryneighbors_period() const;
  virtual const std::vector<Eigen::Vector2i>& getextboundarytypeindex_period() const;

  virtual void generateextNeighbors_period(const BoundaryType& boundarytype);

protected:
  
  const real_t& x1_;
  const real_t& x2_;
  const real_t& y1_;
  const real_t& y2_;
  const int& xDiv_;
  const int& yDiv_;

  int ncell_;
  real_t hx_;
  real_t hy_;
  real_t RefVol1D_ = 1.e0;
  real_t RefVol2D_ = 1.e0;
  Vector CellVol_;
  Vector Jx_;
  Vector Jy_;
  Vector JacobiDet_;
  DiagnalMatrix CellVol_diag_;
  DiagnalMatrix Jx_diag_;
  DiagnalMatrix Jy_diag_;
  DiagnalMatrix JacobiDet_diag_;
  Matrix CellCenter_;
  Matrix P_;
  Eigen::MatrixXi T_;
  BoundaryType boundary_type_;

  int intboundaryNum_;
  std::vector<Eigen::Vector2i> intboundarytypeindex_;
  std::vector<real_t> intboundarylength_;
  std::vector<Eigen::Vector2i> intboundaryneighbors_;
  std::vector<Eigen::Vector2d> intboundarynormal_;

  int extboundaryNum_;
  std::vector<real_t> extboundarylength_;
  std::vector<int> extboundarytypeindex_;
  std::vector<int> extboundaryneighbors_;
  std::vector<Eigen::Vector2d> extboundarynormal_;
  std::vector<Eigen::Vector2d> extboundarycenter_;
  std::vector<QUEST::BoundaryType> extboundarytype_;

  std::vector<Eigen::Vector2i> extboundaryneighbors_period_;
  std::vector<Eigen::Vector2i> extboundarytypeindex_period_;

  virtual void generatencell();
  virtual void generatehx();
  virtual void generatePT();
  virtual void generateCellCenter();
  virtual void generateCellVol();
  virtual void generateJxJy();
  virtual void generateJacobiDet();
  virtual void generateintNeighbors();
  virtual void generateextNeighbors();
};

}; // namespace QUEST

#endif // QUEST_TENSORMESH_HPP

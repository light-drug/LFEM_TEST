#include "hyperbolicproblemsbase.hpp"

namespace QUEST
{

HyperbolicProblems1DBase::HyperbolicProblems1DBase(const fespace1D* fe, const int& t_order)
  : fe_(fe), t_order_(t_order) {
  
  // ********* 传入相关变量 ********** //
  intboundaryNum_ = fe_->getmesh1D()->getintboundaryNum();
  extboundaryNum_ = fe_->getmesh1D()->getextboundaryNum();
  num_equations_ = fe_->getnum_equations();
  polydim_ = fe_->getbasis()->getpolydim();
  ncell_ = fe_->getmesh1D()->getncell();
  dim_ = fe_->getdim();
  numqua_ = fe_->getnumqua();
  hx_ = fe_->getmesh1D()->gethx();
  NTH_ = fe_->getNTH();
  // ******************************** //
  u_modal_.resize(num_equations_, Matrix::Zero(polydim_, ncell_));
  u_nodal_.resize(num_equations_, Matrix::Zero(numqua_, ncell_));
  fu_nodal_.resize(dim_, u_nodal_);

  flux_int_.resize(num_equations_, Matrix::Zero(1, intboundaryNum_));
  flux_ext_.resize(num_equations_, Matrix::Zero(1, extboundaryNum_));

};

const std::vector<Matrix>& HyperbolicProblems1DBase::getumodal() const {
  return u_modal_;
};

const std::vector<Matrix>& HyperbolicProblems1DBase::getunodal() const {
  return u_nodal_;
};

const std::vector<std::vector<Matrix>>& HyperbolicProblems1DBase::getfunodal() const {
  return fu_nodal_;
}



} // namespace QUEST

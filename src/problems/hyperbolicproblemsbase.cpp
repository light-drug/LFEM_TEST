#include "hyperbolicproblemsbase.hpp"

namespace QUEST
{

HyperbolicProblems1DBase::HyperbolicProblems1DBase(const TensorMesh1D* mesh1D,
                          const fespace1D* fe, 
                          const EX_TVDRK* rk_table)
  : mesh1D_(mesh1D), fe_(fe), rk_table_(rk_table) {
  
  // ********* 传入相关变量 ********** //
  num_equations_ = fe_->getnum_equations();
  const int& stages = rk_table_->getstages();
  // ******************************** //
  u_modal_.resize(num_equations_, Matrix::Zero(polydim_, ncell_));
  u_modal_stages_.resize(stages, u_modal);
  Lu_stages_.resize(stages, u_modal);
  fu_nodal_.resize(dim_, u_nodal_);

  flux_int_.resize(num_equations_, Matrix::Zero(1, intboundaryNum_));
  flux_ext_.resize(num_equations_, Matrix::Zero(1, extboundaryNum_));

};

const std::vector<Matrix>& HyperbolicProblems1DBase::getumodal() const {
  return u_modal_;
};

void HyperbolicProblems1DBase::setcfl(const real_t& cfl) 
{
  cfl_ = cfl;
};

} // namespace QUEST

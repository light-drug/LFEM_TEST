#include "hyperbolic.hpp"

namespace QUEST
{
// 基类定义
HyperbolicFormBase::HyperbolicFormBase(const fespace1D* fe) 
    : fe_(fe)
{
  // ********** 传入相关变量 ********** //
  const int& num_equations = fe_->getnum_equations();
  const int& polydim = fe_->getbasis()->getpolydim();
  const int& ncell = fe_->getmesh1D()->getncell();
  // ********************************* //
  num_equations_ = num_equations;
};

const fespace1D* HyperbolicFormBase::getfespace() const {
  return fe_;
};

void HyperbolicFormBase::ComputeFuTest(const std::vector<Matrix>& u_nodal, std::vector<Matrix>& Temp) {
  std::vector<Matrix> Fu;
  ComputeFu(u_nodal, i, Fu);
  fe_->Assemble_F(Fu, 1, &Temp);  
};

void HyperbolicFormBase::ComputeFluxTest(const std::vector<Matrix>& flux_int, 
                      const std::vector<Matrix>& flux_ext, 
                      std::vector<Matrix>& Temp) {
  fe_->Assemble_Flux(flux_int, flux_ext, &Temp);
};

// 迎风通量周期边界条件定义 F(u) = u;
UpwindFluxPeriodicForm::UpwindFluxPeriodicForm(const fespace1D* fe, const Vector Coeff)
: HyperbolicFormBase(fe), Coeff_(Coeff) {
  QUEST_VERIFY((num_equations_ == 1)
    ," Upwind flux needs only one equation ! (num_equations need to be equal to one !)");
};

void UpwindFluxPeriodicForm::ComputeFu(const std::vector<Matrix>& u_nodal, 
                                      const int& dimid, 
                                      std::vector<Matrix>& Fu) {
  QUEST_VERIFY((dimid <= dim_), "dimid > dim !! ");

  Fu = Coeff_(dimid - 1) * u;
};

void UpwindFluxPeriodicForm::ComputeFlux(const std::vector<real_t>& u_int, 
                                        const std::vector<real_t>& u_ext,
                                        const Vector& normal, 
                                        std::vector<real_t>& flux) {
  Vector 
  real_t wind = normal.dot()
}

} // namespace QUEST

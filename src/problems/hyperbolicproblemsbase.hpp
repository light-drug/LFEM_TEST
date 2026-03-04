#ifndef QUEST_PROBLEMSBASE_HPP
#define QUEST_PROBLEMSBASE_HPP

#include "fespace.hpp"
#include "config.hpp"
#include <vector>

namespace QUEST
{

class HyperbolicProblems1DBase
{
public:

  HyperbolicProblems1DBase(const fespace1D* fe, const int& t_order);

  ~HyperbolicProblems1DBase() = default;

  virtual void init() = 0;
  virtual void setdt(real_t* dt) = 0;
  virtual void time_stepping(real_t& dt) = 0;
  virtual void update() = 0;

  const std::vector<Matrix>& getumodal() const;
  const std::vector<Matrix>& getunodal() const;
  const std::vector<std::vector<Matrix>>& getfunodal() const;


protected: 

  const fespace1D* fe_;
  const int& t_order_;

  int num_equations_;
  int intboundaryNum_;
  int extboundaryNum_;
  int dim_;
  int polydim_;
  int numqua_;
  int ncell_;
  int hx_;
  int NTH_;
  std::vector<Matrix> u_modal_;
  std::vector<Matrix> u_nodal_;
  std::vector<std::vector<Matrix>> fu_nodal_;
  std::vector<Matrix> flux_int_;
  std::vector<Matrix> flux_ext_;

};

// class IMEXDriftDiffusion1D: public HyperbolicProblems1DBase
// {
// public:

//   IMEXDriftDiffusion1D(/* args */);

//   ~IMEXDriftDiffusion1D();
// };

} // namespace QUEST



#endif // QUEST_PROBLEMSBASE_HPP

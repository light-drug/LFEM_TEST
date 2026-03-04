#ifndef QUEST_KINETICDRIFTDIFFUSION1D_SL_NEWTON_HPP
#define QUEST_KINETICDRIFTDIFFUSION1D_SL_NEWTON_HPP

#include "config.hpp"
#include "FDmesh.hpp"
#include "error.hpp"
#include "timer.hpp"
#include "Poissonsolver1DFD.hpp"
#include "KineticDriftDiffusion1D_SL.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

namespace QUEST 
{

class KineticDriftDiffusion1D_SL_Newton: public KineticDriftDiffusion1D_SL
{
protected:

  SparseMatrix De_;
  SparseMatrix A_;
  Vector deltarho_;
  Vector deltaE_;
  Vector deltaphi_;
  
public:
  KineticDriftDiffusion1D_SL_Newton(const KineticFDmesh* mesh1D,
                        Poissonsolver1DFD* poisol,
                        const int& x_order,
                        const int& t_order);
  ~KineticDriftDiffusion1D_SL_Newton() override = default;

  void init(const Solver1DTypeFD& soltype) override;
  void update_rho(const real_t& Trun, const real_t& dt) override;
  void generateS(const real_t& Trun, const real_t& dt) override;
  void generatebf(const real_t& Trun, 
                  const real_t& dt,
                  Vector* bf) override;
  void SolveS(const real_t& Trun, const Vector& bf) override;
  void update_E(const real_t& Trun) override;
};


} // namespace QUEST


#endif  // QUEST_KINETICDRIFTDIFFUSION1D_SL_NEWTON_HPP

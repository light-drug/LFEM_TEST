#ifndef QUEST_DISPLAY_HPP
#define QUEST_DISPLAY_HPP

#include "config.hpp"
#include "timer.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace QUEST 
{

void DisplayAccTable1D(const IntVector& xDiv_vec,
                       const Vector& errorL1_vec,
                       const Vector& errorL2_vec, 
                       const Vector& errorLinf_vec,
                       std::ostream& Outfile);


void DisplayAccTable2D(const IntVector& xDiv_vec,
                      const IntVector& yDiv_vec,
                      const Vector& errorL1_vec,
                      const Vector& errorL2_vec, 
                      const Vector& errorLinf_vec,
                      std::ostream& Outfile);

} // namespace QUEST

#endif // QUEST_DISPLAY_HPP 

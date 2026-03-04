#ifndef QUEST_HPP
#define QUEST_HPP
// array
#include "src/array/array.hpp"

// config
#include "src/config/config.hpp"

// general
#include "src/general/error.hpp"
#include "src/general/timer.hpp"
#include "src/general/optparser.hpp"
#include "src/general/display.hpp"
#include "src/general/error.hpp"
#include "src/general/vector_overload.hpp"

// mesh
#include "src/mesh/Tensormesh.hpp"
#include "src/mesh/fespace.hpp"
#include "src/mesh/FDmesh.hpp"

// basis
#include "src/basis/basis.hpp"
#include "src/basis/integralrule.hpp"

// src_test
#include "src/src_test/srctest_1.hpp"

// operator
#include "src/operator/RK_table.hpp"

// prblems
#include "src/problems/Euler1D_acctest.hpp"
#include "src/problems/Poissonsolver1DFD.hpp"
#include "src/problems/DriftDiffusion1D_FD.hpp"
#include "src/problems/KineticDriftDiffusion1D_SL.hpp"
#include "src/problems/KineticLinearDiffusion_SL.hpp"
#include "src/problems/KineticDriftDiffusion1D_SL_Newton.hpp"
#include "src/problems/APDriftDiffusion_SL.hpp"
#include "src/problems/KineticDriftDiffusion1D_SL_lomac.hpp"

// DG problems
#include "src/problems/hyperbolicproblemsbase.hpp"
#include "src/problems/poissonsolver1D.hpp"
#include "src/problems/KineticLD_DG_IMEX_EX.hpp"
#include "src/problems/KineticLD_DG_IMEX_Schur.hpp"
#include "src/problems/KineticLD_DG_IMEX_Schur_Dirichlet.hpp"
#include "src/problems/KineticDriftD_DG_IMEX_Schur.hpp"


// problems 2D
#include "src/problems2D/poissonsolver2D.hpp"
#endif  // QUEST_HPP

#include "QUEST.hpp"
#include <iostream>
#include <cmath>
#include <numbers>
#include <complex>

// enum class Solver1DTypeFD {
//   CG = 0
//   PCG = 1,
//   LDLT = 2,
//   LU = 3,
//   GMRES = 4,
//   BICG  = 5,
// };
int main(int argc, char *argv[]) {

  using cdouble = std::complex<double>;

  Eigen::MatrixXcd A(2,2);
  A << cdouble(1,2), cdouble(3,0),
       cdouble(0,-1), cdouble(2,5);
  
  std::cout << "A = \n" << A << std::endl;
  std::cout << " A(1,1) * A(0,0) = " << A(1,1) * A(0,0) << std::endl;
  PAUSE();

  
  int problem = 2201;
  std::cout << " Please input the problem number to run the test: \n";
  QUEST::introduction();
  std::cin >> problem;
  std::cout << " problem is set as :" << problem << std::endl;
  TIC;
  switch(problem) {
    case 1: QUEST::KineticLinearDiffusion1D_SL_acctest(argc, argv); break;
    case 2: QUEST::KineticDriftDiffusion1D_SL_acctest(argc, argv); break;
    case 3: QUEST::KineticDriftDiffusion1D_SL_unipolar(argc, argv); break;
    case 301: QUEST::KineticDriftDiffusion1D_SL_PNjunction_example(argc, argv); 
      // ./main -m 1 -nx 200 -nv 100 -Tstop 0.1e0 -NTH 10 -knu 1.e-6 -gamma 1.e0 -cf1 1.e0 -output 100
      break;
    case 4: QUEST::KineticDriftDiffusion1D_SL_unipolar_Newton(argc, argv); break;
    case 5: QUEST::KineticLinearDiffusion1D_SL_acctest_twovel(argc, argv); break;
      // ./main -m 5 -nx 20 -nv 2 -Tstop 0.1e0 -ox 1 -ot 1 -NTH 10 -cfl 1.e0 -knu 1.e0
    case 501: QUEST::KineticLinearDiffusion1D_SL_IsotropicBoundary(argc, argv); 
      // ./main -m 1 -nx 200 -nv 32 -Tstop 0.1e0 -ox 1 -ot 1 -NTH 10 -knu 1.e0 -output 100
      // ./main -m 1 -nx 200 -nv 32 -Tstop 0.15e0 -ox 1 -ot 1 -NTH 10 -knu 1.e-6 -output 100
      break;
    case 6: QUEST::KienticDriftDiffusion1D_SL_unipolar_QUASINEUTRAL(argc, argv); break;
    case 7: QUEST::KienticDriftDiffusion1D_SL_acctest_QUASINEUTRAL(argc, argv); break;
    case 8: QUEST::KineticDriftDiffusion1D_SL_acctest_Different_Gamma(argc, argv); break;
    case 9: QUEST::KineticDriftDiffusion1D_SL_PNjunction_QUASINEUTRAL(argc, argv); break;
    case 10: QUEST::KineticLinearDiffusion1D_SL_Implicit_acctest_twovel(argc, argv); break;
    case 11: QUEST::DriftDiffusion1D_FD_unipolar(argc, argv); break;
    case 12: QUEST::DriftDiffusion1D_FD_unipolar_penaltyIter(argc, argv); break;
    case 13: QUEST::DriftDiffusion1D_FD_unipolar_penaltyImplicit(argc, argv); break;
    case 14: QUEST::DriftDiffusion1D_FD_acctest(argc, argv); break;
    case 15: QUEST::DriftDiffusion1D_FD_PNjunction_example(argc, argv); break;
    case 16: QUEST::DriftDiffusion1D_FD_PNjunction_penaltyImplicit_example(argc, argv); break;
    case 17: QUEST::Poisson1D_DG_acctest(argc, argv); 
      // ./main -m 4 -nx 20 -basis 2 -poi 3 -poitol 1.e-10 -c11 1.e0 -c12 0.5
      break;
    case 1701: QUEST::Poisson1D_DG_acctest_period(argc, argv); 
      // ./main -m 4 -nx 20 -basis 2 -poi 3 -poitol 1.e-10 -c11 1.e0 -c12 0.5
      break;
    case 18: QUEST::KineticLinearDiffusion1D_DGIMEX_acctest(argc, argv); 
      // ./main -m 4 -nx 10 -Tstop 1.e0 -basis 2 -ot 3 -NTH 10 -Chy 0.1 -Cdif 0.006 -knu 1.e-6 -output 100
      break;
    case 19: QUEST::KineticLinearDiffusion1D_DGIMEX_Schur_acctest(argc, argv); 
      // ./main -m 6 -nx 10 -Tstop 0.e0 -basis 0 -ot 1 -NTH 1 -knu 1.e-6 -schur 4 -schurtol 1.e-13 -output 100
      break;
    case 20: QUEST::KineticLinearDiffusion1D_DGIMEX_RiemannProblem(argc, argv); 
      // ./main -x1 -1.e0 -x2 1.e0 -nx 100 -Tstop 0.15e0 -basis 2 -ot 3 -NTH 10 -Chy 0.1 -Cdif 0.006 -knu 0.7e0 -output 100 -CR 1.e0
      // ./main -x1 -2.e0 -x2 2.e0 -nx 100 -Tstop 0.15e0 -basis 2 -ot 3 -NTH 10 -Chy 0.1 -Cdif 0.006 -knu 1.e-6 -output 100 -CR 1.e0
      break;
    case 2001: QUEST::KineticLinearDiffusion1D_DGIMEX_IsentropicBC(argc, argv); 
      // ./main -nx 100 -nv 32 -Tstop 0.1e0 -basis 2 -ot 3 -NTH 10 -Chy 0.1 -Cdif 0.006 -knu 1.e0 -output 100 -CR 1.e0
      // ./main -nx 100 -nv 32 -Tstop 0.15e0 -basis 2 -ot 3 -NTH 10 -Chy 0.1 -Cdif 0.006 -knu 1.e-6 -output 100 -CR 1.e0
      break;
    case 21: QUEST::KineticLinearDiffusion1D_DGIMEX_Schur_RiemannProblem(argc, argv); 
      // ./main -x1 -1.e0 -x2 1.e0 -nx 100 -Tstop 0.15e0 -basis 2 -ot 3 -NTH 10 -knu 0.7e0 -schur 3 -schurtol 1.e-13 -output 100 -CR 1.e0
      // ./main -x1 -2.e0 -x2 2.e0 -nx 100 -Tstop 0.15e0 -basis 2 -ot 3 -NTH 10 -knu 1.e-6 -schur 3 -schurtol 1.e-13 -output 100 -CR 1.e0
      break;
    case 2101: QUEST::KineticLinearDiffusion1D_DGIMEX_Schur_IsentropicBC(argc, argv); 
      // ./main -nx 40 -nv 16 -Tstop_vec "0.1e0 0.5e0 1.e0 1.6e0 4.e0" -basis 2 -ot 3 -NTH 10 -knu 1.e0 -schur 4 -schurtol 1.e-13 -output 100 -CR 1.e0
      // ./main -nx 1000 -nv 16 -Tstop_vec "0.15e0 0.25e0 2.e0" -basis 2 -ot 3 -NTH 10 -knu 1.e-8 -schur 4 -schurtol 1.e-13 -output 100 -CR 1.e0
      break;
    case 22: QUEST::KineticDriftDiffusion1D_DGIMEX_Schur_Unipolar(argc, argv); 
      // ./main -m 1 -x1 0.e0 -x2 1.e0 -nx 100 -v1 -10 -v2 10 -nv 100 -Tstop_vec "0.0e0" -basis 2 -ot 3 -NTH 10 -knu 1.e-6 -schur 4 -schurtol 1.e-13 -output 100 -itertol 1.e-13
      break;
    case 2201: QUEST::KineticDriftDiffusion1D_DGIMEX_Schur_acctest(argc, argv); 
      // ./main -m 5 -nx 20 -v1 -10 -v2 10 -nv 40 -Tstop_vec "0.5e0" -basis 0 -ot 1 -NTH 1 -knu 1.e-6 -schur 4 -schurtol 1.e-10 -output 1 -itertol 1.e-9 -gamma 0.8 -plot 0
      break;
    case 2202: QUEST::KineticLinearDiffusion1D_DGIMEX_Schur_withMaxwell_acctest(argc, argv); 
      // ./main -m 4 -nx 20 -v1 -10 -v2 10 -nv 40 -Tstop_vec "0.1e0" -basis 2 -ot 3 -NTH 10 -knu 1.e-6 -schur 4 -schurtol 1.e-12 -output 100 -plot 0
      break;
    case 801: QUEST::KineticDriftDiffusion1D_SL_acctest_Different_Gamma_conservative(argc, argv); break;
    case 802: QUEST::KineticDriftDiffusion1D_SL_Conservative_Lomac_acctest(argc, argv); 
      /* ./main -m 5 -nx 20 -nv 40 -Tstop 1.e0 -ox 1 -ot 1 -gamma 1.e0 -NTH 10 -theta 1.e0 -sigmas 1.e0 -knu 1.e-6 -time_step 1 -dd 3 -poi 3*/
      break;
    case 803: QUEST::KineticDriftDiffusion1D_SL_Conservative_Lomac_different_gamma_acctest(argc, argv); 
      /* ./main -m 5 -nx 20 -nv 40 -Tstop 1.e0 -ox 1 -ot 1 -gamma 1.e0 -NTH 10 -theta 1.e0 -sigmas 1.e0 -knu 1.e-6 -time_step 1 -dd 3 -poi 3*/
      break;
    case 804: QUEST::KineticDriftDiffusion1D_SL_Conservative_Lomac_PNjunction(argc, argv); 
      /* ./main -m 5 -nx 20 -nv 40 -Tstop 1.e0 -ox 1 -ot 1 -gamma 1.e0 -NTH 10 -theta 1.e0 -sigmas 1.e0 -knu 1.e-6 -time_step 1 -dd 3 -poi 2*/
      break;
    default:
      QUEST_ERROR("Invalid problem kind ");
      return 1;
  }
  TOC;

}

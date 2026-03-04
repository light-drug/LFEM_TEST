#include "srctest_1.hpp"

namespace QUEST
{

void introduction() {
  std::cout << Color::blue << " 1 :" << Color::reset
            << " Kinetic Linear Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (accuracy test)" << Color::reset << std::endl;

  std::cout << Color::blue << " 2 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (accuracy test)" << Color::reset << std::endl;

  std::cout << Color::blue << " 3 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian unipolar"
            << Color::red << " (standard)" << Color::reset << std::endl;

  std::cout << Color::blue << " 301 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian PN junction"
            << Color::red << " (standard)" << Color::reset << std::endl;

  std::cout << Color::blue << " 4 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian unipolar"
            << Color::red << " (Newton iteration)" << Color::reset << std::endl;

  std::cout << Color::blue << " 5 :" << Color::reset
            << " Kinetic Linear Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (two-velocity)" << Color::reset << std::endl;
  
  std::cout << Color::blue << " 501 :" << Color::reset
            << " Kinetic Linear Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (isotropic boundary condition)" << Color::reset << std::endl;

  std::cout << Color::blue << " 6 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian unipolar"
            << Color::red << " (Quasineutral)" << Color::reset << std::endl;

  std::cout << Color::blue << " 7 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian acctest"
            << Color::red << " (Quasineutral)" << Color::reset << std::endl;

  std::cout << Color::blue << " 8 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (Different Gamma)" << Color::reset << std::endl;

  std::cout << Color::blue << " 9 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian PN junction"
            << Color::red << " (Quasineutral)" << Color::reset << std::endl;

  std::cout << Color::blue << " 10 :" << Color::reset
            << " Kinetic Linear Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (two-velocity, Implicit in micro equation, different with scheme 5)"
            << Color::reset << std::endl;

  std::cout << Color::blue << " 11 :" << Color::reset
            << " Drift Diffusion 1D Finite Difference unipolar"
            << Color::red << " (Only the Macro solver)" << Color::reset << std::endl;

  std::cout << Color::blue << " 12 :" << Color::reset
            << " Drift Diffusion 1D Finite Difference unipolar with penalty iteration"
            << Color::red << " (Only the Macro solver, maximum principle satisfying)" << Color::reset << std::endl;

  std::cout << Color::blue << " 13 :" << Color::reset
            << " Drift Diffusion 1D Finite Difference unipolar with penalty Implicit"
            << Color::red << " (Only the Macro solver)" << Color::reset << std::endl;

  std::cout << Color::blue << " 14 :" << Color::reset
            << " Drift Diffusion 1D Finite Difference acctest"
            << Color::red << " (Only the Macro solver)" << Color::reset << std::endl;

  std::cout << Color::blue << " 15 :" << Color::reset
            << " Drift Diffusion 1D Finite Difference PNjunction"
            << Color::red << " (Only the Macro solver)" << Color::reset << std::endl;

  std::cout << Color::blue << " 16 :" << Color::reset
            << " Drift Diffusion 1D Finite Difference PNjunction with penalty Implicit"
            << Color::red << " (Only the Macro solver)" << Color::reset << std::endl;

  std::cout << Color::blue << " 17 :" << Color::reset
            << " DG Poisson accuracy test"
            << Color::red << " (Dirichlet condition)" << Color::reset << std::endl;

  std::cout << Color::blue << " 1701 :" << Color::reset
            << " DG Poisson accuracy test"
            << Color::red << " (Period condition)" << Color::reset << std::endl;

  std::cout << Color::blue << " 18 :" << Color::reset
            << " DG-IMEX Linear kinetic accuracy test " 
            << Color::red << " (Mircro-macro decomposition with periodical boundary condition)" << Color::reset << std::endl;

  std::cout << Color::blue << " 19 :" << Color::reset
            << " DG-IMEX-Schur Linear kinetic accuracy test" 
            << Color::red << " (Mircro-macro decomposition with periodical boundary condition)" << Color::reset << std::endl;

  std::cout << Color::blue << " 20 :" << Color::reset
            << " DG-IMEX Linear kinetic " 
            << Color::red << " (Mircro-macro decomposition Riemann problem rho_L = 2, rho_R = 1)" << Color::reset << std::endl;

  std::cout << Color::blue << " 2001 :" << Color::reset
            << " DG-IMEX Linear kinetic " 
            << Color::red << " (Mircro-macro decomposition Isentropic Boundary Condition)" << Color::reset << std::endl;
  
  std::cout << Color::blue << " 21 :" << Color::reset
            << " DG-IMEX-Schur Linear kinetic " 
            << Color::red << " (Mircro-macro decomposition Riemann problem rho_L = 2, rho_R = 1)" << Color::reset << std::endl;

  std::cout << Color::blue << " 2101 :" << Color::reset
            << " DG-IMEX-Schur Linear kinetic " 
            << Color::red << " (Mircro-macro decomposition Isentropic Boundary Condition)" << Color::reset << std::endl;

  std::cout << Color::blue << " 22 :" << Color::reset
            << " DG-IMEX-Schur kinetic Drift Diffusion" 
            << Color::red << " (Mircro-macro decomposition Unipolar rho+rho-rho+)" << Color::reset << std::endl;
  
  std::cout << Color::blue << " 2201 :" << Color::reset
            << " DG-IMEX-Schur kinetic Drift Diffusion" 
            << Color::red << " (Mircro-macro decomposition acctest)" << Color::reset << std::endl;

  std::cout << Color::blue << " 2202 :" << Color::reset
            << " DG-IMEX-Schur kinetic Linear Diffusion with Maxwell" 
            << Color::red << " (Mircro-macro decomposition acctest)" << Color::reset << std::endl;

  std::cout << Color::blue << " 801 :" << Color::reset
            << " Kinetic Drift Diffusion 1D Semi-Lagrangian accuracy test"
            << Color::red << " (Different Gamma and conservative test in 10 - 27)" << Color::reset << std::endl;

  std::cout << Color::blue << " 802 :" << Color::reset
            << " Kinetic Drift Diffusion 1D conservative Semi-Lagrangian accuracy test"
            << Color::red << " (example same with problem 2 !)" << Color::reset << std::endl;
  
  std::cout << Color::blue << " 803 :" << Color::reset
            << " Kinetic Drift Diffusion 1D conservative Semi-Lagrangian accuracy test with different gamma"
            << Color::red << " (example same with problem 2 !)" << Color::reset << std::endl;

  std::cout << Color::blue << " 804 :" << Color::reset
            << " Kinetic Drift Diffusion 1D conservative Semi-Lagrangian PN junction"
            << Color::red << " (example same with problem 15 and 16 !)" << Color::reset << std::endl;
  std::cout << "Please select the solver :" << std::endl;
}

void basis_test() {

  TIC;
  BasisType basistype = QUEST::BasisType::P1;
  QuadratureType quatype = QUEST::QuadratureType::kGaussLegendre;
  BasisFunction1D basis(basistype);
  Vector quax;
  Vector quaw;
  int n_qua = 4;
  QUEST::Internal::initializeGaussLegendre(n_qua, &quax, &quaw);
  Matrix v;
  basis.Map(quax, &v);
  std::cout << " v = \n" << v << std::endl;
  TOC;
  
}

void fespace_test() {
  QUEST::BasisType basistype = QUEST::BasisType::P4;
  QUEST::QuadratureType quatype = QUEST::QuadratureType::kGaussLegendre;
  int qua_order1D = 5;
  real_t x1 = 0.e0;
  real_t x2 = 2.e0;
  int xDiv = 10;
  QUEST::BoundaryType boundarytype = QUEST::BoundaryType::PeriodBoundary;

  QUEST::TensorMesh1D mesh1D(x1, x2, xDiv);
  QUEST::BasisFunction1D basis(basistype);
  QUEST::fespace1D fe(&mesh1D, &basis, qua_order1D, quatype);
  mesh1D.generateextNeighbors_period(boundarytype);
  const Eigen::MatrixXi& Tm =  fe.getTm();
  const int& NTdofs = fe.getNTdofs();
  const std::vector<Matrix>& flux_u_v = fe.getflux_u_v();
  const Matrix& dv_u = fe.getdv_u();
  const Matrix& v_u = fe.getv_u();
  const DiagnalMatrix& v_u_diag = fe.getv_u_diag();
  const DiagnalMatrix& v_u_diaginv = fe.getv_u_diaginv();
  const std::vector<Vector>& boundary_u = fe.getboundary_u();
  const Matrix& test_ref = fe.gettest_ref();
  const Matrix& test_dx_ref = fe.gettest_dx_ref();
  std::cout << v_u << std::endl;
}

void DriftDiffusion1D_FD_unipolar(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_FD_unipolar";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 50;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t gamma = 2.e-3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 0.5e0;
  real_t Tstop = 0.e0;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The Debye length gamma.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of Drift-Diffusion solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration for Drift-Diffusion solver.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "Highly associated with the time step size.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::FDmesh mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_Base poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::DriftDiffusion1D_FD DD_FD(&mesh1D, &poisolver, x_order, t_order);
    DD_FD.setcfl(cfl);
    DD_FD.setgamma(gamma);
    DD_FD.setsparsetol(spartol);
    DD_FD.setiterationtol(ddtol);
    DD_FD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_phi = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_phi" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_phi = OutputDir + "/" + filename_phi;
    std::ofstream outFile_rho(file_rho), outFile_phi(file_phi);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << DD_FD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      DD_FD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      DD_FD.updateAll(Trun, dt);
      Trun += dt;
    };

    DD_FD.update_E(Trun, 0.e0);
    const Vector& rho = DD_FD.getrho();
    const Vector& rhod = DD_FD.getrhod();
    const Vector& phi = DD_FD.getphi();

    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(phi, phi, title, outFile_phi); outFile_phi.close();
  }
}

void DriftDiffusion1D_FD_PNjunction_example(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_FD_PNjunction";  // 当前函数名
  // OutputDir = std::string("result/") + __func__ ;
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 50;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 2;
  int ddsoltype_id = 3;
  real_t gamma = 2.e-3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 0.5e0;
  real_t Tstop = 0.e0;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The Debye length gamma.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of Drift-Diffusion solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration for Drift-Diffusion solver.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "Highly associated with the time step size.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::FDmesh mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_Base poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::DriftDiffusion1D_FD_PNjunction DD_FD(&mesh1D, &poisolver, x_order, t_order);
    DD_FD.setcfl(cfl);
    DD_FD.setgamma(gamma);
    DD_FD.setsparsetol(spartol);
    DD_FD.setiterationtol(ddtol);
    DD_FD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_phi = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_phi" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_phi = OutputDir + "/" + filename_phi;
    std::ofstream outFile_rho(file_rho), outFile_phi(file_phi);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << DD_FD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      DD_FD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      DD_FD.updateAll(Trun, dt);
      Trun += dt;
    };

    DD_FD.update_E(Trun, 0.e0);
    const Vector& rho = DD_FD.getrho();
    const Vector& rhod = DD_FD.getrhod();
    const Vector& phi = DD_FD.getphi();

    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(phi, phi, title, outFile_phi); outFile_phi.close();
  }

};

void DriftDiffusion1D_FD_unipolar_penaltyIter(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_FD_unipolar_penaltyIter";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 50;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 2;
  int ddsoltype_id = 3;
  real_t gamma = 2.e-3;
  real_t mu = 1.e0;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 0.5e0;
  real_t Tstop = 0.e0;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The Debye length gamma.");
  args.AddOption(&mu, "-mu", "--mu",
                  "The penalty term.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of Drift-Diffusion solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration for Drift-Diffusion solver.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "Highly associated with the time step size.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::FDmesh mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_Base poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::DriftDiffusion1D_FD_penaltyIter DD_FD(&mesh1D, &poisolver, x_order, t_order);
    DD_FD.setcfl(cfl);
    DD_FD.setmu(mu);
    DD_FD.setgamma(gamma);
    DD_FD.setsparsetol(spartol);
    DD_FD.setiterationtol(ddtol);
    DD_FD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_phi = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_phi" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_phi = OutputDir + "/" + filename_phi;
    std::ofstream outFile_rho(file_rho), outFile_phi(file_phi);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << DD_FD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      DD_FD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      DD_FD.updateAll(Trun, dt);
      Trun += dt;
    };

    DD_FD.update_E(Trun, 0.e0);
    const Vector& rho = DD_FD.getrho();
    const Vector& rhod = DD_FD.getrhod();
    const Vector& phi = DD_FD.getphi();

    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(phi, phi, title, outFile_phi); outFile_phi.close();
  }
}

void DriftDiffusion1D_FD_unipolar_penaltyImplicit(int& argc, char *argv[]) 
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_FD_unipolar_penaltyImplicit";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 50;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 2;
  int ddsoltype_id = 3;
  real_t gamma = 2.e-3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 0.5e0;
  real_t Tstop = 0.e0;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The Debye length gamma.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of Drift-Diffusion solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration for Drift-Diffusion solver.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "Highly associated with the time step size.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::FDmesh mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_Base poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::DriftDiffusion1D_FD_penaltyImplicit DD_FD(&mesh1D, &poisolver, x_order, t_order);
    DD_FD.setcfl(cfl);
    DD_FD.setgamma(gamma);
    DD_FD.setsparsetol(spartol);
    DD_FD.setiterationtol(ddtol);
    DD_FD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_phi = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_phi" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_phi = OutputDir + "/" + filename_phi;
    std::ofstream outFile_rho(file_rho), outFile_phi(file_phi);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << DD_FD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      DD_FD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      DD_FD.updateAll(Trun, dt);
      Trun += dt;
      std::cout << " Trun = " << Trun << std::endl;
    };

    const Vector& rho = DD_FD.getrho();
    const Vector& rhod = DD_FD.getrhod();
    const Vector& phi = DD_FD.getphi();

    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(phi, phi, title, outFile_phi); outFile_phi.close();
  }
};

void DriftDiffusion1D_FD_PNjunction_penaltyImplicit_example(int& argc, char *argv[]) 
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_FD_PNjunction_penaltyImplicit";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 50;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t gamma = 2.e-3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 0.5e0;
  real_t Tstop = 0.e0;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The Debye length gamma.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of Drift-Diffusion solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration for Drift-Diffusion solver.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "Highly associated with the time step size.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::FDmesh mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_Base poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::DriftDiffusion1D_FD_PNjunction_penaltyImplicit DD_FD(&mesh1D, &poisolver, x_order, t_order);
    DD_FD.setcfl(cfl);
    DD_FD.setgamma(gamma);
    DD_FD.setsparsetol(spartol);
    DD_FD.setiterationtol(ddtol);
    DD_FD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_phi = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_phi" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_phi = OutputDir + "/" + filename_phi;
    std::ofstream outFile_rho(file_rho), outFile_phi(file_phi);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << DD_FD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      DD_FD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      DD_FD.updateAll(Trun, dt);
      Trun += dt;
      std::cout << " Trun = " << Trun << std::endl;
    };

    const Vector& rho = DD_FD.getrho();
    const Vector& rhod = DD_FD.getrhod();
    const Vector& phi = DD_FD.getphi();

    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(phi, phi, title, outFile_phi); outFile_phi.close();
  }
}

void DriftDiffusion1D_FD_acctest(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_FD_acctest";
  
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  int xDiv = 50;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t gamma = 2.e-3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 0.5e0;
  real_t Tstop = 0.e0;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The Debye length gamma.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of Drift-Diffusion solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration for Drift-Diffusion solver.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "Highly associated with the time step size.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector phiL1(m), phiL2(m), phiLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::FDmesh_period mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_Base_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::DriftDiffusion1D_FD_period DD_FD(&mesh1D, &poisolver, x_order, t_order);
    DD_FD.setcfl(cfl);
    DD_FD.setgamma(gamma);
    DD_FD.setsparsetol(spartol);
    DD_FD.setiterationtol(ddtol);
    DD_FD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_phi = "X" + std::to_string(mesh1D.getNx()) +
                              "_T" + std::to_string(Tstop) + 
                              "_phi" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_phi = OutputDir + "/" + filename_phi;
    std::ofstream outFile_rho(file_rho), outFile_phi(file_phi);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << DD_FD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      DD_FD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      DD_FD.updateAll(Trun, dt);
      Trun += dt;
    };

    DD_FD.update_E(Trun, 0.e0);
    const Vector& rho = DD_FD.getrho();
    const Vector& rhod = DD_FD.getrhod();
    const Vector& phi = DD_FD.getphi();
    DD_FD.setrhofinal(Tstop);
    const Vector& rho_exact = DD_FD.getrhofinal();
    const Vector phi_exact = DD_FD.getphifinal(Tstop);

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    mesh1D.computerrorL1(phi, phi_exact, &(phiL1(i)));
    mesh1D.computerrorL2(phi, phi_exact, &(phiL2(i)));
    mesh1D.computerrorLinf(phi, phi_exact, &(phiLinf(i)));
    
    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(phi, phi_exact, title, outFile_phi); outFile_phi.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, phiL1, phiL2, phiLinf, std::cout);
}

void KineticDriftDiffusion1D_SL_acctest(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector EL1(m), EL2(m), ELinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_period kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_E" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieDD.getcfl() << std::endl;
    kieDD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      if (kieDD.getiter() > 1000) {
        break;
      }
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    const Vector& rho = kieDD.getrho();
    const Vector& rho_exact = kieDD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    kieDD.update_E(Tstop);
    Vector E_exact = kieDD.getEfinal(Tstop);
    std::vector<Vector> f_exact = kieDD.getffinal(Tstop);
    mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
    if (kieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, EL1, EL2, ELinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
}

void KineticDriftDiffusion1D_SL_acctest_Different_Gamma(int& argc, char *argv[]) 
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_acctest";
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-9; // 隐式格式每一步误差
  real_t ddtol = 1.e-9; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  // 这里要对不同的gamma进行测试
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma) + "_knu" + std::to_string(knu);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector EL1(m), EL2(m), ELinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_period_Different_Gamma kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_E" + ".dat";
    std::string filename_totaldensity = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_totaldensity" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_totaldensity = OutputDir + "/" + filename_totaldensity;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    std::ofstream outFile_totaldensity(file_totaldensity);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieDD.getcfl() << std::endl;
    kieDD.setrhofinal(Tstop);
    std::vector<real_t> totaldensity;
    std::vector<real_t> timevec;
    totaldensity.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    timevec.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    const Vector& rho = kieDD.getrho();
    totaldensity.push_back(rho.sum() * mesh1D.gethx());
    timevec.push_back(Trun);
    

    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      

      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      totaldensity.push_back(rho.sum() * mesh1D.gethx());
      timevec.push_back(Trun);
      if (kieDD.getiter() > 1000) {
        break;
      }
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    for (int j = 0; j < totaldensity.size(); j++) 
    {
      outFile_totaldensity << timevec[j] << " " <<  totaldensity[j] << std::endl;
    }
    outFile_totaldensity.close();

    const Vector& rho_exact = kieDD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    kieDD.update_E(Tstop);
    Vector E_exact = kieDD.getEfinal(Tstop);
    std::vector<Vector> f_exact = kieDD.getffinal(Tstop);
    mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
    if (kieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, EL1, EL2, ELinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
};

void KineticDriftDiffusion1D_SL_acctest_Different_Gamma_conservative(int& argc, char *argv[]) 
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DDc_acctest";
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-9; // 隐式格式每一步误差
  real_t ddtol = 1.e-9; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  // 这里要对不同的gamma进行测试
  OutputDir = OutputDir + "_gamma" + std::to_string(gamma) + "_knu" + std::to_string(knu);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector EL1(m), EL2(m), ELinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_period_Different_Gamma_conservative 
      kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_E" + ".dat";
    std::string filename_totaldensity = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_totaldensity" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_totaldensity = OutputDir + "/" + filename_totaldensity;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    std::ofstream outFile_totaldensity(file_totaldensity);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieDD.getcfl() << std::endl;
    kieDD.setrhofinal(Tstop);
    std::vector<real_t> totaldensity;
    std::vector<real_t> timevec;
    totaldensity.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    timevec.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    const Vector& rho = kieDD.getrho();
    totaldensity.push_back(rho.sum() * mesh1D.gethx());
    timevec.push_back(Trun);
    

    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      

      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      totaldensity.push_back(rho.sum() * mesh1D.gethx());
      timevec.push_back(Trun);
      if (kieDD.getiter() > 1000) {
        break;
      }
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    for (int j = 0; j < totaldensity.size(); j++) 
    {
      outFile_totaldensity << timevec[j] << " " <<  totaldensity[j] << std::endl;
    }
    outFile_totaldensity.close();

    const Vector& rho_exact = kieDD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    kieDD.update_E(Tstop);
    Vector E_exact = kieDD.getEfinal(Tstop);
    std::vector<Vector> f_exact = kieDD.getffinal(Tstop);
    mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
    if (kieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, EL1, EL2, ELinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
};

void KineticLinearDiffusion1D_SL_acctest(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Linear-Diffusion equation !";
  std::string OutputDir = "result/LD_acctest";
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 1000;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ldsoltype_id = 0;
  real_t poitol = 1.e-9; // 泊松方程误差
  real_t spartol = 1.e-9; // 隐式格式每一步误差
  real_t ldtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int NTH = 10;
  
  int m = 5;
  
  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ldsoltype_id, "-ld", "--ldsoltype",
                  "The macro linear diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  title = title + " T = " + std::to_string(Tstop);
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ldsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ldsoltype_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticLinearDiffusion1D_SL_period kieLD(&mesh1D, &poisolver, x_order, t_order);
    kieLD.settheta(theta);
    kieLD.setcfl(cfl);
    kieLD.setgamma(gamma);
    kieLD.seteps(knu);
    kieLD.setsparsetol(spartol);
    kieLD.setiterationtol(ldtol);
    kieLD.setD(theta);
    kieLD.setNTH(NTH);
    kieLD.init(ldsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho) , outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieLD.getcfl() << std::endl;
    kieLD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = kieLD.getE();
    
    while (Trun < Tstop) {
      kieLD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieLD.updateAll(Trun, dt);
      Trun += dt;
    };

    const Vector& rho = kieLD.getrho();
    const Vector& rho_exact = kieLD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    const std::vector<Vector>& f = kieLD.getf();
    std::vector<Vector> f_exact = kieLD.getffinal(Tstop);
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));
    
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
}

void KineticLinearDiffusion1D_SL_acctest_twovel(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Linear-Diffusion equation !";
  std::string OutputDir = "result/LD_acctest_twovel";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = - std::numbers::pi;
  real_t x2 = std::numbers::pi;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  int xDiv = 100;
  int vDiv = 2;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ldsoltype_id = 3;  // 这个是纯扩散可以用CG
  real_t poitol = 1.e-9; // 泊松方程误差
  real_t spartol = 1.e-9; // 隐式格式每一步误差
  real_t ldtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int type = 1; // 0: 微观方程用半拉格式 1: 微观方程用隐式迎风
  int NTH = 1;
  
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ldsoltype_id, "-ld", "--ldsoltype",
                  "The macro linear diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&type, "-type", "--type",
                  "The type for the evolution of the micro part f.");
  args.ParseCheck(std::cout);
  TOC;
  title = title + " T = " + std::to_string(Tstop);
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ldsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ldsoltype_id);

    QUEST::KineticFDmesh_period_twovel mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticLinearDiffusion1D_SL_period_twovel kieLD(&mesh1D, &poisolver, x_order, t_order);
    kieLD.settheta(theta);
    kieLD.setcfl(cfl);
    kieLD.seteps(knu);
    kieLD.setsparsetol(spartol);
    kieLD.setiterationtol(ldtol);
    kieLD.setD(theta);
    kieLD.setNTH(NTH);
    kieLD.init(ldsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho) , outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieLD.getcfl() << std::endl;
    kieLD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      kieLD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieLD.updateAll(Trun, dt);
      Trun += dt;
    };

    const Vector& rho = kieLD.getrho();
    const Vector& rho_exact = kieLD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    const std::vector<Vector>& f = kieLD.getf();
    std::vector<Vector> f_exact = kieLD.getffinal(Tstop);
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));
    
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
}

void KineticLinearDiffusion1D_SL_IsotropicBoundary(int& argc, char *argv[]) {
  std::string title = "Unipolar diode test for Drift-Diffusion equation !";
  std::string OutputDir = "result/LD_IsotropicBoundary";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = -1.e0;
  real_t v2 = 1.e0;
  int xDiv = 100;
  int vDiv = 32;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 2;
  int ldsoltype_id = 3;
  int outputgap = 10;
  real_t poitol = 1.e-10; // 泊松每一步误差(用不到)
  real_t spartol = 1.e-10; // 隐式格式每一步误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-3;
  real_t Tstop = 0.05e0;
  real_t sigmas = 1.e0;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&ldsoltype_id, "-ld", "--ldsoltype_id",
                  "The implicit solver for the macro equation type value.");
  args.AddOption(&spartol, "-spartol", "--ldsparse_tol",
                  "The tolerance value of the implicit solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The parameter sigmas in the collosion operator.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);
  TOC;

  OutputDir = OutputDir + "_knu" +  std::to_string(knu);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  // Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  // Vector EL1(m), EL2(m), ELinf(m);
  // Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ldsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ldsoltype_id);

    QUEST::KineticFDmesh_Gauss mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticLinearDiffusion1D_SL kieLD(&mesh1D, &poisolver, x_order, t_order);
    kieLD.setcfl(cfl);
    kieLD.seteps(knu);
    kieLD.setsparsetol(spartol);
    kieLD.setsigmas(sigmas);
    kieLD.setNTH(NTH);
    
    kieLD.init(ldsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) +
                              "_T" + std::to_string(Tstop) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_f(file_f);

    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieLD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    const std::vector<Vector>& f = kieLD.getf();
    while (Trun < Tstop) {
      kieLD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      kieLD.updateAll(Trun, dt);

      Trun += dt;
      step_count += 1;
      if (outputgap > 0)
      {
        if (step_count % outputgap == 0) 
        {
          std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
        }
      };
    };

    const Vector& rho = kieLD.getrho();
    Vector rho_init = kieLD.getrhoinit();
    std::vector<Vector> f_init = kieLD.getfinit();

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_init, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(f, f_init, title, outFile_f); outFile_f.close();
    // if (kieDD.getiter() > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
    //   std::abort();
    // }
  }
}

void KineticLinearDiffusion1D_SL_Implicit_acctest_twovel(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Linear-Diffusion equation !";
  std::string OutputDir = "result/LD_acctest_twovel";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = - std::numbers::pi;
  real_t x2 = std::numbers::pi;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  int xDiv = 100;
  int vDiv = 2;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ldsoltype_id = 0;  // 这个是纯扩散可以用CG
  real_t poitol = 1.e-9; // 泊松方程误差
  real_t spartol = 1.e-9; // 隐式格式每一步误差
  real_t ldtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int type = 1; // 0: 微观方程用半拉格式 1: 微观方程用隐式迎风
  int NTH = 1;
  
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ldsoltype_id, "-ld", "--ldsoltype",
                  "The macro linear diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&type, "-type", "--type",
                  "The type for the evolution of the micro part f.");
  args.ParseCheck(std::cout);
  TOC;
  title = title + " T = " + std::to_string(Tstop);
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ldsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ldsoltype_id);

    QUEST::KineticFDmesh_period_twovel mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticLinearDiffusion1D_SL_period_twovel_v2 kieLD(&mesh1D, &poisolver, x_order, t_order);
    kieLD.settheta(theta);
    kieLD.setcfl(cfl);
    kieLD.seteps(knu);
    kieLD.setsparsetol(spartol);
    kieLD.setiterationtol(ldtol);
    kieLD.setD(theta);
    kieLD.setNTH(NTH);
    kieLD.init(ldsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho) , outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieLD.getcfl() << std::endl;
    kieLD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;

    while (Trun < Tstop) {
      kieLD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieLD.updateAll(Trun, dt);
      Trun += dt;
    };

    const Vector& rho = kieLD.getrho();
    const Vector& rho_exact = kieLD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    const std::vector<Vector>& f = kieLD.getf();
    std::vector<Vector> f_exact = kieLD.getffinal(Tstop);
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));
    
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
}

} // namespace QUEST







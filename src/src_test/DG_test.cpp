#include "srctest_1.hpp"

namespace QUEST
{

void Poisson1D_DG_acctest(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Poisson equation using DG method !";
  std::string OutputDir = "result_DG/Poisson1D_DG_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  int xDiv = 100;
  int basis_id = 2;
  int poisoltype_id = 3;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  real_t poitol = 1.e-9; // 泊松方程误差
  int m = 3;
  real_t C11 = 1.e0;
  real_t C12 = 0.5e0;
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&C11, "-c11", "--c11",
                  "The C11 parameter in the LDG method for Poisson.");
  args.AddOption(&C12, "-c12", "--c12",
                  "The C12 parameter in the LDG method for Poisson."); 
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;

  IntVector xDiv_vec(m);
  Vector uL1(m), uL2(m), uLinf(m);
  Vector qL1(m), qL2(m), qLinf(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DType poisoltype = 
          static_cast<QUEST::PoissonSolver1DType>(poisoltype_id);
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::TensorMesh1D mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.init();
    PoissonSolver1DParameter poi_param;
    poi_param.C11 = C11;
    poi_param.C12 = C12;
    QUEST::Poisson_acctest_1D poi(&fe, poi_param);
    poi.init(poisoltype, poitol);

    std::string filename_u = "X" + std::to_string(mesh1D.getncell()) + 
                              "_o" + std::to_string(basis.getpolydim()) +
                              "_u" + ".dat";
    std::string filename_q = "X" + std::to_string(mesh1D.getncell()) +
                              "_o" + std::to_string(basis.getpolydim()) +
                              "_q" + ".dat";
    std::string file_u = OutputDir + "/" + filename_u;
    std::string file_q = OutputDir + "/" + filename_q;
    std::ofstream outFile_u(file_u), outFile_q(file_q);
    
    Matrix Dirichlet_u(1,2);
    Dirichlet_u(0,0) = poi.u_bc(x1);
    Dirichlet_u(0,1) = poi.u_bc(x2);

    Matrix RHS_nodal;
    fe.Interpolate_Initial(
      [&](const Matrix& x) { return poi.RHS(x);}, &RHS_nodal);

    Vector bg, bf;
    poi.Assemble_bg(Dirichlet_u, &bg);
    poi.Assemble_bf(Dirichlet_u, RHS_nodal, &bf);

    Matrix u_real_modal, u_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x) { return poi.u_real(x);}, &u_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x) { return poi.u_real(x);}, &u_real_nodal);

    Matrix q_real_modal, q_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x) { return poi.u_real_dx(x);}, &q_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x) { return poi.u_real_dx(x);}, &q_real_nodal);

    std::vector<Matrix> u_numerical_modal, u_numerical_nodal;
    poi.solveall(bg, bf, &u_numerical_modal);
    fe.modal_to_nodal1D(u_numerical_modal, &u_numerical_nodal);

    // 计算误差 //
    fe.computerrorL1(u_numerical_nodal[0], u_real_nodal, &(uL1(i)));
    // std::cout << "uL1(" << i << ") = " << uL1(i) << std::endl;
    // PAUSE(); PAUSE();
    fe.computerrorL2(u_numerical_nodal[0], u_real_nodal, &(uL2(i)));
    fe.computerrorLinf(u_numerical_nodal[0], u_real_nodal, &(uLinf(i)));

    fe.computerrorL1(u_numerical_nodal[1], q_real_nodal, &(qL1(i)));
    fe.computerrorL2(u_numerical_nodal[1], q_real_nodal, &(qL2(i)));
    fe.computerrorLinf(u_numerical_nodal[1], q_real_nodal, &(qLinf(i)));

    Vector u_numerical_aver = u_numerical_modal[0].row(0);
    Vector u_exact_aver = u_real_modal.row(0);
    Vector q_numerical_aver = u_numerical_modal[1].row(0);
    Vector q_exact_aver = q_real_modal.row(0);
    fe.DisplayResult(u_numerical_aver, 
                        u_exact_aver, title, outFile_u); outFile_u.close();
    fe.DisplayResult(q_numerical_aver, 
                        q_exact_aver, title, outFile_q); outFile_q.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, uL1, uL2, uLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, qL1, qL2, qLinf, std::cout);
};

void Poisson1D_DG_acctest_period(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Poisson equation using DG method !";
  std::string OutputDir = "result/DG/Poisson1D_DG_acctest_period";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  int xDiv = 100;
  int basis_id = 2;
  int poisoltype_id = 3;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  real_t poitol = 1.e-9; // 泊松方程误差
  int m = 3;
  real_t C11 = 1.e0;
  real_t C12 = 0.5e0;
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&C11, "-c11", "--c11",
                  "The C11 parameter in the LDG method for Poisson.");
  args.AddOption(&C12, "-c12", "--c12",
                  "The C12 parameter in the LDG method for Poisson."); 
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;

  IntVector xDiv_vec(m);
  Vector uL1(m), uL2(m), uLinf(m);
  Vector qL1(m), qL2(m), qLinf(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::PoissonSolver1DType poisoltype = 
          static_cast<QUEST::PoissonSolver1DType>(poisoltype_id);
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::TensorMesh1D mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.init();
    PoissonSolver1DParameter poi_param;
    poi_param.C11 = C11;
    poi_param.C12 = C12;
    QUEST::Poisson_acctest_1D_period poi(&fe, poi_param);
    poi.init(poisoltype, poitol);

    std::string filename_u = "X" + std::to_string(mesh1D.getncell()) + 
                              "_o" + std::to_string(basis.getpolydim()) +
                              "_u" + ".dat";
    std::string filename_q = "X" + std::to_string(mesh1D.getncell()) +
                              "_o" + std::to_string(basis.getpolydim()) +
                              "_q" + ".dat";
    std::string file_u = OutputDir + "/" + filename_u;
    std::string file_q = OutputDir + "/" + filename_q;
    std::ofstream outFile_u(file_u), outFile_q(file_q);
    
    Matrix Dirichlet_u(1,2);
    Dirichlet_u(0,0) = poi.u_bc(x1);
    Dirichlet_u(0,1) = poi.u_bc(x2);

    Matrix RHS_nodal;
    fe.Interpolate_Initial(
      [&](const Matrix& x) { return poi.RHS(x);}, &RHS_nodal);

    Vector bg, bf;
    poi.Assemble_bg(Dirichlet_u, &bg);
    poi.Assemble_bf(Dirichlet_u, RHS_nodal, &bf);

    Matrix u_real_modal, u_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x) { return poi.u_real(x);}, &u_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x) { return poi.u_real(x);}, &u_real_nodal);

    Matrix q_real_modal, q_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x) { return poi.u_real_dx(x);}, &q_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x) { return poi.u_real_dx(x);}, &q_real_nodal);

    std::vector<Matrix> u_numerical_modal, u_numerical_nodal;
    poi.solveall(bg, bf, &u_numerical_modal);
    fe.modal_to_nodal1D(u_numerical_modal, &u_numerical_nodal);

    // 计算误差 //
    fe.computerrorL1(u_numerical_nodal[0], u_real_nodal, &(uL1(i)));
    // std::cout << "uL1(" << i << ") = " << uL1(i) << std::endl;
    // PAUSE(); PAUSE();
    fe.computerrorL2(u_numerical_nodal[0], u_real_nodal, &(uL2(i)));
    fe.computerrorLinf(u_numerical_nodal[0], u_real_nodal, &(uLinf(i)));

    fe.computerrorL1(u_numerical_nodal[1], q_real_nodal, &(qL1(i)));
    fe.computerrorL2(u_numerical_nodal[1], q_real_nodal, &(qL2(i)));
    fe.computerrorLinf(u_numerical_nodal[1], q_real_nodal, &(qLinf(i)));

    Vector u_numerical_aver = u_numerical_modal[0].row(0);
    Vector u_exact_aver = u_real_modal.row(0);
    Vector q_numerical_aver = u_numerical_modal[1].row(0);
    Vector q_exact_aver = q_real_modal.row(0);
    fe.DisplayResult(u_numerical_aver, 
                        u_exact_aver, title, outFile_u); outFile_u.close();
    fe.DisplayResult(q_numerical_aver, 
                        q_exact_aver, title, outFile_q); outFile_q.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, uL1, uL2, uLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, qL1, qL2, qLinf, std::cout);
};

void KineticLinearDiffusion1D_DGIMEX_RiemannProblem(int& argc, char *argv[])
{
  std::string title = "Riemann Problem for Poisson equation using DG-IMEX method !";
  std::string OutputDir = "result/DG/LD_DG_IMEX_RiemannProblem";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = - 1.e0;
  real_t x2 = 1.e0;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  real_t knu = 0.5e0;
  real_t Chy = 0.2e0;
  real_t Cdif = 0.1e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  real_t Tstop = 0.e0;  // 终止时刻
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 1;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&Chy, "-Chy", "--Chy",
                  "The cfl number for hyperbolic restriction.");
  args.AddOption(&Cdif, "-Cdif", "--Cdif",
                  "The cfl number for diffusion restriction.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::KineticTensorMesh1D_twovel mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::Kinetic1D_LD_DG_IMEX_EX kieLD(&mesh1D, &fe, &rk_table);
    kieLD.seteps(knu);
    kieLD.setChy(Chy);
    kieLD.setCdif(Cdif);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
    // rk_table.printall();
    // PAUSE();
    
    std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_knu" + std::to_string(knu) + 
                              "_rho" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::ofstream outFile_rho(file_rho);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;

    while (Trun < Tstop) 
    {
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
          std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
        }
      };
    };
    
    const Matrix& rho_modal = kieLD.getrho_modal();
    const std::vector<Matrix>& g_modal = kieLD.getg_modal();
    Matrix rho_numerical_nodal;
    fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
    std::vector<Matrix> g_numerical_nodal;
    fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

    Vector rho_numerical_aver = rho_modal.row(0);
    Matrix rho_init_modal;
    fe.Project_Initial(
      [&](const Matrix& x) { return kieLD.rho_init(x); }, &(rho_init_modal));
    Vector rho_init_aver = rho_init_modal.row(0);
    fe.DisplayResult(rho_numerical_aver, 
                        rho_init_aver, title, outFile_rho); outFile_rho.close();
  };
};

void KineticLinearDiffusion1D_DGIMEX_IsentropicBC(int& argc, char *argv[])
{
  std::string title = "Riemann Problem for Poisson equation using DG-IMEX method !";
  std::string OutputDir = "result/DG/LD_DG_IMEX_IsentropicBC";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  real_t knu = 0.5e0;
  real_t Chy = 0.2e0;
  real_t Cdif = 0.1e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  real_t Tstop = 0.e0;  // 终止时刻
  int xDiv = 100;
  int vDiv = 64;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 1;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&Chy, "-Chy", "--Chy",
                  "The cfl number for hyperbolic restriction.");
  args.AddOption(&Cdif, "-Cdif", "--Cdif",
                  "The cfl number for diffusion restriction.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::KineticTensorMesh1D_Gauss mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::Kinetic1D_LD_DG_IMEX_EX_IsentropicBC kieLD(&mesh1D, &fe, &rk_table);
    kieLD.seteps(knu);
    kieLD.setChy(Chy);
    kieLD.setCdif(Cdif);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
    // rk_table.printall();
    // PAUSE();
    
    std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_knu" + std::to_string(knu) + 
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::ofstream outFile_rho(file_rho);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;

    while (Trun < Tstop) 
    {
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
          std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
        }
      };
    };
    
    const Matrix& rho_modal = kieLD.getrho_modal();
    const std::vector<Matrix>& g_modal = kieLD.getg_modal();
    Matrix rho_numerical_nodal;
    fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
    std::vector<Matrix> g_numerical_nodal;
    fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

    Vector rho_numerical_aver = rho_modal.row(0);
    Matrix rho_init_modal;
    fe.Project_Initial(
      [&](const Matrix& x) { return kieLD.rho_init(x); }, &(rho_init_modal));
    Vector rho_init_aver = rho_init_modal.row(0);
    fe.DisplayResult(rho_numerical_aver, 
                        rho_init_aver, title, outFile_rho); outFile_rho.close();
  };
};

void KineticLinearDiffusion1D_DGIMEX_acctest(int& argc, char *argv[])
{
  std::string title = "Riemann Problem for Poisson equation using DG-IMEX method !";
  std::string OutputDir = "result/DG/LD_DG_IMEX_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = - std::numbers::pi;
  real_t x2 = std::numbers::pi;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  real_t knu = 0.5e0;
  real_t Chy = 0.2e0;
  real_t Cdif = 0.1e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  real_t Tstop = 0.e0;  // 终止时刻
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 3;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&Chy, "-Chy", "--Chy",
                  "The cfl number for hyperbolic restriction.");
  args.AddOption(&Cdif, "-Cdif", "--Cdif",
                  "The cfl number for diffusion restriction.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  std::vector<Vector> gL1(2), gL2(2), gLinf(2);
  for (int j = 0; j < vDiv; j++)
  {
    gL1[j].resize(m); 
    gL2[j].resize(m);
    gLinf[j].resize(m);
  };

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::KineticTensorMesh1D_twovel mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::Kinetic1D_LD_DG_IMEX_EX_twovel_period kieLD(&mesh1D, &fe, &rk_table);
    kieLD.seteps(knu);
    kieLD.setChy(Chy);
    kieLD.setCdif(Cdif);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
    // rk_table.printall();
    // PAUSE();
    
    std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_rho" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::ofstream outFile_rho(file_rho);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;

    while (Trun < Tstop) 
    {
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
          std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
        }
      };
    };
    
    const Matrix& rho_modal = kieLD.getrho_modal();
    const std::vector<Matrix>& g_modal = kieLD.getg_modal();
    Matrix rho_numerical_nodal;
    fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
    std::vector<Matrix> g_numerical_nodal;
    fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

    Matrix rho_real_nodal;
    std::vector<Matrix> g_real_nodal;
    kieLD.getrho_real_nodal(Tstop, &rho_real_nodal);
    kieLD.getg_real_nodal(Tstop, &g_real_nodal);

    // 计算误差 //
    fe.computerrorL1(rho_numerical_nodal, rho_real_nodal, &(rhoL1(i)));
    fe.computerrorL2(rho_numerical_nodal, rho_real_nodal, &(rhoL2(i)));
    fe.computerrorLinf(rho_numerical_nodal, rho_real_nodal, &(rhoLinf(i)));

    const int& Nv = mesh1D.getNv();
    for (int j = 0; j < Nv; j++)
    {
      fe.computerrorL1(g_numerical_nodal[j], g_real_nodal[j], &(gL1[j](i)));
      fe.computerrorL2(g_numerical_nodal[j], g_real_nodal[j], &(gL2[j](i)));
      fe.computerrorLinf(g_numerical_nodal[j], g_real_nodal[j], &(gLinf[j](i)));
    };
    
    Vector rho_numerical_aver = rho_modal.row(0);
    Matrix rho_real_modal;
    kieLD.getrho_real_modal(Tstop, &rho_real_modal);
    Vector rho_exact_aver = rho_real_modal.row(0);
    fe.DisplayResult(rho_numerical_aver, 
                        rho_exact_aver, title, outFile_rho); outFile_rho.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  for (int j = 0; j < vDiv; j++)
  {
    QUEST::DisplayAccTable1D(xDiv_vec, gL1[j], gL2[j], gLinf[j], std::cout);
  };
};

void KineticLinearDiffusion1D_DGIMEX_Schur_acctest(int& argc, char *argv[])
{
  std::string title = "Riemann Problem for Poisson equation using DG-IMEX method !";
  std::string OutputDir = "result/DG/LD_DG_IMEX_Schur_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = - std::numbers::pi;
  real_t x2 = std::numbers::pi;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  real_t knu = 0.5e0;
  // real_t Chy = 0.2e0;
  // real_t Cdif = 0.1e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  real_t Tstop = 0.e0;  // 终止时刻
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 3;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  int schur_solver_type_id = 3;
  real_t schur_tol = 1.e-10;

  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  // args.AddOption(&Chy, "-Chy", "--Chy",
  //                 "The cfl number for hyperbolic restriction.");
  // args.AddOption(&Cdif, "-Cdif", "--Cdif",
  //                 "The cfl number for diffusion restriction.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.AddOption(&schur_solver_type_id, "-schur", "--schur_solver_type",
                  "The solver type for the Schur complement system.");
  args.AddOption(&schur_tol, "-schurtol", "--schurtol",
                  "The sparse solver tolerance.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  std::vector<Vector> gL1(2), gL2(2), gLinf(2);
  for (int j = 0; j < vDiv; j++)
  {
    gL1[j].resize(m); 
    gL2[j].resize(m);
    gLinf[j].resize(m);
  };

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);
    QUEST::Solver1DType schur_solver_type = 
          static_cast<QUEST::Solver1DType>(schur_solver_type_id);

    QUEST::KineticTensorMesh1D_twovel mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::KineticLD_DG_IMEX_IM_Schur_twovel_period kieLD(&mesh1D, &fe, &rk_table, schur_solver_type);
    kieLD.seteps(knu);
    // kieLD.setChy(Chy);
    // kieLD.setCdif(Cdif);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    kieLD.setsparsetol(schur_tol);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
    // rk_table.printall();
    // PAUSE();
    
    std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_rho" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::ofstream outFile_rho(file_rho);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;

    while (Trun < Tstop) 
    {
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
    
    const Matrix& rho_modal = kieLD.getrho_modal();
    const std::vector<Matrix>& g_modal = kieLD.getg_modal();
    Matrix rho_numerical_nodal;
    fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
    std::vector<Matrix> g_numerical_nodal;
    fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

    Matrix rho_real_nodal;
    std::vector<Matrix> g_real_nodal;
    kieLD.getrho_real_nodal(Tstop, &rho_real_nodal);
    kieLD.getg_real_nodal(Tstop, &g_real_nodal);

    // 计算误差 //
    fe.computerrorL1(rho_numerical_nodal, rho_real_nodal, &(rhoL1(i)));
    fe.computerrorL2(rho_numerical_nodal, rho_real_nodal, &(rhoL2(i)));
    fe.computerrorLinf(rho_numerical_nodal, rho_real_nodal, &(rhoLinf(i)));

    const int& Nv = mesh1D.getNv();
    for (int j = 0; j < Nv; j++)
    {
      fe.computerrorL1(g_numerical_nodal[j], g_real_nodal[j], &(gL1[j](i)));
      fe.computerrorL2(g_numerical_nodal[j], g_real_nodal[j], &(gL2[j](i)));
      fe.computerrorLinf(g_numerical_nodal[j], g_real_nodal[j], &(gLinf[j](i)));
    };
    
    Vector rho_numerical_aver = rho_modal.row(0);
    Matrix rho_real_modal;
    kieLD.getrho_real_modal(Tstop, &rho_real_modal);
    Vector rho_exact_aver = rho_real_modal.row(0);
    fe.DisplayResult(rho_numerical_aver, 
                        rho_exact_aver, title, outFile_rho); outFile_rho.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  for (int j = 0; j < vDiv; j++)
  {
    QUEST::DisplayAccTable1D(xDiv_vec, gL1[j], gL2[j], gLinf[j], std::cout);
  };
};

void KineticLinearDiffusion1D_DGIMEX_Schur_RiemannProblem(int& argc, char *argv[])
{
  std::string title = "Riemann Problem for Poisson equation using DG-IMEX-Schur method !";
  std::string OutputDir = "result/DG/LD_DG_IMEX_Schur_RiemannProblem";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = - 1.e0;
  real_t x2 = 1.e0;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  real_t knu = 0.5e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  real_t Tstop = 0.e0;  // 终止时刻
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 1;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  int schur_solver_type_id = 3;
  real_t schur_tol = 1.e-10;
  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.AddOption(&schur_solver_type_id, "-schur", "--schur_solver_type",
                  "The solver type for the Schur complement system.");
  args.AddOption(&schur_tol, "-schurtol", "--schurtol",
                  "The sparse solver tolerance.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);
    QUEST::Solver1DType schur_solver_type = 
          static_cast<QUEST::Solver1DType>(schur_solver_type_id);

    QUEST::KineticTensorMesh1D_twovel mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::KineticLD_DG_IMEX_IM_Schur_New kieLD(&mesh1D, &fe, &rk_table, schur_solver_type);
    kieLD.seteps(knu);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    kieLD.setsparsetol(schur_tol);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
    // rk_table.printall();
    // PAUSE();
    
    std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_knu" + std::to_string(knu) + 
                              "_rho" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::ofstream outFile_rho(file_rho);

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;

    while (Trun < Tstop) 
    {
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
          std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
        }
      };
    };
    
    const Matrix& rho_modal = kieLD.getrho_modal();
    const std::vector<Matrix>& g_modal = kieLD.getg_modal();
    Matrix rho_numerical_nodal;
    fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
    std::vector<Matrix> g_numerical_nodal;
    fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

    Vector rho_numerical_aver = rho_modal.row(0);
    Matrix rho_init_modal;
    fe.Project_Initial(
      [&](const Matrix& x) { return kieLD.rho_init(x); }, &(rho_init_modal));
    Vector rho_init_aver = rho_init_modal.row(0);
    fe.DisplayResult(rho_numerical_aver, 
                        rho_init_aver, title, outFile_rho); outFile_rho.close();

    std::cout << "============================================" << std::endl;
    std::cout << "Simulation ending" << std::endl;
    std::cout << "============================================" << std::endl;
  };
};

void KineticLinearDiffusion1D_DGIMEX_Schur_IsentropicBC(int& argc, char *argv[])
{
  std::string title = "Riemann Problem for Linear Kinetic equation using DG-IMEX-Schur method !";
  std::string OutputDir = "result/DG/LD_DG_IMEX_Schur_IsentropicBC";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = - 1.e0;
  real_t v2 = 1.e0;
  real_t knu = 0.5e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  // real_t Tstop = 0.e0;  // 终止时刻
  Vector Tstop_vec(1);
  Tstop_vec << 0.e0;
  int xDiv = 100;
  int vDiv = 16;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 1;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  int schur_solver_type_id = 3;
  real_t schur_tol = 1.e-10;
  
  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  // args.AddOption(&Tstop, "-Tstop", "--Tstop",
  //                 "The final time value.");
  args.AddOption(&Tstop_vec, "-Tstop_vec", "--Tstop_vec",
                  "A set of final time value.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.AddOption(&schur_solver_type_id, "-schur", "--schur_solver_type",
                  "The solver type for the Schur complement system.");
  args.AddOption(&schur_tol, "-schurtol", "--schurtol",
                  "The sparse solver tolerance.");
  args.ParseCheck(std::cout);
  IntVector xDiv_vec(m), vDiv_vec(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);
    QUEST::Solver1DType schur_solver_type = 
          static_cast<QUEST::Solver1DType>(schur_solver_type_id);

    QUEST::KineticTensorMesh1D_Gauss mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::KineticLD_DG_IMEX_IM_Schur_New_IsentropicBC kieLD(&mesh1D, &fe, &rk_table, schur_solver_type);
    kieLD.seteps(knu);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    kieLD.setsparsetol(schur_tol);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
    // rk_table.printall();
    // PAUSE();
    
    
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    real_t Tstop = 0.e0;
    for (int k = 0; k < Tstop_vec.size(); k++)
    {
      if (k != 0)
      {
        QUEST_VERIFY((Tstop_vec(k) > Tstop_vec(k-1)), 
                    " The Tstop_vec must be an increasing sequence !");
      }
      
      Tstop = Tstop_vec(k);
      std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_T" + std::to_string(Tstop) +
                              "_knu" + std::to_string(knu) + 
                              "_rho" + ".dat";
      std::string file_rho = OutputDir + "/" + filename_rho;
      std::ofstream outFile_rho(file_rho);

      while (Trun < Tstop) 
      {
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
            std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
          }
        };
      };

      const Matrix& rho_modal = kieLD.getrho_modal();
      const std::vector<Matrix>& g_modal = kieLD.getg_modal();
      Matrix rho_numerical_nodal;
      fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
      std::vector<Matrix> g_numerical_nodal;
      fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

      Vector rho_numerical_aver = rho_modal.row(0);
      Matrix rho_init_modal;
      fe.Project_Initial(
        [&](const Matrix& x) { return kieLD.rho_init(x); }, &(rho_init_modal));
      Vector rho_init_aver = rho_init_modal.row(0);
      fe.DisplayResult(rho_numerical_aver, 
                          rho_init_aver, title, outFile_rho); outFile_rho.close();

    };

    std::cout << "============================================" << std::endl;
    std::cout << "Simulation ending" << std::endl;
    std::cout << "============================================" << std::endl;
  };
};

void KineticDriftDiffusion1D_DGIMEX_Schur_acctest(int& argc, char *argv[])
{
  std::string title = "Acctest Problem for Kinetic Drift-Diffusion-Poisson equation using DG_IMEX_Schur!";
  std::string OutputDir = "result/DG/DD_DG_IMEX_Schur_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = - 10.e0;
  real_t v2 = 10.e0;
  real_t knu = 0.5e0;
  real_t theta = 1.e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t gamma = 0.8e0;
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  // real_t Tstop = 0.e0;  // 终止时刻
  Vector Tstop_vec(1);
  Tstop_vec << 0.e0;
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 3;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  int plot = 0;
  int schur_solver_type_id = 3;
  int poi_solver_type_id = 3;
  real_t C11 = 1.e0;
  real_t C12 = 0.5e0;
  real_t schur_tol = 1.e-10;
  real_t poi_tol = 1.e-10;
  real_t iter_tol = 1.e-10;

  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  // args.AddOption(&Tstop, "-Tstop", "--Tstop",
  //                 "The final time value.");
  args.AddOption(&Tstop_vec, "-Tstop_vec", "--Tstop_vec",
                  "A set of final time value.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The gamma in the collision term Poisson equation.");
  args.AddOption(&theta, "-theta", "--theta",
                  "The theta parameter in the Maxwell distribution.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&C11, "-C11", "--C11",
                  "The penalty term in the Poisson equation.");
  args.AddOption(&C12, "-C12", "--C12",
                  "The alternating flux in the Poisson equation.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.AddOption(&plot, "-plot", "--if_plot",
                  "The plot bool variable.");
  args.AddOption(&schur_solver_type_id, "-schur", "--schur_solver_type",
                  "The solver type for the Schur complement system.");
  args.AddOption(&poi_solver_type_id, "-poi", "--poi_solver_type",
                  "The solver type for the Poisson equation.");
  args.AddOption(&schur_tol, "-schurtol", "--schurtol",
                  "The sparse solver tolerance.");
  args.AddOption(&iter_tol, "-itertol", "--itertol",
                  "The iteration tolerance for the drift diffusion iteration ");
  args.AddOption(&poi_tol, "-poitol", "--poitol",
                  "The Poisson tolerance for the drift diffusion equation ");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);
    QUEST::Solver1DType schur_solver_type = 
          static_cast<QUEST::Solver1DType>(schur_solver_type_id);
    QUEST::PoissonSolver1DType poi_solver_type = 
          static_cast<QUEST::PoissonSolver1DType>(poi_solver_type_id);
    
    PoissonSolver1DParameter pa;
    pa.C11 = C11;
    pa.C12 = C12;

    QUEST::KineticTensorMesh1D mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::PoissonSolver1D_period poi_solver(&fe, pa);
    poi_solver.init(poi_solver_type, poi_tol);
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    rk_table.printall();
    PAUSE();
    QUEST::KineticDriftD_DG_IMEX_IM_Schur_period kieDD(&mesh1D, &fe, &rk_table, 
                                &poi_solver, schur_solver_type);
    kieDD.seteps(knu);
    kieDD.setgamma(gamma);
    kieDD.setCR(CR);
    kieDD.setbeta1(beta);
    kieDD.setsigmas(sigmas);
    kieDD.settheta(theta);
    kieDD.setsparsetol(schur_tol);
    kieDD.setitertol(iter_tol);
    kieDD.setNTH(NTH);
    std::cout << "============================================" << std::endl;
    kieDD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
        
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    real_t Tstop = 0.e0;
    std::ofstream outFile_rho;
    for (int k = 0; k < Tstop_vec.size(); k++)
    {
      if (k != 0)
      {
        QUEST_VERIFY((Tstop_vec(k) > Tstop_vec(k-1)), 
                    " The Tstop_vec must be an increasing sequence !");
      }
      
      Tstop = Tstop_vec(k);
      std::cout << " Tstop = " << Tstop << std::endl;
      if (plot > 0)
      {
        std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                                "_ox" + std::to_string(basis.getpolydim()) +
                                "_ot" + std::to_string(t_order) +
                                "_knu" + std::to_string(knu) +
                                "_T" + std::to_string(Tstop) +
                                "_rho" + ".dat";
        std::string file_rho = OutputDir + "/" + filename_rho;
        outFile_rho.open(file_rho);
      };
      while (Trun < Tstop) 
      {
        kieDD.setdt(&dt);
        if (Trun + dt > Tstop) {
            dt = Tstop - Trun;
        };
        if (outputgap > 0)
        {
          if (step_count % outputgap == 0) 
          {
            std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
          }
        };        
        kieDD.updateAll(Trun, dt);
        Trun += dt;
        step_count += 1;
      };
      
      const Matrix& rho_modal = kieDD.getrho_modal();
      const std::vector<Matrix>& g_modal = kieDD.getg_modal();
      Matrix rho_numerical_nodal;
      fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
      std::vector<Matrix> g_numerical_nodal;
      fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);
      
      Matrix rho_real_nodal;
      std::vector<Matrix> g_real_nodal;
      kieDD.getrho_real_nodal(Tstop, &rho_real_nodal);
      kieDD.getg_real_nodal(Tstop, &g_real_nodal);
      
      // 计算误差 //
      fe.computerrorL1(rho_numerical_nodal, rho_real_nodal, &(rhoL1(i)));
      fe.computerrorL2(rho_numerical_nodal, rho_real_nodal, &(rhoL2(i)));
      fe.computerrorLinf(rho_numerical_nodal, rho_real_nodal, &(rhoLinf(i)));

      if (plot > 0) 
      {
        Vector rho_numerical_aver = rho_modal.row(0);
        Matrix rho_real_modal;
        kieDD.getrho_real_modal(Tstop, &rho_real_modal);
        Vector rho_real_aver = rho_real_modal.row(0);
        fe.DisplayResult(rho_numerical_aver, 
                            rho_real_aver, title, outFile_rho); outFile_rho.close();
      };
    }; // for

    std::cout << "============================================" << std::endl;
    std::cout << "Simulation ending for " << mesh1D.getxDiv() << std::endl;
    std::cout << "============================================" << std::endl;

  };
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
};

void KineticDriftDiffusion1D_DGIMEX_Schur_Unipolar(int& argc, char *argv[])
{
  std::string title = "Unipolar Problem for Kinetic Drift-Diffusion-Poisson equation using DG_IMEX_Schur!";
  std::string OutputDir = "result/DG/DD_DG_IMEX_Schur_Unipolar";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = - 10.e0;
  real_t v2 = 10.e0;
  real_t knu = 0.5e0;
  real_t theta = 1.e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t gamma = 2.e-3;
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  // real_t Tstop = 0.e0;  // 终止时刻
  Vector Tstop_vec(1);
  Tstop_vec << 0.e0;
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 3;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  int schur_solver_type_id = 3;
  int poi_solver_type_id = 0;
  real_t C11 = 1.e0;
  real_t C12 = 0.5e0;
  real_t schur_tol = 1.e-10;
  real_t poi_tol = 1.e-10;
  real_t iter_tol = 1.e-10;
  int plot = 0;

  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  // args.AddOption(&Tstop, "-Tstop", "--Tstop",
  //                 "The final time value.");
  args.AddOption(&Tstop_vec, "-Tstop_vec", "--Tstop_vec",
                  "A set of final time value.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The gamma in the collision term Poisson equation.");
  args.AddOption(&theta, "-theta", "--theta",
                  "The theta parameter in the Maxwell distribution.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&C11, "-C11", "--C11",
                  "The penalty term in the Poisson equation.");
  args.AddOption(&C12, "-C12", "--C12",
                  "The alternating flux in the Poisson equation.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.AddOption(&schur_solver_type_id, "-schur", "--schur_solver_type",
                  "The solver type for the Schur complement system.");
  args.AddOption(&poi_solver_type_id, "-poi", "--poi_solver_type",
                  "The solver type for the Poisson equation.");
  args.AddOption(&schur_tol, "-schurtol", "--schurtol",
                  "The sparse solver tolerance.");
  args.AddOption(&iter_tol, "-itertol", "--itertol",
                  "The iteration tolerance for the drift diffusion iteration ");
  args.AddOption(&poi_tol, "-poitol", "--poitol",
                  "The Poisson tolerance for the drift diffusion equation ");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);
  // Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  // std::vector<Vector> gL1(2), gL2(2), gLinf(2);
  // for (int j = 0; j < vDiv; j++)
  // {
  //   gL1[j].resize(m); 
  //   gL2[j].resize(m);
  //   gLinf[j].resize(m);
  // };

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);
    QUEST::Solver1DType schur_solver_type = 
          static_cast<QUEST::Solver1DType>(schur_solver_type_id);
    QUEST::PoissonSolver1DType poi_solver_type = 
          static_cast<QUEST::PoissonSolver1DType>(poi_solver_type_id);
    
    PoissonSolver1DParameter pa;
    pa.C11 = C11;
    pa.C12 = C12;

    QUEST::KineticTensorMesh1D mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::PoissonSolver1D poi_solver(&fe, pa);
    poi_solver.init(poi_solver_type, poi_tol);
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    QUEST::KineticDriftD_DG_IMEX_IM_Schur kieDD(&mesh1D, &fe, &rk_table, 
                                &poi_solver, schur_solver_type);
    kieDD.seteps(knu);
    kieDD.setgamma(gamma);
    kieDD.setCR(CR);
    kieDD.setbeta1(beta);
    kieDD.setsigmas(sigmas);
    kieDD.settheta(theta);
    kieDD.setsparsetol(schur_tol);
    kieDD.setitertol(iter_tol);
    kieDD.setNTH(NTH);
    std::cout << "============================================" << std::endl;
    kieDD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
        
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    real_t Tstop = 0.e0;

    for (int k = 0; k < Tstop_vec.size(); k++)
    {
      if (k != 0)
      {
        QUEST_VERIFY((Tstop_vec(k) > Tstop_vec(k-1)), 
                    " The Tstop_vec must be an increasing sequence !");
      }
      
      Tstop = Tstop_vec(k);
      std::cout << " Tstop = " << Tstop << std::endl;
      std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                              "_ox" + std::to_string(basis.getpolydim()) +
                              "_ot" + std::to_string(t_order) +
                              "_knu" + std::to_string(knu) +
                              "_T" + std::to_string(Tstop) +
                              "_rho" + ".dat";
      std::string file_rho = OutputDir + "/" + filename_rho;
      std::ofstream outFile_rho(file_rho);

      while (Trun < Tstop) 
      {
        kieDD.setdt(&dt);
        if (Trun + dt > Tstop) {
            dt = Tstop - Trun;
        };
        if (outputgap > 0)
        {
          if (step_count % outputgap == 0) 
          {
            std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
          }
        };
        PAUSE();
        kieDD.updateAll(Trun, dt);
        Trun += dt;
        step_count += 1;
      };
    
      const Matrix& rho_modal = kieDD.getrho_modal();
      const std::vector<Matrix>& g_modal = kieDD.getg_modal();
      Matrix rho_numerical_nodal;
      fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
      std::vector<Matrix> g_numerical_nodal;
      fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);

      
      Vector rho_numerical_aver = rho_modal.row(0);
      const Matrix& rho_d_modal = kieDD.getrho_d_modal();
      Vector rho_d_aver = rho_d_modal.row(0);
      fe.DisplayResult(rho_numerical_aver, 
                          rho_d_aver, title, outFile_rho); outFile_rho.close();
    }; // for

    std::cout << "============================================" << std::endl;
    std::cout << "Simulation ending for " << mesh1D.getxDiv() << std::endl;
    std::cout << "============================================" << std::endl;

  };
  // QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  // for (int j = 0; j < vDiv; j++)
  // {
  //   QUEST::DisplayAccTable1D(xDiv_vec, gL1[j], gL2[j], gLinf[j], std::cout);
  // };
};

void KineticLinearDiffusion1D_DGIMEX_Schur_withMaxwell_acctest(int& argc, char *argv[])
{
  std::string title = "Acctest Problem for Kinetic Linear-Diffusion-Poisson equation using DG_IMEX_Schur!";
  std::string OutputDir = "result/DG/LD_DG_IMEX_Schur_withMaxwell_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = - 10.e0;
  real_t v2 = 10.e0;
  real_t knu = 0.5e0;
  real_t theta = 1.e0;
  real_t beta = 0.5e0;  // 交替通量取+-0.5e0
  real_t gamma = 0.8e0;
  real_t sigmas = 1.e0; // 碰撞项大小
  real_t CR = 1.e0;     // 边界加罚大小
  // real_t Tstop = 0.e0;  // 终止时刻
  Vector Tstop_vec(1);
  Tstop_vec << 0.e0;
  int xDiv = 100;
  int vDiv = 2;
  int basis_id = 1;
  int qua_order = 4;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  int m = 3;
  int t_order = 2;
  int NTH = 10;
  int outputgap = 10;
  int schur_solver_type_id = 3;
  int poi_solver_type_id = 3;
  real_t C11 = 1.e0;
  real_t C12 = 0.5e0;
  real_t schur_tol = 1.e-10;
  real_t poi_tol = 1.e-10;
  int plot = 0;

  
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
                  "The test first vDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The order in the space direction.");  
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  // args.AddOption(&Tstop, "-Tstop", "--Tstop",
  //                 "The final time value.");
  args.AddOption(&Tstop_vec, "-Tstop_vec", "--Tstop_vec",
                  "A set of final time value.");
  args.AddOption(&CR, "-CR", "--CR",
                  "The CR number for boundary condition penalty.");
  args.AddOption(&beta, "-beta", "--beta",
                  "The alternative flux choice in the LDG scheme.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The sigmas in the collision term Q(f).");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The gamma in the collision term Poisson equation.");
  args.AddOption(&theta, "-theta", "--theta",
                  "The theta parameter in the Maxwell distribution.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&C11, "-C11", "--C11",
                  "The penalty term in the Poisson equation.");
  args.AddOption(&C12, "-C12", "--C12",
                  "The alternating flux in the Poisson equation.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.AddOption(&plot, "-plot", "--if_plot",
                  "The plot bool variable.");
  args.AddOption(&schur_solver_type_id, "-schur", "--schur_solver_type",
                  "The solver type for the Schur complement system.");
  args.AddOption(&poi_solver_type_id, "-poi", "--poi_solver_type",
                  "The solver type for the Poisson equation.");
  args.AddOption(&schur_tol, "-schurtol", "--schurtol",
                  "The sparse solver tolerance.");
  args.AddOption(&poi_tol, "-poitol", "--poitol",
                  "The Poisson tolerance for the drift diffusion equation ");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);
    QUEST::Solver1DType schur_solver_type = 
          static_cast<QUEST::Solver1DType>(schur_solver_type_id);
    QUEST::PoissonSolver1DType poi_solver_type = 
          static_cast<QUEST::PoissonSolver1DType>(poi_solver_type_id);
    
    PoissonSolver1DParameter pa;
    pa.C11 = C11;
    pa.C12 = C12;

    QUEST::KineticTensorMesh1D mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);
    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype);
    fe.setNTH(NTH);
    fe.init();
    QUEST::PoissonSolver1D_period poi_solver(&fe, pa);
    poi_solver.init(poi_solver_type, poi_tol);
    QUEST::IMEX_RK rk_table(t_order);
    rk_table.init();
    rk_table.printall();
    QUEST::KineticLinearD_DG_IMEX_IM_Schur_period kieLD(&mesh1D, &fe, &rk_table, 
                                &poi_solver, schur_solver_type);
    kieLD.seteps(knu);
    kieLD.setgamma(gamma);
    kieLD.setCR(CR);
    kieLD.setbeta1(beta);
    kieLD.setsigmas(sigmas);
    kieLD.settheta(theta);
    kieLD.setsparsetol(schur_tol);
    kieLD.setNTH(NTH);
    std::cout << "============================================" << std::endl;
    kieLD.init();
    std::cout << "initialize ending" << std::endl;
    std::cout << "============================================" << std::endl;
        
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    real_t Tstop = 0.e0;
    std::ofstream outFile_rho;

    for (int k = 0; k < Tstop_vec.size(); k++)
    {
      if (k != 0)
      {
        QUEST_VERIFY((Tstop_vec(k) > Tstop_vec(k-1)), 
                    " The Tstop_vec must be an increasing sequence !");
      }
      
      Tstop = Tstop_vec(k);
      std::cout << " Tstop = " << Tstop << std::endl;
      if (plot > 0)
      {
        std::string filename_rho = "X" + std::to_string(mesh1D.getncell()) + 
                                  "_ox" + std::to_string(basis.getpolydim()) +
                                  "_ot" + std::to_string(t_order) +
                                  "_knu" + std::to_string(knu) +
                                  "_T" + std::to_string(Tstop) +
                                  "_rho" + ".dat";
        std::string file_rho = OutputDir + "/" + filename_rho;
        outFile_rho.open(file_rho);
      };

      while (Trun < Tstop) 
      {
        kieLD.setdt(&dt);
        if (Trun + dt > Tstop) {
            dt = Tstop - Trun;
        };
        if (outputgap > 0)
        {
          if (step_count % outputgap == 0) 
          {
            std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
          }
        };        
        kieLD.updateAll(Trun, dt);
        Trun += dt;
        step_count += 1;
      };
      
      const Matrix& rho_modal = kieLD.getrho_modal();
      const std::vector<Matrix>& g_modal = kieLD.getg_modal();
      Matrix rho_numerical_nodal;
      fe.modal_to_nodal1D(rho_modal, &rho_numerical_nodal);
      std::vector<Matrix> g_numerical_nodal;
      fe.modal_to_nodal1D(g_modal, &g_numerical_nodal);
      
      Matrix rho_real_nodal;
      std::vector<Matrix> g_real_nodal;
      kieLD.getrho_real_nodal(Tstop, &rho_real_nodal);
      kieLD.getg_real_nodal(Tstop, &g_real_nodal);
      
      // 计算误差 //
      fe.computerrorL1(rho_numerical_nodal, rho_real_nodal, &(rhoL1(i)));
      fe.computerrorL2(rho_numerical_nodal, rho_real_nodal, &(rhoL2(i)));
      fe.computerrorLinf(rho_numerical_nodal, rho_real_nodal, &(rhoLinf(i)));

      if (plot > 0)
      {
        Vector rho_numerical_aver = rho_modal.row(0);
        Matrix rho_real_modal;
        kieLD.getrho_real_modal(Tstop, &rho_real_modal);
        Vector rho_real_aver = rho_real_modal.row(0);
        fe.DisplayResult(rho_numerical_aver, 
                            rho_real_aver, title, outFile_rho); outFile_rho.close();
      };
    }; // for

    std::cout << "============================================" << std::endl;
    std::cout << "Simulation ending for " << mesh1D.getxDiv() << std::endl;
    std::cout << "============================================" << std::endl;

  };
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
};

} // namespace QUEST

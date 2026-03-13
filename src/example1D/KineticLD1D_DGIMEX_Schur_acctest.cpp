#include "QUEST.hpp"

/*
  ./bin/KineticLD1D_DGIMEX_Schur_acctest -m 6 -nx 10 -Tstop 0.e0 \
  -basis 2 -ot 3 -NTH 1 -knu 1.e-6 \
  -schur 4 -schurtol 1.e-13 -output 100
*/

int main(int argc, char *argv[])
{
  std::string title = "Two velocity telegraph equation using DG-IMEX-Schur method !";
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
  return 0;
};

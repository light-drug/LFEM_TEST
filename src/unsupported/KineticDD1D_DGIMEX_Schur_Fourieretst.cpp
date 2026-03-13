#include "QUEST.hpp"

int main(int argc, char *argv[])
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
  real_t timeratio = 0.8e0;
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
  args.AddOption(&timeratio, "-timeratio", "timeratio",
                  "The time ratio for the sake of stability.");
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




}
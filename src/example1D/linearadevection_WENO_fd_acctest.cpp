#include "QUEST.hpp"

#include <filesystem>

/*
  ./bin/linearadevection_WENO_fd_acctest -m 4 -nx 40 -ox 5 -ot 3 -period 1 -Tstop 0.2
*/

int main(int argc, char *argv[])
{
  std::string title = "Accuracy test for 1D linear advection using WENO-FD";
  std::string output_dir = "result/WENO_FD/LinearAdvection1D";
  std::filesystem::create_directories(output_dir);

  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 40;
  int x_order = 5;
  int t_order = 3;
  real_t a = 1.e0;
  real_t cfl = 0.4e0;
  real_t Tstop = 0.1e0;
  int m = 4;
  int period = 0;
  int plot = 0;

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&x1, "-x1", "--x1", "Left boundary of x.");
  args.AddOption(&x2, "-x2", "--x2", "Right boundary of x.");
  args.AddOption(&xDiv, "-nx", "--xDiv", "Base number of cells.");
  args.AddOption(&x_order, "-ox", "--x_order", "WENO order in space (3 or 5).");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order in time (1~3).");
  args.AddOption(&a, "-a", "--advection_speed", "Advection speed.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&m, "-m", "--max", "Number of mesh refinements.");
  args.AddOption(&period, "-period", "--period", "Use periodic BC or not.");
  args.AddOption(&plot, "-plot", "--if_plot", "Output pointwise result if >0.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m);
  Vector uL1(m), uL2(m), uLinf(m);

  for (int i = 0; i < m; ++i) 
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;

    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();
    // rk_table.printall();
    // PAUSE();

    Vector u_num;
    Vector u_ex;

    if (period > 0) {
      QUEST::FDmesh_period mesh1D(x1, x2, xDiv_vec(i));
      mesh1D.init();
      QUEST::LinearAdvection_WENO_FD_period solver(&mesh1D, &rk_table, a, x_order);
      solver.setCFL(cfl);
      solver.init();

      real_t Trun = 0.e0;
      real_t dt = 0.e0;
      while (Trun < Tstop) 
      {
        solver.setdt(&dt);
        if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
        }
        solver.updateAll(Trun, dt);
        Trun += dt;
      }
      const Vector& xc = mesh1D.getCellCenter_vec();
      u_num = solver.getu();
      u_ex = solver.u_real(xc, Tstop);

      mesh1D.computerrorL1(u_num, u_ex, &uL1(i));
      mesh1D.computerrorL2(u_num, u_ex, &uL2(i));
      mesh1D.computerrorLinf(u_num, u_ex, &uLinf(i));

      if (plot > 0) {
        std::string filename = "WENO_period_X" + std::to_string(mesh1D.getxDiv()) +
                               "_ox" + std::to_string(x_order) + 
                               "_ot" + std::to_string(t_order) + ".dat";
        std::ofstream outFile(output_dir + "/" + filename);
        mesh1D.DisplayResult_withX(xc, u_num, u_ex, title, outFile);
      }
    } else 
    {
      QUEST::FDmesh mesh1D(x1, x2, xDiv_vec(i));
      mesh1D.init();
      QUEST::LinearAdvection_WENO_FD solver(&mesh1D, &rk_table, a, x_order);
      solver.setCFL(cfl);
      solver.init();

      real_t Trun = 0.e0;
      real_t dt = 0.e0;
      while (Trun < Tstop) {
        solver.setdt(&dt);
        if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
        }
        solver.updateAll(Trun, dt);
        Trun += dt;
      }
      const Vector& xc = mesh1D.getCellCenter_vec();
      u_num = solver.getu();
      u_ex = solver.u_real(xc, Tstop);

      mesh1D.computerrorL1(u_num, u_ex, &uL1(i));
      mesh1D.computerrorL2(u_num, u_ex, &uL2(i));
      mesh1D.computerrorLinf(u_num, u_ex, &uLinf(i));

      if (plot > 0) {
        std::string filename = "WENO_outflow_X" + std::to_string(mesh1D.getxDiv()) +
                               "_ox" + std::to_string(x_order) + 
                               "_ot" + std::to_string(t_order) + ".dat";
        std::ofstream outFile(output_dir + "/" + filename);
        mesh1D.DisplayResult_withX(xc, u_num, u_ex, title, outFile);
      }
    }
  }

  QUEST::DisplayAccTable1D(xDiv_vec, uL1, uL2, uLinf, std::cout);
  return 0;
}

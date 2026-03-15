#include "QUEST.hpp"

#include <cmath>

/*
  ./bin/Euler2D_WENO_fd_period_acctest \
  -m 5 -nx 20 -ny 20 -ox 5 -ot 3\
  -Tstop 0.1 -output 10 -cfl 0.8e0
*/

int main(int argc, char* argv[])
{
  real_t x1 = 0.e0, x2 = 2.e0 * 3.14159265358979323846;
  real_t y1 = 0.e0, y2 = 2.e0 * 3.14159265358979323846;
  int xDiv = 20, yDiv = 20;
  int x_order = 5;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.8e0;
  real_t Tstop = 0.1e0;
  int outputgap = 1;
  int m = 3;

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&xDiv, "-nx", "--xDiv", "Base x cells.");
  args.AddOption(&yDiv, "-ny", "--yDiv", "Base y cells.");
  args.AddOption(&x_order, "-ox", "--x_order", "WENO order (3 or 5).");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order.");
  args.AddOption(&gamma, "-gamma", "--gamma", "Gamma.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&m, "-m", "--max", "Refinement levels.");
  args.AddOption(&outputgap, "-output", "--outputgap", "Output every outputgap time steps if > 0.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), yDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int lv = 0; lv < m; ++lv) {

    xDiv_vec(lv) = std::pow(2, lv) * xDiv;
    yDiv_vec(lv) = std::pow(2, lv) * yDiv;

    QUEST::FDmesh2D mesh2D(x1, x2, y1, y2, xDiv_vec(lv), yDiv_vec(lv));
    mesh2D.init();

    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    QUEST::Euler2D_WENO_FD_period solver(&mesh2D, &rk_table, x_order);
    solver.setgamma(gamma);
    solver.setCFL(cfl);
    solver.init();

    TIC;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    while (Trun < Tstop) 
    {
      solver.setdt(&dt);
      if (Trun + dt > Tstop) {
        dt = Tstop - Trun;
      }
      solver.updateAll(Trun, dt);
      Trun += dt;
      ++step_count;

      if (outputgap > 0 && step_count % outputgap == 0)
      {
        std::cout << "Trun = " << Trun << ", dt = " << dt << std::endl;
      }
    }
    std::cout << " Tstop = " << Tstop 
              << " for the mesh Nx = " << mesh2D.getxDiv() 
              << ", Ny = " << mesh2D.getyDiv()
              << std::endl;
    TOC;

    const std::vector<Matrix>& cc = mesh2D.getCellCenter();
    Matrix rho_ex = solver.rho_real(cc[0], cc[1], Tstop);
    mesh2D.computerrorL1(solver.getrho(), rho_ex, &rhoL1(lv));
    mesh2D.computerrorL2(solver.getrho(), rho_ex, &rhoL2(lv));
    mesh2D.computerrorLinf(solver.getrho(), rho_ex, &rhoLinf(lv));
  }

  QUEST::DisplayAccTable2D(xDiv_vec, yDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  return 0;
}

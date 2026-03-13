#include "QUEST.hpp"

#include <cmath>

/*
  ./bin/Euler2D_WENO_fd_period_acctest -m 3 -nx 20 -ny 20 -ox 5 -ot 3 -Tstop 0.1
*/

int main(int argc, char* argv[])
{
  real_t x1 = 0.e0, x2 = 2.e0 * 3.14159265358979323846;
  real_t y1 = 0.e0, y2 = 2.e0 * 3.14159265358979323846;
  int xDiv = 20, yDiv = 20;
  int x_order = 5;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.4e0;
  real_t Tstop = 0.1e0;
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
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m), yDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int lv = 0; lv < m; ++lv) {
    xDiv_vec(lv) = std::pow(2, lv) * xDiv;
    yDiv_vec(lv) = std::pow(2, lv) * yDiv;

    QUEST::TensorMesh2D mesh2D(x1, x2, y1, y2, xDiv_vec(lv), yDiv_vec(lv));
    mesh2D.init();
    mesh2D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);

    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    QUEST::Euler2D_WENO_FD_period solver(&mesh2D, &rk_table, x_order);
    solver.setgamma(gamma);
    solver.setCFL(cfl);
    solver.init();

    real_t Trun = 0.e0, dt = 0.e0;
    while (Trun < Tstop) {
      solver.setdt(&dt);
      if (Trun + dt > Tstop) dt = Tstop - Trun;
      solver.updateAll(Trun, dt);
      Trun += dt;
    }

    const Matrix& cc = mesh2D.getCellCenter();
    Matrix xmat(yDiv_vec(lv), xDiv_vec(lv));
    Matrix ymat(yDiv_vec(lv), xDiv_vec(lv));
    for (int j = 0; j < yDiv_vec(lv); ++j) {
      for (int i = 0; i < xDiv_vec(lv); ++i) {
        const int id = mesh2D.CellIndex(i, j);
        xmat(j, i) = cc(0, id);
        ymat(j, i) = cc(1, id);
      }
    }
    Matrix rho_ex = solver.rho_real(xmat, ymat, Tstop);
    Matrix err = (solver.getrho() - rho_ex).cwiseAbs();
    const real_t area = mesh2D.gethx() * mesh2D.gethy();
    rhoL1(lv) = err.sum() * area;
    rhoL2(lv) = std::sqrt((err.array().square().sum()) * area);
    rhoLinf(lv) = err.maxCoeff();
  }

  QUEST::DisplayAccTable2D(xDiv_vec, yDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  return 0;
}

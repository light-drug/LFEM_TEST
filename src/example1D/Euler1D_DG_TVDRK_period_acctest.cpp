#include "QUEST.hpp"

#include <cmath>

/*
  ./bin/Euler1D_DG_TVDRK_period_acctest -m 4 -nx 20 -basis 2 -ot 3 -cfl 0.8 -Tstop 0.4 -output 10
*/

int main(int argc, char* argv[])
{
  real_t x1 = 0.e0;
  real_t x2 = 2.e0;
  int xDiv = 20;
  int basis_id = 1;
  int qua_order = 4;
  int quatype_id = 1;    // 积分类型
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.2e0;
  real_t Tstop = 0.1e0;
  int outputgap = 20;
  int m = 4;

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&x1, "-x1", "--x1", "Left boundary of x.");
  args.AddOption(&x2, "-x2", "--x2", "Right boundary of x.");
  args.AddOption(&xDiv, "-nx", "--xDiv", "Base number of cells.");
  args.AddOption(&basis_id, "-basis", "--basis_id", "DG basis order (0~4). Useful: 1 for P1 + limiter.");
  args.AddOption(&qua_order, "-q", "--qua_order", "Quadrature order.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id",
                  "The quadrature order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order in time (1~3).");
  args.AddOption(&gamma, "-gamma", "--gamma", "Gas constant gamma.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&m, "-m", "--max", "Number of mesh refinements.");
  args.AddOption(&outputgap, "-output", "--outputgap", "Output every outputgap time steps if >0.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; ++i)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::TensorMesh1D mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    mesh1D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);

    QUEST::BasisFunction1D basis(basistype);
    QUEST::fespace1D fe(&mesh1D, &basis, qua_order, quatype, 3);
    fe.init();

    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    QUEST::Euler1D_DG_TVDRK_period solver(&mesh1D, &fe, &rk_table);
    solver.setgamma(gamma);
    solver.setcfl(cfl);
    solver.setlimiter(false);
    solver.init();

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

    std::vector<Matrix> u_nodal;
    fe.modal_to_nodal1D(solver.getumodal(), &u_nodal);
    Matrix rho_ex;
    fe.Interpolate_Final([&solver](const Matrix& x, const real_t t) { return solver.rho_real(x, t); }, Tstop, &rho_ex);

    fe.computerrorL1(u_nodal[0], rho_ex, &rhoL1(i));
    fe.computerrorL2(u_nodal[0], rho_ex, &rhoL2(i));
    fe.computerrorLinf(u_nodal[0], rho_ex, &rhoLinf(i));
  }

  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  return 0;
}

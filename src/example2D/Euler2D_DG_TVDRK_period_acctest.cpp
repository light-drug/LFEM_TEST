#include "QUEST.hpp"

#include <cmath>

/*
  ./bin/Euler2D_DG_TVDRK_period_acctest -m 3 -nx 10 -ny 10 \
  -basis 1 -ot 3 -cfl 0.8 \
  -Tstop 0.1 -output 10
*/

int main(int argc, char* argv[])
{
  real_t x1 = 0.e0, x2 = 2.e0 * 3.14159265358979323846;
  real_t y1 = 0.e0, y2 = 2.e0 * 3.14159265358979323846;
  int xDiv = 10, yDiv = 10;
  int basis_id = 1;
  int qua_order = 4;
  int quatype_id = 1;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.8e0;
  real_t Tstop = 0.1e0;
  int m = 3;
  int outputgap = 0;

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&xDiv, "-nx", "--xDiv", "Base x cells.");
  args.AddOption(&yDiv, "-ny", "--yDiv", "Base y cells.");
  args.AddOption(&basis_id, "-basis", "--basis_id", "DG basis order.");
  args.AddOption(&qua_order, "-q", "--qua_order", "Quadrature order.");
  args.AddOption(&quatype_id, "-quatype", "--quatype_id", "Quadrature type.");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order.");
  args.AddOption(&gamma, "-gamma", "--gamma", "Gamma.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&m, "-m", "--max", "Refinement levels.");
  args.AddOption(&outputgap, "-output", "--outputgap", "Output every outputgap time steps if > 0.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m);
  IntVector yDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  int num_equations = 4;

  for (int lv = 0; lv < m; ++lv) {
    xDiv_vec(lv) = std::pow(2, lv) * xDiv;
    yDiv_vec(lv) = std::pow(2, lv) * yDiv;

    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::TensorMesh2D mesh2D(x1, x2, y1, y2, xDiv_vec(lv), yDiv_vec(lv));
    mesh2D.init();
    mesh2D.generateextNeighbors_period(QUEST::BoundaryType::PeriodBoundary);

    QUEST::BasisFunction2D basis(static_cast<QUEST::BasisType>(basis_id));
    QUEST::fespace2D fe(&mesh2D, &basis, qua_order, quatype, num_equations);
    fe.init();

    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    QUEST::Euler2D_DG_TVDRK_period solver(&mesh2D, &fe, &rk_table);
    solver.setgamma(gamma);
    solver.setcfl(cfl);
    solver.init();

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    while (Trun < Tstop) 
    {
      solver.setdt(&dt);
      if (Trun + dt > Tstop) 
      {
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
    fe.modal_to_nodal2D(solver.getumodal(), &u_nodal);
    Matrix rho_ex;
    fe.Interpolate_Final([&solver](const Matrix& x, const Matrix& y, const real_t t) {
      return solver.rho_real(x, y, t);
    }, Tstop, &rho_ex);

    fe.computerrorL1(u_nodal[0], rho_ex, &rhoL1(lv));
    fe.computerrorL2(u_nodal[0], rho_ex, &rhoL2(lv));
    fe.computerrorLinf(u_nodal[0], rho_ex, &rhoLinf(lv));
  }

  QUEST::DisplayAccTable2D(xDiv_vec, yDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  return 0;
}

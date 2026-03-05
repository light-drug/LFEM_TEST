#include "QUEST.hpp"
#include "Euler1D_WENO_fd.hpp"

#include <cmath>

class EulerPeriodicExact final : public QUEST::Euler1D_WENO_FD_period
{
public:
  using QUEST::Euler1D_WENO_FD_period::Euler1D_WENO_FD_period;

  Vector rho_exact(const Vector& x, const real_t& t) const
  {
    const real_t x1 = mesh1D_period_->getx1();
    const real_t x2 = mesh1D_period_->getx2();
    const real_t L = x2 - x1;
    Vector y = x.array() - t;
    y = (y.array() - x1 - ((y.array() - x1) / L).floor() * L + x1).matrix();
    return rho_init(y);
  }

  Vector vx_exact(const Vector& x, const real_t&) const
  {
    return Vector::Ones(x.size());
  }

  Vector pre_exact(const Vector& x, const real_t&) const
  {
    return Vector::Ones(x.size());
  }
};

int main(int argc, char* argv[])
{
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 40;
  int x_order = 5;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.4e0;
  real_t Tstop = 0.1e0;
  int m = 4;
  int outputgap = 0;

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&x1, "-x1", "--x1", "Left boundary of x.");
  args.AddOption(&x2, "-x2", "--x2", "Right boundary of x.");
  args.AddOption(&xDiv, "-nx", "--xDiv", "Base number of cells.");
  args.AddOption(&x_order, "-ox", "--x_order", "WENO order in space (3 or 5).");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order in time (1~3).");
  args.AddOption(&gamma, "-gamma", "--gamma", "Gas constant gamma.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&m, "-m", "--max", "Number of mesh refinements.");
  args.AddOption(&outputgap, "-output", "--outputgap", "Output every outputgap time steps if >0.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);

  for (int i = 0; i < m; ++i) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;

    QUEST::FDmesh_period mesh1D(x1, x2, xDiv_vec(i));
    mesh1D.init();
    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    EulerPeriodicExact solver(&mesh1D, &rk_table, x_order);
    solver.setgamma(gamma);
    solver.setCFL(cfl);
    solver.init();

    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    while (Trun < Tstop) {
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

    const Vector& xc = mesh1D.getCellCenter_vec();
    Vector rho_ex = solver.rho_exact(xc, Tstop);
    mesh1D.computerrorL1(solver.getrho(), rho_ex, &rhoL1(i));
    mesh1D.computerrorL2(solver.getrho(), rho_ex, &rhoL2(i));
    mesh1D.computerrorLinf(solver.getrho(), rho_ex, &rhoLinf(i));
  }

  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  return 0;
}

#include "QUEST.hpp"
#include "Euler1D_WENO_fd.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>

class SodEulerWENO final : public QUEST::Euler1D_WENO_FD
{
public:
  using QUEST::Euler1D_WENO_FD::Euler1D_WENO_FD;

  real_t rho_init(const real_t& x) const override { return (x < 0.5e0) ? 1.e0 : 0.125e0; }
  real_t vx_init(const real_t&) const override { return 0.e0; }
  real_t pre_init(const real_t& x) const override { return (x < 0.5e0) ? 1.e0 : 0.1e0; }

  real_t rho_bc(const real_t& x, const real_t&) const override { return rho_init(x); }
  real_t vx_bc(const real_t& x, const real_t&) const override { return vx_init(x); }
  real_t pre_bc(const real_t& x, const real_t&) const override { return pre_init(x); }
};

int main(int argc, char *argv[])
{
  std::string output_dir = "result/WENO_FD/Euler1D";
  std::filesystem::create_directories(output_dir);

  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  int xDiv = 400;
  int x_order = 5;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.4e0;
  real_t Tstop = 0.2e0;
  int outputgap = 20;
  int m = 1;

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&x1, "-x1", "--x1", "Left boundary of x.");
  args.AddOption(&x2, "-x2", "--x2", "Right boundary of x.");
  args.AddOption(&xDiv, "-nx", "--xDiv", "Number of cells.");
  args.AddOption(&x_order, "-ox", "--x_order", "WENO order in space (3 or 5).");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order in time (1~3).");
  args.AddOption(&gamma, "-gamma", "--gamma", "Gas constant gamma.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&outputgap, "-output", "--outputgap", "Output every outputgap time steps if >0.");
  args.AddOption(&m, "-m", "--max", "Number of mesh refinements.");
  args.ParseCheck(std::cout);

  for (int i = 0; i < m; ++i) {
    const int xDiv_i = std::pow(2, i) * xDiv;

    QUEST::FDmesh mesh1D(x1, x2, xDiv_i);
    mesh1D.init();
    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    SodEulerWENO solver(&mesh1D, &rk_table, x_order);
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

    std::ofstream outFile(output_dir + "/Euler1D_Riemann_final_X" + std::to_string(xDiv_i) + ".dat");
    const Vector& xc = mesh1D.getCellCenter_vec();
    const Vector& rho = solver.getrho();
    const Vector& rhovx = solver.getrhovx();
    const Vector& E = solver.getE();
    Vector vx = rhovx.array() / rho.array();
    Vector pre = solver.getpre(solver.getu());
    outFile << "TITLE = \"Euler 1D Sod Riemann\"\n";
    outFile << "VARIABLES = \"X\", \"Rho\", \"Vx\", \"P\", \"E\"\n";
    outFile << "ZONE I = " << mesh1D.getNx() << " F = POINT\n";
    for (int j = 0; j < mesh1D.getNx(); ++j) {
      outFile << xc(j) << " " << rho(j) << " " << vx(j) << " " << pre(j) << " " << E(j) << "\n";
    }
  }
  return 0;
}

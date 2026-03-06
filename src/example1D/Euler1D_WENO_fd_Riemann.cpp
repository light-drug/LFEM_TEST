#include "QUEST.hpp"
#include "Euler1D_WENO_fd.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>

/*
  ./bin/Euler1D_WENO_fd_Rieman -m 1 -nx 200 -ox 5 -ot 3 -problem 1 -Tstop 1.3
*/
namespace QUEST
{

struct RiemannInitialCondition {
  real_t x0;
  real_t rhoL, rhoR;
  real_t vxL,  vxR;
  real_t pL,   pR;
};

class Euler1D_WENO_FD_Riemann final : public Euler1D_WENO_FD {
public:
  Euler1D_WENO_FD_Riemann(const FDmesh* mesh1D,
                          const EX_TVDRK* rk_table,
                          const int x_order,
                          const RiemannInitialCondition& ic)
  : Euler1D_WENO_FD(mesh1D, rk_table, x_order), ic_(ic) {};
  ~Euler1D_WENO_FD_Riemann() override = default;

  real_t rho_init(const real_t& x) const override { return (x < ic_.x0) ? ic_.rhoL : ic_.rhoR; }
  Vector rho_init(const Vector& x) const override 
  { 
    const int Nx = x.size(); 
    Vector rho(Nx); 
    for (int i = 0; i < Nx; ++i) 
    {
      rho(i) = rho_init(x(i)); 
    } 
    return rho; 
  };
  real_t rho_bc(const real_t& x, const real_t&) const override { return rho_init(x); }

  real_t vx_init(const real_t& x) const override { return (x < ic_.x0) ? ic_.vxL  : ic_.vxR;  }
  Vector vx_init(const Vector& x) const override 
  { 
    const int Nx = x.size(); 
    Vector vx(Nx); 
    for (int i = 0; i < Nx; ++i) 
    {
      vx(i) = vx_init(x(i)); 
    } 
    return vx; 
  };
  real_t vx_bc(const real_t& x, const real_t&) const override { return vx_init(x); }
  
  real_t pre_init(const real_t& x) const override { return (x < ic_.x0) ? ic_.pL   : ic_.pR;   }
  Vector pre_init(const Vector& x) const override 
  { 
    const int Nx = x.size(); 
    Vector pre(Nx); 
    for (int i = 0; i < Nx; ++i) 
    {
      pre(i) = pre_init(x(i)); 
    } 
    return pre; 
  };
  real_t pre_bc(const real_t& x, const real_t&) const override { return pre_init(x); }
  
private:
  RiemannInitialCondition ic_;
};

}; // namespace QUEST

int main(int argc, char *argv[])
{
  std::string output_dir = "result/WENO_FD/Euler1D";
  std::filesystem::create_directories(output_dir);

  real_t x1 = -5.e0;
  real_t x2 = 5.e0;
  int xDiv = 400;
  int x_order = 5;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.4e0;
  real_t Tstop = 0.2e0;
  int outputgap = 20;
  int problem = 0; // 0: Sod, 1: Lax
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
  args.AddOption(&problem, "-problem", "--problem", "Set the type of Riemann problem.");
  args.AddOption(&m, "-m", "--max", "Number of mesh refinements.");
  args.ParseCheck(std::cout);

  IntVector xDiv_vec(m);
  for (int i = 0; i < m; ++i) 
  {
    const int xDiv_i = std::pow(2, i) * xDiv;
    xDiv_vec(i) = xDiv_i;

    QUEST::FDmesh mesh1D(x1, x2, xDiv_i);
    mesh1D.init();
    QUEST::EX_TVDRK rk_table(t_order);
    rk_table.init();

    QUEST::RiemannInitialCondition ic;
    switch(problem){
      case 0: ic = {0.e0, 1.0, 0.125, 0.0, 0.0, 1.0, 0.1}; break;
      case 1: ic = {0.e0, 0.445, 0.5, 0.698 / 0.445, 0.0, 3.528, 0.571}; break;
    }
    QUEST::Euler1D_WENO_FD_Riemann solver(&mesh1D, &rk_table, x_order, ic);
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
    std::cout << " Tstop = " << Tstop 
              << " for the mesh Nx = " << mesh1D.getxDiv() << std::endl;
    TOC;

    std::string filename = output_dir + "/Euler1D_Riemann_X" + std::to_string(xDiv_i) 
                                      + "_ox" + std::to_string(x_order) + 
                                      + "_ot" + std::to_string(t_order) + ".dat";
    std::ofstream outFile(filename);
    const Vector& xc = mesh1D.getCellCenter_vec();
    const Vector& rho = solver.getrho();
    const Vector& rhovx = solver.getrhovx();
    const Vector& E = solver.getE();
    Vector vx = rhovx.array() / rho.array();
    Vector pre = solver.getpre(solver.getu());
    outFile << "TITLE = \"Euler 1D Sod Riemann\"\n";
    outFile << "VARIABLES = \"X\", \"Rho\", \"Vx\", \"P\", \"E\"\n";
    outFile << "ZONE I = " << mesh1D.getxDiv() << " F = POINT\n";
    for (int j = 0; j < mesh1D.getxDiv(); ++j) {
      outFile << xc(j) << " " << rho(j) << " " << vx(j) << " " << pre(j) << " " << E(j) << "\n";
    }
  }
  return 0;
}

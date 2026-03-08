#include "QUEST.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>

namespace QUEST
{

struct DGRiemannInitialCondition {
  real_t x0;
  real_t rhoL, rhoR;
  real_t vxL, vxR;
  real_t pL, pR;
};

class Euler1D_DG_TVDRK_Riemann final : public Euler1D_DG_TVDRK {
public:
  Euler1D_DG_TVDRK_Riemann(const TensorMesh1D* mesh1D,
                          const fespace1D* fe,
                          const EX_TVDRK* rk_table,
                          const DGRiemannInitialCondition& ic)
    : Euler1D_DG_TVDRK(mesh1D, fe, rk_table), ic_(ic) {}

  real_t rho_init(const real_t& x) override { return (x < ic_.x0) ? ic_.rhoL : ic_.rhoR; }
  Matrix rho_init(const Matrix& x) override {
    Matrix out = Matrix::Zero(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i) {
      out(i) = rho_init(x(i));
    }
    return out;
  }
  real_t rho_bc(const real_t& x, const real_t& t) override { return rho_init(x); }

  real_t vx_init(const real_t& x) override { return (x < ic_.x0) ? ic_.vxL : ic_.vxR; }
  Matrix vx_init(const Matrix& x) override {
    Matrix out = Matrix::Zero(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i) {
      out(i) = vx_init(x(i));
    }
    return out;
  }
  real_t vx_bc(const real_t& x, const real_t& t) override { return vx_init(x); }

  real_t pre_init(const real_t& x) override { return (x < ic_.x0) ? ic_.pL : ic_.pR; }
  Matrix pre_init(const Matrix& x) override {
    Matrix out = Matrix::Zero(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i) {
      out(i) = pre_init(x(i));
    }
    return out;
  }
  real_t pre_bc(const real_t& x, const real_t& t) override { return pre_init(x); }

private:
  DGRiemannInitialCondition ic_;
};

} // namespace QUEST

/*
  ./bin/Euler1D_DG_TVDRK_Riemann -nx 400 -basis 1 -problem 0 -Tstop 0.2
*/
int main(int argc, char* argv[])
{
  std::string output_dir = "result/DG/Euler1D";
  std::filesystem::create_directories(output_dir);

  real_t x1 = -5.e0;
  real_t x2 = 5.e0;
  int xDiv = 400;
  int basis_order = 1;
  int qua_order = 4;
  int t_order = 3;
  real_t gamma = 1.4e0;
  real_t cfl = 0.15e0;
  real_t Tstop = 0.2e0;
  int outputgap = 20;
  int problem = 0; // 0: Sod, 1: Lax

  QUEST::OptionsParser args(argc, argv);
  args.AddOption(&x1, "-x1", "--x1", "Left boundary of x.");
  args.AddOption(&x2, "-x2", "--x2", "Right boundary of x.");
  args.AddOption(&xDiv, "-nx", "--xDiv", "Number of cells.");
  args.AddOption(&basis_order, "-basis", "--basis_order", "DG basis order.");
  args.AddOption(&qua_order, "-q", "--qua_order", "Quadrature order.");
  args.AddOption(&t_order, "-ot", "--t_order", "TVD-RK order in time (1~3).");
  args.AddOption(&gamma, "-gamma", "--gamma", "Gas constant gamma.");
  args.AddOption(&cfl, "-cfl", "--cfl", "CFL number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop", "Final time.");
  args.AddOption(&outputgap, "-output", "--outputgap", "Output every outputgap time steps if >0.");
  args.AddOption(&problem, "-problem", "--problem", "0: Sod, 1: Lax.");
  args.ParseCheck(std::cout);

  QUEST::TensorMesh1D mesh1D(x1, x2, xDiv);
  mesh1D.init();

  QUEST::BasisFunction1D basis(static_cast<QUEST::BasisType>(basis_order));
  QUEST::fespace1D fe(&mesh1D, &basis, qua_order, QUEST::QuadratureType::GaussQuadrature, 3);
  fe.init();

  QUEST::EX_TVDRK rk_table(t_order);
  rk_table.init();

  QUEST::DGRiemannInitialCondition ic;
  switch(problem) {
    case 0: ic = {0.e0, 1.0, 0.125, 0.0, 0.0, 1.0, 0.1}; break;
    case 1: ic = {0.e0, 0.445, 0.5, 0.698 / 0.445, 0.0, 3.528, 0.571}; break;
    default: QUEST_ERROR("Unsupported Riemann problem type.");
  }

  QUEST::Euler1D_DG_TVDRK_Riemann solver(&mesh1D, &fe, &rk_table, ic);
  solver.setgamma(gamma);
  solver.setcfl(cfl);
  solver.setlimiter(true, 1.5e0);
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
  Matrix rho = u_nodal[0];
  Matrix rhovx = u_nodal[1];
  Matrix E = u_nodal[2];
  Matrix vx = rhovx.array() / rho.array();
  Matrix p = (gamma - 1.e0) * (E.array() - 0.5e0 * rhovx.array().square() / rho.array());

  std::string filename = output_dir + "/Euler1D_DG_Riemann_N" + std::to_string(xDiv)
                                  + "_basis" + std::to_string(basis_order)
                                  + "_problem" + std::to_string(problem) + ".dat";
  std::ofstream outFile(filename);
  const Matrix& x = mesh1D.getCellCenter();
  outFile << "TITLE = \"Euler 1D DG Riemann\"\n";
  outFile << "VARIABLES = \"X\", \"Rho\", \"Vx\", \"P\", \"E\"\n";
  outFile << "ZONE I = " << xDiv << " F = POINT\n";
  for (int j = 0; j < xDiv; ++j) {
    outFile << x(0, j) << " " << rho(0, j) << " " << vx(0, j) << " " << p(0, j) << " " << E(0, j) << "\n";
  }

  return 0;
}

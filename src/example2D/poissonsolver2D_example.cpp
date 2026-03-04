#include "QUEST.hpp"

/*
  ./bin/poissonsolver2D_example -m 4 -nx 20 -ny 20 -basis 2 -plot 0 -poi 0
*/

int main(int argc, char *argv[])
{ 
  std::string title = "Accuracy test for Poisson equation using DG method !";
  std::string OutputDir = "result/DG/Poisson2D_DG_acctest";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t y1 = 0.e0;
  real_t y2 = 2.e0 * std::numbers::pi;
  int xDiv = 100;
  int yDiv = 100;
  int basis_id = 2;
  int poisoltype_id = 3;
  int qua_order = 5;     // 这个是积分点个数
  int quatype_id = 1;    // 积分类型
  real_t poitol = 1.e-9; // 泊松方程误差
  int m = 3;
  real_t C11 = 1.e0;
  Eigen::Vector2d C12;
  C12(0) = 0.5e0;
  C12(1) = 0.5e0;
  int plot = 0;
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&y1, "-y1", "--y1",
                  "The test y1 value.");
  args.AddOption(&y2, "-y2", "--y2",
                  "The test y2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&yDiv, "-ny", "--yDiv",
                  "The test first yDiv value.");
  args.AddOption(&basis_id, "-basis", "--basis_id",
                  "The order in the space direction.");
  args.AddOption(&qua_order, "-qua_order", "--qua_order",
                  "The quadrature order in the space direction.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&C11, "-c11", "--c11",
                  "The C11 parameter in the LDG method for Poisson.");
  args.AddOption(&(C12(0)), "-c12_0", "--c12_0",
                  "The C12 parameter in the LDG method for Poisson."); 
  args.AddOption(&(C12(1)), "-c12_1", "--c12_1",
                  "The C12 parameter in the LDG method for Poisson."); 
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&plot, "-plot", "--if_plot",
                  "The plot bool variable.");
  args.ParseCheck(std::cout);
  TOC;

  IntVector xDiv_vec(m);
  IntVector yDiv_vec(m);
  Vector uL1(m), uL2(m), uLinf(m);
  Vector qxL1(m), qxL2(m), qxLinf(m);
  Vector qyL1(m), qyL2(m), qyLinf(m);

  for (int i = 0; i < m; i++)
  {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    yDiv_vec(i) = std::pow(2, i) * yDiv;
    QUEST::PoissonSolver2DType poisoltype = 
          static_cast<QUEST::PoissonSolver2DType>(poisoltype_id);
    QUEST::BasisType basistype = 
          static_cast<QUEST::BasisType>(basis_id);
    QUEST::QuadratureType quatype = 
          static_cast<QUEST::QuadratureType>(quatype_id);

    QUEST::TensorMesh2D mesh2D(x1, x2, y1, y2, xDiv_vec(i), yDiv_vec(i));
    mesh2D.init();
    QUEST::BasisFunction2D basis(basistype);
    QUEST::fespace2D fe(&mesh2D, &basis, qua_order, quatype);
    fe.init();
    QUEST::PoissonSolver2DParameter poi_param;
    poi_param.C11 = C11;
    poi_param.C12 = C12;
    QUEST::Poisson_acctest_2D poi(&fe, poi_param);
    poi.init(poisoltype, poitol);

    // std::cout << " hx = " << mesh2D.gethx() << std::endl;
    // std::cout << " hy = " << mesh2D.gethy() << std::endl;
    // PAUSE();

    std::string filename_u = "X" + std::to_string(mesh2D.getxDiv()) + 
                              "_Y" + std::to_string(mesh2D.getyDiv()) +
                              "_o" + std::to_string(basis.getk1D()) +
                              "_u" + ".dat";
    std::string filename_qx = "X" + std::to_string(mesh2D.getxDiv()) +
                              "_Y" + std::to_string(mesh2D.getyDiv()) +
                              "_o" + std::to_string(basis.getk1D()) +
                              "_qx" + ".dat";
    std::string filename_qy = "X" + std::to_string(mesh2D.getxDiv()) +
                              "_Y" + std::to_string(mesh2D.getyDiv()) +
                              "_o" + std::to_string(basis.getk1D()) +
                              "_qy" + ".dat";
    std::string file_u = OutputDir + "/" + filename_u;
    std::string file_qx = OutputDir + "/" + filename_qx;
    std::string file_qy = OutputDir + "/" + filename_qy;
    std::ofstream outFile_u, outFile_qx, outFile_qy;
    if (plot > 0)
    {
      outFile_u.open(file_u);
      outFile_qx.open(file_qx);
      outFile_qy.open(file_qy);
    };
    
    Matrix Dirichlet_u;
    fe.Interpolate_Initial_Bdr(
      [&](const real_t& x, const real_t& y) { return poi.u_real(x, y);}, &Dirichlet_u);
    // const std::vector<Matrix>& CoorBdrRef = fe.getCoorBdrRef();
    // const std::vector<Eigen::Vector2d> extboundarycenter = 
    //   mesh2D.getextboundarycenter();
    // const int& extboundaryNum = mesh2D.getextboundaryNum();
    // const std::vector<int>& ExtBTypeIndex = mesh2D.getextboundarytypeindex();
    // const std::vector<int>& extNei = mesh2D.getextboundaryneighbors();
    // const Matrix& Cellcenter = mesh2D.getCellCenter();
    // const real_t& hx = mesh2D.gethx();
    // const real_t& hy = mesh2D.gethy();
    // Dirichlet_u.resize(qua_order, extboundaryNum);
    // Dirichlet_u.setZero();
    // Matrix temp;
    // real_t x0, y0;
    // int cell_index;
    // for (int i = 0; i < extboundaryNum; i++)
    // {
    //   temp = CoorBdrRef[ExtBTypeIndex[i]];
    //   cell_index = extNei[i];
    //   for (int q = 0; q < qua_order; q++)
    //   {
    //     x0 = Cellcenter(0, cell_index) + temp(0, q) * hx;
    //     y0 = Cellcenter(1, cell_index) + temp(1, q) * hy;
    //     std::cout << " x0 = " << x0 << " y0 = " << y0 << std::endl;
    //     Dirichlet_u(q, i) = poi.u_real(x0, y0);
    //   }
    // }
    // PAUSE();

    Matrix RHS_nodal;
    fe.Interpolate_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.RHS(x, y);}, &RHS_nodal);

    Vector bg, bf;
    poi.Assemble_bg(Dirichlet_u, &bg);
    poi.Assemble_bf(Dirichlet_u, RHS_nodal, &bf);

    Matrix u_real_modal, u_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.u_real(x, y);}, &u_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.u_real(x, y);}, &u_real_nodal);

    Matrix qx_real_modal, qx_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.u_real_dx(x, y);}, &qx_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.u_real_dx(x, y);}, &qx_real_nodal);

    Matrix qy_real_modal, qy_real_nodal;
    fe.Project_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.u_real_dy(x, y);}, &qy_real_modal);
    fe.Interpolate_Initial(
      [&](const Matrix& x, const Matrix& y) { return poi.u_real_dy(x, y);}, &qy_real_nodal);

    std::vector<Matrix> u_numerical_modal, u_numerical_nodal;
    poi.solveall(bg, bf, &u_numerical_modal);
    fe.modal_to_nodal2D(u_numerical_modal, &u_numerical_nodal);

    // 计算误差 //
    fe.computerrorL1(u_numerical_nodal[0], u_real_nodal, &(uL1(i)));
    fe.computerrorL2(u_numerical_nodal[0], u_real_nodal, &(uL2(i)));
    fe.computerrorLinf(u_numerical_nodal[0], u_real_nodal, &(uLinf(i)));

    fe.computerrorL1(u_numerical_nodal[1], qx_real_nodal, &(qxL1(i)));
    fe.computerrorL2(u_numerical_nodal[1], qx_real_nodal, &(qxL2(i)));
    fe.computerrorLinf(u_numerical_nodal[1], qx_real_nodal, &(qxLinf(i)));

    fe.computerrorL1(u_numerical_nodal[2], qy_real_nodal, &(qyL1(i)));
    fe.computerrorL2(u_numerical_nodal[2], qy_real_nodal, &(qyL2(i)));
    fe.computerrorLinf(u_numerical_nodal[2], qy_real_nodal, &(qyLinf(i)));

    // 试一下投影误差
    // Matrix u_project_modal, u_project_nodal;
    // fe.Project_Initial(
    //   [&](const Matrix& x, const Matrix& y) { return poi.u_real(x, y);}, &u_project_modal);
    // fe.modal_to_nodal2D(u_project_modal, &u_project_nodal);

    // Matrix qx_project_modal, qx_project_nodal;
    // fe.Project_Initial(
    //   [&](const Matrix& x, const Matrix& y) { return poi.u_real_dx(x, y);}, &qx_project_modal);
    // fe.modal_to_nodal2D(qx_project_modal, &qx_project_nodal);

    // Matrix qy_project_modal, qy_project_nodal;
    // fe.Project_Initial(
    //   [&](const Matrix& x, const Matrix& y) { return poi.u_real_dy(x, y);}, &qy_project_modal);
    // fe.modal_to_nodal2D(qy_project_modal, &qy_project_nodal);
    
    // fe.computerrorL1(u_project_nodal, u_real_nodal, &(uL1(i)));
    // fe.computerrorL2(u_project_nodal, u_real_nodal, &(uL2(i)));
    // fe.computerrorLinf(u_project_nodal, u_real_nodal, &(uLinf(i)));

    // fe.computerrorL1(qx_project_nodal, qx_real_nodal, &(qxL1(i)));
    // fe.computerrorL2(qx_project_nodal, qx_real_nodal, &(qxL2(i)));
    // fe.computerrorLinf(qx_project_nodal, qx_real_nodal, &(qxLinf(i)));

    // fe.computerrorL1(qy_project_nodal, qy_real_nodal, &(qyL1(i)));
    // fe.computerrorL2(qy_project_nodal, qy_real_nodal, &(qyL2(i)));
    // fe.computerrorLinf(qy_project_nodal, qy_real_nodal, &(qyLinf(i)));

    if (plot > 0) 
    {
      Vector u_numerical_aver = u_numerical_modal[0].row(0);
      Vector u_exact_aver = u_real_modal.row(0);
      Vector qx_numerical_aver = u_numerical_modal[1].row(0);
      Vector qx_exact_aver = qx_real_modal.row(0);
      Vector qy_numerical_aver = u_numerical_modal[2].row(0);
      Vector qy_exact_aver = qy_real_modal.row(0);
      fe.DisplayResult(u_numerical_aver, 
                          u_exact_aver, title, outFile_u); outFile_u.close();
      fe.DisplayResult(qx_numerical_aver, 
                          qx_exact_aver, title, outFile_qx); outFile_qx.close();
      fe.DisplayResult(qy_numerical_aver, 
                          qy_exact_aver, title, outFile_qy); outFile_qy.close();
    };
  }
  QUEST::DisplayAccTable2D(xDiv_vec, yDiv_vec, uL1, uL2, uLinf, std::cout);
  QUEST::DisplayAccTable2D(xDiv_vec, yDiv_vec, qxL1, qxL2, qxLinf, std::cout);
  QUEST::DisplayAccTable2D(xDiv_vec, yDiv_vec, qyL1, qyL2, qyLinf, std::cout);
  return 0;
}

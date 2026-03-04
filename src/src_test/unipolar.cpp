#include "srctest_1.hpp"

namespace QUEST
{

void KineticDriftDiffusion1D_SL_unipolar(int& argc, char *argv[]) {
  std::string title = "Unipolar diode test for Drift-Diffusion equation !";
  std::string OutputDir = "result/Unipolartest";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 2;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-6;
  real_t spartol = 1.e-6; // 隐式格式每一步误差
  real_t ddtol = 1.e-6; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-3;
  real_t Tstop = 0.05e0;
  real_t theta = 1.e0;
  real_t gamma = 0.002;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;

  OutputDir = OutputDir + "_gamma" +  std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  // Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  // Vector EL1(m), EL2(m), ELinf(m);
  // Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_E" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) +
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieDD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      kieDD.updateAll(Trun, dt);
      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      // if (kieDD.getiter() > 1000) {
      //   break;
      // }
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    const Vector& rho = kieDD.getrho();
    Vector rho_init = kieDD.getrhoinit();
    const Vector& rhod = kieDD.getrhod();
    std::vector<Vector> f_init = kieDD.getfinit();
    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    kieDD.update_E(Tstop);
    // Vector E_exact = kieDD.getEfinal(Tstop);
    // mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    // mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    // mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    // mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    // mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    // mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult_ex(E, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_init, title, outFile_f); outFile_f.close();
    // if (kieDD.getiter() > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
    //   std::abort();
    // }
  }
}

void KineticDriftDiffusion1D_SL_PNjunction_example(int& argc, char *argv[]) {
  std::string title = "PN junction diode test for semiconductor Boltzmann equation !";
  std::string OutputDir = "result/DD_SL_PNjunction";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 200;
  int x_order = 1;
  int t_order = 1;
  int poisoltype_id = 2;
  int ddsoltype_id = 3;
  int outputgap = 10;
  real_t poitol = 1.e-10; // 泊松每一步误差(用不到)
  real_t sparsetol = 1.e-10; // 隐式格式每一步误差
  real_t ddtol = 1.e-10;
  real_t cfl = 1.e0;
  real_t knu = 1.e-3;
  real_t Tstop = 0.05e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int NTH = 10;
  int m = 1;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&ddsoltype_id, "-ld", "--ddsoltype_id",
                  "The implicit solver for the macro equation type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&sparsetol, "-sparsetol", "--ldsparse_tol",
                  "The tolerance value of the implicit solver type value.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value iteration solver in the macro equation.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&knu, "-gamma", "--gamma",
                  "The gamma Debye length in the drift diffusion equation.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "The cfl condition for the time step size.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);
  TOC;

  OutputDir = OutputDir + "_knu" +  std::to_string(knu);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  // Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  // Vector EL1(m), EL2(m), ELinf(m);
  // Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_PNjunction kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(sparsetol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_T" + std::to_string(Tstop) + 
                              "_rho" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) +
                              "_T" + std::to_string(Tstop) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_f(file_f);

    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieDD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    const std::vector<Vector>& f = kieDD.getf();
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      kieDD.updateAll(Trun, dt);

      Trun += dt;
      step_count += 1;
      if (outputgap > 0)
      {
        if (step_count % outputgap == 0) 
        {
          std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
        }
      };
    };

    const Vector& rho = kieDD.getrho();
    Vector rho_init = kieDD.getrhoinit();
    std::vector<Vector> f_init = kieDD.getfinit();

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_init, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(f, f_init, title, outFile_f); outFile_f.close();
    // if (kieDD.getiter() > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
    //   std::abort();
    // }
  }
}

void KineticDriftDiffusion1D_SL_unipolar_Newton(int& argc, char *argv[]) {
  std::string title = "Unipolar diode test for Drift-Diffusion equation with Newton iteration !";
  std::string OutputDir = "result/Unipolartest_Newton";
  std::filesystem::create_directories(OutputDir);
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 2;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 0.002;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  // Vector EL1(m), EL2(m), ELinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_Newton kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_E" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << kieDD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      // if (kieDD.getiter() > 1000) {
      //   break;
      // }
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    const Vector& rho = kieDD.getrho();
    Vector rho_init = kieDD.getrhoinit();
    std::vector<Vector> f_init = kieDD.getfinit();
    // mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    // mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    // mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));
    kieDD.update_E(Tstop);
    // Vector E_exact = kieDD.getEfinal(Tstop);
    // mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    // mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    // mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    // mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    // mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    // mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_init, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult_ex(E, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_init, title, outFile_f); outFile_f.close();
    // if (kieDD.getiter() > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
    //   std::abort();
    // }
  }
}

void KienticDriftDiffusion1D_SL_unipolar_QUASINEUTRAL(int& argc, char *argv[]) {
  std::string title = "AP for Quasi-Neutral and Diffusion-limit Unipolar diode test for Drift-Diffusion equation !";
  std::string OutputDir = "result/Unipolartest_AP";
  
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 2;
  int quasi_poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-8;
  real_t spartol = 1.e-8; // 隐式格式每一步误差
  real_t ddtol = 1.e-8; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-3;
  real_t Tstop = 0.05e0;
  real_t theta = 1.e0;
  real_t gamma = 0.002;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&quasi_poisoltype_id, "-qpoi", "--quasi_poisoltype",
                  "The test quasi-neutral poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of update procedure for density.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the implicit iteration for drift-diffusion.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  
  OutputDir = OutputDir + "_gamma" +  std::to_string(gamma);
  std::filesystem::create_directories(OutputDir);
  
  IntVector xDiv_vec(m), vDiv_vec(m);
  // Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  // Vector EL1(m), EL2(m), ELinf(m);
  // Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD quasi_poisoltype = 
          static_cast<QUEST::Solver1DTypeFD>(quasi_poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::APDriftDiffusion_SL APkieDD(&mesh1D, &poisolver, x_order, t_order);
    APkieDD.settheta(theta);
    APkieDD.setcfl(cfl);
    APkieDD.setgamma(gamma);
    APkieDD.seteps(knu);
    APkieDD.setsparsetol(spartol);
    APkieDD.setiterationtol(ddtol);
    APkieDD.setsoltype_poi(quasi_poisoltype);
    APkieDD.setD(theta);
    APkieDD.setquasipoi_tol(poitol);
    APkieDD.setNTH(NTH);
    
    APkieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_E" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) +
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << APkieDD.getcfl() << std::endl;
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = APkieDD.getE();
    const std::vector<Vector>& f = APkieDD.getf();
    while (Trun < Tstop) {
      APkieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      APkieDD.updateAll(Trun, dt);
      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      // if (APkieDD.getiter() > 1000) {
      //   break;
      // }
      // std::vector<Vector> f_exact = APkieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    
    const Vector& rho = APkieDD.getrho();
    Vector rho_init = APkieDD.getrhoinit();
    const Vector& rhod = APkieDD.getrhod();
    std::vector<Vector> f_init = APkieDD.getfinit();
    PAUSE();
    // APkieDD.update_E_final(Tstop, 0.e0);
    
    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult_ex(E, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_init, title, outFile_f); outFile_f.close();
    // if (APkieDD.getiter() > 1000) {
    //   std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
    //   std::abort();
    // }
  }
}

void KienticDriftDiffusion1D_SL_acctest_QUASINEUTRAL(int& argc, char *argv[]) {
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DDAP_acctest";
  real_t x1 = 0.e0;
  real_t x2 = 2.e0 * std::numbers::pi;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式线性方程组误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  int NTH = 10;
  int m = 5;
  
  TIC;
  QUEST::OptionsParser args(argc, argv);

  args.AddOption(&x1, "-x1", "--x1",
                  "The test x1 value.");
  args.AddOption(&x2, "-x2", "--x2",
                  "The test x1 value.");
  args.AddOption(&v1, "-v1", "--v1",
                  "The test v1 value.");
  args.AddOption(&v2, "-v2", "--v2",
                  "The test v2 value.");
  args.AddOption(&xDiv, "-nx", "--xDiv",
                  "The test first xDiv value.");
  args.AddOption(&vDiv, "-nv", "--vDiv",
                  "The test vDiv value.");
  args.AddOption(&x_order, "-ox", "--x_order",
                  "The order in the space direction.");
  args.AddOption(&t_order, "-ot", "--t_order",
                  "The order in the time direction.");
  args.AddOption(&is, "-is", "--is",
                  "The test is value.");
  args.AddOption(&poisoltype_id, "-poi", "--poisoltype",
                  "The test poisson solver type value.");
  args.AddOption(&ddsoltype_id, "-dd", "--ddsoltype",
                  "The test drift diffusion solver type value.");
  args.AddOption(&poitol, "-poitol", "--poitol",
                  "The tolerance value of Poisson solver type value.");
  args.AddOption(&spartol, "-spartol", "--spartol",
                  "The tolerance value of the linear system for the density.");
  args.AddOption(&ddtol, "-ddtol", "--ddtol",
                  "The tolerance value of the iteration between electric potential and density.");
  args.AddOption(&knu, "-knu", "--Knudsen",
                  "The Knudsen dimensionless number.");
  args.AddOption(&Tstop, "-Tstop", "--Tstop",
                  "The final time value.");
  args.AddOption(&cfl, "-cfl", "--cfl",
                  "Highly associated with the time step size.");
  args.AddOption(&gamma, "-gamma", "--gamma",
                  "The parameter gamma in the Poisson equation.");
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.ParseCheck(std::cout);
  TOC;
  
  OutputDir = OutputDir + "_gamma" +  std::to_string(gamma) + "_knu" + std::to_string(knu);
  std::filesystem::create_directories(OutputDir);

  IntVector xDiv_vec(m), vDiv_vec(m);
  Vector rhoL1(m), rhoL2(m), rhoLinf(m);
  Vector EL1(m), EL2(m), ELinf(m);
  Vector fL1(m), fL2(m), fLinf(m);

  for (int i = 0; i < m; i++) {
    xDiv_vec(i) = std::pow(2, i) * xDiv;
    vDiv_vec(i) = std::pow(2, i) * vDiv;
    QUEST::PoissonSolver1DTypeFD poisoltype = 
          static_cast<QUEST::PoissonSolver1DTypeFD>(poisoltype_id);
    QUEST::Solver1DTypeFD ddsoltype = 
          static_cast<QUEST::Solver1DTypeFD>(ddsoltype_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::APDriftDiffusion1D_SL_period APkieDD(&mesh1D, &poisolver, x_order, t_order);
    APkieDD.settheta(theta);
    APkieDD.setcfl(cfl);
    APkieDD.setgamma(gamma);
    APkieDD.seteps(knu);
    APkieDD.setsparsetol(spartol);
    APkieDD.setiterationtol(ddtol);
    APkieDD.setD(theta);
    APkieDD.setNTH(NTH);
    APkieDD.setquasipoi_tol(poitol);
    APkieDD.init(ddsoltype);

    std::string filename_rho = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) +
                              "_cfl" + std::to_string(cfl) + 
                              "_rho" + ".dat";
    std::string filename_E = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_E" + ".dat";
    std::string filename_f = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_cfl" + std::to_string(cfl) + 
                              "_f" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    // , outFile_E(file_E)
    std::cout << " hx = " << mesh1D.gethx() << std::endl;
    std::cout << " Nx = " << mesh1D.getNx() << std::endl;
    std::cout << " cfl = " << APkieDD.getcfl() << std::endl;
    APkieDD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    const Vector& E = APkieDD.getE();
    const std::vector<Vector>& f = APkieDD.getf();
    
    while (Trun < Tstop) {
      APkieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      APkieDD.updateAll(Trun, dt);
      std::cout << " Trun = " << Trun << ", dt = " << dt << std::endl;
      Trun += dt;
      if (APkieDD.getiter() > 1000) {
        break;
      }
      // std::vector<Vector> f_exact = APkieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };

    const Vector& rho = APkieDD.getrho();
    const Vector& rho_exact = APkieDD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    APkieDD.update_E_final(Tstop);
    Vector E_exact = APkieDD.getEfinal(Tstop);
    std::vector<Vector> f_exact = APkieDD.getffinal(Tstop);
    mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
    if (APkieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, EL1, EL2, ELinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
};

} // namespace QUEST

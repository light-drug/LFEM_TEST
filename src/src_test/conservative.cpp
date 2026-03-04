#include "srctest_1.hpp"

namespace QUEST
{

void KineticDriftDiffusion1D_SL_Conservative_Lomac_acctest(int& argc, char *argv[])
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_acctest_lomac";
  std::filesystem::create_directories(OutputDir);
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
  int time_stepping_id = 1;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  real_t sigmas = 1.e0;
  int NTH = 10;
  int m = 5;
  int outputgap = 10;
  
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
  args.AddOption(&time_stepping_id, "-time_step", "--time_stepping_method",
                  "The test time stepping method for micro equation. 0 for Exp and 1 for RK.");
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
  args.AddOption(&theta, "-theta", "--theta",
                  "The parameter theta.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The parameter sigmas in the collision operator.");                
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);
  TOC;
  
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
    QUEST::TimeSteppingType time_stepping = 
          static_cast<QUEST::TimeSteppingType>(time_stepping_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_lomac_period kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.setsigmas(sigmas);
    kieDD.setTimeSteppingType(time_stepping);
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
    std::string filename_totaldensity = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_totaldensity" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::string file_totaldensity = OutputDir + "/" + filename_totaldensity;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    std::ofstream outFile_totaldensity(file_totaldensity);
    // , outFile_E(file_E)
    kieDD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    const Vector& rho = kieDD.getrho();
    const real_t& hx = mesh1D.gethx();

    std::vector<real_t> totaldensity;
    std::vector<real_t> timevec;
    totaldensity.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    timevec.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    totaldensity.push_back(rho.sum() * hx);
    timevec.push_back(Trun);
    
    
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      Trun += dt;
      step_count += 1;
      totaldensity.push_back(rho.sum() * hx);
      timevec.push_back(Trun);
      if (kieDD.getiter() > 1000) {
        break;
      };
      if (outputgap > 0)
      {
        if (step_count % outputgap == 0) 
        {
          std::cout << " Trun = " << Trun << ", dt = " << dt 
                    << ", Total density = " << rho.sum() * hx
                    << std::endl;
        }
      };
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };
    
    const Vector& rho_exact = kieDD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    kieDD.update_E(Tstop);
    Vector E_exact = kieDD.getEfinal(Tstop);
    std::vector<Vector> f_exact = kieDD.getffinal(Tstop);
    mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
    if (kieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }

    for (int j = 0; j < totaldensity.size(); j++) 
    {
      outFile_totaldensity << timevec[j] << " " <<  totaldensity[j] << std::endl;
    }
    outFile_totaldensity.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, EL1, EL2, ELinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
};

void KineticDriftDiffusion1D_SL_Conservative_Lomac_different_gamma_acctest(int& argc, char *argv[])
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_acctest_lomac";
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
  int time_stepping_id = 1;
  real_t poitol = 1.e-9;
  real_t spartol = 1.e-7; // 隐式格式每一步误差
  real_t ddtol = 1.e-7; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  real_t sigmas = 1.e0;
  int NTH = 10;
  int m = 5;
  int outputgap = 10;
  
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
  args.AddOption(&time_stepping_id, "-time_step", "--time_stepping_method",
                  "The test time stepping method for micro equation. 0 for Exp and 1 for RK.");
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
  args.AddOption(&theta, "-theta", "--theta",
                  "The parameter theta.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The parameter sigmas in the collision operator.");                
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);
  TOC;

  OutputDir = OutputDir + "_gamma" +  std::to_string(gamma);
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
    QUEST::TimeSteppingType time_stepping = 
          static_cast<QUEST::TimeSteppingType>(time_stepping_id);

    QUEST::KineticFDmesh_period mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD_period poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_lomac_period_Different_Gamma kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.setsigmas(sigmas);
    kieDD.setTimeSteppingType(time_stepping);
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
    std::string filename_totaldensity = "X" + std::to_string(mesh1D.getNx()) +
                              "_V" + std::to_string(mesh1D.getNv()) + 
                              "_totaldensity" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    std::string file_totaldensity = OutputDir + "/" + filename_totaldensity;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    std::ofstream outFile_totaldensity(file_totaldensity);
    // , outFile_E(file_E)
    kieDD.setrhofinal(Tstop);
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    const Vector& rho = kieDD.getrho();
    const real_t& hx = mesh1D.gethx();

    std::vector<real_t> totaldensity;
    std::vector<real_t> timevec;
    totaldensity.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    timevec.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    totaldensity.push_back(rho.sum() * hx);
    timevec.push_back(Trun);
    
    
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      Trun += dt;
      step_count += 1;
      totaldensity.push_back(rho.sum() * hx);
      timevec.push_back(Trun);
      if (kieDD.getiter() > 1000) {
        break;
      };
      if (outputgap > 0)
      {
        if (step_count % outputgap == 0) 
        {
          std::cout << " Trun = " << Trun << ", dt = " << dt 
                    << ", Total density = " << rho.sum() * hx
                    << std::endl;
        }
      };
      // std::vector<Vector> f_exact = kieDD.getffinal(Trun);
      // std::string titletemp = "Accuracy test for Drift-Diffusion equation Trun = " + std::to_string(Trun);
      // mesh1D.DisplayResult(f, f_exact, titletemp, outFile_f); 
      // outFile_f.close();
      // PAUSE();
    };
    
    const Vector& rho_exact = kieDD.getrhofinal();

    mesh1D.computerrorL1(rho, rho_exact, &(rhoL1(i)));
    mesh1D.computerrorL2(rho, rho_exact, &(rhoL2(i)));
    mesh1D.computerrorLinf(rho, rho_exact, &(rhoLinf(i)));

    kieDD.update_E(Tstop);
    Vector E_exact = kieDD.getEfinal(Tstop);
    std::vector<Vector> f_exact = kieDD.getffinal(Tstop);
    mesh1D.computerrorL1(E, E_exact, &(EL1(i)));
    mesh1D.computerrorL2(E, E_exact, &(EL2(i)));
    mesh1D.computerrorLinf(E, E_exact, &(ELinf(i)));
    mesh1D.computerrorL1(f, f_exact, &(fL1(i)));
    mesh1D.computerrorL2(f, f_exact, &(fL2(i)));
    mesh1D.computerrorLinf(f, f_exact, &(fLinf(i)));

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rho_exact, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_exact, title, outFile_f); outFile_f.close();
    if (kieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }

    for (int j = 0; j < totaldensity.size(); j++) 
    {
      outFile_totaldensity 
        << std::fixed << std::setprecision(3) << timevec[j] << " "
        << std::scientific << std::setprecision(8) << totaldensity[j]
        << std::endl;
    }
    outFile_totaldensity.close();
  }
  QUEST::DisplayAccTable1D(xDiv_vec, rhoL1, rhoL2, rhoLinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, EL1, EL2, ELinf, std::cout);
  QUEST::DisplayAccTable1D(xDiv_vec, fL1, fL2, fLinf, std::cout);
};

void KineticDriftDiffusion1D_SL_Conservative_Lomac_PNjunction(int& argc, char *argv[])
{
  std::string title = "Accuracy test for Drift-Diffusion equation !";
  std::string OutputDir = "result/DD_lomac_PNjunction";
  real_t x1 = 0.e0;
  real_t x2 = 1.e0;
  real_t v1 = -10.e0;
  real_t v2 = 10.e0;
  int xDiv = 100;
  int vDiv = 20;
  int x_order = 1;
  int t_order = 1;
  int is = 4;
  int poisoltype_id = 3;
  int ddsoltype_id = 3;
  int time_stepping_id = 1;
  real_t poitol = 1.e-10;
  real_t spartol = 1.e-10; // 隐式格式每一步误差
  real_t ddtol = 1.e-10; // 迭代误差
  real_t cfl = 1.e0;
  real_t knu = 1.e-5;
  real_t Tstop = 0.e0;
  real_t theta = 1.e0;
  real_t gamma = 1.e0;
  real_t sigmas = 1.e0;
  int NTH = 10;
  int m = 5;
  int outputgap = 10;
  
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
  args.AddOption(&time_stepping_id, "-time_step", "--time_stepping_method",
                  "The test time stepping method for micro equation. 0 for Exp and 1 for RK.");
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
  args.AddOption(&theta, "-theta", "--theta",
                  "The parameter theta.");
  args.AddOption(&sigmas, "-sigmas", "--sigmas",
                  "The parameter sigmas in the collision operator.");                
  args.AddOption(&NTH, "-NTH", "--NTH",
                  "The number of threads by OPENMP.");
  args.AddOption(&m, "-m", "--max",
                  "Max computational mesh size.");
  args.AddOption(&outputgap, "-output", "--outputgap",
                  "The Trun output gap time steps.");
  args.ParseCheck(std::cout);
  TOC;

  OutputDir = OutputDir + "_gamma" +  std::to_string(gamma);
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
    QUEST::TimeSteppingType time_stepping = 
          static_cast<QUEST::TimeSteppingType>(time_stepping_id);

    QUEST::KineticFDmesh mesh1D(x1, x2, v1, v2, xDiv_vec(i), vDiv_vec(i));
    mesh1D.init();
    QUEST::Poissonsolver1DFD poisolver(&mesh1D, x_order);
    poisolver.init(poisoltype, poitol);
    QUEST::KineticDriftDiffusion1D_SL_lomac_PNjunction kieDD(&mesh1D, &poisolver, x_order, t_order);
    kieDD.settheta(theta);
    kieDD.setcfl(cfl);
    kieDD.setgamma(gamma);
    kieDD.seteps(knu);
    kieDD.setsparsetol(spartol);
    kieDD.setiterationtol(ddtol);
    kieDD.setD(theta);
    kieDD.setNTH(NTH);
    kieDD.setsigmas(sigmas);
    kieDD.setTimeSteppingType(time_stepping);
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
    // std::string filename_totaldensity = "X" + std::to_string(mesh1D.getNx()) +
    //                           "_V" + std::to_string(mesh1D.getNv()) + 
    //                           "_totaldensity" + ".dat";
    std::string file_rho = OutputDir + "/" + filename_rho;
    std::string file_E = OutputDir + "/" + filename_E;
    std::string file_f = OutputDir + "/" + filename_f;
    // std::string file_totaldensity = OutputDir + "/" + filename_totaldensity;
    std::ofstream outFile_rho(file_rho), outFile_E(file_E), outFile_f(file_f);
    // std::ofstream outFile_totaldensity(file_totaldensity);
    // , outFile_E(file_E)
    real_t Trun = 0.e0;
    real_t dt = 0.e0;
    int step_count = 0;
    const Vector& E = kieDD.getE();
    const std::vector<Vector>& f = kieDD.getf();
    const Vector& rho = kieDD.getrho();
    const real_t& hx = mesh1D.gethx();

    // std::vector<real_t> totaldensity;
    // std::vector<real_t> timevec;
    // totaldensity.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    // timevec.reserve(Tstop / (cfl * mesh1D.gethx()) + 1);
    // totaldensity.push_back(rho.sum() * hx);
    // timevec.push_back(Trun);
    
    
    while (Trun < Tstop) {
      kieDD.setdt(&dt);
      if (Trun + dt > Tstop) {
          dt = Tstop - Trun;
      };
      
      kieDD.updateAll(Trun, dt);
      Trun += dt;
      step_count += 1;
      // totaldensity.push_back(rho.sum() * hx);
      // timevec.push_back(Trun);
      if (kieDD.getiter() > 1000) {
        break;
      };
      if (outputgap > 0)
      {
        if (step_count % outputgap == 0) 
        {
          std::cout << " Trun = " << Trun << ", dt = " << dt 
                    // << ", Total density = " << rho.sum() * hx
                    << std::endl;
        }
      };
    };
    
    const Vector& rhod = kieDD.getrhod();
    const Vector& phi = kieDD.getphi();
    std::vector<Vector> f_init = kieDD.getfinit();

    title = title + " T = " + std::to_string(Trun);
    mesh1D.DisplayResult(rho, rhod, title, outFile_rho); outFile_rho.close();
    // mesh1D.DisplayResult(E, E_exact, title, outFile_E); outFile_E.close();
    mesh1D.DisplayResult(f, f_init, title, outFile_f); outFile_f.close();
    if (kieDD.getiter() > 1000) {
      std::cout << "[ERROR] Iteration exceeded 1000 steps. Abort. Write the wrong density value in the document " << std::endl;
      std::abort();
    }
  }
};


}


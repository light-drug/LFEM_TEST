#include "display.hpp"

namespace QUEST
{

void DisplayAccTable1D(const IntVector& xDiv_vec,
                       const Vector& errorL1_vec,
                       const Vector& errorL2_vec, 
                       const Vector& errorLinf_vec,
                       std::ostream& Outfile) {

  // 输出标题行
  Outfile << std::setw(8) << "Mesh"
          << std::setw(20) << "errorL1" << std::setw(10) << "order"
          << std::setw(20) << "errorL2" << std::setw(10) << "order"
          << std::setw(20) << "errorLinf" << std::setw(10) << "order"
          << "\n";

  // 输出数据
  for (int k = 0; k < xDiv_vec.size(); k++) {   
    if (k == 0) {
      Outfile << std::setw(9) << xDiv_vec[k]
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL1_vec[k]
              << std::setw(10) << "" 
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL2_vec[k]
              << std::setw(10) << ""
              << std::setw(20) << std::scientific << std::setprecision(2) << errorLinf_vec[k]
              << std::setw(10) << ""
              << "\n";
    } else {
      Outfile << std::setw(9) << xDiv_vec[k]
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL1_vec[k] 
              << std::setw(10) << std::fixed << std::setprecision(2) 
              << std::log2(errorL1_vec[k-1] / errorL1_vec[k]) / std::log2(2.0)
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL2_vec[k]
              << std::setw(10) << std::fixed << std::setprecision(2) 
              << std::log2(errorL2_vec[k-1] / errorL2_vec[k]) / std::log2(2.0)
              << std::setw(20) << std::scientific << std::setprecision(2) << errorLinf_vec[k]
              << std::setw(10) << std::fixed << std::setprecision(2) 
              << std::log2(errorLinf_vec[k-1] / errorLinf_vec[k]) / std::log2(2.0)
              << "\n";
    };
  };

  Outfile << "\\begin{table}[htbp]\n";
  Outfile << "\\begin{tabular}{c|c c|c c|c c}\n";
  Outfile << "\\hline\n";
  Outfile << "Mesh & $L^1$ error & order & $L^2$ error & order & $L^\\infty$ error & order \\\\ \n";
  Outfile << "\\hline\n";

  for (int k = 0; k < xDiv_vec.size(); ++k) {
    Outfile << xDiv_vec[k] << " & ";
    Outfile << std::scientific << std::setprecision(2) << errorL1_vec[k] << " & ";
    if (k == 0)
      Outfile << " -- ";
    else {
      Outfile << std::fixed << std::setprecision(2)
              << std::log2(errorL1_vec[k - 1] / errorL1_vec[k]);
    };

    Outfile << " & " << std::scientific << std::setprecision(2) << errorL2_vec[k] << " & ";
    if (k == 0)
      Outfile << " -- ";
    else {
      Outfile << std::fixed << std::setprecision(2)
              << std::log2(errorL2_vec[k - 1] / errorL2_vec[k]);
    };

    Outfile << " & " << std::scientific << std::setprecision(2) << errorLinf_vec[k] << " & ";
    if (k == 0)
      Outfile << " -- ";
    else {
      Outfile << std::fixed << std::setprecision(2)
              << std::log2(errorLinf_vec[k - 1] / errorLinf_vec[k]);
    };

    Outfile << " \\\\\n";
  };

  Outfile << "\\hline\n";
  Outfile << "\\end{tabular}\n";
  Outfile << "\\end{table}\n";
}

void DisplayAccTable2D(const IntVector& xDiv_vec,
                      const IntVector& yDiv_vec,
                      const Vector& errorL1_vec,
                      const Vector& errorL2_vec, 
                      const Vector& errorLinf_vec,
                      std::ostream& Outfile) {

  // 输出标题行
  Outfile << std::setw(8) << "Mesh"
          << std::setw(20) << "errorL1" << std::setw(10) << "order"
          << std::setw(20) << "errorL2" << std::setw(10) << "order"
          << std::setw(20) << "errorLinf" << std::setw(10) << "order"
          << "\n";

  // 输出数据
  for (int k = 0; k < xDiv_vec.size(); k++) {   
    if (k == 0) {
      Outfile << std::setw(9) << xDiv_vec[k] << "x" << yDiv_vec[k]
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL1_vec[k]
              << std::setw(10) << "" 
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL2_vec[k]
              << std::setw(10) << ""
              << std::setw(20) << std::scientific << std::setprecision(2) << errorLinf_vec[k]
              << std::setw(10) << ""
              << "\n";
    } else {
      Outfile << std::setw(9) << xDiv_vec[k] << "x" << yDiv_vec[k]
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL1_vec[k] 
              << std::setw(10) << std::fixed << std::setprecision(2) 
              << std::log2(errorL1_vec[k-1] / errorL1_vec[k]) / std::log2(2.0)
              << std::setw(20) << std::scientific << std::setprecision(2) << errorL2_vec[k]
              << std::setw(10) << std::fixed << std::setprecision(2) 
              << std::log2(errorL2_vec[k-1] / errorL2_vec[k]) / std::log2(2.0)
              << std::setw(20) << std::scientific << std::setprecision(2) << errorLinf_vec[k]
              << std::setw(10) << std::fixed << std::setprecision(2) 
              << std::log2(errorLinf_vec[k-1] / errorLinf_vec[k]) / std::log2(2.0)
              << "\n";
    };
  };

  Outfile << "\\begin{table}[htbp]\n";
  Outfile << "\\begin{tabular}{c|c c|c c|c c}\n";
  Outfile << "\\hline\n";
  Outfile << "Mesh & $L^1$ error & order & $L^2$ error & order & $L^\\infty$ error & order \\\\ \n";
  Outfile << "\\hline\n";

  for (int k = 0; k < xDiv_vec.size(); ++k) {
    Outfile << xDiv_vec[k] << "\\times" << yDiv_vec[k] << " & ";
    Outfile << std::scientific << std::setprecision(2) << errorL1_vec[k] << " & ";
    if (k == 0)
      Outfile << " -- ";
    else {
      Outfile << std::fixed << std::setprecision(2)
              << std::log2(errorL1_vec[k - 1] / errorL1_vec[k]);
    };

    Outfile << " & " << std::scientific << std::setprecision(2) << errorL2_vec[k] << " & ";
    if (k == 0)
      Outfile << " -- ";
    else {
      Outfile << std::fixed << std::setprecision(2)
              << std::log2(errorL2_vec[k - 1] / errorL2_vec[k]);
    };

    Outfile << " & " << std::scientific << std::setprecision(2) << errorLinf_vec[k] << " & ";
    if (k == 0)
      Outfile << " -- ";
    else {
      Outfile << std::fixed << std::setprecision(2)
              << std::log2(errorLinf_vec[k - 1] / errorLinf_vec[k]);
    };

    Outfile << " \\\\\n";
  };

  Outfile << "\\hline\n";
  Outfile << "\\end{tabular}\n";
  Outfile << "\\end{table}\n";
}

} // namespace QUEST

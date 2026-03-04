#ifndef QUEST_SRC_TEST_HPP
#define QUEST_SRC_TEST_HPP

#include "basis.hpp"
#include "timer.hpp"
#include "integralrule.hpp"
#include "Tensormesh.hpp"
#include "fespace.hpp"
#include "QUEST.hpp"

#include <filesystem>

namespace QUEST {

void introduction();

void basis_test();

void fespace_test();

void DriftDiffusion1D_FD_unipolar(int& argc, char *argv[]);

// 一维PN单组份简化算例
void DriftDiffusion1D_FD_PNjunction_example(int& argc, char *argv[]);

// 在迭代处加罚保证极大值原理
void DriftDiffusion1D_FD_unipolar_penaltyIter(int& argc, char *argv[]);

// 在隐式部分加罚来保证与gamma无关的时间步长(存疑,这个算例有初始层)以及AP
void DriftDiffusion1D_FD_unipolar_penaltyImplicit(int& argc, char *argv[]);

// PN结算例，隐式加罚得到与gamma无关时间步长，这个PN结算例没有初始层
void DriftDiffusion1D_FD_PNjunction_penaltyImplicit_example(int& argc, char *argv[]);

// DD模型精度算例
void DriftDiffusion1D_FD_acctest(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_acctest(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_acctest_Different_Gamma(int& argc, char *argv[]);

// 尝试的一些守恒方法(差分框架下没啥用)
void KineticDriftDiffusion1D_SL_acctest_Different_Gamma_conservative(int& argc, char *argv[]);

void KineticLinearDiffusion1D_SL_acctest(int& argc, char *argv[]);

// 这个是复现张国梁那篇两速度模型的方法，其中微观方程更新用的半拉
void KineticLinearDiffusion1D_SL_acctest_twovel(int& argc, char *argv[]);

void KineticLinearDiffusion1D_SL_IsotropicBoundary(int& argc, char *argv[]);

// 这个是复现张国梁那篇两速度模型的方法，其中微观方程更新用的隐式
void KineticLinearDiffusion1D_SL_Implicit_acctest_twovel(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_unipolar(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_PNjunction_example(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_unipolar_Newton(int& argc, char *argv[]);

void KienticDriftDiffusion1D_SL_unipolar_QUASINEUTRAL(int& argc, char *argv[]);

void KienticDriftDiffusion1D_SL_acctest_QUASINEUTRAL(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_PNjunction_QUASINEUTRAL(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_Conservative_Lomac_acctest(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_Conservative_Lomac_different_gamma_acctest(int& argc, char *argv[]);

void KineticDriftDiffusion1D_SL_Conservative_Lomac_PNjunction(int& argc, char *argv[]);

// 以下是DG部分的solver
// DG-IMEX Macro-micro decomposition method for Kinetic Drift-Diffusion equation
void Poisson1D_DG_acctest(int& argc, char *argv[]);

void Poisson1D_DG_acctest_period(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_RiemannProblem(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_IsentropicBC(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_acctest(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_Schur_acctest(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_Schur_RiemannProblem(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_Schur_IsentropicBC(int& argc, char *argv[]);

// DG-IMEN solver for the kinetic drift-diffusion equation
void KineticDriftDiffusion1D_DGIMEX_Schur_acctest(int& argc, char *argv[]);

void KineticDriftDiffusion1D_DGIMEX_Schur_Unipolar(int& argc, char *argv[]);

void KineticLinearDiffusion1D_DGIMEX_Schur_withMaxwell_acctest(int& argc, char *argv[]);

} // namespace QUEST 

#endif // QUEST_SRC_TEST_HPP

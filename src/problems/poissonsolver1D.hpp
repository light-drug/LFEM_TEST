#ifndef QUEST_POISSONSOLVER1D_HPP
#define QUEST_POISSONSOLVER1D_HPP

#include <iostream>

#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "hyperbolicproblemsbase.hpp"
#include "config.hpp"
#include "error.hpp"
#include "timer.hpp"

namespace QUEST
{

enum class PoissonSolver1DType {
  CG,
  PCG,
  LDLT,
  LU,
  GMRES,
  BICG,
};

struct PoissonSolver1DParameter {
  real_t C11;
  real_t C12;
};

class PoissonSolver1D
{
protected:

  int num_equations_;
  int intboundaryNum_;
  int extboundaryNum_;
  int dim_;
  int polydim_;
  int numqua_;
  int ncell_;
  int hx_;
  const fespace1D* fe_;
  const PoissonSolver1DParameter& pa_;
  PoissonSolver1DType poitype_;
  SparseMatrix B_;
  SparseMatrix BT_;
  SparseMatrix Kinv_;
  SparseMatrix K_;
  SparseMatrix C_;
  SparseMatrix S_;

  virtual void generateB();
  virtual void generateK();
  virtual void generateKinv();
  virtual void generateC();
  virtual void generateS();
#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> chol_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> chol_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;

public:
  PoissonSolver1D(const fespace1D* fe, const PoissonSolver1DParameter& pa);

  virtual ~PoissonSolver1D() = default;
  virtual void init(const PoissonSolver1DType& poitype,
      const real_t& tol);
  virtual void solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const;
  virtual void solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Matrix>* All) const;
  virtual void solveall(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Matrix* E,
                Matrix* Phi) const;
  
  virtual void Assemble_bg(const Matrix& Dirichlet,
                Vector* bg) const;
  virtual void Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const;

};

class PoissonSolver1D_period
 : public PoissonSolver1D
{
protected:

#ifndef QUEST_USE_MKL
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif
  Eigen::GMRES<SparseMatrix> gmres_;

public:
  PoissonSolver1D_period(const fespace1D* fe, const PoissonSolver1DParameter& pa);
  ~PoissonSolver1D_period() override = default;
  void init(const PoissonSolver1DType& poitype, const real_t& tol) override;

  void generateB() override;
  void generateC() override;
  void generateS() override;
  virtual void generateBT();

  void solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const override;
  void solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Matrix>* All) const override;
  void solveall(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Matrix* E,
                Matrix* Phi) const override;

  void Assemble_bg(const Matrix& Dirichlet,
                  Vector* bg) const override;
  void Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const override;

};

class Poisson_acctest_1D : public PoissonSolver1D
{
public:
  Poisson_acctest_1D(const fespace1D* fe, const PoissonSolver1DParameter& pa);
  ~Poisson_acctest_1D() override = default;

  virtual Matrix RHS(const Matrix& x);
  virtual real_t u_bc(const real_t& x);
  virtual Matrix u_real(const Matrix& x);
  virtual real_t u_real(const real_t& x);
  virtual real_t u_real_dx(const real_t& x);
  virtual Matrix u_real_dx(const Matrix& x);

};

class Poisson_acctest_1D_period : public PoissonSolver1D_period
{
public:
  Poisson_acctest_1D_period(const fespace1D* fe, const PoissonSolver1DParameter& pa);
  ~Poisson_acctest_1D_period() override = default;

  virtual Matrix RHS(const Matrix& x);
  virtual real_t u_bc(const real_t& x);
  virtual Matrix u_real(const Matrix& x);
  virtual real_t u_real(const real_t& x);
  virtual real_t u_real_dx(const real_t& x);
  virtual Matrix u_real_dx(const Matrix& x);

};




} // namespace QUEST


#endif // QUEST_POISSONSOLVER1D_HPP

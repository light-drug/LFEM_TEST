#ifndef QUEST_POISSONSOLVER2D_HPP
#define QUEST_POISSONSOLVER2D_HPP

#include <iostream>

#include "fespace.hpp"
#include "Tensormesh.hpp"
#include "hyperbolicproblemsbase.hpp"
#include "config.hpp"
#include "error.hpp"
#include "timer.hpp"

namespace QUEST
{

enum class PoissonSolver2DType {
  CG,
  PCG,
  LDLT,
  LU,
  GMRES,
  BICG,
};

struct PoissonSolver2DParameter {
  real_t C11;
  Eigen::Vector2d C12;
};

class PoissonSolver2D
{
protected:

  int num_equations_;
  int intboundaryNum_;
  int extboundaryNum_;
  int dim_;
  int polydim_;
  int numqua_;
  int ncell_;
  const fespace2D* fe_;
  const PoissonSolver2DParameter& pa_;
  PoissonSolver2DType poitype_;
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
  PoissonSolver2D(const fespace2D* fe, const PoissonSolver2DParameter& pa);

  virtual ~PoissonSolver2D() = default;
  virtual void init(const PoissonSolver2DType& poitype,
      const real_t& tol);
  virtual void solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Vector>* All) const;
  virtual void solveall(const Vector& bg,
                const Vector& bf,
                std::vector<Matrix>* All) const;
  virtual void solveall(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                std::vector<Matrix>* E,
                Matrix* Phi) const;
  
  virtual void Assemble_bg(const Matrix& Dirichlet,
                Vector* bg) const;
  virtual void Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const;

};

class PoissonSolver2D_period
 : public PoissonSolver2D
{
protected:

#ifndef QUEST_USE_MKL
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif
  Eigen::GMRES<SparseMatrix> gmres_;

public:
  PoissonSolver2D_period(const fespace2D* fe, const PoissonSolver2DParameter& pa);
  ~PoissonSolver2D_period() override = default;
  void init(const PoissonSolver2DType& poitype, const real_t& tol) override;

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
                std::vector<Matrix>* E,
                Matrix* Phi) const override;

  void Assemble_bg(const Matrix& Dirichlet,
                  Vector* bg) const override;
  void Assemble_bf(const Matrix& Dirichlet,
                const Matrix& f_nodal,
                Vector* bf) const override;

};

class Poisson_acctest_2D : public PoissonSolver2D
{
public:
  Poisson_acctest_2D(const fespace2D* fe, const PoissonSolver2DParameter& pa);
  ~Poisson_acctest_2D() override = default;

  virtual Matrix RHS(const Matrix& x, const Matrix& y);
  virtual real_t u_bc(const real_t& x, const real_t& y);
  virtual Matrix u_real(const Matrix& x, const Matrix& y);
  virtual real_t u_real(const real_t& x, const real_t& y);
  virtual real_t u_real_dx(const real_t& x, const real_t& y);
  virtual Matrix u_real_dx(const Matrix& x, const Matrix& y);
  virtual real_t u_real_dy(const real_t& x, const real_t& y);
  virtual Matrix u_real_dy(const Matrix& x, const Matrix& y);

};

class Poisson_acctest_2D_period : public PoissonSolver2D_period
{
public:
  Poisson_acctest_2D_period(const fespace2D* fe, const PoissonSolver2DParameter& pa);
  ~Poisson_acctest_2D_period() override = default;

  virtual Matrix RHS(const Matrix& x, const Matrix& y);
  virtual real_t u_bc(const real_t& x, const real_t& y);
  virtual Matrix u_real(const Matrix& x, const Matrix& y);
  virtual real_t u_real(const real_t& x, const real_t& y);
  virtual real_t u_real_dx(const real_t& x, const real_t& y);
  virtual Matrix u_real_dx(const Matrix& x, const Matrix& y);
  virtual real_t u_real_dy(const real_t& x, const real_t& y);
  virtual Matrix u_real_dy(const Matrix& x, const Matrix& y);

};




} // namespace QUEST


#endif // QUEST_POISSONSOLVER2D_HPP

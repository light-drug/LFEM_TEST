#ifndef QUEST_POISSONSOLVER1DFD_HPP
#define QUEST_POISSONSOLVER1DFD_HPP

#include "config.hpp"
#include "FDmesh.hpp"
#include "error.hpp"
#include <vector>
#include <cmath>
#include "timer.hpp"

namespace QUEST
{

enum class PoissonSolver1DTypeFD {
  CG,
  PCG,
  LDLT,
  LU,
  GMRES,
  BICG,
};

class Poissonsolver1DFD_Base
{

protected:
  const FDmesh* mesh1D_;
  const int& x_order_;

  PoissonSolver1DTypeFD poitype_;
  SparseMatrix S_;
  // SparseMatrix Sinv_;

private:
#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> chol_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> chol_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;

public:

  Poissonsolver1DFD_Base(const FDmesh* mesh1D,
                    const int& x_order);

  virtual ~Poissonsolver1DFD_Base() = default;

  virtual void generateS();
  virtual const SparseMatrix& getS() const;
  // virtual void initSinv();
  // virtual const SparseMatrix getSinv() const;
  virtual void init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol);
  virtual void Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf);
  virtual void Solve(const Vector& bf,
            Vector* solution);
  virtual void SolveAll(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* solution,
                Vector* solution_dx);
};

class Poissonsolver1DFD_Base_period : public Poissonsolver1DFD_Base
{
private:

#ifndef QUEST_USE_MKL
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif 
  Eigen::BiCGSTAB<SparseMatrix> bicg_;


public:

  Poissonsolver1DFD_Base_period(const FDmesh_period* mesh1D, 
                    const int& x_order);
  ~Poissonsolver1DFD_Base_period() override = default;

  void generateS() override;
  void init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol) override;
  void Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf) override;
  void Solve(const Vector& bf,
            Vector* solution) override;
  void SolveAll(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* solution,
                Vector* solution_dx) override;
                    
};



class Poissonsolver1DFD
{

protected:
  const KineticFDmesh* mesh1D_;
  const int& x_order_;

  PoissonSolver1DTypeFD poitype_;
  SparseMatrix S_;
  // SparseMatrix Sinv_;

private:
#ifndef QUEST_USE_MKL
  Eigen::SimplicialLDLT<SparseMatrix> chol_;
#else 
  Eigen::PardisoLDLT<SparseMatrix> chol_;
#endif 
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower> cg_;

public:

  Poissonsolver1DFD(const KineticFDmesh* mesh1D,
                    const int& x_order);

  virtual ~Poissonsolver1DFD() = default;

  virtual void generateS();
  virtual const SparseMatrix& getS() const;
  // virtual void initSinv();
  // virtual const SparseMatrix getSinv() const;
  virtual void init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol);
  virtual void Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf);
  virtual void Solve(const Vector& bf,
            Vector* solution);
  virtual void SolveAll(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* solution,
                Vector* solution_dx);
};


class Poissonsolver1DFD_period : public Poissonsolver1DFD
{

private:

#ifndef QUEST_USE_MKL
  Eigen::SparseLU<SparseMatrix> lu_;
#else 
  Eigen::PardisoLU<SparseMatrix> lu_;
#endif 
  Eigen::BiCGSTAB<SparseMatrix> bicg_;


public:

  Poissonsolver1DFD_period(const KineticFDmesh_period* mesh1D, 
                    const int& x_order);
  ~Poissonsolver1DFD_period() override = default;

  void generateS() override;
  void init(const PoissonSolver1DTypeFD& poitype,
            const real_t& tol) override;
  void Assemble_bf(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* bf) override;
  void Solve(const Vector& bf,
            Vector* solution) override;
  void SolveAll(const real_t& bc_L,
                const real_t& bc_R,
                const Vector& f,
                Vector* solution,
                Vector* solution_dx) override;
                    
};

std::string name(const QUEST::PoissonSolver1DTypeFD& poitype);

} // namespace QUEST

#endif // QUEST_POISSONSOLVER1DFD_HPP

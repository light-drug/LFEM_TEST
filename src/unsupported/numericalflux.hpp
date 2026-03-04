#ifndef QUEST_NUMERICALFLUX_HPP
#define QUEST_NUMERICALFLUX_HPP

#include "config.hpp"
#include "error.hpp"

namespace QUEST
{
  
class FluxFunction
{
protected:

  const int num_equations_;

  const int dim_;

public:

  FluxFunction(const int num_equations, const int dim): 
      num_equations_(num_equations), dim_(dim) { };

  // 这个dimid从1开始 //
  virtual void ComputeFx(const Vector& u, const int& dimid, Vector& Fu) const = 0;

  virtual void ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju) const = 0;

  virtual void ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju, 
        Matrix& Lu, DiagonalMatrix& lambdau, Matrix& Ru) const = 0;
      
  virtual void ComputeMaxSpeed(const Vector& u, const int& dimid, real_t& alpha) const = 0;
  
  ~FluxFunction() = default;

};

// U_t + a(0)U_x + a(1)U_y + a(2)U_z = 0
class TransportFlux: public FluxFunction
{
private:

  const Vector a_;

public:

  TransportFlux(const int num_equations, const int dim, const Vector a);

  void ComputeFx(const Vector& u, const int& dimid, Vector& Fu) const override;

  void ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju) const override;

  void ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju, 
        Matrix& Lu, DiagonalMatrix& lambda, Matrix& Ru) const override;
      
  void ComputeMaxSpeed(const Vector& u, const int& dimid, real_t& alpha) const override;
  
  ~TransportFlux() = default;

};

class EulerFlux: public FluxFunction
{
private:

  const real_t gasgamma_;

public:

  EulerFlux(const int dim, const real_t gasgamma);

  void ComputeFx(const Vector& u, const int& dimid, Vector& Fu) const override;

  void ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju) const override;

  void ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju, 
        Matrix& Lu, DiagonalMatrix& lambda, Matrix& Ru) const override;
      
  void ComputeMaxSpeed(const Vector& u, const int& dimid, real_t& alpha) const override;
  
  ~EulerFlux() = default;

};

class EulerianNumericalFlux
{
protected:

  const int num_equations_;

  const int dim_;

public:

  NumericalFlux(const int num_equations, const int dim);

  ~NumericalFlux();

};

class UpwindFlux: public EulerianNumericalFlux
{
private:
  /* data */
public:
  UpwindFlux(/* args */);
  ~UpwindFlux();
};

class DownwindFlux: public EulerianNumericalFlux
{
private:
  /* data */
public:
  DownwindFlux(/* args */);
  ~DownwindFlux();
};

class AverageFlux: public NumericalFlux
{
private:
  /* data */
public:
  AverageFlux(/* args */);
  ~AverageFlux();
};

class LaxFriedchesFlux: public NumericalFlux
{
private:
  /* data */
public:
  LaxFriedchesFlux(/* args */);
  ~LaxFriedchesFlux();
};

} // namespace QUEST





#endif  // QUEST_NUMERICALFLUX_HPP 

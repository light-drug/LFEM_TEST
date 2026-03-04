#include "numericalflux.hpp"

namespace QUEST
{

TransportFlux::TransportFlux(const int num_equations, const int dim, const Vector a)
  : FluxFunction(num_equations, dim), a_(a) {
    if (num_equations != 1)
    {
      QUEST_ABORT("Linear Flux need num_equation equal to one ! ")
    }
  };

void TransportFlux::ComputeFx(const Vector& u, const int& dimid, Vector& Fu) const {
  Fu = a_(dimid) * u;
};

void TransportFlux::ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju) const {
  Ju.setZero();
  Ju.resize(num_equations_, num_equations_);
  Ju(0,0) = a_(dimid);
};

void TransportFlux::ComputeJacobian(const Vector& u, const int& dimid, Matrix& Ju, 
        Matrix& Lu, DiagonalMatrix& lambda, Matrix& Ru) const {
  Lu.resize(num_equations_, num_equations_);
  Lu(0,0) = 1;
  Ru.resize(num_equations_, num_equations_);
  Ru(0,0) = 1;
  lambda.resize(num_equations_);
  lambda.diagonal()[0] = a_(dimid);
};

EulerFlux::EulerFlux(const int dim, const real_t gasgamma) 
    : FluxFunction(dim + 2, dim), gasgamma_(gasgamma) { };


void EulerFlux::ComputeFx(const Vector& u, const int& dimid, Vector& Fu) const {
  real_t density = u(0);
  Vector momentum(dim_);
  real_t energy = u(dim_ + 1);
  real_t pre;
  for (int d = 0; d < dim_; d++)
  {
    momentum(d) = u(d + 1);
  }
  real_t kinetic_energy = 0.5e0 * (momentum.dot(momentum)) / rho;
  pre = (gasgamma_ - 1) * (energy - kinetic_energy);

  Fu.resize(num_equations_);
  Fu(0) = momentum(dimid - 1);
  real_t speed = Fu(0) / density;
  for (int d = 0; d < dim_; d++)
  {
    Fu(d + 1) = momentum(d) * Fu(0) / density;
  };
  Fu(dim_ + 1) = 
  
  

};

} // namespace QUEST


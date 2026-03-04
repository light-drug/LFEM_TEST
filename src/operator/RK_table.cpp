#include "RK_table.hpp"

namespace QUEST
{

RK_table::RK_table(const int& t_order)
  : t_order_(t_order) {};

const Matrix RK_table::getA() const
{
  return A_;
};

const Vector RK_table::getb() const
{
  return b_;
};

const Vector RK_table::getc() const
{
  return c_;
};

const int RK_table::getstages() const
{
  return stages_;
};

const int RK_table::gett_order() const
{
  return t_order_;
};

void RK_table::printall() const
{ 
  std::cout << " RK method of order " << t_order_ << " with " 
            << stages_ << " stages: " << std::endl;
  std::cout << "\n";
  std::cout << " A = \n" << A_ << std::endl;
  std::cout << "\n";
  std::cout << " b = \n" << b_.transpose() << std::endl;
  std::cout << "\n";
  std::cout << " c = \n" << c_.transpose() << std::endl;
  std::cout << "\n";
};

EX_TVDRK::EX_TVDRK(const int& t_order)
  : RK_table(t_order) {};

void EX_TVDRK::init() {
  switch (t_order_)
  {
  case 1: stages_ = 2; break;
  case 2: stages_ = 3; break;
  case 3: stages_ = 4; break;
  default:
    QUEST_ERROR(" EX_TVDRK only support order <= 3 and >= 1! ");
    break;
  }
  Matrix zero_mat = Matrix::Zero(stages_, stages_);
  Vector zero_vec = Vector::Zero(stages_);
  A_ = zero_mat;
  b_ = zero_vec;
  c_ = zero_vec;
  Be_ = zero_mat;

  
  if (t_order_ == 1) {
    c_(0) = 0.e0;
    
    c_(1) = 1.e0;
    A_(1, 0) = 1.e0;
    Be_(1, 0) = 1.e0;

  } else if (t_order_ == 2) {
    c_(0) = 0.e0;
    
    c_(1) = 1.e0;
    A_(1, 0) = 1.e0;
    Be_(1, 0) = 1.e0;

    c_(2) = 1.e0;
    A_(2, 0) = 0.5e0;
    A_(2, 1) = 0.5e0;
    Be_(2, 1) = 0.5e0;
  } else if (t_order_ == 3) {
    c_(0) = 1.e0;
    
    c_(1) = 1.e0;
    A_(1, 0) = 1.0e0;
    Be_(1, 0) = 1.e0;

    c_(2) = 0.5e0;
    A_(2, 0) = 0.75e0;
    A_(2, 1) = 0.25e0;
    Be_(2, 1) = 0.25e0;

    c_(3) = 1.e0;
    A_(3, 0) = 1.e0 / 3.e0;
    A_(3, 2) = 1.e0 / 3.e0;
    Be_(3, 2) = 2.e0 / 3.e0;
  }
  
};

const Vector EX_TVDRK::getb() const
{
  QUEST_ERROR(" EX_TVDRK has no use of vector b! ");
  return b_;
};

const Matrix EX_TVDRK::getBe() const
{
  return Be_;
};

void EX_TVDRK::printall() const
{ 
  std::cout << " TVD-RK (SSP-RK) method of order " << t_order_ << " with " 
            << stages_ << " stages: " << std::endl;
  std::cout << " A = \n" << A_ << std::endl;
  std::cout << " Be = \n" << Be_ << std::endl;
  std::cout << " c = \n" << c_.transpose() << std::endl;
};

IMEX_RK::IMEX_RK(const int& t_order)
  : RK_table(t_order) {};

void IMEX_RK::init() {
  switch (t_order_)
  {
  case 1: stages_ = 2; break;
  case 2: stages_ = 3; break;
  case 3: stages_ = 5; break;
  default:
    QUEST_ERROR(" EX_TVDRK only support order <= 3 and >= 1! ");
    break;
  }
  Matrix zero_mat = Matrix::Zero(stages_, stages_);
  Vector zero_vec = Vector::Zero(stages_);
  A_ = zero_mat;
  b_ = zero_vec;
  c_ = zero_vec;

  Ai_ = zero_mat;
  bi_ = zero_vec;
  ci_ = zero_vec;

  
  if (t_order_ == 1) {
    // ARS(1,1,1)
    c_(0) = 0.e0;
    
    c_(1) = 1.e0;
    A_(1, 0) = 1.e0;
    Ai_(1, 1) = 1.e0;
    b_ = A_.row(1).transpose();
    bi_ = Ai_.row(1).transpose();
    
    ci_ = c_;
  } else if (t_order_ == 2) {
    // ARS(2,2,2)
    real_t ga = 1.0e0 - 1.0e0 / std::sqrt(2.0e0);
    real_t de = 1.e0 - 1.0e0 / (2.0e0 * ga);
    c_(0) = 0.e0;
    
    c_(1) = ga;
    A_(1, 0) = ga;
    Ai_(1, 1) = ga;
    
    c_(2) = 1.e0;
    A_(2, 0) = de;
    A_(2, 1) = 1.e0 - de;
    b_ = A_.row(2).transpose();
    Ai_(2, 1) = 1.e0 - ga;
    Ai_(2, 2) = ga;
    bi_ = Ai_.row(2).transpose();

    ci_ = c_; 
  } else if (t_order_ == 3) {
    // ARS(3,4,3)
    c_(0) = 0.e0;

    c_(1) = 0.5e0;
    A_(1, 0) = 0.5e0;
    Ai_(1, 1) = 0.5e0;

    c_(2) = 2.e0 / 3.e0;
    A_(2, 0) = 11.e0 / 18.e0; A_(2, 1) = 1.e0 / 18.e0;
    Ai_(2, 1) = 1.e0 / 6.e0; Ai_(2, 2) = 0.5e0;

    c_(3) = 0.5e0;
    A_(3, 0) = 5.e0 / 6.e0; A_(3, 1) = - 5.e0 / 6.e0; A_(3, 2) = 0.5e0;
    Ai_(3, 1) = - 0.5e0; Ai_(3, 2) = 0.5e0; Ai_(3, 3) = 0.5e0;

    c_(4) = 1.e0;
    A_(4, 0) = 1.e0 / 4.e0; A_(4, 1) = 7.e0 / 4.e0; 
    A_(4, 2) = 3.e0 / 4.e0; A_(4, 3) = - 7.e0 / 4.e0;
    b_ = A_.row(4).transpose();
    Ai_(4, 1) = 3.e0 / 2.e0; Ai_(4, 2) = - 3.e0 / 2.e0; 
    Ai_(4, 3) = 0.5e0; Ai_(4, 4) = 0.5e0;
    bi_ = Ai_.row(4).transpose();

    ci_ = c_;
  }
};

const Matrix IMEX_RK::getAi() const
{
  return Ai_;
};

const Vector IMEX_RK::getbi() const
{
  return bi_;
};

const Vector IMEX_RK::getci() const
{
  return ci_;
};

void IMEX_RK::printall() const
{ 
  std::cout << " IMEX-RK method of order " << t_order_ << " with " 
            << stages_ << " stages: " << std::endl;
  std::cout << "\n";
  std::cout << " Explicit A = \n" << A_ << std::endl;
  std::cout << "\n";
  std::cout << " Explicit b = \n" << b_.transpose() << std::endl;
  std::cout << "\n";
  std::cout << " Explicit c = \n" << c_.transpose() << std::endl;
  std::cout << "\n";
  std::cout << " Implicit Ai = \n" << Ai_ << std::endl;
  std::cout << "\n";
  std::cout << " Implicit bi = \n" << bi_.transpose() << std::endl;
  std::cout << "\n";
  std::cout << " Implicit ci = \n" << ci_.transpose() << std::endl;
};

} // namespace QUEST


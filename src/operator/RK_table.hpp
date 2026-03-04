#ifndef QUEST_RK_TABLE_HPP
#define QUEST_RK_TABLE_HPP

#include <iostream>

#include "config.hpp"
#include "error.hpp"

namespace QUEST
{

class RK_table
{
protected:

  Matrix A_;
  Vector b_;
  Vector c_;

  int t_order_;
  int stages_;
public:
  RK_table(const int& t_order);
  virtual ~RK_table() = default;

  virtual void init() = 0;
  virtual const Matrix getA() const;
  virtual const Vector getb() const;
  virtual const Vector getc() const;

  virtual const int getstages() const;
  virtual const int gett_order() const;
  virtual void printall() const;
};


class EX_TVDRK : public RK_table
{
protected:
  Matrix Be_;
public:
  EX_TVDRK(const int& t_order);
  ~EX_TVDRK() override = default;

  void init() override;

  const Vector getb() const override;

  virtual const Matrix getBe() const;
  void printall() const override;

};


class IMEX_RK : public RK_table
{
protected:

  Matrix Ai_;
  Vector bi_;
  Vector ci_;
  
public:
  IMEX_RK(const int& t_order);
  ~IMEX_RK() override = default;

  void init() override;

  virtual const Matrix getAi() const;
  virtual const Vector getbi() const;
  virtual const Vector getci() const;
  void printall() const override;

};

} // namespace QUEST

#endif // QUEST_RK_TABLE_HPP

#ifndef QUEST_FDMESH_HPP
#define QUEST_FDMESH_HPP

#include "Tensormesh.hpp"
#include "integralrule.hpp"
#include "config.hpp"
#include "error.hpp"
#include "vector_overload.hpp"

namespace QUEST
{

class FDmesh
{
public:
  FDmesh(const real_t& x1, const real_t& x2,
        const int& xDiv);
  virtual ~FDmesh() = default;
  virtual void init();

  virtual const real_t& getx1() const;
  virtual const real_t& getx2() const;
  virtual const int& getxDiv() const;
  virtual const real_t& gethx() const;
  virtual const int& getNx() const;
  virtual const Vector& getxpoints() const;
  virtual const Vector& getX() const;
  virtual const Vector& getCellCenter_vec() const;
  virtual const Matrix& getCellCenter() const;

  virtual void DisplayResult_withX(const Vector& xtemp,
                    const Vector& numerical, const Vector& exact,  
                    const std::string& title, std::ostream& Outfile) const;
  virtual void DisplayResult(const Vector& numerical, const Vector& exact,  
                    const std::string& title, std::ostream& Outfile) const;
  virtual void DisplayResult_ex(const Vector& numerical,  
                    const std::string& title, std::ostream& Outfile) const;

  virtual void computerrorL1(const Vector& numerical, const Vector& exact, 
                          real_t* error);
  virtual void computerrorL2(const Vector& numerical, const Vector& exact, 
                          real_t* error);
  virtual void computerrorLinf(const Vector& numerical, const Vector& exact, 
                          real_t* error);

protected:
  const real_t& x1_;
  const real_t& x2_;
  const int& xDiv_;

  real_t hx_;
  int Nx_;
  Vector xpoints_;
  Vector X_;
  Matrix CellCenter_;
  Vector CellCenter_vec_;

  virtual void generatehx();
  virtual void generatexpoints();
  virtual void generateCellCenter();
};

class FDmesh_period : public FDmesh
{
public: 

  FDmesh_period(const real_t& x1, const real_t& x2,
        const int& xDiv);
  ~FDmesh_period() override = default;

protected:

  void generatexpoints() override;
  
};

class KineticFDmesh : public FDmesh
{
public:
  KineticFDmesh(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);
  ~KineticFDmesh() override = default;
  void init() override;

  virtual const real_t& getv1() const;
  virtual const real_t& getv2() const;
  virtual const int& getvDiv() const;
  virtual const real_t& gethv() const;
  virtual const int& getNv() const;
  virtual const Vector& getvpoints() const;
  virtual const Vector& getV() const;
  virtual const Vector& getvweights() const;
  
  using FDmesh::DisplayResult; 
  using FDmesh::computerrorL1;
  using FDmesh::computerrorL2;
  using FDmesh::computerrorLinf;
  virtual void DisplayResult(const std::vector<Vector>& numerical, 
                    const std::vector<Vector>& exact,  
                    const std::string& title, std::ostream& Outfile) const;

  virtual void computerrorL1(const std::vector<Vector>& numerical, 
                          const std::vector<Vector>& exact, 
                          real_t* error);
  virtual void computerrorL2(const std::vector<Vector>& numerical, 
                          const std::vector<Vector>& exact, 
                          real_t* error);
  virtual void computerrorLinf(const std::vector<Vector>& numerical, 
                          const std::vector<Vector>& exact, 
                          real_t* error);

protected:

  const real_t& v1_;
  const real_t& v2_;
  const int& vDiv_;

  real_t hv_;
  int Nv_;
  Vector vpoints_;
  Vector V_;
  Vector vweights_;

  virtual void generatehv();
  virtual void generatevpoints();
  virtual void generatevweights();
};

// class KineticFVmesh : public KineticFDmesh
// {
// public:
//   KineticFVmesh(const real_t& x1, const real_t& x2,
//         const real_t& v1, const real_t& v2,
//         const int& xDiv, const int& vDiv,
//         const int& qua_num, const QuadratureType& qua_type);
//   ~KineticFVmesh() override = default;

//   void init() override;

// protected:
//   const int& qua_num_;
//   const QuadratureType& qua_type_;

//   Matrix xPb_;
//   Matrix vPb_;
//   Vector wqua_;
//   Vector xqua_;

//   virtual const Matrix& getxPb() const;
//   virtual const Matrix& getvPb() const;
//   virtual const Vector& getwqua() const;
//   virtual const Vector& getxqua() const;
//   virtual const int& getquanum() const;
  
//   void generatexpoints() override;
//   void generatevpoints() override;
//   void generatevweights() override;
//   virtual void generatexPb() override;
//   virtual void generatevPb() override;
  
//   template<typename Func>
//   void average_initial(Func& func, Vector* M) const
//   {
//     Vector point(Nx_);
//     Vector f_quavalue;
//     M->resize(Nx);
//     M->setZero();
//     for (int i = 0; i < qua_num_ i++)
//     {
//       point = xPb_.row(i).transpose();
//       f_quavalue = std::forward<Func>(func)(point);
//       *M = *M + wqua_(i) * f_quavalue * hx_;
//     };
//   };

//   template<typename Func>
//   void average_final(Func& func, const real_t& Tstop, Vector* M) const
//   {
//     Vector point(Nx_);
//     Vector f_quavalue;
//     M->resize(Nx);
//     M->setZero();
//     for (int i = 0; i < qua_num_ i++)
//     {
//       point = xPb_.row(i).transpose();
//       f_quavalue = std::forward<Func>(func)(point, Tstop);
//       *M = *M + wqua_(i) * f_quavalue * hx_;
//     };
//   };

//   template<typename Func>
//   void average_initial(Func& func, std::vector<Vector>* M) const
//   {
//     Vector point(Nx_);
//     Vector f_quavalue;
//     M->resize(Nv_);
//     M->setZero();
//     Vector temp;
//     temp.resize(Nx_); 
//     temp.setZero();
//     for (int j = 0; j < Nv_; j++)
//     {
//       M->at(j).resize(Nx_);
//       M->at(j).setZero();
//       for (int iv = 0; iv < qua_num_; iv++)
//       {
//         for (int i = 0; i < qua_num_ i++)
//         {
//           point = xPb_.row(i).transpose();
//           f_quavalue = std::forward<Func>(func)(point, vPb_(iv, j));
//           temp = temp + wqua_(i) * wqua_(iv) * f_quavalue * hx_ * hv_;
//         };
//       };
//       M->at(j) == temp;
//     };
//   };

//   template<typename Func>
//   void average_final(Func& func, const real_t& Tstop, std::vector<Vector>* M) const
//   {
//     Vector point(Nx_);
//     Vector f_quavalue;
//     M->resize(Nv_);
//     M->setZero();
//     Vector temp;
//     temp.resize(Nx_); 
//     temp.setZero();
//     for (int j = 0; j < Nv_; j++)
//     {
//       M->at(j).resize(Nx_);
//       M->at(j).setZero();
//       for (int iv = 0; iv < qua_num_; iv++)
//       {
//         for (int i = 0; i < qua_num_ i++)
//         {
//           point = xPb_.row(i).transpose();
//           f_quavalue = std::forward<Func>(func)(point, vPb_(iv, j), Tstop);
//           temp = temp + wqua_(i) * wqua_(iv) * f_quavalue * hx_ * hv_;
//         };
//       };
//       M->at(j) == temp;
//     };
//   };
// };

class KineticFDmesh_Gauss : 
  public KineticFDmesh
{
public: 
  KineticFDmesh_Gauss(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);

  ~KineticFDmesh_Gauss() override = default;

protected:
  void generatevpoints() override;
  void generatevweights() override; 
};

class KineticFDmesh_period : 
  public KineticFDmesh
{
public: 
  KineticFDmesh_period(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);

  ~KineticFDmesh_period() override = default;

protected:

  void generatexpoints() override;
  
};

class KineticFDmesh_period_Gauss : 
  public KineticFDmesh_period
{
public: 
  KineticFDmesh_period_Gauss(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);

  ~KineticFDmesh_period_Gauss() override = default;

protected:

  void generatevpoints() override;
  void generatevweights() override; 
  
};

class KineticFDmesh_twovel : public KineticFDmesh
{
public: 
  KineticFDmesh_twovel(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);

  ~KineticFDmesh_twovel() override = default;

protected:

  void generatevpoints() override;
  void generatevweights() override;
  
};


class KineticFDmesh_period_twovel : public KineticFDmesh_period
{
public: 
  KineticFDmesh_period_twovel(const real_t& x1, const real_t& x2,
        const real_t& v1, const real_t& v2,
        const int& xDiv, const int& vDiv);

  ~KineticFDmesh_period_twovel() override = default;

protected:

  void generatevpoints() override;
  void generatevweights() override;
  
};

} // namespace QUEST

#endif // QUEST_FDMESH_HPP

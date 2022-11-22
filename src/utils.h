#ifndef __UTILS_H__
#define __UTILS_H__

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <functional>
#include <random>
#include <cstdio>
#include <locale>
#include <unistd.h>
#include <sys/stat.h>
#include "libmesh/point.h"
#include "libmesh/getpot.h"
#include "libmesh/parallel.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/time_solver.h"
#include "libmesh/diff_solver.h"
#include "libmesh/fe.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/boundary_info.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/fem_system.h"
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/auto_ptr.h"

using namespace libMesh;

// Solaris Studio has no NAN
#ifdef __SUNPRO_CC
#define NAN (1.0/0.0)
#endif

//-------------------------------------------------------------------------------------------------
inline Real pow2 (Real v) { return v*v; }
inline Real pow3 (Real v) { return v*pow2(v); }
inline Real pow4 (Real v) { return pow2(pow2(v)); }
inline Real pow5 (Real v) { return v*pow4(v); }
inline Real pow6 (Real v) { return pow2(pow3(v)); }
inline Real pow7 (Real v) { return v*pow6(v); }
inline Real pow8 (Real v) { return pow2(pow4(v)); }
inline Real pow9 (Real v) { return pow3(pow3(v)); }
//-------------------------------------------------------------------------------------------------
inline
Real degrees_to_radians (const Real d) { return d*(libMesh::pi/180.0); }
inline
Real radians_to_degrees (const Real r) { return r*(180.0/libMesh::pi); }
//-------------------------------------------------------------------------------------------------
inline
Real apply_lbound (const Real& L, const Real& X) { return ( X < L ? L : X ); }
inline
Real apply_ubound (const Real& X, const Real& U) { return ( X > U ? U : X ); }
inline
Real apply_bounds (const Real& L, const Real& X, const Real& U) { return ( X < L ? L : ( X > U ? U : X ) ); }
//-------------------------------------------------------------------------------------------------
inline
int sign (Real r, Real tol =0.0)
{
  if (r >  tol) return  1;
  if (r < -tol) return -1;
  return 0;
}
//-------------------------------------------------------------------------------------------------
inline // Pi or rectangular function
Real Pi_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  0.0;
  else if (C < c1) return  cM;
  else             return  0.0;
}
//-------------------------------------------------------------------------------------------------
inline // step-decay function
Real SD_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  cM;
  else if (C < c1) return  cM*(c1-C)/(c1-c0);
  else             return  0.0;
}
inline // step-decay function - derivative
Real deriv_SD_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  0.0;
  else if (C < c1) return -cM/(c1-c0);
  else             return  0.0;
}
//-------------------------------------------------------------------------------------------------
inline // step-growth function
Real SG_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  cM;
  else if (C < c1) return  cM*(C-c0)/(c1-c0);
  else             return  0.0;
}
inline // step-growth function - derivative
Real deriv_SG_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  0.0;
  else if (C < c1) return  cM/(c1-c0);
  else             return  0.0;
}
//-------------------------------------------------------------------------------------------------
inline // trapezoidal function
Real Tr_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  const Real & c2 = p_[3];
  const Real & c3 = p_[4];
  if      (C < c0) return  0.0;
  else if (C < c1) return  cM*(C-c0)/(c1-c0);
  else if (C < c2) return  cM;
  else if (C < c3) return  cM*(c3-C)/(c3-c2);
  else             return  0.0;
}
inline // trapezoidal function - derivative
Real deriv_Tr_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  const Real & c2 = p_[3];
  const Real & c3 = p_[4];
  if      (C < c0) return  0.0;
  else if (C < c1) return  cM/(c1-c0);
  else if (C < c2) return  0.0;
  else if (C < c3) return -cM/(c3-c2);
  else             return  0.0;
}
//-------------------------------------------------------------------------------------------------
inline // Boltzmann (sigmoidal) increase function
Real Bsi_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return G/(1.0+G);
}
inline // Boltzmann (sigmoidal) increase function - derivative
Real deriv_Bsi_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return G/(dC*(1.0+G)*(1.0+G));
}
//-------------------------------------------------------------------------------------------------
inline // Boltzmann (sigmoidal) dencrease function
Real Bsd_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return 1.0/(1.0+G);
}
inline // Boltzmann (sigmoidal) dencrease function - derivative
Real deriv_Bsd_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0>=cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return -G/(dC*(1.0+G)*(1.0+G));
}
//-------------------------------------------------------------------------------------------------
inline
double uniform_rand (const double from, const double to) {
    std::random_device rd;
    std::default_random_engine rgen(rd());
    std::uniform_real_distribution<double> dist(from, to);
    return dist(rgen);
}
inline
double normal_rand (const double mean, const double stdev) {
    std::random_device rd;
    std::mt19937 rgen(rd());
    std::normal_distribution<double> dist(mean, stdev);
    return dist(rgen);
}
//-------------------------------------------------------------------------------------------------
inline
std::ostream& format_date_time (std::ostream& out, const std::tm& t, const char* fmt) {
  const std::time_put<char>& dateWriter =
    std::use_facet< std::time_put<char> >(out.getloc());
  int n = strlen(fmt);
  if (dateWriter.put(out, out, ' ', &t, fmt, fmt+n).failed())
    throw std::runtime_error("failure to format date time");
  return out;
}
inline
std::string date_time_to_string (const std::tm& t, const char* format) {
  std::stringstream s;
  format_date_time(s, t, format);
  return s.str();
}
inline
std::tm date_now () {
  std::time_t now = std::time(0);
  return *std::localtime(&now);
}
//-------------------------------------------------------------------------------------------------
inline
std::set<int> export_integers (const std::string& s)
{
  std::set<int> numbers;
  // create a copy of the original string
  std::stringstream ss;
  ss << s;
  // iterate till the end of the stream
  std::string tmp;
  while (!ss.eof())
    {
      // extract word by word from the stream
      ss >> tmp;
      // check the given word is integer or not
      int n;
      if (std::stringstream(tmp) >> n)
          numbers.insert(n);
      // to save from space at the end of string
      tmp = "";
    }
  return numbers;
}
//-------------------------------------------------------------------------------------------------
template <typename T>
inline
std::string number2string (const T& Number) {
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}
template <typename T>
inline
T string2number (const std::string& Text) {
    std::istringstream ss(Text);
    T result;
    return ( ss >> result ? result : 0 );
}
//-------------------------------------------------------------------------------------------------
inline
Point rotate (const Point& v, const Real theta_x, const Real theta_y, const Real theta_z) {
    const Real v_x = v(0),
               v_y = v(1),
               v_z = v(2);
    const Real S_x = sin(theta_x), C_x = cos(theta_x),
               S_y = sin(theta_y), C_y = cos(theta_y),
               S_z = sin(theta_z), C_z = cos(theta_z);
    Point r;
    r(0) = v_z*(S_x*S_z + C_x*C_z*S_y) - v_y*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*v_x;
    r(1) = v_y*(C_x*C_z + S_x*S_y*S_z) - v_z*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*v_x;
    r(2) = C_x*C_y*v_z - S_y*v_x + C_y*S_x*v_y;
    return r;
}
//-------------------------------------------------------------------------------------------------
inline
RealTensorValue tensor (const RealVectorValue& a, const RealVectorValue& b)
{
  RealTensorValue a_b;
  a_b(0,0) = a(0)*b(0); a_b(0,1) = a(0)*b(1); a_b(0,2) = a(0)*b(2);
  a_b(1,0) = a(1)*b(0); a_b(1,1) = a(1)*b(1); a_b(1,2) = a(1)*b(2);
  a_b(2,0) = a(2)*b(0); a_b(2,1) = a(2)*b(1); a_b(2,2) = a(2)*b(2);
  return a_b;
}
inline
RealTensorValue tensor (const RealVectorValue& a)
{
  RealTensorValue a_a;
  a_a(0,0) = a(0)*a(0); a_a(0,1) = a(0)*a(1); a_a(0,2) = a(0)*a(2);
  a_a(1,0) = a_a(0,1);  a_a(1,1) = a(1)*a(1); a_a(1,2) = a(1)*a(2);
  a_a(2,0) = a_a(0,2);  a_a(2,1) = a_a(1,2);  a_a(2,2) = a(2)*a(2);
  return a_a;
}
//-------------------------------------------------------------------------------------------------
// Symmetric matrix "A":
// eigenvectors in columns of "eVec" that correspond to eigenvalues in vector "eVal"
void eigen_decomposition(double A[3][3], double eVec[3][3], double eVal[3]);
//-------------------------------------------------------------------------------------------------

#include "./ida.h"
#include "./paraview.h"

#endif // __UTILS_H__

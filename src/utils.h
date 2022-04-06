#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include "libmesh/point.h"
#include "libmesh/getpot.h"
#include "libmesh/parallel.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"

using namespace libMesh;

//-------------------------------------------------------------------------------------------------
inline
Real apply_lbound (const Real& L, const Real& X) { return ( X < L ? L : X ); }
inline
Real apply_ubound (const Real& X, const Real& U) { return ( X > U ? U : X ); }
inline
Real apply_bounds (const Real& L, const Real& X, const Real& U) { return ( X < L ? L : ( X > U ? U : X ) ); }
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

#endif // __UTILS_H__

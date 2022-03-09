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

#endif // __UTILS_H__

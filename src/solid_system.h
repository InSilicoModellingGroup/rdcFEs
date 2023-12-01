// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Source has been modified based on the original code produced by Robert Weidlich (2012)


#ifndef __SOLID_SYSTEM_H__
#define __SOLID_SYSTEM_H__

#include "./hyperelastic.h"

using namespace libMesh;

//-------------------------------------------------------------------------------------------------
class SolidSystem: public FEMSystem
{
public:
  // Constructor
  SolidSystem (EquationSystems& es, std::string name, unsigned int number) :
    FEMSystem(es, name, number)
  {
    // Add a time solver. We are just looking at a steady state problem here.
    this->time_solver = libmesh_make_unique<SteadySolver>(*this);
  }

  // system initialization
  virtual
  void init_data ();

  // context initialization
  virtual
  void init_context (DiffContext & );

  // finite element residual and jacobian calculations
  virtual
  bool element_time_derivative (bool , DiffContext & );

  // contributions for adding boundary conditions
  virtual
  bool side_time_derivative (bool , DiffContext & );

  virtual
  bool eulerian_residual (bool , DiffContext & ) { return false; }

  virtual
  std::string system_type () const { return "SolidSystem"; }

  // override method to update mesh also
  virtual void update ();

  // save the undeformed mesh to an auxiliary system
  void save_initial_mesh ();

  // calculate gradient fields after system solution
  void post_process ();

  // update the auxiliary system data after system solution
  void update_auxiliary ();

public:
  // variable numbers of primary variables in the current system
  unsigned int var[3];
  // variable numbers of primary for an auxiliary system (used to store
  // the reference / undeformed configuration data)
  unsigned int undefo_var[3];
};
//-------------------------------------------------------------------------------------------------

#endif // __SOLID_SYSTEM_H__

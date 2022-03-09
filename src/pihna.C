#include "./utils.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/getpot.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

static void input (GetPot & , EquationSystems & );
static void initial_structure (EquationSystems & , const std::string & );
static void initial_pihna (EquationSystems & , const std::string & );
static void assemble_pihna (EquationSystems & , const std::string & );

extern PerfLog plog;

void pihna (LibMeshInit & init)
{

  // ...done
}

void input (GetPot & in, EquationSystems & es)
{
  std::string name;

  // ...done
}

void initial_structure (EquationSystems & es,
                        const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "uStructure");

  // ...done
}

void initial_pihna (EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "PIHNA");

  // ...done
}

void assemble_pihna (EquationSystems & es,
                     const std::string & system_name)
{
  libmesh_ignore(es, system_name);
  libmesh_assert_equal_to(system_name, "PIHNA");
  libmesh_assert_equal_to(system.n_vars(), 5);

  // ...done
}

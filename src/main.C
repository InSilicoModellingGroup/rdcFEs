#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/perf_log.h"

using namespace libMesh;

PerfLog plog("rdcFEs");

int main (int argc, char* argv[])
{
  LibMeshInit init(argc, argv);

  libMesh::out << std::endl << plog.get_log() << std::endl;
  // ...done
  return 0;
}

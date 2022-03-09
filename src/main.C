#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/perf_log.h"

using namespace libMesh;

PerfLog plog("rdcFEs");

extern void adpm (LibMeshInit & );

int main (int argc, char* argv[])
{
  LibMeshInit init(argc, argv);
  GetPot command_line(argc, argv);

  std::string model;
  if (command_line.search(1, "-m"))
    model = command_line.next(model);
  else return 1;

  // Alzheimer's disease progression model
  if ("adpm"==model) adpm(init);

  libMesh::out << std::endl << plog.get_log() << std::endl;
  // ...done
  return 0;
}

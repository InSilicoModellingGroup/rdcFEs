#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/perf_log.h"

using namespace libMesh;

PerfLog plog("rdcFEs");

extern void adpm (LibMeshInit & );
extern void pihna (LibMeshInit & );
extern void ripf (LibMeshInit & );

int main (int argc, char* argv[])
{
  LibMeshInit init(argc, argv);
  GetPot command_line(argc, argv);

  std::string s;
  if (command_line.search(1, "-m"))
    {
      s = command_line.next(s);
      // Alzheimer's disease progression model
      if ("adpm"==s) adpm(init);
      // PIHNA cancer model
      else if ("pihna"==s) pihna(init);
      // radiation-induced pulmonary fibrosis model
      else if ("ripf"==s) ripf(init);
    }
  else return 1;

  libMesh::out << std::endl << plog.get_log() << std::endl;
  // ...done
  return 0;
}

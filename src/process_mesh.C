#include "./utils.h"

void write_mesh (Mesh& mesh, const char* fname)
{
  // Build a list of (elem, side, bc) tuples.
  auto boundary_element_triples = mesh.get_boundary_info().build_side_list();
  std::cout << boundary_element_triples.size() << std::endl;

  int n_nodes = mesh.n_nodes();
  int n_cells = 0;
  for (const auto & elem : mesh.active_element_ptr_range())
    {
      for (auto s : elem->side_index_range())
        if (nullptr == elem->neighbor_ptr(s))
          ++n_cells;
      ++n_cells;
    }
  int index;

  std::ofstream fmsh(fname);

  fmsh << "$MeshFormat" << std::endl;
  fmsh << "2.2 0 8" << std::endl;
  fmsh << "$EndMeshFormat" << std::endl;
  fmsh << "$Nodes" << std::endl;
  fmsh << n_nodes << std::endl;
  index = 1;
  for (const auto & node : mesh.node_ptr_range())
    {
      if (1+node->id()!=index) libmesh_error();
      const Point coord = (*node);
      fmsh << index
           << ' ' << coord(0)
           << ' ' << coord(1)
           << ' ' << coord(2)
           << std::endl;
      //
      ++index;
    }
  fmsh << "$EndNodes" << std::endl;
  fmsh << "$Elements" << std::endl;
  fmsh << n_cells << std::endl;
  index = 1;
  for (const auto & triplet : boundary_element_triples)
    {
      const Elem & elem = mesh.elem_ref(std::get<0>(triplet));
      //
      std::unique_ptr<const Elem> side =
        elem.build_side_ptr(std::get<1>(triplet));

      int gmsh_type;
      if      (TRI3 ==side->type()) gmsh_type=2;
      else if (QUAD4==side->type()) gmsh_type=3;
      else libmesh_error();
      //
      fmsh << index
           << ' ' << gmsh_type
           << ' ' << 2 << ' ' << std::get<2>(triplet) << ' ' << 0;
      for (unsigned int i=0; i<side->n_nodes(); i++)
        fmsh << ' ' << (1+side->node_id(i));
      fmsh << std::endl;
      //
      ++index;
    }
  for (const auto & elem : mesh.active_element_ptr_range())
    {
      int gmsh_type;
      if      (TET4    ==elem->type()) gmsh_type=4;
      else if (HEX8    ==elem->type()) gmsh_type=5;
      else if (PRISM6  ==elem->type()) gmsh_type=6;
      else if (PYRAMID5==elem->type()) gmsh_type=7;
      else libmesh_error();
      //
      fmsh << index
           << ' ' << gmsh_type
           << ' ' << 2 << ' ' << elem->subdomain_id() << ' ' << 0;
      for (unsigned int i=0; i<elem->n_nodes(); i++)
        fmsh << ' ' << (1+elem->node_id(i));
      fmsh << std::endl;
      //
      ++index;
    }
  fmsh << "$EndElements" << std::endl;
}

void process_mesh (LibMeshInit & init)
{
  std::string in_str;
  // create the mesh object
  Mesh msh(init.comm(), 3);
  //
  std::cout << "Give name of the Gmsh-formatted input file: " << std::flush;
  std::cin >> in_str;
  const std::string input_Gmsh_file = in_str;
  //
  std::cout << std::endl << "FE mesh is now loading... " << std::flush;
  GmshIO(msh).read(input_Gmsh_file.c_str());
  std::cout << " ok" << std::endl;
  //
  std::cout << "Give value to scale nodes' coordinates: " << std::flush;
  std::cin >> in_str;
  const Real scale = string2number<Real>(in_str);
  //
  std::cout << "Give value for mesh translation (X-axis): " << std::flush;
  std::cin >> in_str;
  const Real translate_X = string2number<Real>(in_str);
  std::cout << "Give value for mesh translation (Y-axis): " << std::flush;
  std::cin >> in_str;
  const Real translate_Y = string2number<Real>(in_str);
  std::cout << "Give value for mesh translation (Z-axis): " << std::flush;
  std::cin >> in_str;
  const Real translate_Z = string2number<Real>(in_str);
  //
  std::cout << "Give value for mesh rotation (X-axis) in degrees: " << std::flush;
  std::cin >> in_str;
  const Real rotate_X = degrees_to_radians(string2number<Real>(in_str));
  std::cout << "Give value for mesh rotation (Y-axis) in degrees: " << std::flush;
  std::cin >> in_str;
  const Real rotate_Y = degrees_to_radians(string2number<Real>(in_str));
  std::cout << "Give value for mesh rotation (Z-axis) in degrees: " << std::flush;
  std::cin >> in_str;
  const Real rotate_Z = degrees_to_radians(string2number<Real>(in_str));
  //
  std::cout << "Skip node renumbering? True or false? Insert '1' or '0' respectively: " << std::flush;
  std::cin >> in_str;
  std::cout << "FE mesh is under preparation... " << std::flush;
  msh.prepare_for_use(string2number<int>(in_str));
  std::cout << " ok" << std::endl;
  // print out some info about the FE mesh
  std::cout << std::endl;
  msh.print_info();
  //msh.get_boundary_info().print_info();
  std::cout << std::endl;
  // process the space vector of all nodes of the FE mesh
  std::cout << "FE mesh is now under processing... " << std::flush;
  for (MeshBase::node_iterator ni=msh.nodes_begin(); ni!=msh.nodes_end(); ni++)
    {
      Node& node = *(*ni);
      // scale (up|down) the node coordinates
      Point xyz = scale * node;
      // translate the space vector
      xyz += Point(translate_X, translate_Y, translate_Z);
      // rotate the space vector
      xyz = rotate(xyz, rotate_X, rotate_Y, rotate_Z);
      // update the coordinates of the node
      for (unsigned int i=0; i<3; i++)
        node(i) = xyz(i);
    }
  std::cout << " ok" << std::endl;
  //
  std::cout << "Give name of the output files: " << std::flush;
  std::cin >> in_str;
  const std::string output_filename = in_str;
  //
  std::cout << "Mesh and configuration data is now saving... " << std::flush;
  write_mesh(msh, (output_filename+".msh").c_str());
  {
    std::ofstream fout(output_filename+".config");
    fout << "Gmsh input: " << input_Gmsh_file << std::endl;
    fout << "mesh scaling: " << scale << std::endl;
    fout << "translation (X-axis): " << translate_X << std::endl;
    fout << "translation (Y-axis): " << translate_Y << std::endl;
    fout << "translation (Z-axis): " << translate_Z << std::endl;
    fout << "rotation (X-axis) in degrees: " << radians_to_degrees(rotate_X) << std::endl;
    fout << "rotation (Y-axis) in degrees: " << radians_to_degrees(rotate_Y) << std::endl;
    fout << "rotation (Z-axis) in degrees: " << radians_to_degrees(rotate_Z) << std::endl;
    fout << "output file name: " << output_filename << std::endl;
  }
  ExodusII_IO(msh).write((output_filename+".ex2").c_str());
  std::cout << " ok" << std::endl;
  // ...done
}

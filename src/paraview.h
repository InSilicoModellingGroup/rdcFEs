#ifndef __PARAVIEW_H__
#define __PARAVIEW_H__

//-------------------------------------------------------------------------------------------------
#ifdef SMALLEST_NUMBER
#error "macro parameter \"SMALLEST_NUMBER\" has been defined elsewhere!"
#endif
#define SMALLEST_NUMBER   1.0e-24
//-------------------------------------------------------------------------------------------------
class Paraview_IO : public MeshInput<MeshBase>, public MeshOutput<MeshBase> {
public:
  Paraview_IO (MeshBase & mesh) : MeshInput<MeshBase>(mesh), MeshOutput<MeshBase>(mesh) {}
  ~Paraview_IO () { this->close_pvd(); }
  //
  virtual
  void write_nodal_data ( const std::string & fn, EquationSystems & es,
                          const std::set<std::string> * system_names = nullptr )
    {
      std::vector<std::string> names;
      es.build_variable_names(names, nullptr, system_names);

      std::vector<Number> soln;
      es.build_solution_vector(soln);

      this->write_nodal_data(fn, soln, names);
    }
  virtual
  void write_nodal_data ( const std::string & fn, const std::vector<Number> & soln,
                          const std::vector<std::string> & names )
    {
      if (libMesh::global_processor_id()) return;

      const MeshBase& msh = MeshOutput<MeshBase>::mesh();
      MeshBase::const_node_iterator cni_begin = msh.nodes_begin();
      MeshBase::const_node_iterator cni_end   = msh.nodes_end();
      MeshBase::const_element_iterator cei_begin = msh.active_elements_begin();
      MeshBase::const_element_iterator cei_end   = msh.active_elements_end();

      std::vector<bool> nodes_bool(msh.n_nodes(), false);
      unsigned int n_cell = 0;
      for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
        {
          for (unsigned int i=0; i<(*cei)->n_nodes(); i++)
            nodes_bool[(*cei)->node_id(i)] = true;
          ++n_cell;
        }

      std::vector<int> nodes_vtk_id(msh.n_nodes(), -1);
      unsigned int n_point = 0;
      for (MeshBase::const_node_iterator cni=cni_begin; cni!=cni_end; cni++)
        {
          if ( ! nodes_bool[(*cni)->id()] )
            continue;
          nodes_vtk_id[(*cni)->id()] = n_point;
          ++n_point;
        }

      std::ofstream vtu(fn.c_str());
      // now start writing...
      if (vtu.is_open())
        {
          vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
          vtu << "  <UnstructuredGrid>" << std::endl;
          vtu << "    <Piece  NumberOfPoints=\"" << n_point << "\" NumberOfCells=\"" << n_cell << "\">" << std::endl;
          vtu << "      <Points>" << std::endl;
          vtu << "        <DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_node_iterator cni=cni_begin; cni!=cni_end; cni++)
            {
              if ( ! nodes_bool[(*cni)->id()] )
                continue;
              for (unsigned int i=0; i<3; i++)
                vtu << ' ' << (*cni)->operator()(i) << std::flush;
            }
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "      </Points>" << std::endl;
          vtu << "      <PointData>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"node_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_node_iterator cni=cni_begin; cni!=cni_end; cni++) {
              if ( ! nodes_bool[(*cni)->id()] )
                continue;
              vtu << ' ' << ((*cni)->id()+1) << std::flush;
          }
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          for (unsigned int j=0; j<names.size(); j++)
            {
              vtu << "        <DataArray type=\"Float64\" Name=\"" << names[j] << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
              for (MeshBase::const_node_iterator cni=cni_begin; cni!=cni_end; cni++)
                {
                  if ( ! nodes_bool[(*cni)->id()] )
                    continue;
                  const Real r = soln[(*cni)->id()*names.size()+j];
                  vtu << ' ' << (fabs(r)<=SMALLEST_NUMBER ? 0.0 : r) << std::flush;
                }
              vtu << std::endl;
              vtu << "        </DataArray>" << std::endl;
            }
          vtu << "      </PointData>" << std::endl;
          vtu << "      <CellData>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
            vtu << ' ' << 1+(*cei)->id() << std::flush;
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"region_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
            vtu << ' ' << (*cei)->subdomain_id() << std::flush;
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"processor_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
            vtu << ' ' << (*cei)->processor_id() << std::flush;
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "      </CellData>" << std::endl;
          vtu << "      <Cells>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
            {
              std::vector<unsigned int> conn;
              this->_connectivity(*cei, conn);
              for (auto c : conn)
                vtu << ' ' << nodes_vtk_id[c] << std::flush;
            }
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          unsigned int offset = 0;
          for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
            {
              offset += this->_offset(*cei);
              vtu << ' ' << offset << std::flush;
            }
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "        <DataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (MeshBase::const_element_iterator cei=cei_begin; cei!=cei_end; cei++)
            vtu << ' ' << this->_type(*cei) << std::flush;
          vtu << std::endl;
          vtu << "        </DataArray>" << std::endl;
          vtu << "      </Cells>" << std::endl;
          vtu << "    </Piece>" << std::endl;
          vtu << "  </UnstructuredGrid>" << std::endl;
          vtu << "</VTKFile>" << std::endl;
          // ...done
        }
    }
  //
  virtual
  void read  (const std::string& name) { libmesh_error(); }
  virtual
  void write (const std::string& name) { libmesh_error(); }
  //
  inline
  void open_pvd (std::string fn)
    {
      // sanity check
      if (this->_pvd.is_open()) libmesh_error();

      this->_filename = fn;

      this->_pvd.open(this->_filename + ".pvd");
      if (!this->_pvd.good()) libmesh_error();
      if (!this->_pvd.is_open()) libmesh_error();
      //
      this->_pvd << "<?xml version=\"1.0\"?>\n"
                 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                 << "  <Collection>\n"
                 << std::flush;
    }
  inline
  void update_pvd (EquationSystems & es, unsigned int t = 0)
    {
      // sanity check
      if (!this->_pvd.is_open()) libmesh_error();

      const std::string vtu = this->_filename + "-" + std::to_string(t) + ".vtu";

      this->write_equation_systems(vtu, es);

      this->_pvd << "    <DataSet timestep=\"" << t << "\" group=\"\" part=\"0\" file=\"../" << vtu << "\"/>\n"
                 << std::flush;
    }
  inline
  void close_pvd ()
    {
      if (!this->_pvd.is_open()) return;

      this->_pvd << "  </Collection>\n"
                 << "</VTKFile>\n"
                 << std::flush;
      //
      this->_pvd.close();
    }

private:
  // obtain the connectivity of the element (cell)
  inline
  void _connectivity (const Elem * elem, std::vector<unsigned int> & c) const
    {
      c.clear();
      elem->connectivity(0, IOPackage::VTK, c);
    }
  // obtain the type of the element (cell)
  inline
  unsigned int _type (const Elem * elem) const
    {
      unsigned int t(0);
      switch ( elem->type() )
        {
          case ElemType::EDGE2:    t = 3;  break;// VTK_LINE
          case ElemType::EDGE3:    t = 21; break;// VTK_QUADRATIC_EDGE
          case ElemType::TRI3:     t = 5;  break;// VTK_TRIANGLE
          case ElemType::TRI6:     t = 22; break;// VTK_QUADRATIC_TRIANGLE
          case ElemType::QUAD4:    t = 9;  break;// VTK_QUAD
          case ElemType::QUAD8:    t = 23; break;// VTK_QUADRATIC_QUAD
          case ElemType::QUAD9:    t = 28; break;// VTK_BIQUADRATIC_QUAD
          case ElemType::TET4:     t = 10; break;// VTK_TETRA
          case ElemType::TET10:    t = 24; break;// VTK_QUADRATIC_TETRA
          case ElemType::HEX8:     t = 12; break;// VTK_HEXAHEDRON
          case ElemType::HEX20:    t = 25; break;// VTK_QUADRATIC_HEXAHEDRON
          case ElemType::HEX27:    t = 33; break;// VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON
          case ElemType::PRISM6:   t = 13; break;// VTK_WEDGE
          case ElemType::PRISM15:  t = 26; break;// VTK_QUADRATIC_WEDGE
          case ElemType::PRISM18:  t = 32; break;// VTK_BIQUADRATIC_QUADRATIC_WEDGE
          case ElemType::PYRAMID5: t = 14; break;// VTK_PYRAMID
          default : libmesh_error();
        }
      return t;
    }
  // obtain the element (cell) offset based on the number of nodes
  inline
  unsigned int _offset (const Elem * elem) const
    {
      std::vector<unsigned int> c;
      this->_connectivity(elem, c);
      return c.size();
    }

private:
  // a PVD file (and filename) to store time-dependent simulation results for Paraview
  std::ofstream _pvd;
  std::string _filename;
};
//-------------------------------------------------------------------------------------------------
#undef SMALLEST_NUMBER
//-------------------------------------------------------------------------------------------------

#endif // __PARAVIEW_H__

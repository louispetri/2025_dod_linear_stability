#include "refinedvtkwriter.hh"
#include "binaryio.hh"

#include <algorithm>
#include <unistd.h>
#include <iterator>

namespace {

  template<int MODE, class GV, class T>
  struct ModeTraits;

  template<class GV, class T>
  struct ModeTraits<Dune::writeVolume,GV,T>
  {
    enum { dim = GV::dimension, vCodim = dim };
    typedef typename GV::Grid::ctype ctype;

    typedef typename Dune::SubTriangulation::EntityPartInterface<GV> EntityPartInterface;

    typedef typename Dune::SubTriangulation::DefaultEntityPartGeometry<GV> EntityPartGeometry;
    typedef typename Dune::SubTriangulation::DefaultEntityPart<EntityPartGeometry> EntityPart;
    typedef EntityPartInterface* EntityPartPointer;

    typedef std::vector<EntityPartPointer> Container;
    typedef typename Container::iterator Iterator;
    typedef typename EntityPartInterface::Traits::Geometry Geometry;
    typedef typename GV::template Codim <0>::Entity Entity;

    static void fillContainer(Container & geometries, const Entity & e, const T & t)
    {
      geometries.clear();
      t.create_entity_parts(e, geometries);
    }
    static const Geometry & geometry(const Iterator & it)
    {
      return (*it)->geometry();
    }
    static const auto refElem(const Iterator & it)
    {
      return Dune::ReferenceElements<ctype,dim>::general(geometry(it).type());
    }
  };

  template<class GV, class T>
  struct ModeTraits<Dune::writeBoundary,GV,T>
  {
    enum { dim = GV::dimension, vCodim = dim-1 };
    typedef typename GV::Grid::ctype ctype;
    typedef typename GV::template Codim <0>::Entity Entity;
    typedef Dune::SubTriangulation::IntersectionPartInterface<GV> IntersectionPart;
    typedef std::vector<IntersectionPart*> Container;
    typedef typename Container::iterator Iterator;
    typedef typename IntersectionPart::Traits::FaceGeometry Geometry;
    typedef typename IntersectionPart::Traits::FaceGeometry LocalGeometry;

    //! \todo The following method is not general enough for a multi layer
    //!       sub triangulation like MarchingCube33SubTriangulation which
    //!       has domain boundaries on different layers. The method calls
    //!       (*it)->boundary() which means that boundary() of the
    //!       automatically set default layer is used. Thus only
    //!       intersection parts are used for VTK output which have
    //!       outsideDomainIndex_ != insideDomainIndex_ on the default
    //!       layer!
    static void fillContainer(Container & geometries, const Entity & e, const T & t)
    {
      geometries.clear();
      t.create_intersection_parts(e, geometries);
      typename Container::const_iterator end = geometries.end();
      typename Container::iterator it = geometries.begin();
      while (it!=end)
      {
        if (!(*it)->boundary())
          it = geometries.erase(it);
        else
          ++it;
      }
    }

    static const Geometry & geometry(const Iterator & it)
    {
      return (*it)->geometry();
    }

    static const LocalGeometry & localGeometry(const Iterator & it)
    {
      return (*it)->geometryInInside();
    }

    static const auto refElem(const Iterator & it)
    {
      return Dune::ReferenceElements<ctype,dim-1>::general(geometry(it).type());
    }
  };

  template<class GV, class T>
  struct ModeTraits<Dune::writeFaces,GV,T> :
    public ModeTraits<Dune::writeBoundary,GV,T>
  {
    typedef ModeTraits<Dune::writeBoundary,GV,T> Base;
    typedef typename Base::Container Container;
    typedef typename Base::Entity Entity;
    static void fillContainer(Container & geometries, const Entity & e, const T & t)
    {
      geometries.clear();
      t.create_intersection_parts(e, geometries);
    }
  };

  template<class D>
  struct VTKDataPrint
  {
    static void print (std::ostream & s, D & d)
    {
      s << d << " ";
    }
    static void append (std::ofstream & s, D d)
    {
      char* p = reinterpret_cast<char*>(&d);
      s.write(p,sizeof(D));
    }
    static int ncomps(D & a)
    {
      return 1;
    }
    static int sizeoftype(D & a)
    {
      return sizeof(D);
    }
  };

  template<>
  struct VTKDataPrint<unsigned char>
  {
    typedef unsigned char D;
    static void print (std::ostream & s, D & d)
    {
      s << (int)d << " ";
    }
    static void append (std::ostream & s, D d)
    {
      char* p = reinterpret_cast<char*>(&d);
      s.write(p,sizeof(D));
    }
    static int ncomps(D & a)
    {
      return 1;
    }
    static int sizeoftype(D & a)
    {
      return sizeof(D);
    }
  };

  template<class K, int n>
  struct VTKDataPrint< Dune::FieldVector<K,n> >
  {
    static void print (std::ostream & s, Dune::FieldVector<K,n> & x)
    {
      static_assert(n <= 3, "only up to 3 dimensions supported");
      for (int i=0; i<n; i++) s << x[i] << " ";
      for (int i=n; i<3; i++) s << 0 << " ";
    }
    static void append (std::ostream & s, Dune::FieldVector<K,n> & x)
    {
      static_assert(n <= 3, "only up to 3 dimensions supported");

      Dune::FieldVector<K,3> X(0);

      for (int i=0; i<n; i++) X[i] = x[i];

      char* p = reinterpret_cast<char*>(&X);
      s.write(p,sizeof(K)*3);
    }
    static int ncomps(Dune::FieldVector<K,n> & a)
    {
      return 3;
    }
    static int sizeoftype(Dune::FieldVector<K,n> & a)
    {
      return 3*sizeof(K);
    }
  };

  template<>
  struct VTKDataPrint< std::vector<float> >
  {
    static void print (std::ostream & s, std::vector<float> & x)
    {
      assert(x.size() <= 3);

      for (size_t i=0; i<x.size(); i++) s << x[i] << " ";
      for (int i=x.size(); i<ncomps(x); i++) s << 0 << " ";
    }
    static void append (std::ostream & s, std::vector<float> & x)
    {
      assert(x.size() <= 3);

      Dune::FieldVector<float,3> X(0);

      for (size_t i=0; i<x.size(); i++) X[i] = x[i];

      char* p = reinterpret_cast<char*>(&X);
      s.write(p,sizeof(float)*ncomps(x));
    }
    static int ncomps(std::vector<float> & a)
    {
      return a.size() > 1 ? 3 : 1;
    }
    static int sizeoftype(std::vector<float> & a)
    {
      return ncomps(a)*sizeof(float);
    }
  };

}

namespace Dune {

  template < class T, class RT>
  template < class D >
  void
  RefinedVtkWriter<T,RT>::vtkDataArrayEntry(std::ostream & s, std::string t, std::string n, int size, int ncomp, D d)
  {
    // TODO: switch auf VTKOptions::OutputType
    std::string mode = getFormatString();
    if (mode == "ascii")
    {
      s << indent << "<DataArray type=\"" << t << "\" Name=\"" << n
        << "\" NumberOfComponents=\"" << ncomp
        << "\" format=\"" << getFormatString()
        << "\">" << "\n";
      //       for (typename C::iterator it=array.begin(); it!=array.end(); it++)
      //       {
      //         VTKDataPrint<D>::print(s,*it);
      //       }
      assert(false);
      s << "</DataArray>" << "\n";
      return;
    }
    if (mode == "appended")
    {
      s << indent << "<DataArray type=\"" << t << "\" Name=\"" << n
        << "\" NumberOfComponents=\"" << ncomp
        << "\" format=\"" << getFormatString()
        << "\" offset=\"" << offset_ << "\" />" << "\n";
      offset_+=4; // header
      offset_+=size * sizeof(d); // data
      return;
    }
    DUNE_THROW(Dune::IOError, "wrong data mode, must be ascii or appended");
  }

  template < class T, class RT>
  template < class D >
  void
  RefinedVtkWriter<T,RT>::vtkDataArrayAppend(std::ostream & s, std::string name, D d)
  {
    std::string mode = getFormatString();
    if (mode == "ascii") { return; }
    if (mode == "appended")
    {
      // offset
      int offset = datasize_[name] * sizeof(d);
      s << binary(offset); // datasize[name] * sizeof(d));
      // copy old data
      std::ifstream in;
      in.open(tempFileName(name).c_str(), std::ios::binary);
      if (!in.good())
        DUNE_THROW(Dune::IOError, "could not read temporary file "+tempFileName(name));
      s << in.rdbuf();
      in.close();
      return;
    }
    DUNE_THROW(Dune::IOError, "wrong data mode, must be ascii or appended");
  }

  template < class T, class RT>
  template<int MODE>
  void RefinedVtkWriter<T,RT>::countEntities(int &nvertices, int &ncells, int &ncorners)
  {
    /* --- "register" data ostreams --- */
    for (FunctionIterator i=this->vertexdata.begin(); i!=this->vertexdata.end(); ++i)
      datastream_[i->name()];
    datastream_["coordinates"];
    datastream_["connections"];
    datastream_["offsets"];
    datastream_["vtktypes"];

    /* Determine whether all functions provide unfitted vtk function interface */
    bool unfittedVTKFunctionSwitch(true);
    for (std::size_t i = 0; i<vertexDataFunctions_.size(); ++i)
      unfittedVTKFunctionSwitch &=
        bool(dynamic_cast< const Dune::SubTriangulation::UnfittedVTKFunction<T> *>(&(*(vertexDataFunctions_[i]))));


    /* --- create data ostreams --- */
    for (typename std::map<std::string, std::ofstream*>::iterator it = datastream_.begin();
         it != datastream_.end(); ++it)
    {
      std::string name = it->first;
      datastream_[name] = new std::ofstream;
      datastream_[name]->open(tempFileName(name).c_str(), std::ios::binary);
    }

    /* --- Collect Data --- */
    ncells = 0;
    nvertices = 0;
    offset_ = 0;
    std::vector< Dune::FieldVector<float,dim> > coordinates;
    std::vector< int > connections;
    //    std::vector< int > offsets;
    std::vector< unsigned char > vtktypes;
    ElementIterator eit = gridView_.template begin<0, ptype>();
    ElementIterator endit = gridView_.template end<0, ptype>();
    for (; eit != endit; ++eit) {
      typedef ModeTraits<MODE,GridView,Triangulation> ModeTraits;
      static const int vCodim = ModeTraits::vCodim;
      typedef typename ModeTraits::Container Container;
      Container geometries;
      ModeTraits::fillContainer(geometries, *eit, triangulation_);

      ncells += geometries.size();

      // loop over all geometries
      typename Container::iterator geomend = geometries.end();
      typename Container::iterator geomit = geometries.begin();
      for (; geomit != geomend; ++geomit)
      {

        // Determine whether vtk function should be evaluated on this
        // entity part (depending on provided interface)
        bool evaluateHere = true;
        if(unfittedVTKFunctionSwitch) {
          for (std::size_t i = 0; i<vertexDataFunctions_.size(); ++i)
              evaluateHere &= static_cast< const Dune::SubTriangulation::UnfittedVTKFunction<T> *>(&(*(vertexDataFunctions_[i])))
                            ->evaluateOn(*geomit);
        }
        else{
          DUNE_THROW(Dune::Exception, "Called with object which does not inherit from Dune::SubTriangulation::UnfittedVTKFunction");
        }

        // Skip this entity part
        if(!evaluateHere) {
          ncells--;
          continue;
        }

        // get reference element
        auto refElem = ModeTraits::refElem(geomit);

        for (int c=0; c<refElem.size(vCodim); c++)
        {
          Dune::FieldVector<ctype, dim> global =
            ModeTraits::geometry(geomit).
            global(refElem.position(c,vCodim));
          Dune::FieldVector<ctype, dim> local = eit->geometry().local(global);
          // vtkfunctions
          //! \todo The next line has been removed as RefinedVtkWriter also seems
          //!       to work for SubTriangulationVTKWriteMode == writeBoundary. It can be
          //!       completely removed if this impression proves to be right.
          // if (MODE == writeVolume || MODE == writeFaces)
          {
            for (std::size_t i = 0; i<vertexDataFunctions_.size(); ++i)
            {
              int n = 0;
              for (; n < vertexDataFunctions_[i]->ncomps() && n < 3; n++)
              {
                float data = vertexDataFunctions_[i]->evaluate(n, *eit, local);
                // write data to file
                VTKDataPrint<float>::append(*datastream_[vertexDataFunctions_[i]->name()],data);
                datasize_[vertexDataFunctions_[i]->name()]++;
              }
              if(n > 1)
                for(; n < 3; ++n)
                {
                  // write padding
                  VTKDataPrint<float>::append(*datastream_[vertexDataFunctions_[i]->name()],float(0));
                  datasize_[vertexDataFunctions_[i]->name()]++;
                }
            }
          }
          // global position
          for (int d=0; d<dim; d++)
            VTKDataPrint<float>::append(*datastream_["coordinates"], (float)global[d]);
          if(dim<3) for (int d=3; d>dim; --d) VTKDataPrint<float>::append(*datastream_["coordinates"], 0.0);
          datasize_["coordinates"] += 3;
          // connections
          VTKDataPrint<int>::append(*datastream_["connections"],
                                    nvertices + vtkRenumber(ModeTraits::geometry(geomit).type(), c));
          datasize_["connections"]++;
        }
        nvertices += refElem.size(vCodim);

        // offsets
        VTKDataPrint<int>::append(*datastream_["offsets"], nvertices);
        datasize_["offsets"]++;
        // vtk geometry types
        VTKDataPrint<unsigned char>::append(*datastream_["vtktypes"], vtkType(ModeTraits::geometry(geomit).type()));
        datasize_["vtktypes"]++;
      }
    }

    /* --- close data streams --- */
    for (typename std::map<std::string, std::ofstream*>::iterator it = datastream_.begin();
         it != datastream_.end(); ++it)
    {
      std::string name = it->first;
      datastream_[name]->close();
      delete datastream_[name];
    }

    /* --- vertices == corners --- */
    ncorners = nvertices;
  }

  template < class T, class RT>
  void RefinedVtkWriter<T,RT>::writeCellData (VTK::VTUWriter& writer)
  {
    return;
  }

  template < class T, class RT>
  void RefinedVtkWriter<T,RT>::writeVertexData (VTK::VTUWriter& writer)
  {
    if(writer.phase == VTK::VTUWriter::appended) {
      writeAppendedData(writer);
      return;
    }

    // --- Store Point Data ---
    std::string scalars = "";

    for (std::size_t i = 0; i<cellDataFunctions_.size(); ++i)
      if (cellDataFunctions_[i]->ncomps()==1)
      {
        scalars = cellDataFunctions_[i]->name();
        break;
      }
    std::string vectors = "";
    for (std::size_t i = 0; i<cellDataFunctions_.size(); ++i)
      if (cellDataFunctions_[i]->ncomps()>1)
      {
        vectors = cellDataFunctions_[i]->name();
        break;
      }

    writer.beginPointData(scalars, vectors);

    ++indent;
    for (std::size_t i = 0; i<vertexDataFunctions_.size(); ++i)
    {
      vtkDataArrayEntry(writer.stream,
                        "Float32", vertexDataFunctions_[i]->name(), datasize_[vertexDataFunctions_[i]->name()], (vertexDataFunctions_[i]->ncomps()>1 ? 3 : 1), float());
    }

    --indent;

    writer.endPointData();
  }

  template < class T, class RT>
  void RefinedVtkWriter<T,RT>::writeGridPoints (VTK::VTUWriter& writer)
  {
    if(writer.phase == VTK::VTUWriter::appended) {
      return;
    }

    std::ostream & s = writer.stream;
    // --- Store Coordinates ---
    writer.beginPoints();
    ++indent;
    vtkDataArrayEntry(s, "Float32", "Coordinates", datasize_["coordinates"], 3, float());
    --indent;
    writer.endPoints();
  }

  template < class T, class RT>
  void RefinedVtkWriter<T,RT>::writeGridCells (VTK::VTUWriter& writer)
  {
    if(writer.phase == VTK::VTUWriter::appended) {
      return;
    }

    std::ostream & s = writer.stream;
    // --- Store Grid Cells ---
    writer.beginCells();
    ++indent;
    vtkDataArrayEntry(s, "Int32", "connectivity", datasize_["connections"], 1, int());
    vtkDataArrayEntry(s, "Int32", "offsets", datasize_["offsets"], 1, int());
    vtkDataArrayEntry(s, "UInt8", "types", datasize_["vtktypes"], 1, (unsigned char)(0));
    --indent;
    writer.endCells();
  }

  template < class T, class RT>
  void RefinedVtkWriter<T,RT>::writeAppendedData (VTK::VTUWriter& writer)
  {
    std::ostream & s = writer.stream;
    /* --- copy appended data from tmp files --- */
    //s << indent << "<AppendedData encoding=\"raw\">" << "\n";
    // ++indent;
    // s << indent << "_"; // indicates start of binary data
    for (FunctionIterator i=this->vertexdata.begin(); i!=this->vertexdata.end(); ++i)
    {
      vtkDataArrayAppend(s, i->name(), float());
    }
    vtkDataArrayAppend(s, "coordinates", float());
    vtkDataArrayAppend(s, "connections", int());
    vtkDataArrayAppend(s, "offsets", int());
    vtkDataArrayAppend(s, "vtktypes", (unsigned char)(0));
    s << "\n";
    --indent;
    //s << indent << "</AppendedData>" << "\n";

    /* --- delete tmp ostream files --- */
    for (typename std::map<std::string, std::ofstream*>::iterator it = datastream_.begin();
         it != datastream_.end(); ++it)
    {
      std::string name = it->first;
      unlink(tempFileName(name).c_str());
    }
  }

} // end namespace Dune

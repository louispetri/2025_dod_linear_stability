// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef REFINEDVTKWRITER_HH
#define REFINEDVTKWRITER_HH

#include <vector>
#include <string>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/subtriangulation/common/entitypart.hh>
#include <dune/subtriangulation/common/intersectionpart.hh>
#include <dune/subtriangulation/io/vtkfunction.hh>

namespace Dune {

  enum UDGVTKWriteMode {
    writeVolume   = 1,
    writeFaces    = 2,
    writeBoundary = 4
  };


  template<class T, class RT>
  class RefinedVtkWriter : public VTKWriter<typename T::LevelSetGridView>
  {
  public:
    using GridView = typename T::LevelSetGridView;
  private:
    typedef VTKWriter<GridView> Base;
    typedef typename Base::VTKFunction VTKFunction;
    typedef typename Base::FunctionIterator FunctionIterator;

    Indent indent;
    // using Base::indent;
    // using Base::indentDown;
    // using Base::indentUp;
    using Base::getFormatString;
  public:
    typedef typename std::shared_ptr<typename Base::VTKFunction> VTKFunctionPtr;
    typedef T Triangulation;
    typedef typename GridView::Grid Grid;
    enum
    { dimension = Grid::dimension /* !<dimension of the grid we are using */  };
    enum
    { dimensionworld = Grid::dimensionworld /* !<dimension of the world */  };

    RefinedVtkWriter (const GridView& g, const Triangulation& t) :
      VTKWriter<GridView>(g),
      gridView_ (g), triangulation_ (t),
      domainTag_(0), writeMode_(Dune::writeVolume)
    {};

    //! destructor
    virtual ~RefinedVtkWriter () {};

    /**
     * @brief write output; interface might change later
     * @param name The name of the file to write to.
     * @param ot The output type for the file.
     */
    template <int MODE>
    std::string write(std::string name, std::string mode = "ascii", const int domain_tag=255)
    {
      VTK::OutputType type = VTK::ascii;
      if (mode == "appended") type = VTK::appendedraw;
      return write( name, type, domain_tag, UDGVTKWriteMode(MODE) );
    }

    std::string write(std::string name, std::string mode = "ascii", const int domain_tag=255)
    {
      return write<Dune::writeVolume>( name, mode, domain_tag );
    }

    std::string write (std::string name, VTK::OutputType type = VTK::ascii,
                       const int domain_tag=255, UDGVTKWriteMode mode = Dune::writeVolume)
    {
      // strip trailing ".vtu" if necessary
      if (name.size() > 4 && name.substr(name.size()-4,name.size()-1) == ".vtu")
        name = name.substr(0, name.size()-4);
      writeMode_ = mode;
      baseName_ = name;
      domainTag_ = domain_tag;
      return VTKWriter<GridView>::write( name, type );
    }

    std::string pwrite (std::string name,  const char* path, const char* extendpath,
                        VTK::OutputType type = VTK::ascii,
                        const int domain_tag=255, UDGVTKWriteMode mode = Dune::writeVolume)
    {
      // strip trailing ".vtu" if necessary
      if (name.substr(name.size()-4,name.size()-1) == ".vtu")
        name = name.substr(0, name.size()-4);
      writeMode_ = mode;
      baseName_ = name;
      domainTag_ = domain_tag;
      return VTKWriter<GridView>::pwrite( name.c_str(), path, extendpath, type );
    }

    void addCellData (VTKFunction* p)
    {
      VTKFunctionPtr ptr(p);
      Base::addCellData(ptr);
      cellDataFunctions_.push_back(ptr);
    }

    void addVertexData (VTKFunction* p)
    {
      VTKFunctionPtr ptr(p);
      Base::addVertexData(ptr);
      vertexDataFunctions_.push_back(ptr);
    }

    void addCellData (VTKFunctionPtr p)
    {
      Base::addCellData(p);
      cellDataFunctions_.push_back(p);
    }

    void addVertexData (VTKFunctionPtr p)
    {
      Base::addVertexData(p);
      vertexDataFunctions_.push_back(p);
    }

  protected:
    //! count the vertices, cells and corners
    virtual void countEntities(int &nvertices, int &ncells, int &ncorners)
    {
      switch (writeMode_) {
      case writeVolume :
        return countEntities<writeVolume>(nvertices, ncells, ncorners);
      case writeFaces :
        return countEntities<writeFaces>(nvertices, ncells, ncorners);
      case writeBoundary :
        return countEntities<writeBoundary>(nvertices, ncells, ncorners);
      };
    }

    //! write cell data
    virtual void writeCellData (VTK::VTUWriter& writer);

    //! write vertex data
    virtual void writeVertexData (VTK::VTUWriter& writer);

    //! write the positions of vertices
    virtual void writeGridPoints (VTK::VTUWriter& writer);

    //! write the connectivity array
    virtual void writeGridCells (VTK::VTUWriter& writer);

    //! write the appended data sections
    virtual void writeAppendedData (VTK::VTUWriter& writer);

  private:

    // dimensions
    enum
    { dim = Grid::dimension /* !< dimension of the grid we are using */  };
    enum
    { dimw = Grid::dimensionworld /* !< dimension of the world */  };
    static const Dune::PartitionIteratorType ptype = Dune::Interior_Partition;
    //
    typedef typename Grid::ctype ctype;
    // Iterators
    typedef typename GridView::template Codim < 0 >::template Partition<ptype>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    // Geometries
    typedef Dune::Cuboid < dim, ctype > Geometry;
    typedef typename Triangulation::EntityPart EntityPart;
    typedef typename Triangulation::IntersectionPart IntersectionPart;

    enum VTKGeometryType
    {
      vtkLine = 3,
      vtkTriangle = 5,
      vtkQuadrilateral = 9,
      vtkTetrahedron = 10,
      vtkHexahedron = 12,
      vtkPrism = 13,
      vtkPyramid = 14,
      vtkQTriangle = 22,
      vtkQQuadrilateral = 23,
      vtkQTetrahedron = 24,
      vtkQHexahedron = 25
    };

    template <class D>
    void vtkDataArrayEntry(std::ostream & s, std::string t, std::string n, int size, int ncomp, D);

    template <class D>
    void vtkDataArrayAppend(std::ostream & s, std::string name, D);

    //! count the vertices, cells and corners
    //! and prepare data of type MODE
    template<int MODE>
    void countEntities(int &nvertices, int &ncells, int &ncorners);

    std::string tempFileName(std::string name) const
    {
      char fullname[ 8192 ];
      sprintf(fullname,"%s:s%04d:p%04d:%s.tmp",
              baseName_.c_str(), gridView_.comm().size(), gridView_.comm().rank(), name.c_str());
      // std::cout << fullname << std::endl;
      return fullname;
    };

    //! mapping from GeometryType to VTKGeometryType
    static VTKGeometryType vtkType(const Dune::GeometryType & t)
    {
      if (t.isLine())
        return vtkLine;
      if (t.isTriangle())
        return vtkTriangle;
      if (t.isQuadrilateral())
        return vtkQuadrilateral;
      if (t.isTetrahedron())
        return vtkTetrahedron;
      if (t.isPyramid())
        return vtkPyramid;
      if (t.isPrism())
        return vtkPrism;
      if (t.isHexahedron())
        return vtkHexahedron;
      DUNE_THROW(Dune::IOError,"VTKWriter: unsupported GeometryType "
                 << t <<std::endl);
    }

    // renumber VTK -> Dune
    static int vtkRenumber (const Dune::GeometryType & t, int i)
    {
      static const int quadRenumbering[4] = {0,1,3,2};
      static const int cubeRenumbering[8] = {0,1,3,2,4,5,7,6};
      static const int prismRenumbering[6] = {0,2,1,3,5,4};
      //std::cout << "(" << vtkType(t) << " " << i << ")";
      switch (vtkType(t))
      {
      case vtkQuadrilateral :
        return quadRenumbering[i];
      case vtkHexahedron :
        return cubeRenumbering[i];
      case vtkPrism :
        return prismRenumbering[i];
      default :
        return i;
      }
    }

    /* TODO:
       - make filename accessible from VTKFunction
       - make gridView accessible from VTKFunction
     */

    GridView gridView_;
    const Triangulation & triangulation_;

    //! which subdomains shall be written?
    mutable int domainTag_;
    //! which parts of the domain shall be written?
    mutable UDGVTKWriteMode writeMode_;
    //! baseName of the vtu file
    mutable std::string baseName_;
    //! size of data
    mutable std::map<std::string, size_t> datasize_;
    //! streams for temporary data
    mutable std::map<std::string, std::ofstream*> datastream_;
    //! offsetcounter
    mutable int offset_;

    // store ptrs to the vtkfunctions in order to access the udg interface
    std::vector<VTKFunctionPtr> vertexDataFunctions_;
    std::vector<VTKFunctionPtr> cellDataFunctions_;

  };

}

#include "refinedvtkwriter.cc"

#endif // REFINEDVTKWRITER_HH

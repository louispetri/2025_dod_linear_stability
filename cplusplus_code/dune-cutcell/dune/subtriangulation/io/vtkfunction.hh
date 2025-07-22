// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SubTriangulation_IO_VTKFUNCTION_HH
#define DUNE_SubTriangulation_IO_VTKFUNCTION_HH

#include <dune/grid/io/file/vtk.hh>

#include <dune/subtriangulation/common/entitypart.hh>
#include <dune/subtriangulation/common/intersectionpart.hh>

namespace Dune {

  namespace SubTriangulation {

    // Base type for VTK functions for refined VTK writer
    template<class ST_>
    class UnfittedVTKFunction : public Dune::VTKFunction<typename ST_::LevelSetGridView>
    {
    public:
      using ST = ST_;

      using IntersectionPartPointer = typename ST::IntersectionPartPointer;
      using EntityPartPointer =  typename ST::EntityPartPointer;

      virtual bool evaluateOn(const EntityPartPointer& part) const = 0;
      virtual bool evaluateOn(const IntersectionPartPointer& part) const = 0;
    };

    template<typename GV>
    class HostCellIndexUnfittedVTKGridFunction
      : public UnfittedVTKFunction<GV>
    {
      typedef typename GV::ctype DF;
      enum {n=GV::dimension};
      typedef typename GV::template Codim<0>::Entity Entity;
      typedef UnfittedVTKFunction<GV> Base;
    public:

      HostCellIndexUnfittedVTKGridFunction(const GV & _gv, const int layer = 0)
        : gv_(_gv), layer_(layer) {}

      virtual int ncomps () const
      {
        return 1;
      }

      virtual bool evaluateOn( const typename Base::EntityPartPointer & part) const
      {
        (part.get())->setLayer(layer_);
        hostIndex = gv_.indexSet().index(part->boundingBox().entity());
        domainIndex = part->domainIndex();
        (part.get())->resetLayer();
        return true;
      }

      virtual bool evaluateOn( const typename Base::IntersectionPartPointer & part) const
      {
        return false;
      }

      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DF,n>& xi) const
      {
        return hostIndex;
      }

      virtual std::string name () const
      {
        std::stringstream sname;
        sname << "HostIndex";

        if(layer_ >=0)
          sname << "Layer" << layer_;
        return sname.str();
      }

    private:
      GV gv_;
      const int layer_;
      mutable int hostIndex;
      mutable int domainIndex;
    };

    template<typename GV>
    class DomainIndexUnfittedVTKGridFunction
      : public UnfittedVTKFunction<GV>
    {
      typedef typename GV::ctype DF;
      enum {n=GV::dimension};
      typedef typename GV::template Codim<0>::Entity Entity;
      typedef UnfittedVTKFunction<GV> Base;
    public:

      DomainIndexUnfittedVTKGridFunction
        (const GV & _gv, const bool insideDomain = true, const int layer = 0)
        : gv_(_gv), insideDomain_(insideDomain), layer_(layer) {}

      virtual int ncomps () const
      {
        return 1;
      }

      virtual bool evaluateOn( const typename Base::EntityPartPointer & part) const
      {
        (part.get())->setLayer(layer_);
        domainIndex = part->domainIndex();
        (part.get())->resetLayer();
        return true;
      }

      virtual bool evaluateOn( const typename Base::IntersectionPartPointer & part) const
      {
        (part.get())->setLayer(layer_);
        if(insideDomain_)
          domainIndex = part->insideDomainIndex();
        else
          domainIndex = part->outsideDomainIndex();
        (part.get())->resetLayer();
        return true;
      }

      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DF,n>& xi) const
      {
        return domainIndex;
      }

      virtual std::string name () const
      {
        std::stringstream sname;
        if(insideDomain_)
          sname << "InsideDomainIndex";
        else
          sname << "OutsideDomainIndex";

        if(layer_ >=0)
          sname << "Layer" << layer_;
        return sname.str();
      }

    private:
      const GV & gv_;
      bool insideDomain_;
      const int layer_;
      mutable int domainIndex;
    };

  }
}
#endif

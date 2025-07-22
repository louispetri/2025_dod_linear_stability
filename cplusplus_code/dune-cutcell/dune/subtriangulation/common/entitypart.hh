// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef ENTITY_PART_HH
#define ENTITY_PART_HH

#include <cstdarg>
#include <list>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/utility/typefromvertexcount.hh>
#include <dune/subtriangulation/common/cuboid.hh>

namespace Dune {
  namespace SubTriangulation {

    /** \brief Traits class for entity parts */
    template< typename GV >
    struct EntityPartTraits
    {

      //! Grid view containing this part
      typedef GV GridView;
      //! Field type of domain coordinates
      typedef typename GridView::ctype DomainField;
      //! Local coordinate
      typedef Dune::FieldVector<DomainField,GV::dimension>
      LocalCoordinate;
      //! Global coordinate
      typedef Dune::FieldVector<DomainField,GV::dimension>
      GlobalCoordinate;

      //! Dimension of computational domain
      enum { dimension=GV::dimension };

      //! Geometry for entity parts
      using Geometry = Dune::CachedMultiLinearGeometry<DomainField,dimension,dimension>;

      //! Bounding box for shape functions
      typedef Dune::GridCuboid<typename GridView::Grid> BoundingBox;

      typedef typename GV::template Codim<0>::Entity Entity;
    };

    /** \brief Interface for the geometry engine of an entity part. */
    template <typename GV>
    class EntityPartGeometryInterface
    {
    public:
      //! Traits class
      typedef EntityPartTraits<GV> Traits;

      //! Access the global geometry
      virtual const typename Traits::Geometry& geometry () const = 0;

      //! Access the entity part's home entity in level set mesh or in
      //! fundamental mesh if the entity part was constructed from a whole
      //! fundamental mesh cell
      virtual const typename Traits::Entity & entity() const = 0;

      //! Access the bounding box
      virtual const typename Traits::BoundingBox & boundingBox() const = 0;
    };

    /** \brief Interface for an entity part. */
    template <typename GV>
    class EntityPartInterface : public virtual EntityPartGeometryInterface<GV>
    {
    public:
      //! Traits class
      typedef EntityPartTraits<GV> Traits;

      //! Shared pointer on entity part
      typedef std::shared_ptr<EntityPartInterface> EntityPartPointer;
      //! List for entity parts
      typedef std::list< EntityPartPointer > EntityPartList;

      virtual ~EntityPartInterface() {}

      /** \brief Domain index of this entity part */
      virtual int domainIndex () const = 0;
    };


    /** \brief Base class for implementations of
        IntersectionPartGeometryInterface .

        This call represents intersections with the generic DUNE
        geometries.

        \tparam GV The grid view for which this part was constructed
     */
    template <class GV>
    class DefaultEntityPartGeometry :
      virtual public EntityPartGeometryInterface<GV>
    {
    public:

      //! Traits class
      typedef EntityPartTraits<GV> Traits;

      //! Constructor
      DefaultEntityPartGeometry
        (const typename Traits::Geometry & geometry,
        const typename Traits::Entity & entity,
        const typename Traits::BoundingBox & bbox) :
        geometry_(geometry),
        entity_(entity),
        bbox_(bbox)
      {}

      //! Access the global geometry
      const typename Traits::Geometry& geometry () const
      { return geometry_; }

      //! Access the entity part's home entity in level set mesh or in
      //! fundamental mesh if the entity part was constructed from a whole
      //! fundamental mesh cell
      const typename Traits::Entity & entity() const
      { return entity_; }

      //! Access the bounding box
      const typename Traits::BoundingBox & boundingBox() const
      { return bbox_; }

    private:
      const typename Traits::Geometry geometry_;
      const typename Traits::Entity entity_;
      const typename Traits::BoundingBox bbox_;
    };

    /** \brief Default implementation for entity part.

        \tparam EPG Implementation of EntityPartGeometryInterface.
     */
    template <class EPG>
    class DefaultEntityPart :
      public EPG,
      public EntityPartInterface<typename EPG::Traits::GridView>
    {
    public:
      //! Geometry engine
      typedef EPG EntityPartGeometry;
      //! Traits class
      typedef typename EPG::Traits Traits;

      //! Constructor
      DefaultEntityPart
        (const EPG & geometryEngine,
        const int domainIndex) :
        EPG(geometryEngine),
        domainIndex_(domainIndex)
      {}

      /** \brief Domain index of this entity part*/
      int domainIndex () const {
        return domainIndex_;
      }

    private:
      int domainIndex_;
    };


  }
}

#endif

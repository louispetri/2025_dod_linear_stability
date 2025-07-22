// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef INTERSECTION_PART_HH
#define INTERSECTION_PART_HH

#include <list>
#include <stdarg.h>

#include <dune/geometry/utility/typefromvertexcount.hh>
#include <dune/subtriangulation/common/cuboid.hh>

namespace Dune {
  namespace SubTriangulation {

    /** \brief Traits class for intersection parts */
    template< typename GV >
    struct IntersectionPartTraits
    {

      //! Grid view containing this part
      typedef GV GridView;
      //! Field type of domain coordinates
      typedef typename GridView::ctype DomainField;
      //! Local coordinate
      typedef Dune::FieldVector<DomainField,GV::dimension-1>
      LocalCoordinate;
      //! Global coordinate
      typedef Dune::FieldVector<DomainField,GV::dimension>
      GlobalCoordinate;

      //! Dimension of computational domain
      enum { dimension=GV::dimension };

      //! Geometry for intersection parts
      typedef CachedMultiLinearGeometry<DomainField,dimension-1,dimension> FaceGeometry;

      //! Geometry for entity parts
      typedef CachedMultiLinearGeometry<DomainField,dimension,dimension> Geometry;

      //! Bounding box for shape functions
      typedef Dune::GridCuboid<typename GridView::Grid> BoundingBox;

      //! Entity for cells in grid view
      typedef typename GridView::template Codim<0>::Entity Entity;
    };


    /** \brief Interface for the geometry engine of an intersection
        part.
     */
    template <typename GV>
    class IntersectionPartGeometryInterface
    {
    public:
      //! Traits class
      typedef IntersectionPartTraits<GV> Traits;

      //! Access the host entity on the outside of intersection part
      //! \note This is not comparable to method entity() of entity part
      virtual const typename Traits::Entity & outside() const = 0;
      //! Access the host entity on the inside of intersection part
      //! \note This is not comparable to method entity() of entity part
      virtual const typename Traits::Entity & inside() const = 0;

      /** \brief Access the geometry mapping from the local
         intersection domain to the local domain of the bounding box on
         inside of this intersection part. */
      virtual const typename Traits::FaceGeometry& geometryInInside () const = 0;
      /** \brief Access the geometry mapping from the local
         intersection domain to the local domain of the bounding box on
         outside of this intersection part. */
      virtual const typename Traits::FaceGeometry& geometryInOutside () const = 0;
      //! Access the global geometry
      virtual const typename Traits::FaceGeometry& geometry () const = 0;

      //! Access the bounding box on the inside
      virtual const typename Traits::BoundingBox & boundingBoxInside() const = 0;
      //! Access the bounding box on the outside
      virtual const typename Traits::BoundingBox & boundingBoxOutside() const = 0;

      //! Return unit outer normal at local face coordinate
      virtual typename Traits::GlobalCoordinate unitOuterNormal
        (const typename Traits::LocalCoordinate x) const = 0;
    };

    template <typename GV>
    class IntersectionPartInterface :
      public virtual IntersectionPartGeometryInterface<GV>
    {
    public:
      //! Traits class
      typedef IntersectionPartTraits<GV> Traits;

      //! Shared pointer on intersection part
      typedef std::shared_ptr<IntersectionPartInterface>
      IntersectionPartPointer;
      //! List for intersection parts
      typedef std::list< IntersectionPartPointer >
      IntersectionPartList;

      virtual ~IntersectionPartInterface() {}

      /** \brief Domain index on the inside of this intersection part */
      virtual int insideDomainIndex () const = 0;
      /** \brief Domain index on the outside of this intersection part */
      virtual int outsideDomainIndex () const = 0;

      /** \brief Indicates whether this intersection part is on a
          domain boundary. */
      virtual bool boundary() const = 0;

      /** \brief Indicates whether this intersection part has a domain
          cell patch neighbor on the outside. */
      virtual bool neighbor() const = 0;
    };


    /** \brief Base class for implementations of
        IntersectionPartGeometryInterface .

        This call represents intersections with the generic DUNE
        geometries.

        \tparam GV The grid view for which this part was constructed
     */
    template <class GV>
    class IntersectionPartGeometryBase :
      virtual public IntersectionPartGeometryInterface<GV>
    {
    public:

      //! Traits class
      typedef IntersectionPartTraits<GV> Traits;

      //! Constructor
      IntersectionPartGeometryBase
        (const typename Traits::FaceGeometry & faceGeometry,
        const typename Traits::BoundingBox & insideBBox,
        const typename Traits::BoundingBox & outsideBBox) :
        geometryInInside_(computeFaceGeometryInBoundingBox(faceGeometry,insideBBox)),
        geometryInOutside_(computeFaceGeometryInBoundingBox(faceGeometry,outsideBBox)),
        geometry_(faceGeometry),
        insideBBox_(insideBBox), outsideBBox_(outsideBBox)
      {}

      //! Access the host entity on the outside of intersection part
      //! \note This is not comparable to method entity() of entity part
      const typename Traits::Entity & outside() const
      { return outsideBBox_.entity(); }
      //! Access the host entity on the inside of intersection part
      //! \note This is not comparable to method entity() of entity part
      const typename Traits::Entity & inside() const
      { return insideBBox_.entity(); }

      /** \brief Access the geometry mapping from the local
         intersection domain to the local domain of the bounding box on
         inside of this intersection part. */
      const typename Traits::FaceGeometry& geometryInInside () const
      { return geometryInInside_; }
      /** \brief Access the geometry mapping from the local
         intersection domain to the local domain of the bounding box on
         outside of this intersection part. */
      const typename Traits::FaceGeometry& geometryInOutside () const
      { return geometryInOutside_; }
      //! Access the global geometry
      const typename Traits::FaceGeometry& geometry () const
      { return geometry_; }

      //! Access the bounding box on the inside
      const typename Traits::BoundingBox & boundingBoxInside() const
      { return insideBBox_; }
      //! Access the bounding box on the outside
      const typename Traits::BoundingBox & boundingBoxOutside() const
      { return outsideBBox_; }

    private:

      // Return a mapping that maps the local domain of the global
      // face geometry into the local domain of the bounding box.
      static typename Traits::FaceGeometry computeFaceGeometryInBoundingBox
        (const typename Traits::FaceGeometry & geometry,
        const typename Traits::BoundingBox & bbox)
      {
        std::vector<typename Traits::GlobalCoordinate>
        cornersInBBox(geometry.corners());
        for(int c=0; c<geometry.corners(); ++c)
          cornersInBBox[c] = bbox.local(geometry.corner(c));
        return typename Traits::FaceGeometry(geometry.type(),cornersInBBox);
      }

      const typename Traits::FaceGeometry geometryInInside_;
      const typename Traits::FaceGeometry geometryInOutside_;
      const typename Traits::FaceGeometry geometry_;
      const typename Traits::BoundingBox insideBBox_;
      const typename Traits::BoundingBox outsideBBox_;
    };


    /** \brief Implementation of IntersectionPartGeometryInterface
        interface.

        Intersection is represented by generic DUNE geoemtries and the
        unit outer normal vector of the intersection is assumed to be
        constant.

        \tparam GV The grid view for which this part was constructed
     */
    template <class GV>
    class AffineIntersectionPartGeometry :
      public IntersectionPartGeometryBase<GV>
    {
      // Base class
      typedef IntersectionPartGeometryBase<GV> Base;
    public:

      //! Traits class
      typedef typename Base::Traits Traits;

      //! Constructor
      AffineIntersectionPartGeometry
        (const typename Traits::FaceGeometry & faceGeometry,
        const typename Traits::BoundingBox & insideBBox,
        const typename Traits::BoundingBox & outsideBBox,
        const typename Traits::GlobalCoordinate & unitOuterNormal) :
        Base(faceGeometry,insideBBox,outsideBBox),
        unitOuterNormal_(unitOuterNormal)
      {}

      //! Return unit outer normal at local face coordinate
      typename Traits::GlobalCoordinate unitOuterNormal
        (const typename Traits::LocalCoordinate) const
      { return unitOuterNormal_; }

    private:
      const typename Traits::GlobalCoordinate unitOuterNormal_;
    };


    /** \brief Default implementation for intersection part.

        \tparam IPG Implementation of IntersectionPartGeometryInterface.
     */
    template <class IPG>
    class DefaultIntersectionPart :
      public IPG,
      public IntersectionPartInterface<typename IPG::Traits::GridView>
    {
    public:
      //! Geometry engine
      typedef IPG IntersectionPartGeometry;
      //! Traits class
      typedef typename IPG::Traits Traits;

      //! Constructor
      DefaultIntersectionPart
        (const IPG & geometryEngine,
        const int insideDomainIndex, const int outsideDomainIndex) :
        IPG(geometryEngine),
        insideDomainIndex_(insideDomainIndex), outsideDomainIndex_(outsideDomainIndex)
      {}

      /** \brief Domain index on the inside of this intersection part */
      int insideDomainIndex () const {
        return insideDomainIndex_;
      }
      /** \brief Domain index on the outside of this intersection part */
      int outsideDomainIndex () const {
        return outsideDomainIndex_;
      }

      /** \brief Indicates whether this intersection part is on a
          domain boundary. */
      bool boundary() const {
        return outsideDomainIndex_ != insideDomainIndex_;
      }

      /** \brief Indicates whether this intersection part has a domain
          cell patch neighbor on the outside. */
      bool neighbor() const {
        return outsideDomainIndex_ != -1;
      }

    private:
      int insideDomainIndex_;
      int outsideDomainIndex_;
    };

  }
}
#endif

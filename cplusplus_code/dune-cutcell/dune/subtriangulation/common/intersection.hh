// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_SUBTRIANGULATION_COMMON_INTERSECTION_HH
#define DUNE_SUBTRIANGULATION_COMMON_INTERSECTION_HH

#include <dune/subtriangulation/common/intersectionpart.hh>

namespace Dune {
  namespace SubTriangulation {

    template <typename GridImp>
    class Intersection
    {
    public:
      //! Domain tag
      using DomainTag = std::vector<int>;
      //! Field type of domain coordinates
      using ctype = typename GridImp::ctype;
      //! List of intersections
      using IntersectionList = std::list<Intersection>;
      //! List of intersection parts
      using IntersectionPartList = std::vector<IntersectionPartInterface<GridImp>*>;

      //! Dimension of the coordinates
      enum { coorddimension=GridImp::dimension };

      //! Codim-0-Entity for this grid
      using Entity = typename GridImp::template Codim<0>::Entity;


      //! Constructor
      Intersection(const IntersectionPartList& intersectionParts) :
        intersectionParts_(intersectionParts)
      {
        // verify that the intersection part list is not empty
        if (intersectionParts.empty())
          DUNE_THROW(Dune::Exception,"List of intersection parts should not be empty!");

        // compute the volume of this intersection by summing up the volumina
        // of the intersection parts belonging to it
        volume_ = 0.;
        for (const auto& ip : intersectionParts_)
        {
          volume_ += ip->geometry().volume();
        }
      }


      //! Access the list of intersection parts that form this intersection
      const IntersectionPartList& intersectionParts() const
      {
        return intersectionParts_;
      }

      /*! \brief Sign vector indicating the placement of this
       *         intersection relative to the interfaces separating
       *         the sub triangulation's domains.
       *
       * For each interface of the sub triangulation, this sign
       * vector has one entry holding an element of {-1,0,1},
       * where 0 indicates that this intersection is directly
       * on the corresponding interface.
       */
      const DomainTag& newDomainTag() const
      {
        return (*intersectionParts_.front()).domainTag();
      }

      //! Domain index on the inside of this intersection
      int insideDomainIndex() const
      {
        return (*intersectionParts_.front()).insideDomainIndex();
      }

      //! Domain index on the outside of this intersection
      int outsideDomainIndex() const
      {
        return (*intersectionParts_.front()).outsideDomainIndex();
      }

      /*! \brief Indicates whether this intersection is on a
       *         domain boundary
       */
      bool boundary() const
      {
        return outsideDomainIndex() != insideDomainIndex();
      }

      /*! \brief Indicates whether this intersection has a domain
       *         cell patch neighbor on the outside
       */
      bool neighbor() const
      {
        return outsideDomainIndex() != -1;
      }

      //! \brief Returns the volume of this intersection
      ctype volume() const
      {
        return volume_;
      }

      /*! \brief Returns the element on the inside of this intersection
       *         (determined by the first intersection part belonging to this intersection)
       */
      const Entity& inside() const
      {
        return (*intersectionParts_.front()).inside();
      }

      /*! \brief Returns the element on the outside of this intersection
       *         (determined by the first intersection part belonging to this intersection)
       */
      const Entity& outside() const
      {
        return (*intersectionParts_.front()).outside();
      }

    private:
      IntersectionPartList intersectionParts_;
      ctype volume_;
    };

  } //namespace SubTriangulation
} //namespace Dune

#endif

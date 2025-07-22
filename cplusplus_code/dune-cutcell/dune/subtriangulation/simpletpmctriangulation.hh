#ifndef DUNE_SUBTRIANGULATION_SIMPLETPMCTRIANGULATION_HH
#define DUNE_SUBTRIANGULATION_SIMPLETPMCTRIANGULATION_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/utility/typefromvertexcount.hh>
#include <dune/subtriangulation/common/entitypart.hh>
#include <dune/subtriangulation/common/intersection.hh>
#include <dune/subtriangulation/common/intersectionpart.hh>
#include <dune/subtriangulation/simpletpmctriangulation/domainconfiguration.hh>
#include <dune/subtriangulation/simpletpmctriangulation/cutcellinformation.hh>
#include <dune/subtriangulation/simpletpmctriangulation/localsubtriangulation.hh>

namespace Dune
{
  namespace SubTriangulation
  {
    namespace SimpleTpmcTriangulationDetail
    {
      template <class GOut, class G>
      GOut createIdentityGeometry(const G& geometry)
      {
        using ctype = typename G::ctype;
        const auto& referenceElement
            = ReferenceElements<ctype, G::mydimension>::general(geometry.type());
        typename ReservedStorageMultiLinearGeometryTraits<ctype>::template CornerStorage<GOut::mydimension, GOut::coorddimension>::Type corners;

        for (int i = 0; i < referenceElement.size(G::mydimension); ++i) {
          corners.push_back(referenceElement.position(i, G::mydimension));
        }
        return GOut(geometry.type(), corners);
      }
    }

    template <class GV, class LGV>
    class SimpleTpmcTriangulation
    {
    public:
      using GridView = GV;
      using Grid = typename GV::Grid;
      enum { dim = GV::dimension };
      using ctype = typename GV::ctype;
      using Entity = typename GV::template Codim<0>::Entity;
      using Geometry = Dune::CachedMultiLinearGeometry<ctype, dim, dim>;
      using EntityPart = Dune::SubTriangulation::EntityPartInterface<GridView>;
      using EntityPartGeometry = Dune::SubTriangulation::DefaultEntityPartGeometry<GV>;
      using EntityPartImpl = Dune::SubTriangulation::DefaultEntityPart<EntityPartGeometry>;
      using EntityPartInterface = Dune::SubTriangulation::EntityPartInterface<GV>;
      using EntityPartPointer = EntityPartInterface*;
      using EntityPartList = std::vector<EntityPartPointer>;
      using IntersectionPart = Dune::SubTriangulation::IntersectionPartInterface<GridView>;
      using IntersectionPartGeometry = AffineIntersectionPartGeometry<GV>;
      using IntersectionPartImpl = DefaultIntersectionPart<IntersectionPartGeometry>;
      using IntersectionPartInterface = Dune::SubTriangulation::IntersectionPartInterface<GV>;
      using IntersectionPartPointer = IntersectionPartInterface*;
      using IntersectionPartList = std::vector<IntersectionPartPointer>;

      using BoundingBox = Dune::GridCuboid<Grid>;
      using LevelSetGridView = LGV;

      using LocalSubTriangulation = LocalSubTriangulation<GV, LGV>;

      explicit SimpleTpmcTriangulation(
          const GV& fundamentalGridView, const LGV& levelSetGridView,
          const DomainConfiguration<GV, LGV>& domainConfiguration, bool forceRefinement = false,
          ctype valueTolerance = 1e-12l)
          : fundamentalGridView_(fundamentalGridView)
          , fundamentalElementMapper_(fundamentalGridView_)
          , levelSetGridView_(levelSetGridView)
          , domainConfiguration_(domainConfiguration)
          , forceRefinement_(forceRefinement)
          , valueTolerance_(valueTolerance)
          , cutCellInformation_(fundamentalGridView_, levelSetGridView_, domainConfiguration_, 1)
          , localSubTriangulation_(domainConfiguration_, fundamentalGridView_, levelSetGridView_, cutCellInformation_, valueTolerance_)
      {
        // TODO: Get rid of magic numbers
        entityPartStorage_.reserve(10);
        intersectionPartStorage_.reserve(20);
      }

      void create_entity_parts(const Entity& e, EntityPartList& entityparts) const
      {
        entityPartStorage_.clear();

        if (!forceRefinement_ && cutCellInformation_.filledBySingleCutCell(e)) {
          auto domainIndex = cutCellInformation_.domainOfSingleCutCell(e);
          const auto& info = cutCellInformation_.information(e, domainIndex);
          VolumeSnippet<GV, LGV> snippet(
              0, SimpleTpmcTriangulationDetail::createIdentityGeometry<
                     typename VolumeSnippet<GV, LGV>::GeometryInFundamental>(e.geometry()),
              e);
          entityparts.push_back(createEntityPartFromSnippet(snippet, e,
              domainIndex, info.boundingBox, info.volume));
        } else {
          localSubTriangulation_.bindOnVolume(e);

          for (auto volumeIt = localSubTriangulation_.volumeSnippetsBegin(); volumeIt != localSubTriangulation_.volumeSnippetsEnd(); ++volumeIt) {
            int domainIndex = localSubTriangulation_.snippetDomainIndex(*volumeIt);

            if (localSubTriangulation_.hasCutCell(domainIndex)) {
              const auto& info = cutCellInformation_.information(e, domainIndex);
              entityparts.push_back(createEntityPartFromSnippet(
                *volumeIt, e, domainIndex, info.boundingBox, info.volume));
            }
          }
        }
      }

      void create_intersection_parts(
          const Entity& e, IntersectionPartList& intersectionparts) const
      {
        intersectionPartStorage_.clear();
        DUNE_THROW(Dune::Exception, "SimpleTpmcTriangulation::create_intersection_parts not implemented");
      }

      const BoundingBox boundingBox(const EntityPart& ep) const
      {
        return ep.boundingBox();
      }

      const BoundingBox BBox(const unsigned int domainIndex, const Entity& e) const
      {
        return cutCellInformation_.information(e, domainIndex).boundingBox;
      }

      bool& forceRefinement()
      {
        return forceRefinement_;
      }

      bool forceRefinement() const
      {
        return forceRefinement_;
      }

      bool isInDomain(const Entity& e, const int& domain_index) const
      {
        assert(false);
        return false;
      }
      bool isHostCell(const Entity& hostCandidate, const int& domainIndex) const
      {
        return cutCellInformation_.cutCellsExist(hostCandidate, domainIndex);
      }
      bool isHostCell(const Entity& hostCandidate) const
      {
        for (const auto& domain : domainConfiguration_.domains()) {
          if (isHostCell(hostCandidate, domain.index())) {
            return true;
          }
        }
        return false;
      }

      const GV& gridView() const
      {
        return fundamentalGridView_;
      }

      const LGV& levelSetGridView() const
      {
        return levelSetGridView_;
      }

      const DomainConfiguration<GV,LGV>& domainConfiguration() const
      {
        return domainConfiguration_;
      }

      const CutCellInformation<GV, LGV>& cutCellInformation() const
      {
        return cutCellInformation_;
      }

      CutCellInformation<GV, LGV>& cutCellInformation()
      {
        return cutCellInformation_;
      }

      LocalSubTriangulation& localSubTriangulation()
      {
        return localSubTriangulation_;
      }

    private:
      EntityPartPointer createEntityPartFromSnippet(const VolumeSnippet<GV, LGV>& snippet,
          const Entity& e, int domainIndex, const BoundingBox& boundingBox, ctype volume) const
      {
        // transform coordinates in fundamental cell to global coordinates
        std::vector<FieldVector<ctype, dim> > globalCornerCoordinates;
        for (int i = 0; i < snippet.geometryInFundamental().corners(); ++i) {
          globalCornerCoordinates.push_back(
              e.geometry().global(snippet.geometryInFundamental().corner(i)));
        }

        // create geometry of entity part
        typename EntityPartGeometry::Traits::Geometry geometry(
            Dune::geometryTypeFromVertexCount(dim, globalCornerCoordinates.size()),
            globalCornerCoordinates);

        EntityPartGeometry partGeometry(geometry, snippet.homeElement(), boundingBox);
        entityPartStorage_.emplace_back(partGeometry, domainIndex);

        return &entityPartStorage_.back();
      }

      IntersectionPartPointer createIntersectionPartFromSnippet(
          const InternalInterfaceSnippet<GV, LGV>& snippet, const Entity& e) const
      {
        DUNE_THROW(Dune::Exception, "SimpleTpmcTriangulation::createIntersectionPartFromSnippet not implemented");
      }

      IntersectionPartPointer createIntersectionPartFromSnippet(
          const ExternalInterfaceSnippet<GV, LGV>& snippet, const Entity& e, const Domain& domain) const
      {
        // transform coordinates in fundamental cell to global coordinates
        std::vector<FieldVector<ctype, dim> > globalCornerCoordinates;
        for (int i = 0; i < snippet.globalGeometry().corners(); ++i) {
          globalCornerCoordinates.push_back(snippet.globalGeometry().corner(i));
        }

        // create geometry of intersection part
        typename IntersectionPartGeometry::Traits::FaceGeometry geometry(
            Dune::geometryTypeFromVertexCount(dim - 1, globalCornerCoordinates.size()),
            globalCornerCoordinates);

        const auto& interiorInformation = cutCellInformation_.information(e, domain.index());

        if (snippet.outsideFundamentalElementIndex() == fundamentalElementMapper_.index(e)) {
          IntersectionPartGeometry igeom(geometry, interiorInformation.boundingBox,
              interiorInformation.boundingBox, snippet.globalUnitOuterNormal());

          intersectionPartStorage_.emplace_back(
              igeom, domain.index(), -1);

          return &intersectionPartStorage_.back();
        } else {
          const auto& exteriorInformation = cutCellInformation_.information(
              snippet.outsideFundamentalElementIndex(), domain.index());

          IntersectionPartGeometry igeom(geometry, interiorInformation.boundingBox,
              exteriorInformation.boundingBox, snippet.globalUnitOuterNormal());

          intersectionPartStorage_.emplace_back(
              igeom, domain.index(),
              domain.index());

          return &intersectionPartStorage_.back();
        }
      }

      GV fundamentalGridView_;
      Dune::SingleCodimSingleGeomTypeMapper<GV, 0> fundamentalElementMapper_;
      LGV levelSetGridView_;
      DomainConfiguration<GV, LGV> domainConfiguration_;
      bool forceRefinement_;
      ctype valueTolerance_;
      CutCellInformation<GV, LGV> cutCellInformation_;
      mutable LocalSubTriangulation localSubTriangulation_;

      mutable std::vector<EntityPartImpl> entityPartStorage_;
      mutable std::vector<IntersectionPartImpl> intersectionPartStorage_;
    };
  }
}

#endif // DUNE_SubTriangulation_SIMPLETPMCTRIANGULATION_HH

#ifndef DUNE_SUBTRIANGULATION_LOCALSUBTRIANGULATION_HH
#define DUNE_SUBTRIANGULATION_LOCALSUBTRIANGULATION_HH

#include <algorithm>
#include <bits/ranges_util.h>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/scsgmapper.hh>

#include <dune/subtriangulation/common/cuboid.hh>
#include <dune/subtriangulation/common/entitypart.hh>
#include <dune/subtriangulation/common/intersectionpart.hh>
#include <dune/subtriangulation/simpletpmctriangulation/cutcellinformation.hh>
#include <dune/subtriangulation/simpletpmctriangulation/domainconfiguration.hh>
#include <dune/subtriangulation/simpletpmctriangulation/interface.hh>
#include <dune/subtriangulation/simpletpmctriangulation/interfacesnippet.hh>
#include <dune/subtriangulation/simpletpmctriangulation/volumesnippet.hh>
#include <dune/subtriangulation/tpmc/tpmcrefinement.hh>

namespace Dune::SubTriangulation {

  template <class ctype>
  FieldVector<ctype, 3> crossProduct(
      const FieldVector<ctype, 3>& a, const FieldVector<ctype, 3>& b)
  {
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
  }

  template <class ctype, int size>
  FieldVector<ctype, 3>
  computeUnitOuterNormal(const Dune::ReservedVector<FieldVector<ctype, 3>, size>& corners,
                          const std::vector<ctype>& values)
  {
    // move corner 0 to the origin
    auto a = corners[1];
    auto b = corners[2];
    a -= corners[0];
    b -= corners[0];
    // compute normal to corner 1 and 2
    auto cp = crossProduct(a, b);
    // normalize
    cp /= cp.two_norm();
    return cp;
  }

  template <class ctype, int size>
  FieldVector<ctype, 2>
  computeUnitOuterNormal(const ReservedVector<FieldVector<ctype, 2>, size>& corners,
                         const std::vector<ctype>& values)
  {
    std::array<Dune::FieldVector<ctype, 2>, 4> basisGradients
    {
      { {-1.0, -1.0}, {1.0, -1.0}, {-1.0, 1.0}, {1.0, 1.0} }
    };

    Dune::FieldVector<ctype, 2> gradient;

    for (int i = 0; i < 4; ++i)
    {
      gradient += values[i] * basisGradients[i];
    }

    // move corner 0 to the origin
    auto diff = corners[0];
    diff -= corners[1];
    std::swap(diff[0], diff[1]);
    diff[0] *= -1;

    diff /= diff.two_norm();

    if (diff * gradient < 0.0) {
      diff *= -1.0;
    }

    return diff;
  }

  template <class ctype, int size>
  FieldVector<ctype, 1>
  computeUnitOuterNormal(const Dune::ReservedVector<FieldVector<ctype, 1>, size>& corners,
                          const std::vector<ctype>& values)
  {
    assert(values.size() == 2);
    if (values[0] < values[1])
      return {1};
    else
      return {-1};
  }

  template <class GOut>
  class CreateIdentityGeometry
  {
    using ctype = typename GOut::ctype;

  public:
    template<class GT>
    GOut operator()(const GT& gt)
    {
      const auto& referenceElement = ReferenceElements<ctype, GOut::coorddimension>::general(gt);

      corners.clear();

      for (int i = 0; i < referenceElement.size(GOut::mydimension); ++i) {
        corners.emplace_back(referenceElement.position(i, GOut::coorddimension));
      }

      return GOut(gt, corners);
    }

    private:
      typename ReservedStorageMultiLinearGeometryTraits<ctype>::template CornerStorage<GOut::mydimension, GOut::coorddimension>::Type corners;
  };

  template<class Target>
  class CombineGeometries
  {
  public:
    template<class Guest, class Host>
    constexpr Target operator()(const Guest& guest, const Host& host)
    {
      coords.resize(guest.corners());
  
      for (unsigned int v = 0; v < coords.size(); ++v) {
        coords[v] = host.global(guest.corner(v));
      }

      return Target(guest.type(), coords);
    }

    template<class Guest, class Host>
    constexpr Target operator()(const Guest& guest, const Host& host, std::false_type forward)
    {
      coords.resize(guest.corners());
  
      for (unsigned int v = 0; v < coords.size(); ++v) {
        coords[v] = host.local(guest.corner(v));
      }

      return Target(guest.type(), coords);
    }

  private:
    typename ReservedStorageMultiLinearGeometryTraits<typename Target::ctype>::template CornerStorage<Target::mydimension, Target::coorddimension>::Type coords;
  };

  template <class To>
  class ConvertGeometry
  {
  public:
    template<class From>
    To operator()(const From& from)
    {
      coords.clear();

      for (int i = 0; i < from.corners(); ++i) {
        coords.push_back(from.corner(i));
      }

      return To(from.type(), coords);
    }

  private:
    typename ReservedStorageMultiLinearGeometryTraits<typename To::ctype>::template CornerStorage<To::mydimension, To::coorddimension>::Type coords;
  };

  template <class Geometry, class FundamentalElement, class VolumeSnippetIt, class LST>
  class CutCellGeometry
  {
  public:
    friend LST;

    using ctype = typename Geometry::ctype;
    static const int coorddimension = Geometry::coorddimension;

    using QuadraturePoint = Dune::QuadraturePoint<ctype, coorddimension>;

    void quadratureRule(std::vector<QuadraturePoint>& quadraturePoints, int intorder) const
    {
      quadraturePoints.clear();

      for (auto snippetsIt = snippetsBegin_; snippetsIt != snippetsEnd_; ++snippetsIt) {
        const auto& rule = Dune::QuadratureRules<ctype, coorddimension>::rule(snippetsIt->geometryInFundamental().type(), intorder);

        for (auto qp : rule) {
          auto position = snippetsIt->geometryInFundamental().global(qp.position());
          auto weight = fundamentalGeometry_->integrationElement(position);
          position = fundamentalGeometry_->global(position);
          position = geometry_.local(position);
          weight *= qp.weight() * snippetsIt->geometryInFundamental().integrationElement(qp.position());

          quadraturePoints.emplace_back(position, weight);
        }
      }
    }

    //! \brief Domain index this cut-cell is part of
    int domainIndex() const
    {
      return domainIndex_;
    }

    const Geometry& geometry() const
    {
      return geometry_;
    }

    bool isTriangle() const
    {
      if constexpr (coorddimension == 2) {
        return std::distance(snippetsBegin_, snippetsEnd_) == 1 && snippetsBegin_->geometryInFundamental().corners() == 3;
      }

      return false;
    }

  private:
    Geometry geometry_;
    typename FundamentalElement::Geometry* fundamentalGeometry_;
    VolumeSnippetIt snippetsBegin_;
    VolumeSnippetIt snippetsEnd_;
    int domainIndex_;
  };

  template <class Geometry, class FundamentalElement, class InterfaceSnippetIt, class LST>
  class CutIntersectionGeometry
  {
  public:
    friend LST;

    using ctype = typename Geometry::ctype;
    enum { coorddimension=Geometry::coorddimension };

    using QuadraturePoint = Dune::QuadraturePoint<ctype, coorddimension>;

    void quadratureRule(std::vector<QuadraturePoint>& quadraturePoints, int intorder) const
    {
      quadraturePoints.clear();

      for (auto snippetsIt = snippetsBegin_; snippetsIt != snippetsEnd_; ++snippetsIt) {
        const auto& rule = Dune::QuadratureRules<ctype, coorddimension - 1>::rule(snippetsIt->globalGeometry().type(), intorder);

        for (auto qp : rule) {
          auto position = snippetsIt->globalGeometry().global(qp.position());
          auto weight = snippetsIt->globalGeometry().integrationElement(qp.position()) * qp.weight();

          quadraturePoints.emplace_back(position, weight);
        }
      }
    }

    const Geometry& geometry() const
    {
      return globalGeometry_;
    }

    const Geometry& geometryInInside() const
    {
      return geometryInInside_;
    }

    const Geometry& geometryInOutside() const
    {
      return geometryInOutside_;
    }

    const FundamentalElement& inside() const
    {
      // TODO: We are assuming here that we have a bounding box
      return geometryInInside_.entity();
    }

    const FundamentalElement& outside() const
    {
      // TODO: We are assuming here that we have a bounding box
      return geometryInOutside_.entity();
    }

    int insideDomainIndex() const
    {
      return insideDomainIndex_;
    }

    int outsideDomainIndex() const
    {
      return outsideDomainIndex_;
    }

    bool neighbor() const
    {
      return outsideDomainIndex_ >= 0;
    }

    Dune::FieldVector<ctype, coorddimension> unitOuterNormal(const Dune::FieldVector<ctype, coorddimension>& x) const
    {
      return unitOuterNormal_;
    }

    Dune::FieldVector<ctype, coorddimension> center() const
    {
      return snippetsBegin_->globalGeometry().center();
    }

    bool domainIsOutside(int domainIndex) const
    {
      return domainIndex != insideDomainIndex() && domainIndex == outsideDomainIndex();
    }

    void invertGeometry()
    {
      std::swap(geometryInInside_, geometryInOutside_);
      std::swap(insideDomainIndex_, outsideDomainIndex_);
      unitOuterNormal_ *= -1.0;
    }

  private:
    Geometry globalGeometry_;
    Geometry geometryInInside_;
    Geometry geometryInOutside_;
    Dune::FieldVector<ctype, coorddimension> unitOuterNormal_;
    InterfaceSnippetIt snippetsBegin_;
    InterfaceSnippetIt snippetsEnd_;
    int insideDomainIndex_;
    int outsideDomainIndex_;
  };

  template<class GV, class LGV>
  class LocalSubTriangulation
  {
  public:
    using GridView = GV;
    using LevelSetGridView = LGV;

    using FundamentalElement = typename GridView::template Codim<0>::Entity;
    using FundamentalIntersection = typename GridView::Intersection;
    using LevelSetElement = typename LevelSetGridView::template Codim<0>::Entity;
    using ctype = typename GridView::ctype;
    enum { dim = GridView::dimension };

    using EntityPartGeometry = Dune::SubTriangulation::DefaultEntityPartGeometry<GV>;
    using EntityPart = DefaultEntityPart<EntityPartGeometry>;
    using IntersectionPartGeometry = AffineIntersectionPartGeometry<GV>;
    using IntersectionPart = DefaultIntersectionPart<IntersectionPartGeometry>;

    using EntityPartList = std::vector<EntityPart*>;
    using IntersectionPartList = std::vector<IntersectionPart*>;

    using VolumeSnippet = VolumeSnippet<GridView, LevelSetGridView>;
    using InternalInterfaceSnippet = InternalInterfaceSnippet<GridView, LevelSetGridView>;
    using ExternalInterfaceSnippet = ExternalInterfaceSnippet<GridView, LevelSetGridView>;

    using VolumeSnippetIdentityGeometry = CreateIdentityGeometry<typename VolumeSnippet::GeometryInFundamental>;
    using VolumeSnippetCombineGeometries = CombineGeometries<typename VolumeSnippet::GeometryInFundamental>;

    using InternalInterfaceCombineGeometries = CombineGeometries<typename InternalInterfaceSnippet::GlobalGeometry>;
    using ExternalInterfaceIdentityGeometry = CreateIdentityGeometry<typename ExternalInterfaceSnippet::GlobalGeometry>;
    using ExternalInterfaceCombineGeometries = CombineGeometries<typename ExternalInterfaceSnippet::GlobalGeometry>;
    using ExternalInterfaceConvertGeometry = ConvertGeometry<typename ExternalInterfaceSnippet::GlobalGeometry>;

    using VolumeRefinement = TpmcRefinement<ctype, dim>;
    using FaceRefinement = TpmcRefinement<ctype, dim - 1>;

    using BoundingBox = Dune::GridCuboid<typename GridView::Grid>;
    using CutCell = CutCellGeometry<BoundingBox, FundamentalElement, typename std::vector<VolumeSnippet>::iterator, LocalSubTriangulation<GridView, LevelSetGridView>>;
    using CutIntersection = CutIntersectionGeometry<BoundingBox, FundamentalElement, typename std::vector<InternalInterfaceSnippet>::iterator, LocalSubTriangulation<GridView, LevelSetGridView>>;

    LocalSubTriangulation(const DomainConfiguration<GridView, LevelSetGridView>& domainConfiguration,
                          const GridView& fundamentalGridView, const LevelSetGridView& levelSetGridView,
                          CutCellInformation<GridView, LevelSetGridView>& cutCellInformation,
                          bool forceRefinement = false, ctype valueTolerance = 1e-12l)
      : domainConfiguration_(domainConfiguration)
      , fundamentalGridView_(fundamentalGridView)
      , levelSetGridView_(levelSetGridView)
      , elementMapper_(fundamentalGridView)
      , fundamentalGeometry_(typename GridView::Grid::Traits::template Codim<0>::GeometryImpl(Dune::FieldVector<ctype, dim>(0.0), Dune::FieldVector<ctype, dim>(1.0)))
      , cutCellInformation_(cutCellInformation)
      , valueTolerance_(valueTolerance)
      , forceRefinement_(forceRefinement)
    {
      cutCells_.resize(domainConfiguration.size());
      // TODO: We should be able to estimate the proper number in advance
      // same for everything below
      cutIntersections_.resize(domainConfiguration.numberOfInterfaces() * 10);

      volumeSnippets_.resize(10, VolumeSnippet(0, volumeSnippetIdentityGeometry_(Dune::GeometryTypes::cube(dim)), FundamentalElement()));
      newVolumeSnippets_.resize(10, VolumeSnippet(0, volumeSnippetIdentityGeometry_(Dune::GeometryTypes::cube(dim)), FundamentalElement()));

      internalInterfaceSnippets_.reserve(20);
      internalIntersectionToInterface_.resize(30);

      externalInterfaceSnippets_.resize(30, InternalInterfaceSnippet(
              externalInterfaceIdentityGeometry_(Dune::GeometryTypes::cube(dim-1)),
              Dune::FieldVector<ctype, dim>(),
              FundamentalElement(), -1));

      newExternalInterfaceSnippets_.resize(30, InternalInterfaceSnippet(
              externalInterfaceIdentityGeometry_(Dune::GeometryTypes::cube(dim-1)),
              Dune::FieldVector<ctype, dim>(),
              FundamentalElement(), -1));

      externalInterfaceSnippetEnd_ = externalInterfaceSnippets_.begin();
      newExternalInterfaceSnippetEnd_ = newExternalInterfaceSnippets_.begin();

      domainToVolume_.resize(domainConfiguration.size());
      domainToCenterOfMass_.resize(domainConfiguration.size());
      domainToBoundingBox_.resize(domainConfiguration.size());

      cutsPerIniternalInterface_ = 1 << (domainConfiguration_.numberOfInterfaces() - 1);

      for (std::size_t i = 0; i < domainConfiguration_.numberOfInterfaces(); ++i) {
        std::size_t interfaceOffset = i * cutsPerIniternalInterface_;
        internalInterfaceRelativePositions_.push_back(std::vector<InterfaceRelativePosition>());
        std::size_t currentSize = 1;

        for (std::size_t j = 0; j < domainConfiguration_.numberOfInterfaces(); ++j) {
          if (i != j) {
            for (std::size_t k = 0; k < currentSize; ++k) {
              internalInterfaceRelativePositions_.push_back(internalInterfaceRelativePositions_[interfaceOffset + k]);
              internalInterfaceRelativePositions_[interfaceOffset + k].push_back(InterfaceRelativePosition::interior);
              internalInterfaceRelativePositions_.back().push_back(InterfaceRelativePosition::exterior);
            }

            currentSize *= 2;
          } else {
            for (std::size_t k = 0; k < currentSize; ++k) {
              internalInterfaceRelativePositions_[interfaceOffset + k].push_back(InterfaceRelativePosition::interface);
            }
          }
        }
      }

      externalInterfaceRelativePositions_.push_back(std::vector<InterfaceRelativePosition>());

      for (std::size_t i = 0; i < domainConfiguration_.numberOfInterfaces(); ++i) {
        std::size_t currentSize = externalInterfaceRelativePositions_.size();

        for (std::size_t j = 0; j < currentSize; ++j) {
          externalInterfaceRelativePositions_.push_back(externalInterfaceRelativePositions_[j]);
          externalInterfaceRelativePositions_[j].push_back(InterfaceRelativePosition::interior);
          externalInterfaceRelativePositions_.back().push_back(InterfaceRelativePosition::exterior);
        }
      }

      volumeRelativePositions_.push_back(std::vector<InterfaceRelativePosition>());

      for (std::size_t i = 0; i < domainConfiguration_.numberOfInterfaces(); ++i) {
        std::size_t currentSize = volumeRelativePositions_.size();

        for (std::size_t j = 0; j < currentSize; ++j) {
          volumeRelativePositions_.push_back(volumeRelativePositions_[j]);
          volumeRelativePositions_[j].push_back(InterfaceRelativePosition::interior);
          volumeRelativePositions_.back().push_back(InterfaceRelativePosition::exterior);
        }
      }

      cutCellInformation_.init(*this, domainConfiguration_);
      cutCellInformation_.printStatistics();
    }

    LocalSubTriangulation(const LocalSubTriangulation&) = delete;
    LocalSubTriangulation& operator=(const LocalSubTriangulation&) = delete;

    void bind(const FundamentalElement& fundamentalElement)
    {
      fundamentalElement_ = fundamentalElement;
      fundamentalGeometry_ = fundamentalElement_.geometry();

      constructVolumeSnippets();
      constructInternalInterfaceSnippets();
      constructExternalInterfaceSnippets();
    }

    void bindOnVolume(const FundamentalElement& fundamentalElement)
    {
      fundamentalElement_ = fundamentalElement;
      fundamentalGeometry_ = fundamentalElement_.geometry();

      constructVolumeSnippets();
    }

    void setupInitialVolumeSnippets()
    {
      auto maxLevel = levelSetGridView_.grid().maxLevel();

      if (fundamentalElement_.isLeaf()) {
        *volumeSnippetEnd_ = VolumeSnippet(0, volumeSnippetIdentityGeometry_(fundamentalGeometry_.type()), fundamentalElement_);
        ++volumeSnippetEnd_;
      } else {
        // loop through all child elements in the level set mesh and create snippets
        for (auto it = fundamentalElement_.hbegin(maxLevel);
              it != fundamentalElement_.hend(maxLevel); ++it) {
          if (!levelSetGridView_.indexSet().contains(*it)) {
            continue;
          }

          *volumeSnippetEnd_ = VolumeSnippet(0, volumeSnippetCombineGeometries_(it->geometry(), fundamentalGeometry_, std::false_type()), *it);
          ++volumeSnippetEnd_;
        }
      }
    }

    void setupInitialInternalSnippets(const Interface<LevelSetGridView>& interface, const LevelSetElement& levelSetElement)
    {
      const auto& lgeo = levelSetElement.geometry();
      auto localInterfaceFunction = localFunction(interface.function());
      localInterfaceFunction.bind(levelSetElement);
      const auto& lref = Dune::ReferenceElements<ctype, dim>::general(lgeo.type());
      functionValues_.resize(lgeo.corners(), 0.0);

      for (int i = 0; i < lgeo.corners(); ++i) {
        functionValues_[i] = localInterfaceFunction(lref.position(i, dim));
      }

      std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val > 0.0 && val < valueTolerance_;}, valueTolerance_);
      std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val < 0.0 && val > -valueTolerance_;}, -valueTolerance_);

      volumeRefinement_.bind(functionValues_.begin(), functionValues_.end());
      volumeRefinement_.reconstructInterface();

      for (const auto& i : volumeRefinement_.interface()) {
        typename ReservedStorageMultiLinearGeometryTraits<ctype>::template CornerStorage<std::max(dim - 1, 0), dim>::Type globalCornerCoordinates;

        for (int c = 0; c < i.corners(); ++c) {
          globalCornerCoordinates.push_back(lgeo.global(i.corner(c)));
        }

        internalInterfaceSnippets_.emplace_back(InternalInterfaceSnippet(
            i.type(), globalCornerCoordinates,
            computeUnitOuterNormal(globalCornerCoordinates, functionValues_),
            levelSetElement, interface.index() * cutsPerIniternalInterface_));
      }

      internalIntersectionToInterface_[interface.index() * cutsPerIniternalInterface_] = interface.index();
    }

    void setupInitialExternalSnippets()
    {
      if (fundamentalElement_.isLeaf()) {
        for (const auto& intersection : intersections(fundamentalGridView_, fundamentalElement_)) {
          const auto& intersectionGeometry = intersection.geometry();

          // TODO: This way the complete external boundary of a cell will form the initial intersection
          // This is the old implemented behavior which might not have been intended.
          // Instead we could form seperate intersections from the fundamental intersections
          *externalInterfaceSnippetEnd_ = InternalInterfaceSnippet(
              externalInterfaceConvertGeometry_(intersectionGeometry),
              intersection.unitOuterNormal(
                  ReferenceElements<ctype, dim - 1>::general(intersectionGeometry.type())
                      .position(0, 0)),
              fundamentalElement_, 0, elementMapper_.index(intersection.neighbor() ? intersection.outside() : fundamentalElement_));
          ++externalInterfaceSnippetEnd_;
        }
      } else {
        auto maxLevel = levelSetGridView_.grid().maxLevel();

        // loop through all child elements in the level set mesh
        for (const auto& lelem : Dune::descendantElements(fundamentalElement_, maxLevel)) {
          if (!levelSetGridView_.indexSet().contains(lelem)) {
            continue;
          }

          for (const auto& intersection : intersections(levelSetGridView_, lelem)) {
            FundamentalElement outsideFundamental = fundamentalElement_;
            // check if the intersection is internal to this fundamental element
            if (intersection.neighbor()) {
              outsideFundamental = intersection.outside();
              while (!fundamentalGridView_.indexSet().contains(outsideFundamental)) {
                outsideFundamental = outsideFundamental.father();
              }
              if (outsideFundamental == fundamentalElement_) {
                continue;
              }
            }

            const auto& intersectionGeometry = intersection.geometry();
            FundamentalIntersection fundamentalIntersection;

            // find intersection in fundamental element
            if (outsideFundamental != fundamentalElement_) {
              for (const auto& fintersection : intersections(fundamentalGridView_, fundamentalElement_)) {
                if (fintersection.neighbor() && fintersection.outside() == outsideFundamental) {
                  fundamentalIntersection = fintersection;
                  break;
                }
              }
            } else { // meaning, no neighbor
              // find the intersection of the fundamental element which contains all corners of
              // the levelset intersection
              for (const auto& fintersection : intersections(fundamentalGridView_, fundamentalElement_)) {
                const auto& fgeo = fintersection.geometry();
                bool inside = true;

                for (int i = 0; i < intersectionGeometry.corners(); ++i) {
                  auto diff = intersectionGeometry.corner(i);
                  diff -= fgeo.global(fgeo.local(intersectionGeometry.corner(i)));

                  if (diff.two_norm() > 1e-12) {
                    inside = false;
                    break;
                  }
                }

                if (inside) {
                  fundamentalIntersection = fintersection;
                  break;
                }
              }
            }

            *externalInterfaceSnippetEnd_ = InternalInterfaceSnippet(
                externalInterfaceConvertGeometry_(intersectionGeometry),
                intersection.unitOuterNormal(
                    ReferenceElements<ctype, dim - 1>::general(intersectionGeometry.type())
                        .position(0, 0)),
                lelem, 0, elementMapper_.index(outsideFundamental));
            ++externalInterfaceSnippetEnd_;
          }
        }
      }
    }

    void constructVolumeSnippets()
    {
      volumeSnippetEnd_ = volumeSnippets_.begin();
      newVolumeSnippetEnd_ = newVolumeSnippets_.begin();

      volumeRanges_.clear();

      std::fill(domainToVolume_.begin(), domainToVolume_.end(), 0.0);
      std::fill(domainToCenterOfMass_.begin(), domainToCenterOfMass_.end(), Dune::FieldVector<ctype, dim>(0.0));
      std::fill(domainToBoundingBox_.begin(), domainToBoundingBox_.end(), RawBoundingBox<ctype, dim>());

      setupInitialVolumeSnippets();

      volumeIndexOffset_ = 1;

      // cut all interfaces into cut cells recursively
      for (const Interface<LGV>& interface : domainConfiguration_.interfaces()) {
        newVolumeSnippetEnd_ = newVolumeSnippets_.begin();
        cutVolume(interface, volumeSnippets_.begin(), volumeSnippetEnd_);
        volumeSnippets_.swap(newVolumeSnippets_);
        volumeSnippetEnd_ = newVolumeSnippetEnd_;

        volumeIndexOffset_ *= 2;
      }

      snippetIndexToDomain_.resize(volumeRelativePositions_.size());

      auto fillInDomainIndex = [this] (const auto& snippet)
      {
        auto domainIt = domainConfiguration_.findDomain(volumeRelativePositions_[snippet.index()]);

        if (domainIt != domainConfiguration_.domainsEnd()) {
          snippetIndexToDomain_[snippet.index()] = domainIt->index();
        } else {
          snippetIndexToDomain_[snippet.index()] = -1;
        }
      };

      std::for_each(volumeSnippets_.begin(), volumeSnippetEnd_, fillInDomainIndex);

      auto cmp = [this] (const auto& a, const auto& b) { return snippetIndexToDomain_[a.index()] < snippetIndexToDomain_[b.index()]; };
      std::sort(volumeSnippets_.begin(), volumeSnippetEnd_, cmp);

      auto it = volumeSnippets_.begin();

      while (it != volumeSnippetEnd_) {
        it = std::find_if(it, volumeSnippetEnd_, [=](const auto& b) { return cmp(*it, b); });
        volumeRanges_.push_back(it);
      }

      auto volumeScaling = fundamentalGeometry_.integrationElement(
        ReferenceElements<ctype, dim>::general(fundamentalGeometry_.type()).position(0, 0));

      // accumulate volumes of all cut cells belonging to the same domain
      it = volumeSnippets_.begin();

      for (auto volumeEnd = volumeRanges_.begin(); volumeEnd != volumeRanges_.end(); ++volumeEnd) {
        int domainIndex = snippetIndexToDomain_[it->index()];

        if (domainIndex > -1) {
          ctype volumeInFundamental = 0.0;

          for (auto partIt = it; partIt != *volumeEnd; ++partIt) {
            for (int i = 0; i < partIt->geometryInFundamental().corners(); ++i) {
              domainToBoundingBox_[domainIndex].update(fundamentalGeometry_.global(partIt->geometryInFundamental().corner(i)));
            }

            const auto& rule = Dune::QuadratureRules<ctype, dim>::rule(
              partIt->geometryInFundamental().type(), 2 * dim + 1);
            ctype snippetVolume = 0.;

            for (const auto& qp : rule) {
              snippetVolume += qp.weight() * partIt->geometryInFundamental().integrationElement(qp.position());
            }

            volumeInFundamental += snippetVolume;
          }

          domainToVolume_[domainIndex] += volumeInFundamental * volumeScaling;
        }

        it = *volumeEnd;
      }

      // rescale center of mass by complete domain volume
      for (std::size_t i = 0; i < domainConfiguration_.numberOfDomains(); ++i) {
        if (domainToVolume_[i] > 0.0) {
          domainToCenterOfMass_[i] /= domainToVolume_[i];
        }
      }
    }

    void constructInternalInterfaceSnippets()
    {
      internalInterfaceSnippets_.clear();
      newInternalInterfaceSnippets_.clear();
      internalInterfaceRanges_.clear();

      if (fundamentalElement_.isLeaf()) {
        for (const Interface<LevelSetGridView>& interface : domainConfiguration_.interfaces()) {
          setupInitialInternalSnippets(interface, fundamentalElement_);
        }
      } else {
        auto maxLevel = levelSetGridView_.grid().maxLevel();
        // loop through all child elements in the level set mesh
        for (auto it = fundamentalElement_.hbegin(maxLevel);
              it != fundamentalElement_.hend(maxLevel); ++it) {
          if (!levelSetGridView_.indexSet().contains(*it)) {
            continue;
          }

          for (const Interface<LGV>& interface : domainConfiguration_.interfaces()) {
            setupInitialInternalSnippets(interface, *it);
          }
        }
      }

      internalIntersectionsIndexOffsets_.clear();
      internalIntersectionsIndexOffsets_.resize(domainConfiguration_.numberOfInterfaces(), 1);

      // cut all other interfaces into each internalinterface
      for (const Interface<LevelSetGridView>& interface : domainConfiguration_.interfaces()) {
        newInternalInterfaceSnippets_.clear();
        cutInternalInterface(interface, internalInterfaceSnippets_.begin(), internalInterfaceSnippets_.end());
        internalInterfaceSnippets_.swap(newInternalInterfaceSnippets_);

        for (const Interface<LevelSetGridView>& interfaceToAdvance : domainConfiguration_.interfaces()) {
          if (interfaceToAdvance.index() != interface.index()) {
            internalIntersectionsIndexOffsets_[interfaceToAdvance.index()] *= 2;
          }
        }
      }
    }

    void constructExternalInterfaceSnippets()
    {
      externalInterfaceSnippetEnd_ = externalInterfaceSnippets_.begin();
      newExternalInterfaceSnippetEnd_ = newExternalInterfaceSnippets_.begin();
      externalInterfaceRanges_.clear();

      setupInitialExternalSnippets();

      externalIntersectionIndexOffset_ = 1;

      for (const Interface<LGV>& interface : domainConfiguration_.interfaces()) {
        newExternalInterfaceSnippetEnd_ = newExternalInterfaceSnippets_.begin();
        cutExternalInterface(interface, externalInterfaceSnippets_.begin(), externalInterfaceSnippetEnd_);
        externalInterfaceSnippets_.swap(newExternalInterfaceSnippets_);
        externalInterfaceSnippetEnd_ = newExternalInterfaceSnippetEnd_;

        externalIntersectionIndexOffset_ *= 2;
      }

      auto cmp = [] (const auto& a, const auto& b) { return a.index() < b.index(); };

      std::sort(externalInterfaceSnippets_.begin(), externalInterfaceSnippetEnd_, cmp);

      auto it = externalInterfaceSnippets_.begin();

      while (it != externalInterfaceSnippetEnd_) {
        it = std::find_if(it, externalInterfaceSnippetEnd_, [=] (const auto& a) { return cmp(*it, a); });
        externalInterfaceRanges_.push_back(it);
      }
    }

    void cutVolume(const Interface<LevelSetGridView>& interface, typename std::vector<VolumeSnippet>::iterator it, typename std::vector<VolumeSnippet>::iterator end)
    {
      auto localInterfaceFunction = localFunction(interface.function());

      for ( ; it != end; ++it) {
        localInterfaceFunction.bind(it->homeElement());
        const auto& hgeo = it->homeElement().geometry();
        functionValues_.resize(it->geometryInFundamental().corners());

        for (int i = 0; i < it->geometryInFundamental().corners(); ++i) {
          functionValues_[i] = localInterfaceFunction(hgeo.local(fundamentalGeometry_.global(it->geometryInFundamental().corner(i))));
        }

        std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val > 0.0 && val < valueTolerance_;}, valueTolerance_);
        std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val < 0.0 && val > -valueTolerance_;}, -valueTolerance_);
        bool allInside = std::all_of(functionValues_.begin(), functionValues_.end(), [] (ctype val) { return val < 0.0;});
        bool allOutside = std::all_of(functionValues_.begin(), functionValues_.end(), [] (ctype val) { return val > 0.0;});

        if (allInside) {
          *newVolumeSnippetEnd_ = *it;
          ++newVolumeSnippetEnd_;
        } else if (allOutside) {
          *newVolumeSnippetEnd_ = *it;
          newVolumeSnippetEnd_->addIndexOffset(volumeIndexOffset_);
          ++newVolumeSnippetEnd_;
        } else {
          // retrieve interior and exterior elements using tpmc
          volumeRefinement_.bind(functionValues_.begin(), functionValues_.end());
          volumeRefinement_.reconstructVolume(tpmc::InteriorDomain);
          volumeRefinement_.reconstructVolume(tpmc::ExteriorDomain);

          for (const auto& part : volumeRefinement_.volume(tpmc::InteriorDomain)) {
            if (part.volume() == 0.0) {
              DUNE_THROW(Dune::Exception, "got an empty part");
            }

            *newVolumeSnippetEnd_ = VolumeSnippet(it->index(),
                volumeSnippetCombineGeometries_(part, it->geometryInFundamental()),
                it->homeElement());
            ++newVolumeSnippetEnd_;
          }

          for (const auto& part : volumeRefinement_.volume(tpmc::ExteriorDomain)) {
            if (part.volume() == 0.0) {
              DUNE_THROW(Dune::Exception, "got an empty part");
            }

            *newVolumeSnippetEnd_ = VolumeSnippet(it->index() + volumeIndexOffset_,
                volumeSnippetCombineGeometries_(part, it->geometryInFundamental()),
                it->homeElement());
            ++newVolumeSnippetEnd_;
          }
        }
      }
    }

    void cutInternalInterface(const Interface<LevelSetGridView>& interface,
                              typename std::vector<InternalInterfaceSnippet>::iterator it,
                              typename std::vector<InternalInterfaceSnippet>::iterator end)
    {
      auto localInterfaceFunction = localFunction(interface.function());

      for ( ; it != end; ++it) {
        if (internalIntersectionToInterface_[it->index()] == interface.index()) {
          newInternalInterfaceSnippets_.emplace_back(*it);
          continue;
        }

        const auto& lgeo = it->levelSetElement().geometry();
        localInterfaceFunction.bind(it->levelSetElement());
        functionValues_.resize(it->globalGeometry().corners(), 0.0);

        for (int i = 0; i < it->globalGeometry().corners(); ++i) {
          functionValues_[i] = localInterfaceFunction(lgeo.local(it->globalGeometry().corner(i)));
        }

        std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val > 0.0 && val < valueTolerance_;}, valueTolerance_);
        std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val < 0.0 && val > -valueTolerance_;}, -valueTolerance_);

        // retrieve interior and exterior elements using tpmc
        faceRefinement_.bind(functionValues_.begin(), functionValues_.end());
        faceRefinement_.reconstructVolume(tpmc::InteriorDomain);
        faceRefinement_.reconstructVolume(tpmc::ExteriorDomain);

        for (auto partIt = faceRefinement_.beginVolume(tpmc::InteriorDomain);
              partIt != faceRefinement_.endVolume(tpmc::InteriorDomain); ++partIt) {
          newInternalInterfaceSnippets_.emplace_back(InternalInterfaceSnippet(
              internalInterfaceCombineGeometries_(*partIt, it->globalGeometry()),
              it->globalUnitOuterNormal(), it->levelSetElement(), it->index()));
        }

        for (auto partIt = faceRefinement_.beginVolume(tpmc::ExteriorDomain);
              partIt != faceRefinement_.endVolume(tpmc::ExteriorDomain); ++partIt) {
          newInternalInterfaceSnippets_.emplace_back(InternalInterfaceSnippet(
              internalInterfaceCombineGeometries_(*partIt, it->globalGeometry()),
              it->globalUnitOuterNormal(), it->levelSetElement(), it->index() + internalIntersectionsIndexOffsets_[internalIntersectionToInterface_[it->index()]]));
          internalIntersectionToInterface_[it->index() + internalIntersectionsIndexOffsets_[internalIntersectionToInterface_[it->index()]]] = internalIntersectionToInterface_[it->index()];
        }
      }
    }

    void cutExternalInterface(const Interface<LevelSetGridView>& interface, 
                              typename std::vector<InternalInterfaceSnippet>::iterator it,
                              typename std::vector<InternalInterfaceSnippet>::iterator end)
    {
      auto localInterfaceFunction = localFunction(interface.function());

      for ( ; it != end; ++it) {
        const auto& lgeo = it->levelSetElement().geometry();
        localInterfaceFunction.bind(it->levelSetElement());
        functionValues_.resize(it->globalGeometry().corners(), 0.0);

        for (int i = 0; i < it->globalGeometry().corners(); ++i) {
          functionValues_[i] = localInterfaceFunction(lgeo.local(it->globalGeometry().corner(i)));
        }

        std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val > 0.0 && val < valueTolerance_;}, valueTolerance_);
        std::replace_if(functionValues_.begin(), functionValues_.end(), [this] (ctype val) { return val < 0.0 && val > -valueTolerance_;}, -valueTolerance_);

        // retrieve interior and exterior elements using tpmc
        faceRefinement_.bind(functionValues_.begin(), functionValues_.end());
        faceRefinement_.reconstructVolume(tpmc::InteriorDomain);
        faceRefinement_.reconstructVolume(tpmc::ExteriorDomain);

        for (auto partIt = faceRefinement_.beginVolume(tpmc::InteriorDomain);
              partIt != faceRefinement_.endVolume(tpmc::InteriorDomain); ++partIt) {
          *newExternalInterfaceSnippetEnd_ = InternalInterfaceSnippet(
              externalInterfaceCombineGeometries_(*partIt, it->globalGeometry()),
              it->globalUnitOuterNormal(), it->levelSetElement(),
              it->index(), it->outsideFundamentalElementIndex());

          ++newExternalInterfaceSnippetEnd_;
        }

        for (auto partIt = faceRefinement_.beginVolume(tpmc::ExteriorDomain);
              partIt != faceRefinement_.endVolume(tpmc::ExteriorDomain); ++partIt) {
          *newExternalInterfaceSnippetEnd_ = InternalInterfaceSnippet(
              externalInterfaceCombineGeometries_(*partIt, it->globalGeometry()),
              it->globalUnitOuterNormal(), it->levelSetElement(),
              it->index() + externalIntersectionIndexOffset_,
              it->outsideFundamentalElementIndex());
          ++newExternalInterfaceSnippetEnd_;
        }
      }
    }

    void createCutCells()
    {
      cutCellsEnd_ = cutCells_.begin();
      auto it = volumeSnippets_.begin();

      assert(cutCells_.size() >= volumeRanges_.size());

      for (auto volumeEnd = volumeRanges_.begin(); volumeEnd != volumeRanges_.end(); ++volumeEnd) {
        int domainIndex = snippetIndexToDomain_[it->index()];

        if (domainIndex > -1) {
          // TODO: Instead of the bounding box we could also insert the fundamental element here
          // or whatever else we want
          cutCellsEnd_->geometry_ = BoundingBox(fundamentalElement_, domainToBoundingBox_[domainIndex].lower,
                                                domainToBoundingBox_[domainIndex].higher, CuboidMode::corners);
          cutCellsEnd_->fundamentalGeometry_ = &fundamentalGeometry_;
          cutCellsEnd_->domainIndex_ = domainIndex;
          cutCellsEnd_->snippetsBegin_ = it;
          cutCellsEnd_->snippetsEnd_ = *volumeEnd;

          ++cutCellsEnd_;
        }

        it = *volumeEnd;
      }
    }

    typename std::vector<CutCell>::iterator cutCellsBegin()
    {
      return cutCells_.begin();
    }

    typename std::vector<CutCell>::iterator cutCellsEnd()
    {
      return cutCellsEnd_;
    }

    void createCutIntersections()
    {
      cutIntersectionsEnd_ = cutIntersections_.begin();

      // TODO: The old code had one snippet per part, which resulted in one part per internal intersection
      // Does it make sense to group them?
      for (auto it = internalInterfaceSnippets_.begin(); it != internalInterfaceSnippets_.end(); ++it) {
        auto interiorDomain
          = domainConfiguration_.findInteriorDomain(internalInterfaceRelativePositions_[it->index()]);
        auto exteriorDomain
          = domainConfiguration_.findExteriorDomain(internalInterfaceRelativePositions_[it->index()]);

        if (interiorDomain == domainConfiguration_.domainsEnd()) {
          if (exteriorDomain != domainConfiguration_.domainsEnd()
              && cutCellInformation_.cutCellsExist(fundamentalElement_, exteriorDomain->index())) {
            // make sure domain with smaller index is exterior
            const auto& exteriorInformation
                = cutCellInformation_.information(fundamentalElement_, exteriorDomain->index());

            cutIntersectionsEnd_->globalGeometry_ = exteriorInformation.boundingBox;
            cutIntersectionsEnd_->geometryInInside_ = exteriorInformation.boundingBox;
            cutIntersectionsEnd_->snippetsBegin_ = it;
            cutIntersectionsEnd_->snippetsEnd_ = it + 1;
            cutIntersectionsEnd_->insideDomainIndex_ = exteriorDomain->index();
            cutIntersectionsEnd_->outsideDomainIndex_ = -1;
            cutIntersectionsEnd_->unitOuterNormal_ = it->globalUnitOuterNormal();
            ++cutIntersectionsEnd_;
          }
        } else if (cutCellInformation_.cutCellsExist(fundamentalElement_, interiorDomain->index())) {
          if (exteriorDomain == domainConfiguration_.domainsEnd()) {
            // make sure domain with smaller index is interior
            const auto& interiorInformation
                = cutCellInformation_.information(fundamentalElement_, interiorDomain->index());

            cutIntersectionsEnd_->globalGeometry_ = interiorInformation.boundingBox;
            cutIntersectionsEnd_->geometryInInside_ = interiorInformation.boundingBox;
            cutIntersectionsEnd_->snippetsBegin_ = it;
            cutIntersectionsEnd_->snippetsEnd_ = it + 1;
            cutIntersectionsEnd_->insideDomainIndex_ = interiorDomain->index();
            cutIntersectionsEnd_->outsideDomainIndex_ = -1;
            cutIntersectionsEnd_->unitOuterNormal_ = it->globalUnitOuterNormal();
            ++cutIntersectionsEnd_;
          } else if (cutCellInformation_.cutCellsExist(fundamentalElement_, exteriorDomain->index())) {
            // make sure domain with smaller index is interior
            if (interiorDomain->index() > exteriorDomain->index()) {
              std::swap(interiorDomain, exteriorDomain);
            }

            const auto& interiorInformation
                = cutCellInformation_.information(fundamentalElement_, interiorDomain->index());
            const auto& exteriorInformation
                = cutCellInformation_.information(fundamentalElement_, exteriorDomain->index());

            cutIntersectionsEnd_->globalGeometry_ = exteriorInformation.boundingBox;
            cutIntersectionsEnd_->geometryInInside_ = interiorInformation.boundingBox;
            cutIntersectionsEnd_->geometryInOutside_ = exteriorInformation.boundingBox;
            cutIntersectionsEnd_->snippetsBegin_ = it;
            cutIntersectionsEnd_->snippetsEnd_ = it + 1;
            cutIntersectionsEnd_->insideDomainIndex_ = interiorDomain->index();
            cutIntersectionsEnd_->outsideDomainIndex_ = exteriorDomain->index();
            cutIntersectionsEnd_->unitOuterNormal_ = it->globalUnitOuterNormal();
            ++cutIntersectionsEnd_;
          }
        }
      }

      auto it = externalInterfaceSnippets_.begin();

      for (auto interfaceEnd = externalInterfaceRanges_.begin(); interfaceEnd != externalInterfaceRanges_.end(); ++interfaceEnd) {
        auto domainIt = domainConfiguration_.findDomain(externalInterfaceRelativePositions_[it->index()]);

        if (domainIt != domainConfiguration_.domainsEnd()) {
          for (auto snippetIt = it; snippetIt != *interfaceEnd; ++snippetIt) {
            const auto& interiorInformation = cutCellInformation_.information(fundamentalElement_, domainIt->index());

            if (snippetIt->outsideFundamentalElementIndex() == elementMapper_.index(fundamentalElement_)) {
              cutIntersectionsEnd_->globalGeometry_ = interiorInformation.boundingBox;
              cutIntersectionsEnd_->geometryInInside_ = interiorInformation.boundingBox;
              cutIntersectionsEnd_->snippetsBegin_ = snippetIt;
              cutIntersectionsEnd_->snippetsEnd_ = snippetIt + 1;
              cutIntersectionsEnd_->insideDomainIndex_ = domainIt->index();
              cutIntersectionsEnd_->outsideDomainIndex_ = -1;
              cutIntersectionsEnd_->unitOuterNormal_ = snippetIt->globalUnitOuterNormal();
              ++cutIntersectionsEnd_;
            } else {
              const auto& exteriorInformation = cutCellInformation_.information(
                  snippetIt->outsideFundamentalElementIndex(), domainIt->index());

              cutIntersectionsEnd_->globalGeometry_ = exteriorInformation.boundingBox;
              cutIntersectionsEnd_->geometryInInside_ = interiorInformation.boundingBox;
              cutIntersectionsEnd_->geometryInOutside_ = exteriorInformation.boundingBox;
              cutIntersectionsEnd_->snippetsBegin_ = snippetIt;
              cutIntersectionsEnd_->snippetsEnd_ = snippetIt + 1;
              cutIntersectionsEnd_->insideDomainIndex_ = domainIt->index();
              cutIntersectionsEnd_->outsideDomainIndex_ = domainIt->index();
              cutIntersectionsEnd_->unitOuterNormal_ = snippetIt->globalUnitOuterNormal();
              ++cutIntersectionsEnd_;
            }
          }
        }

        it = *interfaceEnd;
      }
    }

    typename std::vector<CutIntersection>::iterator cutIntersectionBegin()
    {
      return cutIntersections_.begin();
    }

    typename std::vector<CutIntersection>::iterator cutIntersectionEnd()
    {
      return cutIntersectionsEnd_;
    }

    const std::vector<VolumeSnippet>& volumeSnippets()
    {
      return volumeSnippets_;
    }

    typename std::vector<VolumeSnippet>::iterator volumeSnippetsBegin()
    {
      return volumeSnippets_.begin();
    }

    typename std::vector<VolumeSnippet>::iterator volumeSnippetsEnd()
    {
      return volumeSnippetEnd_;
    }

    int snippetDomainIndex(const VolumeSnippet& snippet)
    {
      return snippetIndexToDomain_[snippet.index()];
    }

    bool hasCutCell(int domainIndex)
    {
      return std::find_if(volumeSnippets_.begin(), volumeSnippetEnd_, [this, domainIndex] (const VolumeSnippet& volumeSnippet) { return snippetIndexToDomain_[volumeSnippet.index()] > -1 && snippetIndexToDomain_[volumeSnippet.index()] == domainIndex; }) !=volumeSnippetEnd_;
    }

    bool fillsFundamentalCell(int domainIndex)
    {
      return std::distance(volumeSnippets_.begin(), volumeSnippetEnd_) == 1 && snippetIndexToDomain_[volumeSnippets_.front().index()] == domainIndex;
    }

    ctype cellVolume(int domainIndex)
    {
      return domainToVolume_[domainIndex];
    }

    Dune::FieldVector<ctype, dim> domainToCenterOfMass(int domainIndex) const
    {
      return domainToCenterOfMass_[domainIndex];
    }

    const RawBoundingBox<ctype, dim>& boundingBox(int domainIndex) const
    {
      return domainToBoundingBox_[domainIndex];
    }

    std::size_t numberOfVolumeSnippets(int domainIndex)
    {
      return std::count_if(volumeSnippets_.begin(), volumeSnippetEnd_, [this, domainIndex] (const auto& volumeSnippet) {
        return snippetIndexToDomain_[volumeSnippet.index()] == domainIndex;
      });
    }

    const GridView& gridView() const
    {
      return fundamentalGridView_;
    }

    const DomainConfiguration<GridView, LevelSetGridView>& domainConfiguration() const
    {
      return domainConfiguration_;
    }
    
    const CutCellInformation<GV, LGV>& cutCellInformation() const
    {
      return cutCellInformation_;
    }

  private:
    const DomainConfiguration<GridView, LevelSetGridView>& domainConfiguration_;

    const GridView& fundamentalGridView_;
    const LevelSetGridView& levelSetGridView_;

    VolumeRefinement volumeRefinement_;
    FaceRefinement faceRefinement_;

    Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> elementMapper_;
    FundamentalElement fundamentalElement_;
    typename FundamentalElement::Geometry fundamentalGeometry_;

    std::vector<CutCell> cutCells_;
    typename std::vector<CutCell>::iterator cutCellsEnd_;

    std::vector<CutIntersection> cutIntersections_;
    typename std::vector<CutIntersection>::iterator cutIntersectionsEnd_;

    std::vector<VolumeSnippet> volumeSnippets_;
    std::vector<InternalInterfaceSnippet> internalInterfaceSnippets_;
    // TODO: We want to use the same class for both snippet types
    // The additional unused data should still be justified
    // if we can have the identical handling for the quadrature rule setup
    // So we need to rename the type here
    // On the other hand, keeping two seperate lists for snippet construction
    // is probably still better
    std::vector<InternalInterfaceSnippet> externalInterfaceSnippets_;

    std::vector<VolumeSnippet> newVolumeSnippets_;
    std::vector<InternalInterfaceSnippet> newInternalInterfaceSnippets_;
    std::vector<InternalInterfaceSnippet> newExternalInterfaceSnippets_;

    typename std::vector<VolumeSnippet>::iterator volumeSnippetEnd_;
    typename std::vector<VolumeSnippet>::iterator newVolumeSnippetEnd_;

    typename std::vector<InternalInterfaceSnippet>::iterator externalInterfaceSnippetEnd_;
    typename std::vector<InternalInterfaceSnippet>::iterator newExternalInterfaceSnippetEnd_;

    std::vector<typename std::vector<VolumeSnippet>::iterator> volumeRanges_;
    std::vector<typename std::vector<InternalInterfaceSnippet>::iterator> internalInterfaceRanges_;
    std::vector<typename std::vector<InternalInterfaceSnippet>::iterator> externalInterfaceRanges_;

    // TODO: Replace this by something more efficient. For the moment it's fine if we preallocate everything.
    std::vector<std::vector<InterfaceRelativePosition>> volumeRelativePositions_;
    std::vector<std::vector<InterfaceRelativePosition>> internalInterfaceRelativePositions_;
    std::vector<std::vector<InterfaceRelativePosition>> externalInterfaceRelativePositions_;

    std::vector<ctype> functionValues_;

    std::vector<Dune::FieldVector<ctype, dim>> domainToCenterOfMass_;
    std::vector<RawBoundingBox<ctype, dim>> domainToBoundingBox_;
    std::vector<ctype> domainToVolume_;
    std::vector<int> snippetIndexToDomain_;
    std::vector<std::size_t> internalIntersectionToInterface_;

    std::vector<FieldVector<ctype, dim> > globalCornerCoordinates_;

    VolumeSnippetIdentityGeometry volumeSnippetIdentityGeometry_;
    VolumeSnippetCombineGeometries volumeSnippetCombineGeometries_;
    InternalInterfaceCombineGeometries internalInterfaceCombineGeometries_;
    ExternalInterfaceIdentityGeometry externalInterfaceIdentityGeometry_;
    ExternalInterfaceCombineGeometries externalInterfaceCombineGeometries_;
    ExternalInterfaceConvertGeometry externalInterfaceConvertGeometry_;

    CutCellInformation<GV, LGV>& cutCellInformation_;

    double valueTolerance_;

    int volumeIndexOffset_;
    std::vector<int> internalIntersectionsIndexOffsets_;
    std::size_t cutsPerIniternalInterface_;
    int externalIntersectionIndexOffset_;

    // TODO: This parameter is currently unused
    bool forceRefinement_;
  };

}

#endif
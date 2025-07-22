#ifndef DUNE_SUBTRIANGULATION_TPMC_TPMCREFINEMENT_HH
#define DUNE_SUBTRIANGULATION_TPMC_TPMCREFINEMENT_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/utility/typefromvertexcount.hh>

#include <dune/subtriangulation/common/geometrydefinitions.hh>
#include <dune/subtriangulation/tpmc/isdegenerated.hh>

#include <tpmc/fieldtraits.hh>
#include <tpmc/marchingcubes.hh>
#include <tpmc/thresholdfunctor.hh>

namespace tpmc {

template <class T, int dim>
struct FieldTraits<Dune::FieldVector<T, dim>>
{
  typedef T field_type;
};

} // namespace tpmc

namespace Dune::SubTriangulation {

template <class ctype, int dim, class TF = tpmc::ThresholdFunctor<ctype>>
class TpmcRefinement
{
  using Coordinate = FieldVector<ctype, dim>;
  using ReferenceElement
      = decltype(ReferenceElements<ctype, dim>::general(GeometryTypes::cube(dim)));

public:
  using VolumeGeometryType = MultiLinearGeometry<
      ctype, dim, dim, Dune::SubTriangulation::ReservedStorageMultiLinearGeometryTraits<ctype>>;
  using InterfaceGeometryType = MultiLinearGeometry<
      ctype, std::max(dim - 1, 0), dim,
      Dune::SubTriangulation::ReservedStorageMultiLinearGeometryTraits<ctype>>;
  using VolumeGeometryContainer    = std::vector<VolumeGeometryType>;
  using InterfaceGeometryContainer = std::vector<InterfaceGeometryType>;
  using const_volume_iterator      = typename VolumeGeometryContainer::const_iterator;
  using const_interface_iterator   = typename InterfaceGeometryContainer::const_iterator;

  TpmcRefinement(TF thresholdFunctor               = TF(),
                 tpmc::AlgorithmType algorithmType = tpmc::AlgorithmType::fullTPMC)
      : marchingCubes(algorithmType, thresholdFunctor)
      , referenceElement(ReferenceElements<ctype, dim>::general(GeometryTypes::cube(dim)))
  {}

  template <class Range>
  void bind(Range&& values)
  {
    bind(std::begin(values), std::end(values));
  }

  template <class I>
  void bind(I valuesBegin, I valuesEnd)
  {
    vertices.clear();
    interiorGeometries.clear();
    exteriorGeometries.clear();
    interfaceGeometries.clear();

    key              = marchingCubes.getKey(valuesBegin, valuesEnd);
    tpmcGeometryType = tpmc::makeGeometryType(dim, std::distance(valuesBegin, valuesEnd));
    marchingCubes.getVertices(valuesBegin, valuesEnd, key, std::back_inserter(vertices));
  }

  const_interface_iterator beginInterface() const { return interfaceGeometries.begin(); }

  const_interface_iterator endInterface() const { return interfaceGeometries.end(); }

  const_volume_iterator beginVolume(tpmc::ReconstructionType type) const
  {
    return type == tpmc::ReconstructionType::InteriorDomain ? interiorGeometries.begin()
                                                            : exteriorGeometries.begin();
  }

  const_volume_iterator endVolume(tpmc::ReconstructionType type) const
  {
    return type == tpmc::ReconstructionType::InteriorDomain ? interiorGeometries.end()
                                                            : exteriorGeometries.end();
  }

  IteratorRange<const_interface_iterator> interface() const
  {
    return IteratorRange<const_interface_iterator>(beginInterface(), endInterface());
  }
  IteratorRange<const_volume_iterator> volume(tpmc::ReconstructionType type) const
  {
    return IteratorRange<const_volume_iterator>(beginVolume(type), endVolume(type));
  }

  void reconstructInterface()
  {
    if (interfaceGeometries.empty()) {
      // retrieve elements
      std::vector<std::vector<int>> elements;
      marchingCubes.getElements(tpmcGeometryType, key, tpmc::ReconstructionType::Interface,
                                std::back_inserter(elements));

      // transformation of an element given as coordinate numbers to an InterfaceGeometry
      for (const std::vector<int>& element : elements) {
        using namespace std::placeholders;

        typename Dune::SubTriangulation::ReservedStorageMultiLinearGeometryTraits<
            ctype>::template CornerStorage<std::max(dim - 1, 0), dim>::Type coordinates;
        coordinates.resize(element.size());
        std::transform(element.begin(), element.end(), coordinates.begin(),
                       std::bind(&TpmcRefinement::transformCoordinate, this, _1));

        if (!IsDegenerated<ctype, dim - 1>::check(coordinates)) {
          GeometryType gt = geometryTypeFromVertexCount(dim - 1, coordinates.size());
          interfaceGeometries.push_back(InterfaceGeometryType(gt, coordinates));
        }
      }
    }
  }

  void reconstructVolume(tpmc::ReconstructionType type)
  {
    VolumeGeometryContainer& container = type == tpmc::ReconstructionType::InteriorDomain
                                             ? interiorGeometries
                                             : exteriorGeometries;

    if (container.empty()) {
      // retrieve elements
      std::vector<std::vector<int>> elements;
      marchingCubes.getElements(tpmcGeometryType, key, type, std::back_inserter(elements));

      for (const std::vector<int>& element : elements) {
        using namespace std::placeholders;

        typename Dune::SubTriangulation::ReservedStorageMultiLinearGeometryTraits<
            ctype>::template CornerStorage<dim, dim>::Type coordinates;
        coordinates.resize(element.size());
        std::transform(element.begin(), element.end(), coordinates.begin(),
                       std::bind(&TpmcRefinement::transformCoordinate, this, _1));

        if (!IsDegenerated<ctype, dim>::check(coordinates)) {
          GeometryType gt = geometryTypeFromVertexCount(dim, coordinates.size());
          container.push_back(VolumeGeometryType(gt, coordinates));
        }
      }
    }
  }

private:
  Coordinate transformCoordinate(int c) const
  {
    return c < 0 ? referenceElement.position(-c - 1, dim) : vertices[c];
  }

  tpmc::MarchingCubes<ctype, dim, Coordinate, TF> marchingCubes;
  std::size_t key;
  tpmc::GeometryType tpmcGeometryType;
  const ReferenceElement referenceElement;
  std::vector<Coordinate> vertices;
  VolumeGeometryContainer interiorGeometries;
  VolumeGeometryContainer exteriorGeometries;
  InterfaceGeometryContainer interfaceGeometries;
};

} // namespace Dune::SubTriangulation

#endif // DUNE_SUBTRIANGULATION_TPMC_TPMCREFINEMENT_HH

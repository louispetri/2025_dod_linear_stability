#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MONOMIALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MONOMIALBASIS_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/monomial.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune
{

namespace Functions
{

template<class K, int dim, int order>
class MonomialLocalFiniteElementCache;

template<class K, int order>
class MonomialLocalFiniteElementCache<K, 2, order>
{
public:
  using LocalFiniteElement = MonomialLocalFiniteElement<K, K, 2, order>;

  MonomialLocalFiniteElementCache() :
    quadrilateralElement(GeometryTypes::cube(2)),
    triangleElement(GeometryTypes::simplex(2))
  {
  }

  const LocalFiniteElement* get(const GeometryType& type) const
  {
    if (type.isQuadrilateral()) {
      return &quadrilateralElement;
    }
    else if (type.isTriangle()) {
      return &triangleElement;
    }
    else {
      DUNE_THROW(Dune::NotImplemented, "Unsupported geometry type");
    }
  }

private:
  LocalFiniteElement quadrilateralElement;
  LocalFiniteElement triangleElement;
};

template<class GV, int k, class K>
class MonomialNode : public LeafBasisNode
{
  static const int dim = GV::dimension;

  using FiniteElementCache = MonomialLocalFiniteElementCache<K, dim, k>;

  static const FiniteElementCache& finiteElementCache()
  {
    static FiniteElementCache instance;
    return instance;
  }


public:
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::LocalFiniteElement;

  MonomialNode() :
    finiteElement_(nullptr),
    element_(nullptr)
  {

  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = finiteElementCache().get(element_->type());
    this->setSize(finiteElement_->size());
  }

  const Element& element() const
  {
    return *element_;
  }

  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

protected:
  const FiniteElement* finiteElement_;
  const Element* element_;
};

template<class GV, int k, class K>
class MonomialPreBasis
{
public:
  static const int dim = GV::dimension;

private:
  using FiniteElementCache = MonomialLocalFiniteElementCache<K, dim, k>;

public:
  using GridView = GV;
  using size_type = std::size_t;
  using LocalFiniteElement = typename FiniteElementCache::LocalFiniteElement;
  using Node = MonomialNode<GridView, k, K>;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  /** \brief Constructor for a given grid view object */
  MonomialPreBasis(const GridView& gv) :
    gv_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gv_;
  }

  void update(const GridView& gv)
  {
    gv_ = gv;
  }

  size_type size() const
  {
    switch (dim)
    {
      case 1:
        return maxNodeSize() * gv_.size(0);
      case 2:
        return maxNodeSize() * (gv_.size(Dune::GeometryTypes::triangle) + gv_.size(Dune::GeometryTypes::quadrilateral));
      case 3:
        return maxNodeSize() * (gv_.size(Dune::GeometryTypes::tetrahedron)
             + gv_.size(Dune::GeometryTypes::pyramid)
             + gv_.size(Dune::GeometryTypes::prism)
             + gv_.size(Dune::GeometryTypes::hexahedron));
    }

    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  size_type dimension() const
  {
    return size();
  }

  static constexpr size_type maxNodeSize()
  {
    return LocalFiniteElement::Traits::LocalBasisType::size();
  }

  /**
  * \brief Create tree node
  */
  Node makeNode() const
  {
    return Node();
  }

  void initializeIndices()
  {
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& gridIndexSet = gridView().indexSet();
    const auto& element = node.element();

    std::size_t offset = gridIndexSet.index(element) * maxNodeSize();

    for (size_type i = 0, end = node.size() ; i < end ; ++i, ++it) {
      *it = {offset + i};
    }

    return it;
  }

protected:
  GridView gv_;
};

namespace BasisFactory {

template<std::size_t k, typename R=double>
auto monomial()
{
  return [](const auto& gridView) {
    return MonomialPreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

} // end namespace BasisFactory

template<typename GV, int k, class K = double>
using MonomialBasis = DefaultGlobalBasis<MonomialPreBasis<GV, k, K> >;

} // end namepsace Functions
} // end namespace Dune

#endif
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBDOMAINBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBDOMAINBASIS_HH

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune
{
namespace Functions
{

template<class ST>
class CutCellIndexMap
{
public:
  using SubTriangulation = ST;
  using Element = typename SubTriangulation::Entity;

  CutCellIndexMap(const SubTriangulation& subTriangulation) :
    subTriangulation_(subTriangulation),
    gridView_(subTriangulation_.gridView()),
    indexSet_(gridView_.indexSet())
  {
    std::vector<std::size_t> currentIndices(subTriangulation.domainConfiguration().numberOfDomains());
    std::fill(currentIndices.begin(), currentIndices.end(), 0);

    indices_.resize(subTriangulation_.gridView().size(0) * subTriangulation.domainConfiguration().numberOfDomains());

    for (int domainIndex = 0; domainIndex < subTriangulation.domainConfiguration().numberOfDomains(); ++domainIndex) {
      for (const auto& element : elements(gridView_)) {
        if (subTriangulation_.isHostCell(element, domainIndex)) {
          indices_[domainIndex * subTriangulation_.gridView().size(0) + indexSet_.index(element)] = currentIndices[domainIndex];
          ++currentIndices[domainIndex];
        }
        else {
          indices_[domainIndex * subTriangulation_.gridView().size(0) + indexSet_.index(element)] = std::numeric_limits<std::size_t>::max();
        }
      }
    }
  }

  bool hasDofs(const Element& element, int domainIndex) const
  {
    return indices_[domainIndex * gridView_.size(0) + indexSet_.index(element)] < std::numeric_limits<std::size_t>::max();
  }

  std::size_t index(const Element& element, int domainIndex) const
  {
    return indices_[domainIndex * gridView_.size(0) + indexSet_.index(element)];
  }

  std::size_t numberOfCutCells(std::size_t domainIndex) const
  {
    return subTriangulation_.cutCellInformation().numberOfCutCells(domainIndex);
  }

private:
  const SubTriangulation& subTriangulation_;
  const typename SubTriangulation::GridView& gridView_;
  const typename SubTriangulation::GridView::IndexSet& indexSet_;
  std::vector<std::size_t> indices_;
};

template<class Node, class CutCellIndexMap>
class SubDomainNode : public Node
{
public:
  using size_type = typename Node::size_type;
  using Element = typename Node::Element;
  using FiniteElement = typename Node::FiniteElement;

  SubDomainNode(const CutCellIndexMap* indexMap, int domainIndex) :
    indexMap_(indexMap),
    domainIndex_(domainIndex)
  {
  }

  void bind(const Element& element)
  {
    if (indexMap_->hasDofs(element, domainIndex_)) {
      Node::bind(element);
    }
    else {
      this->element_ = nullptr;
      this->finiteElement_ = nullptr;
      this->setSize(0);
    }
  }

private:
  const CutCellIndexMap* indexMap_;
  int domainIndex_;
};

template<class PreBasis, class CutCellIndexMap>
class SubDomainPreBasis : public PreBasis
{
public:
  static const int dim = PreBasis::dim;

  using GridView = typename PreBasis::GridView;
  using Node = SubDomainNode<typename PreBasis::Node, CutCellIndexMap>;
  using LocalFiniteElement = typename PreBasis::LocalFiniteElement;
  using size_type = typename PreBasis::size_type;

  static constexpr size_type maxMultiIndexSize = PreBasis::maxMultiIndexSize;
  static constexpr size_type minMultiIndexSize = PreBasis::minMultiIndexSize;
  static constexpr size_type multiIndexBufferSize = PreBasis::multiIndexBufferSize;

  SubDomainPreBasis(const GridView& gv, const CutCellIndexMap* indexMap, int domainIndex) :
    PreBasis(gv),
    indexMap_(indexMap),
    domainIndex_(domainIndex)
  {

  }

  Node makeNode() const
  {
    return Node(indexMap_, domainIndex_);
  }

  size_type size() const
  {
    // Assumes that there are no variable size finite elements
    return indexMap_->numberOfCutCells(domainIndex_) * this->maxNodeSize();
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

  void initializeIndices()
  {
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    if (node.size() == 0) {
      return it;
    }

    const auto& element = node.element();
    std::size_t offset = indexMap_->index(element, domainIndex_) * this->maxNodeSize();

    for (size_type i = 0, end = node.size() ; i < end ; ++i, ++it) {
      *it = {offset + i};
    }

    return it;
  }

private:
  const CutCellIndexMap* indexMap_;
  int domainIndex_;
};

namespace BasisFactory {

template<class PreBasis, class CutCellIndexMap>
auto subdomain(const CutCellIndexMap* cutCellIndexMap, int domainIndex)
{
  return [=](const auto& gridView) {
    return SubDomainPreBasis<PreBasis, CutCellIndexMap>(gridView, cutCellIndexMap, domainIndex);
  };
}

} // end namespace BasisFactory

template<class PreBasis, class CutCellIndexMap>
using SubDomainBasis = DefaultGlobalBasis<SubDomainPreBasis<PreBasis, CutCellIndexMap> >;

} // namespace Functions
} // namespace Dune

#endif
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIDOMAINBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIDOMAINBASIS_HH

#include <dune/common/reservedvector.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/multidomainlocalview.hh>

namespace Dune::Functions
{
  
template<class IMS, class SPB>
class MultiDomainBasis
{
  static const bool isBlocked = std::is_same_v<IMS, BasisFactory::BlockedLexicographic> or std::is_same_v<IMS, BasisFactory::BlockedInterleaved>;

public:

  //! Pre-basis providing the implementation details
  using SubPreBasis = SPB;

  //! The grid view that the FE space is defined on
  using GridView = typename SubPreBasis::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  static constexpr size_type maxMultiIndexSize = SubPreBasis::maxMultiIndexSize + isBlocked;
  static constexpr size_type minMultiIndexSize = SubPreBasis::minMultiIndexSize + isBlocked;
  static constexpr size_type multiIndexBufferSize = SubPreBasis::multiIndexBufferSize + isBlocked;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = MultiDomainLocalView<MultiDomainBasis<IMS, SubPreBasis>>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename LocalView::MultiIndex;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<std::size_t, multiIndexBufferSize>;

  using IndexMergingStrategy = IMS;

  using Node = typename SubPreBasis::Node;

  /**
   * \brief Constructor from a PreBasis factory
   *
   * \param gridView  The GridView this basis is based on
   * \param factory  A factory functor that gets the `gridView` and returns a `PreBasis`
   */
  template<class PreBasisFactory,
    std::enable_if_t<Dune::IsCallable<PreBasisFactory(GridView, int), SubPreBasis>::value, int> = 0>
  MultiDomainBasis(const GridView& gridView, PreBasisFactory&& factory, int numberOfDomains)
  {
    static_assert(models<Concept::PreBasis<GridView>, SubPreBasis>(), "Type passed to DefaultGlobalBasis does not model the PreBasis concept.");

    for (int i = 0; i < numberOfDomains; ++i) {
      subPreBases_.push_back(factory(gridView, i));
      subPreBases_.back().initializeIndices();
    }
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return subPreBases_[0].gridView();
  }

  //! Obtain the pre-basis providing the implementation details
  const SubPreBasis& subPreBasis(int domainIndex) const
  {
    return subPreBases_[domainIndex];
  }

  //! Obtain the pre-basis providing the implementation details
  SubPreBasis& subPreBasis(int domainIndex)
  {
    return subPreBases_[domainIndex];
  }

  /**
   * \brief Update the stored grid view
   *
   * This will update the indexing information of the global basis.
   * It must be called if the grid has changed.
   */
  void update(const GridView & gv)
  {
    for (std::size_t i = 0; i < subPreBases_.size(); ++i) {
      subPreBases_[i].update(gv);
      subPreBases_[i].initializeIndices();
    }
  }

  int numberOfDomains() const
  {
    return subPreBases_.size();
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return std::accumulate(subPreBases_.begin(), subPreBases_.end(), 0,
      [] (size_type dim, const auto& subPreBasis) { return dim + subPreBasis.dimension(); });
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size(SizePrefix{});
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return size(prefix, IndexMergingStrategy{});
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it, int domainIndex) const
  {
    return indices(node, it, domainIndex, IndexMergingStrategy{});
  }

  //! Return local view for basis in subdomain
  LocalView localView(int domainIndex) const
  {
    return LocalView(*this, domainIndex);
  }

  //! Return *this because we are not embedded in a larger basis
  const MultiDomainBasis& rootBasis() const
  {
    return *this;
  }

private:
  template<class MultiIndex>
  static void multiIndexPopFront(MultiIndex& M)
  {
    for(std::size_t i=0; i < M.size()-1; ++i) {
      M[i] = M[i+1];
    }

    M.resize(M.size() - 1);
  }

  size_type size(SizePrefix prefix, BasisFactory::BlockedLexicographic) const
  {
    if (prefix.size() == 0) {
      return subPreBases_.size();
    }

    auto front = prefix.front();
    multiIndexPopFront(prefix);

    return subPreBases_[front].size(prefix);
  }

  size_type size(SizePrefix prefix, BasisFactory::FlatLexicographic) const
  {
    if (prefix.size() == 0) {
      return std::accumulate(subPreBases_.begin(), subPreBases_.end(), 0,
        [](size_type s, const auto& subPreBasis) { return s + subPreBasis.size(); });
    }
    else {
      for (const auto& subPreBasis : subPreBases_) {
        auto firstDigitSize = subPreBasis.size();

        if (prefix[0] < firstDigitSize) {
          return subPreBasis.size(prefix);
        }

        prefix[0] -= firstDigitSize;
      }
    }

    // prefix out of range
    return -1;
  }

  template<class MultiIndex>
  static void multiIndexPushFront(MultiIndex& M, size_type M0)
  {
    M.resize(M.size()+1);
    for(std::size_t i=M.size()-1; i>0; --i) {
      M[i] = M[i-1];
    }

    M[0] = M0;
  }

  template<typename It>
  It indices(const Node& node, It multiIndices, int domainIndex, BasisFactory::BlockedLexicographic) const
  {
    size_type subTreeSize = node.size();
    subPreBases_[domainIndex].indices(node, multiIndices);
  
    for (std::size_t i = 0; i < subTreeSize; ++i) {
      this->multiIndexPushFront(multiIndices[i], domainIndex);
    }
  
    // we are advancing the iterator just in case
    // as of now, this method will never be called from a parent node
    multiIndices += subTreeSize;

    return multiIndices;
  }

  template<typename It>
  It indices(const Node& node, It multiIndices, int domainIndex, BasisFactory::FlatLexicographic) const
  {
    size_type firstComponentOffset = 0;

    for (int i = 0; i < domainIndex; ++i) {
      firstComponentOffset += subPreBases_[i].size();
    }

    subPreBases_[domainIndex].indices(node, multiIndices);

    size_type subTreeSize = node.size();

    for (std::size_t i = 0; i < subTreeSize; ++i) {
      multiIndices[i][0] += firstComponentOffset;
    }

    // we are advancing the iterator just in case
    // as of now, this method will never be called from a parent node
    multiIndices += subTreeSize;
    return multiIndices;
  }

protected:
  std::vector<SubPreBasis> subPreBases_;
};

namespace BasisFactory {

template<class GridView, class ChildPreBasisFactory, class IndexMergingStrategy>
auto makeMultiDomainBasis(const GridView& gridView, ChildPreBasisFactory&& childPreBasisFactory, int numberOfDomains, const IndexMergingStrategy&)
{
  return MultiDomainBasis<IndexMergingStrategy, decltype(childPreBasisFactory(gridView, 0))>(gridView, std::forward<ChildPreBasisFactory>(childPreBasisFactory), numberOfDomains);
};

template<class GridView, class ChildPreBasisFactory>
auto makeMultiDomainBasis(const GridView& gridView, ChildPreBasisFactory&& childPreBasisFactory, int numberOfDomains)
{
  return MultiDomainBasis<BlockedLexicographic, decltype(childPreBasisFactory(gridView, 0))>(gridView, std::forward<ChildPreBasisFactory>(childPreBasisFactory), numberOfDomains);
}

} // end namespace BasisFactory

} // namespace Dune::Functions

#endif
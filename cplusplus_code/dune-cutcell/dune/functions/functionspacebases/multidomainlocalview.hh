#ifndef DUNE_FUNCTIONS_MULTIDOMAINLOCALVIEW_HH
#define DUNE_FUNCTIONS_MULTIDOMAINLOCALVIEW_HH

#include <optional>
#include <type_traits>

#include <dune/common/reservedvector.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/common/multiindex.hh>
#include <dune/functions/common/overflowarray.hh>

namespace Dune::Functions
{

template<class GB>
class MultiDomainLocalView
{
public:

  //! The global FE basis that this is a view on
  using GlobalBasis = GB;

  //! The grid view the global FE basis lives on
  using GridView = typename GlobalBasis::GridView;

  //! Type of the grid element we are bound to
  using Element = typename GridView::template Codim<0>::Entity;

  //! The type used for sizes
  using size_type = std::size_t;

  //! Tree of local finite elements / local shape function sets
  using Tree = typename GlobalBasis::SubPreBasis::Node;

protected:

  using PreBasis = typename GlobalBasis::SubPreBasis;

  // Type used to store the multi indices of the basis vectors.
  // In contrast to MultiIndex this always has dynamic size.
  // It's guaranteed, that you can always cast it to MultiIndex
  using MultiIndexStorage =
      std::conditional_t<(GlobalBasis::minMultiIndexSize == GlobalBasis::maxMultiIndexSize),
        OverflowArray<StaticMultiIndex<size_type, GlobalBasis::maxMultiIndexSize>, GlobalBasis::multiIndexBufferSize>,
        Dune::ReservedVector<size_type, GlobalBasis::multiIndexBufferSize>>;

public:

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex =
      std::conditional_t<(GlobalBasis::minMultiIndexSize == GlobalBasis::maxMultiIndexSize),
        StaticMultiIndex<size_type, GlobalBasis::maxMultiIndexSize>,
        Dune::ReservedVector<size_type, GlobalBasis::multiIndexBufferSize>>;


  /** \brief Construct local view for a given global finite element basis */
  MultiDomainLocalView(const GlobalBasis& globalBasis, int domainIndex) :
    globalBasis_(&globalBasis),
    tree_(globalBasis_->subPreBasis(domainIndex).makeNode()),
    domainIndex_(domainIndex)
  {
    static_assert(models<Concept::BasisTree<GridView>, Tree>(), "Tree type passed to DefaultLocalView does not model the BasisNode concept.");
    initializeTree(tree_);
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = e;
    bindTree(tree_, *element_);
    indices_.resize(size());
    globalBasis_->indices(tree_, indices_.begin(), domainIndex_);
  }

  /** \brief Return if the view is bound to a grid element
   */
  bool bound() const
  {
    return static_cast<bool>(element_);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   */
  void unbind()
  {
    element_.reset();
  }

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Total number of degrees of freedom on this element
   */
  size_type size() const
  {
    return tree_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   */
  size_type maxSize() const
  {
    return globalBasis_->subPreBasis(domainIndex_).maxNodeSize();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex& index(size_type i) const
  {
    return indices_[i];
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

  /**
   * @brief Set the Domain Index object
   *
   * Since setting the view to a different domain invalidates the indices
   * we unbound it
   *
   * @param domainIndex
   */
  void setDomainIndex(int domainIndex)
  {
    domainIndex_ = domainIndex;
    unbind();
  }

protected:
  const GlobalBasis* globalBasis_;
  std::optional<Element> element_;
  Tree tree_;
  std::vector<MultiIndexStorage> indices_;

private:
  int domainIndex_;
};

} // namespace Dune::Functions

#endif
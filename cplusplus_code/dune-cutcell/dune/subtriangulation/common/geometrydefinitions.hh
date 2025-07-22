#ifndef DUNE_SUBTRIANGULATION_SIMPLETPMCTRIANGULATION_GEOMETRYDEFINITIONS_HH
#define DUNE_SUBTRIANGULATION_SIMPLETPMCTRIANGULATION_GEOMETRYDEFINITIONS_HH

#include <dune/common/math.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>

namespace Dune::SubTriangulation {

// ReservedStorageMultiLinearGeometryTraits
// -------------------------

/** \brief Our own variant of the traits class for MultiLinearGeometry
 *
 *  \tparam  ct  coordinate type
 */
template <class ct>
struct ReservedStorageMultiLinearGeometryTraits
{
  /** \brief helper structure containing some matrix routines
   *
   *  This helper allows exchanging the matrix inversion algorithms.
   *  It must provide the following static methods:
   *  \code
   *  template< int m, int n >
   *  static ctype sqrtDetAAT ( const FieldMatrix< ctype, m, n > &A );
   *
   *  template< int m, int n >
   *  static ctype rightInvA ( const FieldMatrix< ctype, m, n > &A,
   *                           FieldMatrix< ctype, n, m > &ret );
   *
   *  template< int m, int n >
   *  static void xTRightInvA ( const FieldMatrix< ctype, m, n > &A,
   *                            const FieldVector< ctype, n > &x,
   *                            FieldVector< ctype, m > &y );
   *  \endcode
   */
  typedef Impl::FieldMatrixHelper<ct> MatrixHelper;

  /** \brief tolerance to numerical algorithms */
  static ct tolerance() { return ct(16) * std::numeric_limits<ct>::epsilon(); }

  /** \brief template specifying the storage for the corners
   *
   *  Internally, the MultiLinearGeometry needs to store the corners of the
   *  geometry.
   *
   *  The corner storage may be chosen depending on geometry dimension and
   *  coordinate dimension. It is required to contain a type named Type, e.g.,
   *  \code
   *  template< int mydim, int cdim >
   *  struct CornerStorage
   *  {
   *    typedef std::vector< FieldVector< ctype, cdim > > Type;
   *  };
   *  \endcode
   *  By default, a std::vector of FieldVector is used.
   *
   *  Apart from being copy constructable and assignable, an \c const corner
   *  storage object \c corners must support the expressions \c
   *  begin(corners), \c end(corners), and subscription \c corners[i].  \c
   *  begin() and \c end() are looked up via ADL and in namespace \c std:
   *  \code
   *  using std::begin;
   *  using std::end;
   *  // it is a const_iterator over the corners in Dune-ordering
   *  auto it = begin(corners);
   *  FieldVector<ctype, cdim> c0 = *it;
   *  auto itend = end(corners);
   *  while(it != itend) {
   *    //...
   *  }
   *
   *  // elements must be accessible by subscription, indexed in
   *  // Dune-ordering
   *  FieldVector<ctype, cdim> c1 = corners[1];
   *  \endcode
   *  This means that all of the following qualify: \c
   *  FieldVector<ctype,cdim>[1<<mydim], \c
   *  std::array<FieldVector<ctype,cdim>,(1<<mydim)>, \c
   *  std::vector<FieldVector<ctype,cdim>>.
   *
   *  \note The expression \c end(corners) isn't actually used by the
   *        implementation currently, but we require it anyway so we can add
   *        runtime checks for the container size when we feel like it.
   *
   *  It is also possible to use a \c std::reference_wrapper of a suitable
   *  container as the type for the corner storage.  The implementation
   *  automatically calls \c corners.get() on internally stored \c
   *  std::reference_wrapper objects before applying \c begin(), \c end(),
   *  or subscription in that case.
   *
   *  \note Using \c std::reference_wrapper of some container as the corner
   *        storage means that the geometry has no control over the lifetime
   *        of or the access to that container.  When the lifetime of the
   *        container ends, or the container itself or its elements are
   *        modified, any geometry object that still references that
   *        container becomes invalid.  The only valid operation on invalid
   *        geometry objects are destruction and assignment from another
   *        geometry.  If invalidation happens concurrently with some
   *        operation (other than destruction or assignment) on the
   *        geometry, that is a race.
   *
   *  \tparam  mydim  geometry dimension
   *  \tparam  cdim   coordinate dimension
   */
  template <int mydim, int cdim>
  struct CornerStorage
  {
    // typedef Dune::ReservedVector< FieldVector< ct, cdim >, power(2, cdim) > Type;
    typedef Dune::ReservedVector<FieldVector<ct, cdim>, power(2, cdim)> Type;
  };

  /** \brief will there be only one geometry type for a dimension?
   *
   *  If there is only a single geometry type for a certain dimension,
   *  <em>hasSingleGeometryType::v</em> can be set to true.
   *  Supporting only one geometry type might yield a gain in performance.
   *
   *  If <em>hasSingleGeometryType::v</em> is set to true, an additional
   *  parameter <em>topologyId</em> is required.
   *  Here's an example:
   *  \code
   *  static const unsigned int topologyId = GeometryTypes::simplex(dim).id();
   *  \endcode
   */
  template <int dim>
  struct hasSingleGeometryType
  {
    static const bool v                  = false;
    static const unsigned int topologyId = ~0u;
  };
};
} // namespace Dune::SubTriangulation

#endif
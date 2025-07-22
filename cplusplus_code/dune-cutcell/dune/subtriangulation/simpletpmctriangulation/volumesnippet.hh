#ifndef DUNE_SUBTRIANGULATION_VOLUMESNIPPET_HH
#define DUNE_SUBTRIANGULATION_VOLUMESNIPPET_HH

#include <dune/subtriangulation/common/geometrydefinitions.hh>
#include <dune/subtriangulation/simpletpmctriangulation/interface.hh>

namespace Dune
{
  namespace SubTriangulation
  {

    /** \brief a snippet is a part of a cut cell
        \tparam GV the type of the fundamental mesh grid view
        \tparam LGV the type of the level set mesh grid view
     */
    template <class GV, class LGV>
    class VolumeSnippet
    {
    public:
      //! \brief type of the fundamental mesh grid view
      using FundamentalGridView = GV;
      using LevelSetGridView = LGV;
      using Element = typename LGV::template Codim<0>::Entity;
      using ctype = typename GV::ctype;
      enum { dim = GV::dimension };
      //! \brief geometry type of this snippet in the fundamental mesh cell
      using GeometryInFundamental = Dune::CachedMultiLinearGeometry<ctype, dim, dim,
          Dune::SubTriangulation::ReservedStorageMultiLinearGeometryTraits<ctype>>;

      VolumeSnippet(std::size_t index, const GeometryInFundamental& geometryInFundamental,
          Element homeElement)
          : index_(index)
          , geometryInFundamental_(geometryInFundamental)
          , homeElement_(homeElement)
      {
      }

      VolumeSnippet(std::size_t index, Dune::GeometryType type,
          const std::vector<FieldVector<ctype, dim> >& corners,
          Element homeElement)
          : index_(index)
          , geometryInFundamental_(type, corners)
          , homeElement_(homeElement)
      {
      }

      std::size_t index() const
      {
        return index_;
      }

      void addIndexOffset(std::size_t offset)
      {
        index_ += offset;
      }

      const GeometryInFundamental& geometryInFundamental() const
      {
        return geometryInFundamental_;
      }

      Element homeElement() const
      {
        return homeElement_;
      }

    private:
      std::size_t index_;
      GeometryInFundamental geometryInFundamental_;
      Element homeElement_;
    };
  }
}

#endif // DUNE_SubTriangulation_VOLUMESNIPPET_HH

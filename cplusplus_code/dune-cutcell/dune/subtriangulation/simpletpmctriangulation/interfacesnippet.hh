#ifndef DUNE_SUBTRIANGULATION_INTERFACESNIPPET_HH
#define DUNE_SUBTRIANGULATION_INTERFACESNIPPET_HH

#include <cstddef>

#include <dune/subtriangulation/common/geometrydefinitions.hh>

namespace Dune
{
  namespace SubTriangulation
  {
    template <class GV, class LGV>
    class InterfaceSnippetBase
    {
    public:
      using FundamentalGridView = GV;
      using LevelSetGridView = LGV;
      using LevelSetElement = typename LGV::template Codim<0>::Entity;
      using ctype = typename GV::ctype;
      enum { dim = GV::dimension };
      //! \brief geometry type of this snippet in the fundamental mesh cell
      using GlobalGeometry = Dune::CachedMultiLinearGeometry<ctype, dim - 1, dim,
          Dune::SubTriangulation::ReservedStorageMultiLinearGeometryTraits<ctype> >;

      InterfaceSnippetBase(const GlobalGeometry& globalGeometry,
          const FieldVector<ctype, dim>& globalUnitOuterNormal,
          const LevelSetElement& levelSetElement, std::size_t index, int outsideFundamentalElementIndex = -1)
          : globalGeometry_(globalGeometry)
          , globalUnitOuterNormal_(globalUnitOuterNormal)
          , levelSetElement_(levelSetElement)
          , index_(index)
          , outsideFundamentalElementIndex_(outsideFundamentalElementIndex)
      {
      }

      InterfaceSnippetBase(Dune::GeometryType type,
          const std::vector<FieldVector<ctype, dim> >& corners,
          const FieldVector<ctype, dim>& globalUnitOuterNormal,
          const LevelSetElement& levelSetElement, std::size_t index, int outsideFundamentalElementIndex = -1)
          : globalGeometry_(type, corners)
          , globalUnitOuterNormal_(globalUnitOuterNormal)
          , levelSetElement_(levelSetElement)
          , index_(index)
          , outsideFundamentalElementIndex_(outsideFundamentalElementIndex)
      {
      }

      const GlobalGeometry& globalGeometry() const
      {
        return globalGeometry_;
      }

      const FieldVector<ctype, dim>& globalUnitOuterNormal() const
      {
        return globalUnitOuterNormal_;
      }

      LevelSetElement levelSetElement() const
      {
        return levelSetElement_;
      }

      std::size_t index() const
      {
        return index_;
      }

      void addIndexOffset(std::size_t offset)
      {
        index_ += offset;
      }

      std::size_t outsideFundamentalElementIndex() const
      {
        return outsideFundamentalElementIndex_;
      }

    private:
      GlobalGeometry globalGeometry_;
      FieldVector<ctype, dim> globalUnitOuterNormal_;
      LevelSetElement levelSetElement_;
      std::size_t index_;
      std::size_t outsideFundamentalElementIndex_;
    };

    template <class GV, class LGV>
    using InternalInterfaceSnippet = InterfaceSnippetBase<GV, LGV>;

    template <class GV, class LGV>
    class ExternalInterfaceSnippet : public InterfaceSnippetBase<GV, LGV>
    {
    public:
      using BaseT = InterfaceSnippetBase<GV, LGV>;
      using FundamentalElement = typename GV::template Codim<0>::Entity;
      using GlobalGeometry = typename BaseT::GlobalGeometry;
      using ctype = typename BaseT::ctype;
      enum { dim = BaseT::dim };
      using LevelSetElement = typename BaseT::LevelSetElement;

      ExternalInterfaceSnippet(const GlobalGeometry& globalGeometry,
          const FieldVector<ctype, dim>& globalUnitOuterNormal,
          std::size_t outsideFundamentalElementIndex, int codim1IndexInInside,
          const LevelSetElement& levelSetElement, std::size_t index)
          : InterfaceSnippetBase<GV, LGV>(globalGeometry, globalUnitOuterNormal, levelSetElement, index, outsideFundamentalElementIndex)
          , codim1IndexInInside_(codim1IndexInInside)
      {
      }

      ExternalInterfaceSnippet(Dune::GeometryType type,
          const std::vector<FieldVector<ctype, dim> >& corners,
          const FieldVector<ctype, dim>& globalUnitOuterNormal,
          std::size_t outsideFundamentalElementIndex, int codim1IndexInInside,
          const LevelSetElement& levelSetElement, std::size_t index)
          : InterfaceSnippetBase<GV, LGV>(type, corners, globalUnitOuterNormal, levelSetElement, index, outsideFundamentalElementIndex)
          , codim1IndexInInside_(codim1IndexInInside)
      {
      }

      int codim1IndexInInside() const
      {
        return codim1IndexInInside_;
      }

    private:
      int codim1IndexInInside_;
    };
  }
}

#endif // DUNE_SubTriangulation_INTERFACESNIPPET_HH

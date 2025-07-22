#ifndef DUNE_SUBTRIANGULATION_INTERFACE_HH
#define DUNE_SUBTRIANGULATION_INTERFACE_HH

#include <sstream>
#include <functional>
#include <dune/common/fvector.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune
{
  namespace SubTriangulation
  {
    /** \brief a position relative to an interface */
    enum class InterfaceRelativePosition { interior, exterior, interface, any };

    inline std::string toString(const std::vector<InterfaceRelativePosition>& pos)
    {
      std::stringstream sstr;
      for (unsigned int i = 0; i < pos.size(); ++i) {
        sstr << (i > 0 ? " " : "");
        switch (pos[i]) {
        case InterfaceRelativePosition::interior:
          sstr << "interior";
          break;
        case InterfaceRelativePosition::exterior:
          sstr << "exterior";
          break;
        case InterfaceRelativePosition::interface:
          sstr << "interface";
          break;
        case InterfaceRelativePosition::any:
          sstr << "any";
          break;
        }
      }
      return sstr.str();
    }

    inline std::vector<int> toDomainTag(const std::vector<InterfaceRelativePosition>& pos)
    {
      std::vector<int> v;
      for (auto p : pos) {
        switch (p) {
        case InterfaceRelativePosition::interior:
          v.push_back(-1);
          break;
        case InterfaceRelativePosition::exterior:
          v.push_back(1);
          break;
        default:
          v.push_back(0);
          break;
        }
      }
      return v;
    }

    template <class LGV>
    class Interface
    {
    public:
      using ctype = typename LGV::ctype;
      enum { dim = LGV::dimension };
      using Domain = Dune::FieldVector<ctype, dim>;
      using Range = ctype;
      using LevelSetFunction = Dune::Functions::GridViewFunction<Range(Domain), LGV>;
      using Element = typename LGV::template Codim<0>::Entity;

      Interface(std::size_t index, LevelSetFunction func) : index_(index), function_(func)
      {
      }

      const LevelSetFunction& function() const
      {
        return function_;
      }

      std::size_t index() const
      {
        return index_;
      }

    private:
      std::size_t index_;
      mutable LevelSetFunction function_;
    };

    template <class LGV>
    bool operator==(const Interface<LGV>& a, const Interface<LGV>& b)
    {
      return a.index() == b.index();
    }

    template <class LGV>
    bool operator!=(const Interface<LGV>& a, const Interface<LGV>& b)
    {
      return !(a == b);
    }
  }
}

#endif

#ifndef DUNE_SUBTRIANGULATION_TPMC_ISDEGENERATED_HH
#define DUNE_SUBTRIANGULATION_TPMC_ISDEGENERATED_HH

#include <dune/common/fvector.hh>

namespace Dune::SubTriangulation {

template <class ctype>
struct IsDegeneratedBase
{
  template <int dim>
  static bool eq(const FieldVector<ctype, dim>& a, const FieldVector<ctype, dim>& b)
  {
    ctype sum = 0;

    for (int d = 0; d < dim; d++) {
      sum += std::abs(a[d] - b[d]);
    }

    return (sum <= eqEpsilon);
  }

  static constexpr ctype eqEpsilon = 1e-8;
};

template <class ctype, int dim>
struct IsDegenerated
{
  typedef FieldVector<ctype, dim> Point;

  template <int d>
  static bool check(std::vector<FieldVector<ctype, d> >& c)
  {
    return false;
  }
};

template <class ctype>
struct IsDegenerated<ctype, 0> : IsDegeneratedBase<ctype>
{
  enum
  {
    dim = 0
  };

  template <class Coordinates>
  static bool check(const Coordinates& c)
  {
    return false;
  }
};

template <class ctype>
struct IsDegenerated<ctype, 1> : IsDegeneratedBase<ctype>
{
  typedef IsDegeneratedBase<ctype> BaseT;

  enum
  {
    dim = 1
  };

  template <class Coordinates>
  static bool check(const Coordinates& c)
  {
    return BaseT::eq(c[0], c[1]);
  }
};

template <class ctype>
struct IsDegenerated<ctype, 2> : IsDegeneratedBase<ctype>
{
  typedef IsDegeneratedBase<ctype> BaseT;

  enum
  {
    dim = 2
  };

  template <class Coordinates>
  static bool check(const Coordinates& c)
  {
    switch (c.size()) {
      case 3:
        return checkTriangle(c[0], c[1], c[2]);
      case 4:
        if (BaseT::eq(c[0], c[1])) {
          return checkTriangle(c[0], c[2], c[3]);
        }
        if (BaseT::eq(c[0], c[2])) {
          return checkTriangle(c[0], c[1], c[3]);
        }
        if (BaseT::eq(c[0], c[3])) {
          return true;
        }
        if (BaseT::eq(c[1], c[2])) {
          return true;
        }
        if (BaseT::eq(c[1], c[3])) {
          return checkTriangle(c[0], c[1], c[2]);
        }
        if (BaseT::eq(c[2], c[3])) {
          return checkTriangle(c[0], c[1], c[2]);
        }
        return ((BaseT::eq(c[0], c[1]) && BaseT::eq(c[2], c[3]))
                || (BaseT::eq(c[0], c[2]) && BaseT::eq(c[1], c[3])));
      default:
        DUNE_THROW(Dune::Exception, "Impossible Geometry. "
                                        << "In 2D there is no known geometry with " << c.size()
                                        << " corners");
    }
  }

  template <class Coordinate>
  static bool checkTriangle(const Coordinate& c0, const Coordinate& c1, const Coordinate& c2)
  {
    return (BaseT::eq(c0, c1) || BaseT::eq(c0, c2) || BaseT::eq(c1, c2));
  }
};

template <class ctype>
struct IsDegenerated<ctype, 3> : IsDegeneratedBase<ctype>
{
  typedef IsDegeneratedBase<ctype> BaseT;
  enum
  {
    dim = 3
  };

  template <class Coordinates>
  static bool check(const Coordinates& c)
  {
    switch (c.size()) {
      case 4:
        return (BaseT::eq(c[0], c[1]) || BaseT::eq(c[0], c[2]) || BaseT::eq(c[0], c[3])
                || BaseT::eq(c[1], c[2]) || BaseT::eq(c[1], c[3]) || BaseT::eq(c[2], c[3]));
      case 5: {
        bool deg = false;
        deg |= BaseT::eq(c[3], c[0]);
        deg |= BaseT::eq(c[2], c[1]);
        int mat(0);
        mat += int(BaseT::eq(c[0], c[1]));
        mat += int(BaseT::eq(c[0], c[2]));
        mat += int(BaseT::eq(c[0], c[4]));
        mat += int(BaseT::eq(c[1], c[3]));
        mat += int(BaseT::eq(c[1], c[4]));
        mat += int(BaseT::eq(c[2], c[3]));
        mat += int(BaseT::eq(c[2], c[4]));
        mat += int(BaseT::eq(c[3], c[4]));
        return deg || mat > 1;
      }
      case 6:
        return ((BaseT::eq(c[0], c[3]) && BaseT::eq(c[1], c[4]) && BaseT::eq(c[2], c[5]))
                || (BaseT::eq(c[0], c[1]) && BaseT::eq(c[3], c[4]))
                || (BaseT::eq(c[0], c[2]) && BaseT::eq(c[3], c[5]))
                || (BaseT::eq(c[1], c[2]) && BaseT::eq(c[4], c[5])));
      case 8:
        return ((BaseT::eq(c[0], c[1]) && BaseT::eq(c[2], c[3]) && BaseT::eq(c[4], c[5])
                 && BaseT::eq(c[6], c[7]))
                || (BaseT::eq(c[0], c[2]) && BaseT::eq(c[1], c[3]) && BaseT::eq(c[4], c[6])
                    && BaseT::eq(c[5], c[7]))
                || (BaseT::eq(c[0], c[4]) && BaseT::eq(c[1], c[5]) && BaseT::eq(c[2], c[6])
                    && BaseT::eq(c[3], c[7])));
      default:
        DUNE_THROW(Dune::Exception, "Impossible Geometry. "
                                        << "In 3D there is no known geometry with " << c.size()
                                        << " corners");
    }
  }
};

} // namespace Dune::SubTriangulation

#endif // DUNE_SUBTRIANGULATION_TPMC_ISDEGENERATED_HH

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef CUBOID_HH
#define CUBOID_HH

#include <assert.h>
#include <iostream>

#include <dune/common/version.hh>
#if !DUNE_VERSION_GTE(DUNE_COMMON,2,4)
#include <dune/common/static_assert.hh>
#endif
#include <dune/grid/yaspgrid.hh>
#include <dune/geometry/multilineargeometry.hh>

namespace Dune
{

  /**
   * \brief Data structure to represent a raw bounding box.
   *
   * \note Freshly initialized objects of this type are in an
   *       invalid state.
   *
   * \todo Avoid raw bounding box objects in an invalid state,
   *       e.g. call update on initialization.
   */
  template <typename DomainField, int dim>
  struct RawBoundingBox
  {
    //! Type of bounding box corners
    typedef Dune::FieldVector<DomainField,dim> Coordinate;

    //! Initializes a new bounding box object
    RawBoundingBox() :
      lower(std::numeric_limits<DomainField>::max()),
      higher(-(std::numeric_limits<DomainField>::max()-1.0))
    {}

    //! Extend bounding box such that it contains a given point
    void update(const Coordinate & x){
      for(int d=0; d<dim; ++d) {
        lower[d] = std::min(x[d],lower[d]);
        higher[d] = std::max(x[d],higher[d]);
      }
    }

    //! Extend bounding box such that it contains a given bounding box
    void update(const RawBoundingBox & b){
      this->update(b.lower);
      this->update(b.higher);
    }

    //! Lower and higher corner defining the bounding box
    Coordinate lower; Coordinate higher;
  };


  /** Cuboid of full world dimension */
  struct CuboidMode
  {
    enum Mode {
      corners, size
    };
  };

  template<int dim, class ct>
  class Cuboid
  {
  public:
    Cuboid() :
      o(0.0), h(1.0) {};
    Cuboid(const Dune::FieldVector<ct, dim> & _o, ct _h) :
      o(_o), h(_h) {};
    Cuboid(const Dune::FieldVector<ct, dim> & _o,
           const Dune::FieldVector<ct, dim> & _h,
           CuboidMode::Mode m = CuboidMode::size) :
      o(_o), h(_h)
    {
      if (m==CuboidMode::corners) {
        for (int i=0; i<dim; i++)
        {
          o[i] = std::min(_o[i],_h[i]);
          h[i] = std::abs(_o[i]-_h[i]);
        }
      }
      jinv = computeJacobianInverseTransposed();
    };

    const Dune::FieldVector<ct, dim> global (Dune::FieldVector<ct, dim> v) const
    {
      for (int i=0; i<dim; ++i) v[i] = v[i] * h[i] + o[i];
      return v;
    }
    const Dune::FieldVector<ct, dim> local (Dune::FieldVector<ct, dim> v) const
    {
      for (int i=0; i<dim; ++i) v[i] = (v[i] - o[i]) / h[i];
      return v;
    }

    Dune::FieldVector<ct, dim> corner (const int c) const
    {
      const auto ref = Dune::ReferenceElements<ct,dim>::cube();
      return global(ref.position(c,dim));
    }

    int corners () const
    {
      const auto ref =
        Dune::ReferenceElements<ct,dim>::cube();
      return ref.size(dim);
    }


    const Dune::FieldVector<ct, dim> origin () const
    {
      return o;
    }

    const Dune::FieldVector<ct, dim> height () const
    {
      return h;
    }

    ct integrationElement(const Dune::FieldVector<ct,dim>& local) const
    {
      ct ie=1;
      for (int i=0; i<dim; i++) ie *= h[i];
      return std::abs( ie );
    }

    ct volume() const
    {
      return integrationElement(Dune::FieldVector<ct,dim>(0));
    }

    bool affine() const { return true; }

    const Dune::FieldMatrix<ct,dim,dim>&
    jacobianInverseTransposed ([[maybe_unused]] const Dune::FieldVector<ct,dim>& local) const
    {
      return jinv;
    }

    const Dune::FieldMatrix<ct,dim,dim>&
    jacobianTransposed (const Dune::FieldVector<ct,dim>& local) const
    {
      return jacobianInverseTransposed(local);
    }

    Dune::GeometryType type() const
    {
      return Dune::GeometryTypes::cube(dim);
    }

    bool operator < (const Cuboid & c) const
    {
      for (int i=0; i<dim; i++)
      {
        if (o[i] != c.o[i])
          return (o[i] < c.o[i]);
        if (h[i] != c.h[i])
          return (h[i] < c.h[i]);
      }
      return false;
    }

    bool operator == (const Cuboid & c) const
    {
      for (int i=0; i<dim; ++i)
        if ( o[i] != c.o[i] || h[i] != c.h[i] )
          return false;
      return true;
    }

    bool operator != (const Cuboid & c) const
    {
      return ! (*this == c);
    }

    void print(std::ostream & s) const
    {
      s << "(" << o << ")/(" << h << ")";
    }

  protected:
    Dune::FieldMatrix<ct,dim,dim> computeJacobianInverseTransposed() {
      Dune::FieldMatrix<ct,dim,dim> Jinv(0);
      for (int i=0; i<dim; ++i) Jinv[i][i] = 1.0 / h[i];
      return Jinv;
    }

    Dune::FieldVector<ct, dim> o;
    Dune::FieldVector<ct, dim> h;
    Dune::FieldMatrix<ct, dim, dim> jinv;
  };

  template <class Grid>
  class GridCuboid :
    public Cuboid<Grid::dimension, typename Grid::ctype> {
  public:
    enum { dim = Grid::dimension };
    enum { dimw = Grid::dimensionworld };

    static const int mydimension = Grid::dimension;
    static const int coorddimension = Grid::dimension;
    typedef typename Grid::ctype ct;
    using ctype = Grid::ctype;
    typedef Cuboid<dim,ct> Base;

    typedef typename Grid::template Codim<0>::Entity Entity;

    GridCuboid() :
      Cuboid<dim,ct>(), entity_(Entity())
    {
    }

    GridCuboid(const Entity & E,
               const Cuboid<dim,ct> & c) :
      Cuboid<dim,ct>(c), entity_(E)
    {
#if DUNE_VERSION_GTE(DUNE_COMMON,2,4)
      static_assert((int)dimw==(int)dim, "dim != dimw not supported");
#else
      dune_static_assert((int)dimw==(int)dim, "dim != dimw not supported");
#endif
    };
    GridCuboid(const Entity & E,
               const Dune::FieldVector<ct, dim> & _o, ct _h) :
      Cuboid<dim,ct>(_o,_h), entity_(E)
    {
#if DUNE_VERSION_GTE(DUNE_COMMON,2,4)
      static_assert((int)dimw==(int)dim, "dim != dimw not supported");
#else
      dune_static_assert((int)dimw==(int)dim, "dim != dimw not supported");
#endif
    };
    GridCuboid(const Entity & E,
               const Dune::FieldVector<ct, dim> & _o,
               const Dune::FieldVector<ct, dim> & _h,
               CuboidMode::Mode m = CuboidMode::size) :
      Cuboid<dim,ct>(_o,_h,m), entity_(E)
    {
#if DUNE_VERSION_GTE(DUNE_COMMON,2,4)
      static_assert((int)dimw==(int)dim, "dim != dimw not supported");
#else
      dune_static_assert((int)dimw==(int)dim, "dim != dimw not supported");
#endif
    };

    GridCuboid(const Entity & E) :
      Cuboid<dim,ct>(origin(E), hsize(E), CuboidMode::corners), entity_(E)
    {
#if DUNE_VERSION_GTE(DUNE_COMMON,2,4)
      static_assert((int)dimw==(int)dim, "Not implemented for dim!=dimw");
#else
      dune_static_assert((int)dimw==(int)dim, "Not implemented for dim!=dimw");
#endif
    };

    const Entity & entity() const { return entity_; }

    Dune::FieldVector<ct, dim> center() const {
      Dune::FieldVector<ct, dim> c(h);
      c*=0.5;
      c+=o;
      return c;
    }

    Dune::FieldVector<ct, dim> coord(int i) const {
      return
        this->global(Dune::ReferenceElements<ct,dim>::general(this->type()).
                     position(i,dim));
    }
    bool isInside (Dune::FieldVector<ct, dim> p) const {
      ct eps = std::numeric_limits<ct>::epsilon();
      bool val = true;
      for (int i=0; i<dim; i++)
        val *= ((o[i]+h[i] - p[i] > -eps)
                && (p[i] - o[i] > -eps));
      return val;
    }
  private:
    static Dune::FieldVector<ct, dim> origin(const Entity & E)
    {
      Dune::FieldVector<ct, dim> curro(9999999);
      Dune::FieldVector<ct, dim> currh(-9999999);
      Dune::GeometryType type = E.type();
      int corners = Dune::ReferenceElements<ct,dim>::general(type).size(dim);;
      for (int i=0; i<corners; i++)
      {
        Dune::FieldVector<ct, dim> pos =
          E.geometry().global(
            Dune::ReferenceElements<ct,dim>::general(type).position(i,dim));
        for (int d=0; d<dim; d++)
        {
          curro[d] = std::min(curro[d],pos[d]);
          currh[d] = std::max(currh[d],pos[d]);
        }
      }
      return curro;
    }
    static Dune::FieldVector<ct, dim> hsize(const Entity & E)
    {
      Dune::FieldVector<ct, dim> curro(9999999);
      Dune::FieldVector<ct, dim> currh(-9999999);
      Dune::GeometryType type = E.type();
      int corners = Dune::ReferenceElements<ct,dim>::general(type).size(dim);;
      for (int i=0; i<corners; i++)
      {
        Dune::FieldVector<ct, dim> pos =
          E.geometry().global(
            Dune::ReferenceElements<ct,dim>::general(type).position(i,dim));
        for (int d=0; d<dim; d++)
        {
          curro[d] = std::min(curro[d],pos[d]);
          currh[d] = std::max(currh[d],pos[d]);
        }
      }
      return currh;
    }
    Entity entity_;
    using Base::o;
    using Base::h;
  };

template <int d, class ct>
inline std::ostream & operator << (std::ostream & s,
                                   const Dune::Cuboid<d,ct> & c)
{
  c.print(s);
  return s;
}

} // end namespace Dune

#endif // CUBOID_HH

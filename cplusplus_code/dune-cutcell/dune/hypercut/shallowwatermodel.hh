#ifndef DUNE_HYPERCUT_SHALLOW_WATER_MODEL_HH
#define DUNE_HYPERCUT_SHALLOW_WATER_MODEL_HH

#include <cstddef>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

namespace Dune::Hypercut
{

template<int m, class K>
class ZeroBoundary
{
public:
ZeroBoundary(const Dune::ParameterTree& paramtree)
{  
}

ZeroBoundary()
{  
}

Dune::FieldVector<K, m> operator()(const Dune::FieldVector<K, m>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
{
    return Dune::FieldVector<K, m>(0.0);
}
};

template<int m, class K>
class OutgoingBoundaryCondition
{
public:
    OutgoingBoundaryCondition(Dune::ParameterTree& paramtree)
    {

    }

    OutgoingBoundaryCondition()
    {

    }

    Dune::FieldVector<K, m> operator()(const Dune::FieldVector<K, m>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
    {
        return u;
    }
};

template<int m, class K>
class ReflectingBoundaryCondition
{
public:
    ReflectingBoundaryCondition(Dune::ParameterTree& paramtree)
    {

    }

    ReflectingBoundaryCondition()
    {

    }

    Dune::FieldVector<K, m> operator()(const Dune::FieldVector<K, m>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
    {
        auto bc = u;

        Dune::FieldVector<K, 2> v(0.0);
        
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = u[i + 1];
        }

        auto normalComponent = (v * n) * n;
        v -= 2.0 * normalComponent;

        for (std::size_t i = 0; i < v.size(); ++i) {
            bc[i + 1] = v[i];
        }

        return bc;
    }
};

template<class K_>
class ShallowWaterModel
{
public:
    static const int dim = 2;
    static const int m = 3;

    using K = K_;

    using Flux = Dune::FieldMatrix<K, m, dim>;
    using Coordinate = Dune::FieldVector<K, dim>;
    using State = Dune::FieldVector<K, m>;

    static constexpr K Epsilon = 1e-10;
    static constexpr K g = 9.80665l;

    ShallowWaterModel(const Dune::ParameterTree& paramtree)
    {
    }

    ShallowWaterModel()
    {
    }

    Flux flux(const Coordinate& x, const State& u) const
    {
        Flux f(0.0);

        K v = u[1] / u[0];
        K w = u[2] / u[0];

        auto p = 0.5 * u[0] * u[0];

        f[0][0] = u[1];
        f[1][0] = u[1] * v + p;
        f[2][0] = u[1] * w;

        f[0][1] = u[2];
        f[1][1] = u[2] * v;
        f[2][1] = u[2] * w + p;

        return f;
    }

    State numericalFlux(const Coordinate& n,
                        const Dune::FieldVector<K, m>& ul,
                        const Dune::FieldVector<K, m>& ur,
                        const Coordinate& x) const
    {
        return LaxFriedrichFlux(n, ul, ur, x);
    }

    State LaxFriedrichFlux(const Coordinate& n,
                           const Dune::FieldVector<K, m>& ul,
                           const Dune::FieldVector<K, m>& ur,
                           const Coordinate& x) const
    {
        auto fl = flux(x, ul).transposed();
        auto fr = flux(x, ur).transposed();

        State fnl = n[0] * fl[0] + n[1] * fl[1];
        State fnr = n[0] * fr[0] + n[1] * fr[1];

        K lam = lambda(ul, ur, n, x);

        return 0.5 * (fnl + fnr) + 0.5 * lambda(ul, ur, n, x) * (ul - ur);
    }

    State centralFlux(const Coordinate& n,
                      const Dune::FieldVector<K, m>& ul,
                      const Dune::FieldVector<K, m>& ur,
                      const Coordinate& x) const
    {
        auto fl = flux(x, ul).transposed();
        auto fr = flux(x, ur).transposed();

        State fnl = n[0] * fl[0] + n[1] * fl[1];
        State fnr = n[0] * fr[0] + n[1] * fr[1];

        return 0.5 * (fnl + fnr);
    }

    State centralFluxHalf(const Coordinate& n,
                          const Dune::FieldVector<K, m>& u,
                          const Coordinate& x) const
    {
        auto f = flux(x, u).transposed();

        State fn = n[0] * f[0] + n[1] * f[1];

        return 0.5 * fn;
    }

    State analyticFlux(const Coordinate& n,
                       const Dune::FieldVector<K, m>& u,
                       const Coordinate& x) const
    {
        auto fu = flux(x, u).transposed();

        State fnu = n[0] * fu[0] + n[1] * fu[1];

        return fnu;
    }

    K lambda(const State& ul, const State& ur, const Coordinate& n, const Coordinate& x) const
    {
        auto eigvalsLeft = eigenvalues(ul, n, x);
        auto eigvalsRight = eigenvalues(ur, n, x);

        return std::max(eigvalsLeft.infinity_norm(), eigvalsRight.infinity_norm());
    }

    K lambda() const
    {
        return 1.0;
    }

    Dune::FieldVector<K, m> eigenvalues(const State& u, const Coordinate& n, const Coordinate& x) const
    {
        Dune::FieldVector<K, m> eigvals(0.0);

        K velocity = (u[1] / u[0]) * n[0] + (u[2] / u[0]) * n[1];
        K c = std::sqrt(u[0]);

        eigvals[0] = velocity - c;
        eigvals[1] = velocity;
        eigvals[2] = velocity + c;

        return eigvals;
    }

    K maxEigenvalue(const State& u, const Coordinate& n, const Coordinate& x) const
    {
        Dune::FieldVector<K, m> eigvals = eigenvalues(u, n, x);
        return eigvals.infinity_norm();
    }
};

} // end namespace Dune::Hypercut

#endif
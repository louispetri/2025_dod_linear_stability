#ifndef DUNE_HYPERCUT_SCALAR_TRANSPORT_MODEL_HH
#define DUNE_HYPERCUT_SCALAR_TRANSPORT_MODEL_HH

#include <array>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parametertree.hh>

namespace Dune {

    namespace Hypercut 
    {
        template<class FT>
        class ConstantVectorfield
        {
        public:
            using K = FT;

            ConstantVectorfield(const Dune::ParameterTree& paramtree) :
                v(paramtree.get<Dune::FieldVector<double, 2>>("direction"))
            {

            }

            ConstantVectorfield(const Dune::FieldVector<K, 2>& v)
                : v(v)
            {

            }

            Dune::FieldVector<K, 2> operator()(const Dune::FieldVector<K, 2>& x) const
            {
                return v;
            }

        private:
            Dune::FieldVector<K, 2> v;
        };

        template<class FT>
        class VaryingRampVectorfield
        {
        public:
            using K = FT;

            VaryingRampVectorfield(const Dune::ParameterTree& paramtree) :
                center(paramtree.get<Dune::FieldVector<double, 2>>("center")),
                normal(paramtree.get<Dune::FieldVector<double, 2>>("normal")),
                angle(paramtree.get<double>("angle"))
            {

            }

            VaryingRampVectorfield(const Dune::FieldVector<double, 2>& center,
                                   const Dune::FieldVector<double, 2>& normal,
                                   K angle) :
                center(center),
                normal(normal),
                angle(angle)
            {

            }

            Dune::FieldVector<K, 2> operator()(const Dune::FieldVector<K, 2>& x) const
            {
                auto levelSet = (x - center) * normal;
                Dune::FieldVector<K, 2> v(1.0);
                v[1] = std::tan(angle);
                v /= 2.0 * v.two_norm();
                v *= (2.0 - levelSet);

                return v;
            }

        private:
            Dune::FieldVector<K, 2> center;
            Dune::FieldVector<K, 2> normal;
            K angle;
        };

        template<class FT>
        class RotatingVectorfield
        {
        public:
            using K = FT;

            RotatingVectorfield(const Dune::ParameterTree& paramtree) :
                center(paramtree.get<Dune::FieldVector<double, 2>>("center"))
            {

            }

            RotatingVectorfield(const Dune::FieldVector<K, 2>& center) :
                center(center)
            {

            }

            Dune::FieldVector<K, 2> operator()(const Dune::FieldVector<K, 2>& x) const {
                Dune::FieldVector<K, 2> v;
                v[0] = 2 * M_PI * (-x[1] + center[1]);
                v[1] = 2 * M_PI * (x[0] - center[0]);

                return v;
            }

        private:
            Dune::FieldVector<K, 2> center;
        };

        template<class Vectorfield>
        class RampBoundary
        {
        public:
            using K = typename Vectorfield::K;

            RampBoundary(const Dune::ParameterTree& paramtree) :
                vectorfield(paramtree),
                center(paramtree.get<Dune::FieldVector<double, 2>>("center")),
                angle(paramtree.get<double>("angle"))
            {

            }

            RampBoundary(const Vectorfield& vectorfield, const Dune::FieldVector<K, 2>& center, K angle) :
                vectorfield(vectorfield),
                center(center),
                angle(angle)
            {

            }

            Dune::FieldVector<K, 1> operator()(const Dune::FieldVector<K, 1>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
            {
                auto vf = vectorfield(x);
                x = x - t * vf;
                x[0] -= center[0];

                Dune::FieldMatrix<K, 2, 2> m;
                m[0][0] = -std::cos(angle + M_PI);
                m[0][1] = std::sin(angle);
                m[1][0] = -std::sin(angle);
                m[1][1] = -std::cos(angle + M_PI);

                Dune::FieldVector<K, 2> y;
                m.mv(x, y);

                Dune::FieldVector<K, 1> res(std::sin(y[0] * M_PI * 2.0 / (std::sqrt(2) * (1 - center[0]))));

                return res;
            }

        private:
            Vectorfield vectorfield;
            Dune::FieldVector<K, 2> center;
            K angle;
        };

        template<class K>
        class CircleBoundary
        {
        public:
            CircleBoundary(const Dune::ParameterTree& paramtree) :
                center(paramtree.get<Dune::FieldVector<double, 2>>("center")),
                radius(paramtree.get<double>("radius"))
            {

            }

            Dune::FieldVector<K, 1> operator()(const Dune::FieldVector<K, 1>& u, const Dune::FieldVector<K, 2>& x, const Dune::FieldVector<K, 2>& n, K t) const
            {
                Dune::FieldVector<K, 2> rotation;
                rotation[0] = -std::cos(2 * M_PI * t);
                rotation[1] = -std::sin(2 * M_PI * t);

                auto c = center + (radius * 0.9) * rotation;
                auto d = (x - c).two_norm();
                auto w = radius / 4;
                auto v = std::exp(-d * d / w / w);

                return Dune::FieldVector<K, 1>(v);
            }
        private:
            Dune::FieldVector<K, 2> center;
            K radius;
        };

        template<class FT>
        class P0WaveBoundary
        {
        public:
            using K = FT;

            P0WaveBoundary(const Dune::ParameterTree& paramtree)
            {

            }

            P0WaveBoundary()
            {

            }

            Dune::FieldVector<K, 1> operator()(const Dune::FieldVector<K, 1>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
            {
                return 1.0;
            }
        };

        template<class Vectorfield>
        class P1WaveBoundary
        {
        public:
            using K = typename Vectorfield::K;

            P1WaveBoundary(const Dune::ParameterTree& paramtree) :
                vectorfield(paramtree)
            {

            }

            P1WaveBoundary(const Vectorfield& vectorfield)
                : vectorfield(vectorfield)
            {

            }

            Dune::FieldVector<K, 1> operator()(const Dune::FieldVector<K, 1>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
            {
                auto vf = vectorfield(x);
                x = x - t * vf;

                return x[0];
            }

        private:
            Vectorfield vectorfield;
        };

        template<class Vectorfield>
        class P2WaveBoundary
        {
        public:
            using K = typename Vectorfield::K;

            P2WaveBoundary(const Dune::ParameterTree& paramtree) :
                vectorfield(paramtree)
            {

            }

            Dune::FieldVector<K, 1> operator()(const Dune::FieldVector<K, 1>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
            {
                auto vf = vectorfield(x);
                x = x - t * vf;

                return x[0] * x[0] + x[1] * x[1];
            }

        private:
            Vectorfield vectorfield;
        };

        template<class VF>
        class ScalarTransportModel
        {
        public:
            static const int dim = 2;
            static const int m = 1;

            using Vectorfield = VF;
            using K = typename Vectorfield::K;

            using Flux = Dune::FieldMatrix<K, m, dim>;
            using RedistributionTerm = Dune::FieldMatrix<K, m, dim>;
            using Coordinate = Dune::FieldVector<K, dim>;
            using State = Dune::FieldVector<K, m>;
            
            static constexpr K Epsilon = 1e-10;

            ScalarTransportModel(const Dune::ParameterTree& paramtree) :
                vectorfield(paramtree)
            {

            }


            ScalarTransportModel(const Vectorfield& vectorfield) :
                vectorfield(vectorfield)
            {

            }

            Flux flux(const Coordinate& x, const Dune::FieldVector<K, 1> u) const
            {
                auto vf = vectorfield(x);
                Flux f;
                f[0][0] = u[0] * vf[0];
                f[0][1] = u[0] * vf[1];

                return f;
            }

            Dune::FieldVector<K, 1> numericalFlux(const Coordinate& n,
                                                       const Dune::FieldVector<K, 1>& ul,
                                                       const Dune::FieldVector<K, 1>& ur,
                                                       const Coordinate& x) const
            {
                auto vf = vectorfield(x);
                auto vfn = vf * n;

                if (std::abs(vfn) < Epsilon)
                {
                    vfn = 0.0;
                }

                return 0.5 * vfn * (ul + ur) + 0.5 * std::abs(vfn) * (ul - ur);
            }

            Dune::FieldMatrix<K, 1, 1> inflowMatrix(const Coordinate& n, const Coordinate& x) const
            {
                auto vf = vectorfield(x);
                auto vfn = vf * n;
                auto inflow = vfn < 0.0 ? vfn : 0.0;

                return inflow;
            }

            Dune::FieldMatrix<K, 1, 1> outflowMatrix(const Coordinate& n, const Coordinate& x) const
            {
                auto vf = vectorfield(x);
                auto vfn = vf * n;
                auto inflow = vfn > 0.0 ? vfn : 0.0;

                return inflow;
            }

            std::array<Dune::FieldMatrix<K, 1, 1>, 2> matrices(const Coordinate& x) const
            {
                auto v = vectorfield(x);

                Dune::FieldMatrix<K, 1, 1> A(v[0]);
                Dune::FieldMatrix<K, 1, 1> B(v[1]);

                return {A, B};
            }

            K seminorm(const Coordinate& n, const Dune::FieldVector<K, m>& u, const Coordinate& x) const
            {
                auto v = vectorfield(x);
                auto v_abs = std::abs(v * n);

                return v_abs * u;
            }

            Dune::FieldVector<K, 1> outflow(const Coordinate& n,
                                            const Dune::FieldVector<K, 1>& u,
                                            const Coordinate& x) const
            {
                auto vf = vectorfield(x);
                auto vfn = vf * n;

                if (vfn < Epsilon)
                {
                    vfn = 0.0;
                }

                return vfn * u;
            }

            Dune::FieldVector<K, 1> inflow(const Coordinate& n,
                                           const Dune::FieldVector<K, 1>& u,
                                           const Coordinate& x) const
            {
                auto vf = vectorfield(x);
                auto vfn = vf * n;

                if (vfn > -Epsilon)
                {
                    vfn = 0.0;
                }

                return vfn * u;
            }

            K lambda() const
            {
                // We added this function tempoarily to allow compiling
                // the wave equation related code but it should never be called
                assert(false);
                return 1.0;
            }

            Dune::FieldVector<K, m> eigenvalues(const State& u, const Coordinate& n, const Coordinate& x) const
            {
                return vectorfield(x) * n;
            }

        private:
            Vectorfield vectorfield;
        };
    }
}

#endif
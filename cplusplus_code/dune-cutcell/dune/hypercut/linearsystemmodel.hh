#ifndef DUNE_HYPERCUT_LINEAR_SYSTEM_MODEL_HH
#define DUNE_HYPERCUT_LINEAR_SYSTEM_MODEL_HH

#include <array>
#include <complex>
#include <dune/common/fmatrix.hh>
#include <dune/common/fmatrixev.hh>
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
  class TransportAlternatingSineCosineBoundaryCondition
  {
  public:
    static constexpr K Epsilon = 1e-13;

    TransportAlternatingSineCosineBoundaryCondition(const Dune::ParameterTree& paramtree)
    {
      eigenvaluesA = paramtree.get<Dune::FieldVector<double, m>>("evA");
      eigenvaluesB = paramtree.get<Dune::FieldVector<double, m>>("evB");
      auto tmpR = paramtree.get<Dune::FieldVector<Dune::FieldVector<double, m>, m>>("matrixR");
  
      for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
          R[i][j] = tmpR[i][j];
        }
      }
    }

    Dune::FieldVector<K, m> operator()(const Dune::FieldVector<K, m>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
    {
      bool output = x[0] + Epsilon > 1.0;

      std::array<Dune::FieldVector<K, 2>, m> y;

      for (int i = 0; i < m; ++i) {
        Dune::FieldVector<K, 2> velocity { eigenvaluesA[i], eigenvaluesB[i] };
        y[i] = x - t * velocity;
      }

      Dune::FieldVector<K, m> w;

      for (int i = 0; i < m; ++i) {
        Dune::FieldVector<K, 2> velocity { eigenvaluesA[i], eigenvaluesB[i] };

        if (i % 2 == 0) {
          w[i] = std::sin((y[i] * velocity) * M_PI * 2.0);
        } 
        else {
          w[i] = std::cos((y[i] * velocity) * M_PI * 2.0);
        }
      }

      Dune::FieldVector<K, m> res;
      R.mv(w, res);

      return res;
    }

  private:
      Dune::FieldMatrix<K, m, m> R;
      Dune::FieldVector<K, m> eigenvaluesA;
      Dune::FieldVector<K, m> eigenvaluesB;
  };

  template<int m, class K>
  class AlternatingSineCosineRampBoundaryCondition
  {
  public:
    static constexpr K Epsilon = 1e-13;

    AlternatingSineCosineRampBoundaryCondition(const Dune::ParameterTree& paramtree)
    {
      angle = paramtree.get<double>("angle");
      center = paramtree.get<Dune::FieldVector<double, 2>>("center");
      auto tmpA = paramtree.get<Dune::FieldVector<Dune::FieldVector<double, m>, m>>("matrixA");
      auto tmpB = paramtree.get<Dune::FieldVector<Dune::FieldVector<double, m>, m>>("matrixB");

      Dune::FieldMatrix<K, m, m> A;
      Dune::FieldMatrix<K, m, m> B;

      for (int i = 0; i < m; ++i) {
          A[i] = tmpA[i];
          B[i] = tmpB[i];
      }

      // matrices should be simultaneously diagonizeable
      assert((A * B - B * A).frobenius_norm() < Epsilon);

      Dune::FieldMatrix<K, m, m> eigenvectorsA;
      Dune::FMatrixHelp::eigenValuesVectors(A, eigenvaluesA, eigenvectorsA);
      Dune::FieldMatrix<K, m, m> transpose;

      for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            transpose[i][j] = eigenvectorsA[j][i];
        }
      }
      

      R = transpose;
      Rinv = R;
      Rinv.invert();

      Dune::FieldMatrix<K, m, m> eigenvalueMatrixB = Rinv * B * R;

      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          if (i == j) {
            eigenvaluesB[i] = eigenvalueMatrixB[i][i];
          }
          else {
            assert(std::abs(eigenvalueMatrixB[i][j]) < Epsilon);
          }
        }
      }
    }

    Dune::FieldVector<K, m> operator()(const Dune::FieldVector<K, m>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
    {
        Dune::FieldMatrix<K, 2, 2> rot;
        rot[0][0] = -std::cos(angle + M_PI);
        rot[0][1] = std::sin(angle);
        rot[1][0] = -std::sin(angle);
        rot[1][1] = -std::cos(angle + M_PI);

        std::array<Dune::FieldVector<K, 2>, m> y;

        for (int i = 0; i < m; ++i)
        {
            Dune::FieldVector<K, 2> velocity { eigenvaluesA[i], eigenvaluesB[i] };
            auto x_shifted = x - t * velocity;
            x_shifted[0] -= center[0];
            rot.mv(x_shifted, y[i]);
        }

        std::array<Dune::FieldVector<K, m>, m> bc;

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                if (j % 2 == 0)
                {
                    bc[i][j] = std::sin(y[i][0] * M_PI * 2.0 / (std::sqrt(2) * (1 - center[0])));
                }
                else
                {
                    bc[i][j] = std::cos(y[i][0] * M_PI * 2.0 / (std::sqrt(2) * (1 - center[0])));
                }
            }
        }

        Dune::FieldVector<K, m> w;

        for (int i = 0; i < m; ++i)
        {
            Dune::FieldVector<K, m> wTmp;
            Rinv.mv(bc[i], wTmp);
            w[i] = wTmp[i];
        }

        Dune::FieldVector<K, m> res;
        R.mv(w, res);

        return res;
    }

  private:
    Dune::FieldMatrix<K, m, m> R;
    Dune::FieldMatrix<K, m, m> Rinv;
    Dune::FieldVector<K, m> eigenvaluesA;
    Dune::FieldVector<K, m> eigenvaluesB;
    Dune::FieldVector<K, 2> center;
    K angle;
};

template<class K>
class AcousticsSineCosineBoundaryCondition
{
public:
    AcousticsSineCosineBoundaryCondition(Dune::ParameterTree& paramtree)
    {
        c = paramtree.get<double>("speedOfSound");
    }

    AcousticsSineCosineBoundaryCondition(K c)
        : c(c)
    {
    }

    Dune::FieldVector<K, 3> operator()(const Dune::FieldVector<K, 3>& u, Dune::FieldVector<K, 2> x, const Dune::FieldVector<K, 2>& n, K t) const
    {
        Dune::FieldVector<K, 3> res(0.0);
        res[0] = -(1 / c) * std::cos(2*M_PI*c*t) * (std::sin(2*M_PI*x[0]) + std::sin(2*M_PI*x[1]));
        res[1] = (1/c) * std::sin(2*M_PI*c*t) * std::cos(2*M_PI*x[0]);
        res[2] = (1/c) * std::sin(2*M_PI*c*t) * std::cos(2*M_PI*x[1]);

        return res;
    }

private:
    K c;
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

template<int m_, class K_, bool UseUpwindFlux = true>
class LinearSystemModel
{
public:
    static const int dim = 2;
    static const int m = m_;

    using K = K_;

    using Flux = Dune::FieldMatrix<K, m, dim>;
    using RedistributionTerm = Dune::FieldMatrix<K, m, dim>;
    using Coordinate = Dune::FieldVector<K, dim>;
    using State = Dune::FieldVector<K, m>;

    static constexpr K Epsilon = 1e-10;

    LinearSystemModel(const Dune::ParameterTree& paramtree)
    {
        auto tmpA = paramtree.get<Dune::FieldVector<Dune::FieldVector<double, m>, m>>("matrixA");
        auto tmpB = paramtree.get<Dune::FieldVector<Dune::FieldVector<double, m>, m>>("matrixB");

        LF_lambda = paramtree.get<double>("LF_lambda");

        for (int i = 0; i < m; ++i)
        {
            A[i] = tmpA[i];
            B[i] = tmpB[i];
        }
    }

    LinearSystemModel(const Dune::FieldMatrix<K, m, m>& A,
                        const Dune::FieldMatrix<K, m, m>& B,
                        K LF_lambda) :
        A(A), B(B), LF_lambda(LF_lambda)
    {
    }

    Flux flux(const Coordinate& x, const State& u) const
    {
        Flux f(0.0);
        State tmp(0.0);

        A.mv(u, tmp);

        for (std::size_t j = 0; j < m; ++j)
        {
            f[j][0] = tmp[j];
        }

        B.mv(u, tmp);

        for (std::size_t j = 0; j < m; ++j)
        {
            f[j][1] = tmp[j];
        }

        return f;
    }

    State numericalFlux(const Coordinate& n,
                        const Dune::FieldVector<K, m>& ul,
                        const Dune::FieldVector<K, m>& ur,
                        const Coordinate& x) const
    {
        if constexpr (UseUpwindFlux) {
            return upwindFlux(n, ul, ur, x);
        } else {
            return LaxFriedrichFlux(n, ul, ur, x);
        }
    }

    State upwindFlux(const Coordinate& n,
                        const Dune::FieldVector<K, m>& ul,
                        const Dune::FieldVector<K, m>& ur,
                        const Coordinate& x) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        Dune::FieldVector<K, m> eigenvalues;
        Dune::FieldMatrix<K, m, m> eigenvectors;

        Dune::FMatrixHelp::eigenValuesVectors(matrix, eigenvalues, eigenvectors);

        Dune::FieldMatrix<K, m, m> transpose;

        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < m; ++j)
            {
                transpose[i][j] = eigenvectors[j][i];
            }
        }

        eigenvectors = transpose;
        Dune::FieldMatrix<K, m, m> eigenvectorsInv = eigenvectors;
        eigenvectorsInv.invert();

        Dune::FieldMatrix<K, m, m> matrixPos(0.0);
        Dune::FieldMatrix<K, m, m> matrixNeg(0.0);

        for (std::size_t i = 0; i < m; ++i)
        {
            if (eigenvalues[i] > Epsilon)
            {
                matrixPos[i][i] = eigenvalues[i];
            }
            else if (eigenvalues[i] < -Epsilon)
            {
                matrixNeg[i][i] = eigenvalues[i];
            }
        }

        auto absMatrix = eigenvectors * (matrixPos - matrixNeg) * eigenvectorsInv;
        State flux(0.0);
        State flux2(0.0);
        matrix.mv((ul + ur), flux);
        absMatrix.mv(ul - ur, flux2);

        return 0.5 * flux + 0.5 * flux2;
    }

    State LaxFriedrichFlux(const Coordinate& n,
                            const Dune::FieldVector<K, m>& ul,
                            const Dune::FieldVector<K, m>& ur,
                            const Coordinate& x) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        State ful(0.0);
        State fur(0.0);

        matrix.mv(ul, ful);
        matrix.mv(ur, fur);

        return 0.5 * (ful + fur) + 0.5 * LF_lambda * (ul - ur);
    }

    State centralFlux(const Coordinate& n,
                      const Dune::FieldVector<K, m>& ul,
                      const Dune::FieldVector<K, m>& ur,
                      const Coordinate& x) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        State ful(0.0);
        State fur(0.0);

        matrix.mv(ul, ful);
        matrix.mv(ur, fur);

        return 0.5 * (ful + fur);
    }

    State centralFluxHalf(const Coordinate& n,
                          const Dune::FieldVector<K, m>& u,
                          const Coordinate& x) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        State fu(0.0);
        matrix.mv(u, fu);

        return 0.5 * fu;
    }

    State analyticFlux(const Coordinate& n,
                       const Dune::FieldVector<K, m>& u,
                       const Coordinate& x) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        State fu(0.0);

        matrix.mv(u, fu);

        return fu;
    }

    Dune::FieldMatrix<K, m, m> inflowMatrix(const Coordinate& n, const Coordinate& x) const
    {
        if constexpr (UseUpwindFlux) {
            return inflowMatrixUpwind(n);
        } else {
            return inflowMatrixLF(n);
        }
    }

    State inflow(const Coordinate& n, const State& u, const Coordinate& x) const
    {
        auto matrix = inflowMatrix(n, x);
        State v(0.0);
        matrix.mv(u, v);
        return v;
    }

    State outflow(const Coordinate& n, const State& u, const Coordinate& x) const
    {
        auto matrix = outflowMatrix(n, x);
        State v(0.0);
        matrix.mv(u, v);
        return v;
    }

    Dune::FieldMatrix<K, m, m> outflowMatrix(const Coordinate& n, const Coordinate& x) const
    {
        if constexpr (UseUpwindFlux) {
            return outflowMatrixUpwind(n);
        } else {
            return outflowMatrixLF(n);
        }
    }

    Dune::FieldMatrix<K, m, m> inflowMatrixUpwind(const Coordinate& n) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        Dune::FieldVector<K, m> eigenvalues;
        Dune::FieldMatrix<K, m, m> eigenvectors;

        Dune::FMatrixHelp::eigenValuesVectors(matrix, eigenvalues, eigenvectors);

        Dune::FieldMatrix<K, m, m> transpose;

        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < m; ++j)
            {
                transpose[i][j] = eigenvectors[j][i];
            }
        }

        eigenvectors = transpose;
        Dune::FieldMatrix<K, m, m> eigenvectorsInv = eigenvectors;
        eigenvectorsInv.invert();

        Dune::FieldMatrix<K, m, m> matrixNeg(0.0);

        for (std::size_t i = 0; i < m; ++i)
        {
            if (eigenvalues[i] < -Epsilon)
            {
                matrixNeg[i][i] = eigenvalues[i];
            }
        }

        return eigenvectors * matrixNeg * eigenvectorsInv;
    }

    Dune::FieldMatrix<K, m, m> inflowMatrixLF(const Coordinate& n) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        for (int i = 0; i < m; ++i) {
            matrix[i][i] -= LF_lambda;
        }

        matrix *= 0.5;
        return matrix;
    }

    Dune::FieldMatrix<K, m, m> outflowMatrixUpwind(const Coordinate& n) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        Dune::FieldVector<K, m> eigenvalues;
        Dune::FieldMatrix<K, m, m> eigenvectors;

        Dune::FMatrixHelp::eigenValuesVectors(matrix, eigenvalues, eigenvectors);

        Dune::FieldMatrix<K, m, m> transpose;

        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < m; ++j)
            {
                transpose[i][j] = eigenvectors[j][i];
            }
        }

        eigenvectors = transpose;
        Dune::FieldMatrix<K, m, m> eigenvectorsInv = eigenvectors;
        eigenvectorsInv.invert();

        Dune::FieldMatrix<K, m, m> matrixPos(0.0);

        for (std::size_t i = 0; i < m; ++i)
        {
            if (eigenvalues[i] > Epsilon)
            {
                matrixPos[i][i] = eigenvalues[i];
            }
        }

        return eigenvectors * matrixPos * eigenvectorsInv;
    }

    Dune::FieldMatrix<K, m, m> outflowMatrixLF(const Coordinate& n) const
    {
        Dune::FieldMatrix<K, m, m> matrix(0.0);

        matrix += n[0] * A;
        matrix += n[1] * B;

        for (int i = 0; i < m; ++i) {
            matrix[i][i] += LF_lambda;
        }

        matrix *= 0.5;
        return matrix;
    }

    std::array<Dune::FieldMatrix<K, m, m>, 2> matrices(const Coordinate& x) const
    {
        return {A, B};
    }

    K lambda(const State& ul, const State& ur, const Coordinate& n, const Coordinate& x) const
    {
        return lambda();
    }

    K lambda() const
    {
        return LF_lambda;
    }

    Dune::FieldVector<K, m> eigenvalues(const State& u, const Coordinate& n, const Coordinate& x) const
    {
        auto matrix = n[0] * A + n[1] * B;
        Dune::FieldVector<K, m> eigvals(0.0);
        Dune::FMatrixHelp::eigenValues(matrix, eigvals);
        return eigvals;
    }

    K maxEigenvalue(const State& u, const Coordinate& n, const Coordinate& x) const
    {
        Dune::FieldVector<K, m> eigvals = eigenvalues(u, n, x);
        return eigvals.infinity_norm();
    }

private:
    Dune::FieldMatrix<K, m, m> A;
    Dune::FieldMatrix<K, m, m> B;
    K LF_lambda;
};

} // end namespace Dune::Hypercut

#endif
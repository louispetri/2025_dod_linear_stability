#ifndef DUNE_HYPERCUT_HELPERFUNCTIONS_HH
#define DUNE_HYPERCUT_HELPERFUNCTIONS_HH

#include <array>
#include <chrono>
#include <cstdint>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune::Hypercut {

template<class Vector, class Matrix, class LocalView, int m>
Vector transformBasisComponents(const Vector& coefficients, const Matrix& transform, const LocalView& targetView, const LocalView& sourceView, std::integral_constant<int, m>)
{
  Vector result(0.0);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < transform.rows; ++j) {
      for (int k = 0; k < transform.cols; ++k) {
          result[targetView.tree().child(i).localIndex(j)] += transform[j][k] * coefficients[sourceView.tree().child(i).localIndex(k)];
      }
    }
  }

  return result;
}

template<class GlobalVector, class LocalView, class LocalVector>
void readLocalVector(const GlobalVector& globalVector, const LocalView& localView, LocalVector& localVector)
{
  for (typename LocalView::size_type i = 0; i < localView.size(); ++i) {
    localVector[i] = globalVector[localView.index(i)];
  }
}

template<class GlobalVector, class LocalView, class LocalVector>
void writeLocalVector(GlobalVector& globalVector, const LocalView& localView, const LocalVector& localVector)
{
  for (typename LocalView::size_type i = 0; i < localView.size(); ++i) {
    globalVector[localView.index(i)] = localVector[i];
  }
}

template<class GlobalVector, class LocalView, class LocalVector>
void addLocalVector(GlobalVector& globalVector, const LocalView& localView, const LocalVector& localVector)
{
  for (typename LocalView::size_type i = 0; i < localView.size(); ++i) {
    globalVector[localView.index(i)] += localVector[i];
  }
}

template<class Element, class Basis, class Matrix>
void computeLocalMassMatrix(const Element& element, const Basis& basis, Matrix& massMatrix)
{
  const auto& rule = Dune::QuadratureRules<typename Basis::Traits::DomainFieldType, Basis::Traits::dimDomain>::rule(element.type(), 2 * basis.order());

  for (auto&& qp : rule) {
    std::vector<typename Basis::Traits::RangeType> base(basis.size());
    basis.evaluateFunction(qp.position(), base); 

    for (std::size_t i = 0; i < base.size(); ++i) {
      for (std::size_t j = 0; j < base.size(); ++j) {
        massMatrix[i][j] += qp.weight() * base[i] * base[j];
      }
    }
  }
}

template<class QuadratureRule, class Basis, class Matrix>
void computeMassMatrix(const QuadratureRule& rule, const Basis& basis, Matrix& massMatrix)
{
  for (auto&& qp : rule) {
    std::vector<typename Basis::Traits::RangeType> base(basis.size());
    basis.evaluateFunction(qp.position(), base);

    for (std::size_t i = 0; i < base.size(); ++i) {
      for (std::size_t j = 0; j < base.size(); ++j) {
        massMatrix[i][j] += qp.weight() * base[i] * base[j];
      }
    }
  }
}

template<class Basis, class Matrix, std::size_t dim>
void computeDifferentiationMatrices(const Basis& basis, std::array<Matrix, dim>& differentationMatrices)
{
  if constexpr (dim != 2) {
    DUNE_THROW(Dune::Exception, "Can only compute differentiation matrices in two dimensions");
  }

  for (int i = 1; i <= basis.order(); ++i) {
    int k = i;
    int l = 0;

    for (int j = 0; j <= i; ++j) {
      int deg = (k - 1) + l;

      if (k > 0) {
        differentationMatrices[0][(deg * (deg + 1)) / 2 + l][((deg + 1) * (deg + 2)) / 2 + l] = k;
      }

      if (l > 0) {
        deg = k + (l - 1);
        differentationMatrices[1][(deg * (deg + 1)) / 2 + (l - 1)][((deg + 1) * (deg + 2)) / 2 + l] = l;
      }

      --k;
      ++l;
    }
  }
}

template<class InterpolationPoints, class Basis, class Matrix>
void computeVandermondeMatrix(const InterpolationPoints& points, const Basis& basis, std::vector<typename Basis::Traits::RangeType>& base, Matrix& vandermondeMatrix)
{
  base.resize(basis.size(), 0.0);

  for (std::size_t i = 0; i < points.size(); ++i) {
    basis.evaluateFunction(points[i], base);

    for (std::size_t j = 0; j < basis.size(); ++j) {
      vandermondeMatrix[i][j] = base[j];
    }
  }
}

template<class LocalView, class Matrix, class Coefficients, class K, int m>
void applyVandermondeMatrix(const LocalView& localView, const Matrix& matrix, const Coefficients& coefficients, std::vector<Dune::FieldVector<K, m>>& interpolationData)
{
  std::fill(interpolationData.begin(), interpolationData.end(), Dune::FieldVector<K, m>(0.0));

  for (std::size_t i = 0; i < interpolationData.size(); ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      for (std::size_t k = 0; k < localView.tree().child(j).size(); ++k) {
        interpolationData[i][j] += matrix[i][k] * coefficients[localView.tree().child(j).localIndex(k)];
      }
    }
  }
}

template<class LocalView, class Matrix, class Coefficients, class K, int m>
void applyTransposedVandermondeMatrix(const LocalView& localView, const Matrix& matrix, const std::vector<Dune::FieldVector<K, m>>& interpolationData, Coefficients& coefficients)
{
  coefficients = 0.0;

  for (std::size_t i = 0; i < interpolationData.size(); ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      for (std::size_t k = 0; k < localView.tree().child(j).size(); ++k) {
        coefficients[localView.tree().child(j).localIndex(k)] += matrix[k][i] * interpolationData[i][j];
      }
    }
  }
}

template<class LocalView, class FunctionValue, class Coefficients, class K, int m>
void fillInFunctionValue(const LocalView& localView, const Coefficients& x, const std::vector<K>& phi, FunctionValue& u, std::integral_constant<int, m>)
{
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < localView.tree().child(i).size(); ++j) {
      u[i] += x[localView.tree().child(i).localIndex(j)] * phi[j];
    }
  }
}

template<class Gradient, class JacobianInverseTransposed>
void transformGradients(std::vector<Gradient>& gradphi, const JacobianInverseTransposed& jacInvTransposed)
{
  Gradient tmp;

  for (int i = 0; i < gradphi.size(); ++i) {
    jacInvTransposed.mv(gradphi[i][0], tmp[0]);
    gradphi[i] = tmp;
  }
}

template<class LocalSubTriangulation, class Function, class Value>
Value cutcellVolumesAccumulate(LocalSubTriangulation& localSubTriangulation, const Function& f, Value current)
{
  for (const auto& element : elements(localSubTriangulation.gridView(), Dune::Partitions::interior)) {
    localSubTriangulation.bind(element);
    localSubTriangulation.createCutCells();

    for (auto cutCellIt = localSubTriangulation.cutCellsBegin();
      cutCellIt != localSubTriangulation.cutCellsEnd(); ++cutCellIt) {
      current = f(*cutCellIt, element, current);
    }
  }

  return current;
}

template<class Basis, class Vector, class Solution, class K, class Model, class LocalVector>
class ComputeInfError
{
public:
  using Value = Dune::FieldVector<K, Model::m>;

  ComputeInfError(const Basis& basis, const Vector& vector, const Solution& solution, K t, int intorderadd) :
    basis_(basis),
    vector_(vector),
    solution_(solution),
    t_(t),
    intorderadd_(intorderadd)
  {

  }

  template<class CutCell, class FundamentalElement>
  Value operator()(const CutCell& cutCell, const FundamentalElement& element, Value& current) const {
    auto localView = basis_.localView(cutCell.domainIndex());
    localView.bind(element);
    auto& finiteElement = localView.tree().child(0).finiteElement();
    LocalVector localValue;
    Dune::Hypercut::readLocalVector(vector_, localView, localValue);
    K error = 0.0;

    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = cutCell.geometry().global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);

      Value u(0.0);
      Dune::Hypercut::fillInFunctionValue(localView, localValue, phi, u, std::integral_constant<int, Model::m>());
      auto exact = solution_(qpg, t_);

      for (int i = 0; i < u.size(); ++i) {
        current[i] = std::max(current[i], std::abs(u[i] - exact[i]));
      }
    }

    return current;
  }

private:
  const Basis& basis_;
  const Vector& vector_;
  const Solution& solution_;
  mutable std::vector<Dune::QuadraturePoint<K, Model::dim>> quadratureRule_;
  const K t_;
  const int intorderadd_;
};

template<class Basis, class Vector, class Solution, class K, class Model, class LocalVector>
class ComputeSquaredL2Error
{
public:
  using Value = Dune::FieldVector<K, Model::m>;

  ComputeSquaredL2Error(const Basis& basis, const Vector& vector, const Solution& solution, K t, int intorderadd) :
    basis_(basis),
    vector_(vector),
    solution_(solution),
    t_(t),
    intorderadd_(intorderadd)
  {

  }

  template<class CutCell, class FundamentalElement>
  Value operator()(const CutCell& cutCell, const FundamentalElement& element, Value& current) const {
    auto localView = basis_.localView(cutCell.domainIndex());
    localView.bind(element);
    auto& finiteElement = localView.tree().child(0).finiteElement();
    LocalVector localValue;
    Dune::Hypercut::readLocalVector(vector_, localView, localValue);
    K error = 0.0;

    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = cutCell.geometry().global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);

      Value u(0.0);
      Dune::Hypercut::fillInFunctionValue(localView, localValue, phi, u, std::integral_constant<int, Model::m>());
      auto exact = solution_(qpg, t_);
      auto error = exact - u;

      for (int i = 0; i < u.size(); ++i) {
        current[i] += error[i] * error[i] * qp.weight();
      }
    }

    return current;
  }

private:
  const Basis& basis_;
  const Vector& vector_;
  const Solution& solution_;
  mutable std::vector<Dune::QuadraturePoint<K, Model::dim>> quadratureRule_;
  const K t_;
  const int intorderadd_;
};

template<class Basis, class Vector, class Solution, class K, class Model, class LocalVector>
class ComputeL1Error
{
public:
  using Value = Dune::FieldVector<K, Model::m>;

  ComputeL1Error(const Basis& basis, const Vector& vector, const Solution& solution, K t, int intorderadd) :
    basis_(basis),
    vector_(vector),
    solution_(solution),
    t_(t),
    intorderadd_(intorderadd)
  {

  }

  template<class CutCell, class FundamentalElement>
  Value operator()(const CutCell& cutCell, const FundamentalElement& element, Value& current) const {
    auto localView = basis_.localView(cutCell.domainIndex());
    localView.bind(element);
    auto& finiteElement = localView.tree().child(0).finiteElement();
    LocalVector localValue;
    Dune::Hypercut::readLocalVector(vector_, localView, localValue);

    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = cutCell.geometry().global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);

      Value u(0.0);
      Dune::Hypercut::fillInFunctionValue(localView, localValue, phi, u, std::integral_constant<int, Model::m>());
      auto exact = solution_(qpg, t_);
      auto error = exact - u;

      for (int i = 0; i < u.size(); ++i) {
        current[i] += std::abs(error[i]) * qp.weight();
      }
    }

    return current;
  }

private:
  const Basis& basis_;
  const Vector& vector_;
  const Solution& solution_;
  mutable std::vector<Dune::QuadraturePoint<K, Model::dim>> quadratureRule_;
  const K t_;
  const int intorderadd_;
};

// This class is based on a presentation by Bryce Adelstein-Lelbach (cppcon 2015)
class HighResolutionTimer
{
public:
  static std::uint64_t takeTimeStamp()
  {
    return std::uint64_t(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }

  HighResolutionTimer() :
    startTime_(takeTimeStamp())
  {

  }

  void restart()
  {
    startTime_ = takeTimeStamp();
  }

  double elapsed() const
  {
    return double(takeTimeStamp() - startTime_) * 1e-9;
  }

  std::uint64_t elapsedNanoseconds() const
  {
    return takeTimeStamp() - startTime_;
  }

private:
  std::uint64_t startTime_;
};

} // end namespace Dune::Hypercut

#endif
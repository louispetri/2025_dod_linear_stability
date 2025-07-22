#ifndef DUNE_HYPERCUT_INTERPOLATION_HH
#define DUNE_HYPERCUT_INTERPOLATION_HH

#include "dune/common/fvector.hh"
#include <dune/hypercut/helper.hh>

namespace Dune::Hypercut {

template<class Element1, class Element2, class Basis1, class Basis2, class Matrix>
void interpolateBasis(const Element1& element1, const Element2& element2,
                      const Basis1& basis1, const Basis2& basis2,
                      const Matrix& invMassMatrix, Matrix& baseChange)
{
  assert(basis1.size() == basis2.size());
  const auto& rule = Dune::QuadratureRules<typename Basis1::Traits::DomainFieldType, Basis1::Traits::dimDomain>::rule(element1.type(), 2 * basis1.order());

  for (auto&& qp : rule) {
    std::vector<typename Basis2::Traits::RangeType> y(basis2.size());
    basis2.evaluateFunction(element2.local(element1.global(qp.position())), y);

    std::vector<typename Basis1::Traits::RangeType> base(basis1.size());
    basis1.evaluateFunction(qp.position(), base);

    for (std::size_t i = 0; i < y.size(); ++i) {
      for (std::size_t j = 0; j < base.size(); ++j) {
        for (std::size_t k = 0; k < base.size(); ++k) {
          baseChange[j][i] += invMassMatrix[j][k] * qp.weight() * y[i] * base[k];
        }
      }
    }
  }
}

template<class Element1, class Element2, class Basis1, class Basis2, class Matrix>
void interpolateBasis(const std::vector<Dune::QuadraturePoint<typename Basis1::Traits::DomainFieldType, Basis1::Traits::dimDomain>>& rule,
                      const Element1& element1, const Element2& element2,
                      const Basis1& basis1, const Basis2& basis2,
                      const Matrix& invMassMatrix, Matrix& baseChange)
{
  assert(basis1.size() == basis2.size());

  for (auto&& qp : rule) {
    std::vector<typename Basis2::Traits::RangeType> y(basis2.size());
    basis2.evaluateFunction(element2.local(element1.global(qp.position())), y);

    std::vector<typename Basis1::Traits::RangeType> base(basis1.size());
    basis1.evaluateFunction(qp.position(), base);
    
    for (std::size_t i = 0; i < y.size(); ++i) {
      for (std::size_t j = 0; j < base.size(); ++j) {
        for (std::size_t k = 0; k < base.size(); ++k) {
          baseChange[j][i] += invMassMatrix[j][k] * qp.weight() * y[i] * base[k];
        }
      }
    }
  }
}

template<class F, class Element, class LocalView, class Matrix, class LocalVector, int m>
void interpolateFunctionLocal(const F& f, const LocalVector& state, const Element& element, const LocalView& localView,
                              const Matrix& invMassMatrix, LocalVector& localVector, int intorderadd, std::integral_constant<int, m>)
{
  const auto& basis = localView.tree().child(0).finiteElement().localBasis();
  const auto& rule = Dune::QuadratureRules<typename std::decay_t<decltype(basis)>::Traits::DomainFieldType, std::decay_t<decltype(basis)>::Traits::dimDomain>::rule(element.type(), 2 * basis.order() + intorderadd);

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());
    basis.evaluateFunction(qp.position(), phi);

    Dune::FieldVector<typename LocalVector::field_type, m> u(0.0);
    fillInFunctionValue(localView, state, phi, u, std::integral_constant<int, m>());

    auto y = f(u, element.global(qp.position()));

    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < phi.size(); ++j) {
        for (std::size_t k = 0; k < phi.size(); ++k) {
          localVector[localView.tree().child(i).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * y[i];
        }
      }
    }
  }
}

template<class F, class Element, class LocalView, class Matrix, class LocalVector, int m>
void interpolateFunctionLocal(const F& f, const LocalVector& state, const Element& element, const LocalView& localView,
                              const Matrix& invMassMatrix, LocalVector& localVector, int intorderadd, std::integral_constant<int, 1>)
{
  const auto& basis = localView.tree().child(0).finiteElement().localBasis();
  const auto& rule = Dune::QuadratureRules<typename std::decay_t<decltype(basis)>::Traits::DomainFieldType, std::decay_t<decltype(basis)>::Traits::dimDomain>::rule(element.type(), 2 * basis.order() + intorderadd);

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());
    basis.evaluateFunction(qp.position(), phi);

    Dune::FieldVector<typename LocalVector::field_type, m> u(0.0);
    fillInFunctionValue(localView, state, phi, u, std::integral_constant<int, 1>());

    auto y = f(u, element.global(qp.position()));

    for (std::size_t j = 0; j < phi.size(); ++j) {
      for (std::size_t k = 0; k < phi.size(); ++k) {
        localVector[localView.tree().child(0).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * y;
      }
    }
  }
}

template<class State, class NormalVector>
void reflectState(State& u, const State& uReflected, const NormalVector& normal)
{
    auto e = 0.5 * (u + uReflected);
    auto o = 0.5 * (u - uReflected);

    u[0] = e[0] - o[0];

    typename NormalVector::field_type dotProduct = 0.0;

    for (std::size_t i = 0; i < normal.size(); ++i) {
      dotProduct += (e[i + 1] - o[i + 1]) * normal[i];
    }

    // we assume that the vector components starts at the second position
    for (std::size_t i = 1; i < 1 + normal.size(); ++i) {
      u[i] = e[i] - o[i];
      u[i] -= 2.0 * dotProduct * normal[i - 1];
    }

    for (std::size_t i = 1 + normal.size(); i < u.size(); ++i) {
      u[i] = e[i] - o[i];
    }
}

template<class Element, class LocalView, class Matrix, class LocalVector, class Reflection, class ReflectionVector, int m>
void interpolateFunctionLocalReflected(const LocalVector& state, const Element& element, const LocalView& localView,
                                       const Matrix& invMassMatrix, const Reflection& reflection, const ReflectionVector& reflectionVector, LocalVector& localVector, int intorderadd, std::integral_constant<int, m>)
{
  assert(m > 1);

  const auto& basis = localView.tree().child(0).finiteElement().localBasis();
  const auto& rule = Dune::QuadratureRules<typename std::decay_t<decltype(basis)>::Traits::DomainFieldType, std::decay_t<decltype(basis)>::Traits::dimDomain>::rule(element.type(), 2 * basis.order() + intorderadd);

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());

    basis.evaluateFunction(element.local(reflection(element.global(qp.position()))), phi);
    Dune::FieldVector<typename LocalVector::field_type, m> uReflected(0.0);
    fillInFunctionValue(localView, state, phi, uReflected, std::integral_constant<int, m>());

    basis.evaluateFunction(qp.position(), phi);
    Dune::FieldVector<typename LocalVector::field_type, m> u(0.0);
    fillInFunctionValue(localView, state, phi, u, std::integral_constant<int, m>());

    auto normal = reflectionVector(qp.position());
    reflectState(u, uReflected, normal);

    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < phi.size(); ++j) {
        for (std::size_t k = 0; k < phi.size(); ++k) {
          localVector[localView.tree().child(i).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * u[i];
        }
      }
    }
  }
}

template<class QuadratureRule, class Element, class LocalView, class Matrix, class LocalVector, class Reflection, class ReflectionVector, int m>
void interpolateFunctionLocalReflected(const QuadratureRule& rule, const LocalVector& state, const Element& element, const LocalView& localView,
                                       const Matrix& invMassMatrix, const Reflection& reflection, const ReflectionVector& reflectionVector, LocalVector& localVector, int intorderadd, std::integral_constant<int, m>)
{
  assert(m > 1);

  const auto& basis = localView.tree().child(0).finiteElement().localBasis();

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());

    basis.evaluateFunction(element.local(reflection(element.global(qp.position()))), phi);
    Dune::FieldVector<typename LocalVector::field_type, m> uReflected(0.0);
    fillInFunctionValue(localView, state, phi, uReflected, std::integral_constant<int, m>());

    basis.evaluateFunction(qp.position(), phi);
    Dune::FieldVector<typename LocalVector::field_type, m> u(0.0);
    fillInFunctionValue(localView, state, phi, u, std::integral_constant<int, m>());

    auto normal = reflectionVector(qp.position());
    reflectState(u, uReflected, normal);

    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < phi.size(); ++j) {
        for (std::size_t k = 0; k < phi.size(); ++k) {
          localVector[localView.tree().child(i).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * u[i];
        }
      }
    }
  }
}

template<class F, class Element, class LocalView, class Matrix, class LocalVector, int m>
void interpolateFunctionLocal(const F& f, const Element& element, const LocalView& localView,
                              const Matrix& invMassMatrix, LocalVector& localVector, int intorderadd, std::integral_constant<int, m>)
{
  const auto& basis = localView.tree().child(0).finiteElement().localBasis();
  const auto& rule = Dune::QuadratureRules<typename std::decay_t<decltype(basis)>::Traits::DomainFieldType, std::decay_t<decltype(basis)>::Traits::dimDomain>::rule(element.type(), 2 * basis.order() + intorderadd);

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());
    basis.evaluateFunction(qp.position(), phi);

    auto y = f(element.global(qp.position()));

    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < phi.size(); ++j) {
        for (std::size_t k = 0; k < phi.size(); ++k) {
          localVector[localView.tree().child(i).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * y[i];
        }
      }
    }
  }
}

template<class F, class Element, class LocalView, class Matrix, class LocalVector>
void interpolateFunctionLocal(const F& f, const Element& element, const LocalView& localView,
                              const Matrix& invMassMatrix, LocalVector& localVector, int intorderadd, std::integral_constant<int, 1>)
{
  const auto& basis = localView.tree().child(0).finiteElement().localBasis();
  const auto& rule = Dune::QuadratureRules<typename std::decay_t<decltype(basis)>::Traits::DomainFieldType, std::decay_t<decltype(basis)>::Traits::dimDomain>::rule(element.type(), 2 * basis.order() + intorderadd);

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());
    basis.evaluateFunction(qp.position(), phi);

    auto y = f(element.global(qp.position()));

    for (std::size_t j = 0; j < phi.size(); ++j) {
      for (std::size_t k = 0; k < phi.size(); ++k) {
        localVector[localView.tree().child(0).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * y;
      }
    }
  }
}

template<class F, class QuadratureRule, class BoundingBox, class LocalView, class Matrix, class LocalVector, int m>
void interpolateFunctionOnCutCell(const F& f, const QuadratureRule& rule, const BoundingBox& boundingBox, const LocalView& localView,
                              const Matrix& invMassMatrix, LocalVector& localVector, int intorderadd, std::integral_constant<int, m>)
{
  const auto& basis = localView.tree().child(0).finiteElement().localBasis();

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());
    basis.evaluateFunction(qp.position(), phi);
    auto y = f(boundingBox.global(qp.position()));

    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < phi.size(); ++j) {
        for (std::size_t k = 0; k < phi.size(); ++k) {
          localVector[localView.tree().child(i).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * y[i];
        }
      }
    }
  }
}

template<class F, class QuadratureRule, class BoundingBox, class LocalView, class Matrix, class LocalVector>
void interpolateFunctionOnCutCell(const F& f, const QuadratureRule& rule, const BoundingBox& boundingBox, const LocalView& localView,
                                  const Matrix& invMassMatrix, LocalVector& localVector, int intorderadd, std::integral_constant<int, 1>)
{
  const auto& basis = localView.tree().child(0).finiteElement().localBasis();

  for (auto&& qp : rule) {
    std::vector<typename std::decay_t<decltype(basis)>::Traits::RangeType> phi(basis.size());
    basis.evaluateFunction(qp.position(), phi);
    auto y = f(boundingBox.global(qp.position()));

    for (std::size_t j = 0; j < phi.size(); ++j) {
      for (std::size_t k = 0; k < phi.size(); ++k) {
        localVector[localView.tree().child(0).localIndex(j)] += qp.weight() * invMassMatrix[j][k] * phi[k] * y;
      }
    }
  }
}

template<class ST, class Basis, class F, class Vector, int m, int ComponentSpaceDimension>
void interpolate(ST& subTriangulation, const Basis& basis, const F& f, Vector& vector,
                 std::integral_constant<int, m>, std::integral_constant<int, ComponentSpaceDimension>,
                 int intorderadd = 3)
{
  std::vector<Dune::QuadraturePoint<typename ST::ctype, ST::dim>> quadratureRule;

  for (const auto& element :  elements(subTriangulation.localSubTriangulation().gridView(), Dune::Partitions::interior)) {
    subTriangulation.localSubTriangulation().bindOnVolume(element);
    subTriangulation.localSubTriangulation().createCutCells();

    for (auto cutCellIt = subTriangulation.localSubTriangulation().cutCellsBegin();
         cutCellIt != subTriangulation.localSubTriangulation().cutCellsEnd(); ++cutCellIt) {
      int currentDomainIndex = cutCellIt->domainIndex();
      auto localView = basis.localView(currentDomainIndex);
      localView.bind(element);

      const auto& boundingBox = subTriangulation.cutCellInformation().information(element, currentDomainIndex).boundingBox;

      Dune::FieldMatrix<typename ST::ctype, ComponentSpaceDimension, ComponentSpaceDimension> invMassMatrix;
      Dune::FieldVector<typename ST::ctype, m * ComponentSpaceDimension> localVector;

      cutCellIt->quadratureRule(
        quadratureRule, 2 * localView.tree().child(0).finiteElement().localBasis().order() + intorderadd);
        computeMassMatrix(quadratureRule, localView.tree().child(0).finiteElement().localBasis(),
                          invMassMatrix);
      invMassMatrix.invert();

      interpolateFunctionOnCutCell(f, quadratureRule, boundingBox, localView, invMassMatrix,
                                   localVector, intorderadd, std::integral_constant<int, m>());

      writeLocalVector(vector, localView, localVector);
    }
  }
}

} // end namespace Dune::Hypercut

#endif
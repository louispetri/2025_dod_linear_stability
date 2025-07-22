#ifndef DUNE_HYPERCUT_LOCALOPERATOR_HH
#define DUNE_HYPERCUT_LOCALOPERATOR_HH

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <vector>

#include <dune/common/fmatrixev.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/hypercut/helper.hh>
#include <dune/hypercut/linalg.hh>

namespace Dune::Hypercut {

template<class Model_, class BoundaryCondition, int ComponentSpaceDimension_, class K_>
class LocalOperator
{
public:
  using K = K_;
  using Model = Model_;
  static const int ComponentSpaceDimension = ComponentSpaceDimension_;
  static const int localDimension = Model::m * ComponentSpaceDimension;

  using LocalVector = FieldVector<K, localDimension>;
  using LocalMatrix = FieldMatrix<K, localDimension, localDimension>;

  using QuadraturePoint = Dune::QuadraturePoint<K, Model::dim>;

  static constexpr K Epsilon = 1e-12;

  LocalOperator(const Model& model, const BoundaryCondition& boundaryCondition, K cfl, K capacityFactor, K tau = 1.0, int intorderadd = 1) :
    model_(model),
    boundaryCondition_(boundaryCondition),
    uWeightFactor_(0.5),
    vWeightFactor_(0.5),
    dt_(-1.0),
    cfl_(cfl),
    tau_(tau),
    capacityFactor_(capacityFactor),
    intorderadd_(intorderadd),
    useMatrixStabilization_(true)
  {

  }

  LocalOperator(const Model& model, const BoundaryCondition& boundaryCondition, K cfl, K capacityFactor, bool useMatrixStabilization) :
    model_(model),
    boundaryCondition_(boundaryCondition),
    uWeightFactor_(0.5),
    vWeightFactor_(0.5),
    dt_(-1.0),
    cfl_(cfl),
    tau_(1.0),
    capacityFactor_(capacityFactor),
    intorderadd_(1),
    useMatrixStabilization_(useMatrixStabilization)
  {

  }

  template<class CutCell, class CutIntersectionIt, class LocalView>
  void computeWeights(const CutCell& cutCell, CutIntersectionIt cutIntersectionBegin, CutIntersectionIt cutIntersectionEnd, const LocalView& localView, int domainIndex, std::vector<LocalMatrix>& corrections) const
  {
    weights_.clear();

    LocalMatrix volumeOperator(0.0);
    LocalMatrix outflowOperator(0.0);

    Dune::FieldMatrix<K, Model::m, Model::m> completeInflowMatrix(0.0);

    computeVolumeOperator(cutCell, localView, volumeOperator);

    for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
      if (cutIntersectionIt->insideDomainIndex() == domainIndex || cutIntersectionIt->outsideDomainIndex() == domainIndex) {
        computeOutflowOperator(*cutIntersectionIt, outflowOperator, localView, cutIntersectionIt->domainIsOutside(domainIndex));
        computeInflowP0Operator(*cutIntersectionIt, completeInflowMatrix, cutIntersectionIt->domainIsOutside(domainIndex));
      }
    }

    auto cellFlowOperator = 0.5 * outflowOperator + volumeOperator;
    pseudoInverse(cellFlowOperator, cellFlowOperator);
    completeInflowMatrix.invert();

    Dune::FieldMatrix<K, Model::m, Model::m> identity(0.0);

    for (int i = 0; i < Model::m; ++i) {
      identity[i][i] = 1.0;
    }

    for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
      if (cutIntersectionIt->insideDomainIndex() == domainIndex || cutIntersectionIt->outsideDomainIndex() == domainIndex) {
        LocalMatrix inflowOperator(0.0);
        computeInflowOperator(*cutIntersectionIt, inflowOperator, localView, cutIntersectionIt->domainIsOutside(domainIndex));

        Dune::FieldMatrix<K, Model::m, Model::m> inflowMatrix(0.0);
        computeInflowP0Operator(*cutIntersectionIt, inflowMatrix, cutIntersectionIt->domainIsOutside(domainIndex));
        inflowMatrix.leftmultiply(completeInflowMatrix);

        LocalMatrix residualOperator(0.0);
        computeResidualOperator(cutCell, localView, inflowMatrix, identity, residualOperator);

        auto volumeCorrection = -0.5 * residualOperator + 0.5 * residualOperator.transposed();

        weights_.emplace_back(cellFlowOperator);
        weights_.back().leftmultiply(-0.5 * inflowOperator + volumeCorrection);
        weights_.back() = weights_.back().transposed();
      }
    }

    std::vector<LocalMatrix> faceOperators;

    for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
      if (cutIntersectionIt->insideDomainIndex() == domainIndex || cutIntersectionIt->outsideDomainIndex() == domainIndex) {
        for (auto weightIt = weights_.begin(); weightIt != weights_.end(); ++weightIt) {
          faceOperators.emplace_back(LocalMatrix(0.0));

          computeOutflowOperator(*cutIntersectionIt, faceOperators.back(), localView, cutIntersectionIt->domainIsOutside(domainIndex));
          faceOperators.back().rightmultiply(*weightIt);
          faceOperators.back() *= 0.5;
          faceOperators.back() = faceOperators.back() + faceOperators.back().transposed();

          corrections.emplace_back(faceOperators.back());
          removeNegativeEigenvalues(corrections.back());
        }
      }
    }
  }

  template<class CutCell, class CutIntersectionIt, class LocalView>
  void computeDoDWaveEquationOperators(const CutCell& cutCell, CutIntersectionIt cutIntersectionBegin, CutIntersectionIt cutIntersectionEnd, const LocalView& localView, const std::vector<LocalVector>& neighborValues, int domainIndex, bool& smallCellIsBoundaryCell) const
  {
    std::vector<typename Model::Coordinate> normals;
    std::vector<typename Model::Coordinate> tangents;
    std::vector<typename Model::Coordinate> reflectedNormals;
    std::vector<typename Model::Coordinate> reflectedTangents;
    std::vector<K> waveSpeeds;

    int intorder = 2 * localView.tree().child(0).finiteElement().localBasis().order() + intorderadd_;
    using RangeType = typename std::decay_t<decltype(localView.tree().child(0).finiteElement().localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi_s(localView.tree().child(0).finiteElement().localBasis().size());

    for (auto inflowIntersectionIt = cutIntersectionBegin; inflowIntersectionIt != cutIntersectionEnd; ++inflowIntersectionIt) {
      inflowIntersectionIt->quadratureRule(quadratureRule_, intorder);
      K maxSpeed = 0.0;

      for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
        std::size_t neighborIndex = std::distance(cutIntersectionBegin, cutIntersectionIt);

        if (cutIntersectionIt->neighbor()) {
          auto cutCellGeometry = cutIntersectionIt->domainIsOutside(domainIndex) ? cutIntersectionIt->geometryInOutside() : cutIntersectionIt->geometryInInside();

          cutIntersectionIt->quadratureRule(quadratureRule_, intorder);

          for (const auto& qp : quadratureRule_) {
            localView.tree().child(0).finiteElement().localBasis().evaluateFunction(cutCellGeometry.local(qp.position()), phi_s);
            Dune::FieldVector<K, Model::m> u(0.0);
            fillInFunctionValue(localView, neighborValues[neighborIndex], phi_s, u, std::integral_constant<int, Model::m>());
            maxSpeed = std::max(model_.maxEigenvalue(u, cutIntersectionIt->unitOuterNormal(qp.position()), qp.position()), maxSpeed);
          }
        }
      }

      waveSpeeds.emplace_back(maxSpeed);
    }

    for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
      normals.emplace_back(cutIntersectionIt->unitOuterNormal(typename Model::Coordinate(0.0)));

      if (cutIntersectionIt->domainIsOutside(domainIndex)) {
        normals.back() *= -1.0;
      }

      tangents.emplace_back(typename Model::Coordinate(0.0));
      tangents.back()[0] = normals.back()[1];
      tangents.back()[1] = -normals.back()[0];

      if (!cutIntersectionIt->neighbor()) {
        assert(smallCellIsBoundaryCell == false);
        smallCellIsBoundaryCell = true;
        reflectionVector_ = normals.back();
        reflectionCenter_ = cutIntersectionIt->center();
      }
    }

    assert(normals.size() == 3);

    Dune::FieldVector<K, 3> LFStabilitySolution(0.0);
    computeLFStabilityWeights(normals, tangents, waveSpeeds, LFStabilitySolution);

    inflowLFVectorfields_.resize(3);
    inflowReflectedLFVectorfields_.resize(3);

    for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
      std::size_t inflowFaceIndex = std::distance(cutIntersectionBegin, cutIntersectionIt);

      auto vf = -normals[inflowFaceIndex] * waveSpeeds[inflowFaceIndex] + LFStabilitySolution[inflowFaceIndex] * tangents[inflowFaceIndex];
      inflowLFVectorfields_[inflowFaceIndex] = vf;
    }

    if (smallCellIsBoundaryCell) {
      reflectedNormals.resize(normals.size());
      std::transform(normals.begin(), normals.end(), reflectedNormals.begin(), [this](const auto& n){ return n - 2.0 * (n * reflectionVector_) * reflectionVector_; });

      reflectedTangents.resize(reflectedNormals.size());

      for (std::size_t i = 0; i < reflectedNormals.size(); ++i) {
        reflectedTangents[i][0] = reflectedNormals[i][1];
        reflectedTangents[i][1] = -reflectedNormals[i][0];
      }

      Dune::FieldVector<K, 3> reflectedLFStabilitySolution(0.0);
      computeLFStabilityWeights(reflectedNormals, reflectedTangents, waveSpeeds, reflectedLFStabilitySolution);

      for (auto cutIntersectionIt = cutIntersectionBegin; cutIntersectionIt != cutIntersectionEnd; ++cutIntersectionIt) {
        std::size_t inflowFaceIndex = std::distance(cutIntersectionBegin, cutIntersectionIt);
        inflowReflectedLFVectorfields_[inflowFaceIndex] = -reflectedNormals[inflowFaceIndex] * waveSpeeds[inflowFaceIndex] + reflectedLFStabilitySolution[inflowFaceIndex] * reflectedTangents[inflowFaceIndex];
      }
    }
  }

  template<class CutCell, class LocalView>
  void dgVolume(const CutCell& cutCell, const LocalView& localView, const LocalVector& x, LocalVector& r) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = cutCell.geometry().global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, cutCell.geometry().jacobianInverseTransposed(qp.position()));

      Dune::FieldVector<K, Model::m> u(0.0);
      fillInFunctionValue(localView, x, phi, u, std::integral_constant<int, Model::m>());
      const auto flux = model_.flux(qpg, u);

      const K factor = qp.weight();

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < gradphi.size(); ++j) {
          r[localView.tree().child(i).localIndex(j)] -= flux[i] * gradphi[j][0] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dgFace(const CutIntersection& cutIntersection,
              const LocalView& insideView, const LocalView& outsideView,
              const LocalVector& x_s, const LocalVector& x_n,
              LocalVector& r_s, LocalVector& r_n) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();
    const auto& outsideFiniteElement = outsideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order();

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());
    std::vector<RangeType> phi_n(outsideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = cutIntersection.geometryInInside();
    const auto& geometryInOutside = cutIntersection.geometryInOutside();

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qppos_s = geometryInInside.local(qp.position());
      const auto qppos_n = geometryInOutside.local(qp.position());
      const auto qpg = qp.position();

      insideFiniteElement.localBasis().evaluateFunction(qppos_s, phi_s);
      outsideFiniteElement.localBasis().evaluateFunction(qppos_n, phi_n);

      Dune::FieldVector<K, Model::m> u_s(0.0);
      Dune::FieldVector<K, Model::m> u_n(0.0);

      fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
      fillInFunctionValue(outsideView, x_n, phi_n, u_n, std::integral_constant<int, Model::m>());

      auto normal = cutIntersection.unitOuterNormal(qp.position());
      const K factor = qp.weight();
      const auto flux = model_.numericalFlux(normal, u_s, u_n, qpg);

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s[j] * flux[i] * factor;
          r_n[outsideView.tree().child(i).localIndex(j)] -= phi_n[j] * flux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dgBoundary(const CutIntersection& cutIntersection, const LocalView& localView,
                  const LocalVector& x, LocalVector& r) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    int intorder = 2 * finiteElement.localBasis().order() + intorderadd_;
    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi(finiteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = cutIntersection.geometryInInside();

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qppos_s = geometryInInside.local(qp.position());
      const auto qpg = qp.position();

      finiteElement.localBasis().evaluateFunction(qppos_s, phi);

      Dune::FieldVector<K, Model::m> u(0.0);
      fillInFunctionValue(localView, x, phi, u, std::integral_constant<int, Model::m>());

      auto normal = cutIntersection.unitOuterNormal(qp.position());
      auto bc = boundaryCondition_(u, qpg, normal, t_);

      const K factor = qp.weight();
      const auto flux = model_.numericalFlux(normal, u, bc, qpg);

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < localView.tree().child(i).size(); ++j) {
          r[localView.tree().child(i).localIndex(j)] += phi[j] * flux[i] * factor;
        }
      }
    }
  }

  template<class CutCell, class LocalView>
  void dodVolume(const CutCell& cutCell, const LocalView& localView, const LocalView& extendedView,
                 const LocalVector& x, const LocalVector& x_extended, std::size_t inflowFaceIndex,
                 const std::vector<LocalVector>& extendedBasisCoefficients, K capacity,
                 LocalVector& r, LocalVector& r_extended) const
  {
    LocalVector x_transformed(0.0);
    LocalVector x_extended_transformed(0.0);

    auto weight = weights_[inflowFaceIndex];

    weight.mv(x, x_transformed);
    weight.mv(x_extended, x_extended_transformed);

    const auto& finiteElement = localView.tree().child(0).finiteElement();

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    LocalVector tmp(0.0);

    const auto& geometry = cutCell.geometry();

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = geometry.global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, geometry.jacobianInverseTransposed(qp.position()));

      Dune::FieldVector<K, Model::m> uiTransformed(0.0);
      Dune::FieldVector<K, Model::m> uTransformed(0.0);

      fillInFunctionValue(localView, x_extended_transformed, phi, uiTransformed, std::integral_constant<int, Model::m>());
      fillInFunctionValue(localView, x_transformed, phi, uTransformed, std::integral_constant<int, Model::m>());

      const K factor = qp.weight();

      const auto flux = model_.flux(qpg, uiTransformed) - model_.flux(qpg, uTransformed);

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < gradphi.size(); ++j) {
          r[localView.tree().child(i).localIndex(j)] -= (1 - capacity) * flux[i] * gradphi[j][0] * factor;
        }
      }

      if (extendedBasisCoefficients.size() > 0) {
        auto extendedBasisCoefficientsIt = extendedBasisCoefficients.begin();

        for (std::size_t i = 0; i < Model::m; ++i) {
            for (std::size_t j = 0; j < extendedView.tree().child(i).size(); j++) {
                Dune::FieldMatrix<K, 1, Model::dim> extendedGradient(0.0);

                for (int k = 0; k < Model::m; ++k) {
                    for (int l = 0; l < localView.tree().child(k).size(); ++l) {
                        extendedGradient += (*extendedBasisCoefficientsIt)[localView.tree().child(k).localIndex(l)] * gradphi[l];
                    }
                }

                r_extended[extendedView.tree().child(i).localIndex(j)] += (1 - capacity) * flux[i] * extendedGradient[0] * factor;
                ++extendedBasisCoefficientsIt;
            }
        }
      }
    }
  }

  template<class CutIntersection, class LocalView, class BaseChange>
  void dodFace(const CutIntersection& cutIntersection,
               const LocalView& insideView, const LocalView& outsideView,
               const LocalVector& x_s, const LocalVector& x_n,
               const LocalVector& x_extended,
               std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
               const BaseChange& baseChange, K capacity,
               LocalVector& r_s, LocalVector& r_n,
               bool smallCutCellIsOutside) const
  {
    LocalVector x_n_transformed(0.0);
    LocalVector x_s_transformed(0.0);
    LocalVector x_extended_transformed(0.0);

    LocalVector tmp = x_n;
    tmp = transformBasisComponents(tmp, baseChange, insideView, outsideView, std::integral_constant<int, Model::m>());

    auto weight = weights_[inflowFaceIndex];

    weight.mv(tmp, x_n_transformed);
    weight.mv(x_s, x_s_transformed);
    weight.mv(x_extended, x_extended_transformed);

    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();
    const auto& outsideFiniteElement = outsideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());
    std::vector<RangeType> phi_n(outsideFiniteElement.localBasis().size());
    std::vector<Dune::FieldVector<K, Model::m>> phi_s_transformed(insideView.tree().size());
    std::vector<Dune::FieldVector<K, Model::m>> phi_n_transformed(outsideView.tree().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();
    const auto& geometryInOutside = smallCutCellIsOutside ? cutIntersection.geometryInInside() : cutIntersection.geometryInOutside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qppos_s = geometryInInside.local(qp.position());
      const auto qppos_n = geometryInOutside.local(qp.position());
      const auto qpg = qp.position();

      insideFiniteElement.localBasis().evaluateFunction(qppos_s, phi_s);
      outsideFiniteElement.localBasis().evaluateFunction(qppos_n, phi_n);

      Dune::FieldVector<K, Model::m> u_s(0.0);
      Dune::FieldVector<K, Model::m> u_s_transformed(0.0);
      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u_extended_transformed(0.0);
      Dune::FieldVector<K, Model::m> u_n_transformed(0.0);
      Dune::FieldVector<K, Model::m> u_n(0.0);

      fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_s_transformed, phi_s, u_s_transformed, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended, phi_s, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended_transformed, phi_s, u_extended_transformed, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_n_transformed, phi_s, u_n_transformed, std::integral_constant<int, Model::m>());
      fillInFunctionValue(outsideView, x_n, phi_n, u_n, std::integral_constant<int, Model::m>());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      const auto factor = qp.weight();

      auto flux = model_.numericalFlux(normal, u_extended_transformed, u_n_transformed, qpg) - model_.numericalFlux(normal, u_s_transformed, u_n_transformed, qpg);
      flux *= (1.0 - capacity) * uWeightFactor_;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s[j] * flux[i] * factor;
          r_n[outsideView.tree().child(i).localIndex(j)] += -phi_n[j] * flux[i] * factor;
        }
      }

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); j++) {
          phi_s_transformed[insideView.tree().child(i).localIndex(j)] = Dune::FieldVector<K, Model::m>(0.0);
          phi_n_transformed[outsideView.tree().child(i).localIndex(j)] = Dune::FieldVector<K, Model::m>(0.0);

          for (std::size_t k = 0; k < Model::m; ++k) {
            for (std::size_t l = 0; l < insideView.tree().child(k).size(); ++l) {
              phi_s_transformed[insideView.tree().child(i).localIndex(j)][k] += weight[insideView.tree().child(k).localIndex(l)][insideView.tree().child(i).localIndex(j)] * phi_s[l];

              for (std::size_t p = 0; p < insideView.tree().child(k).size(); ++p) {
                phi_n_transformed[outsideView.tree().child(i).localIndex(j)][k] += weight[insideView.tree().child(k).localIndex(l)][insideView.tree().child(i).localIndex(p)] * baseChange[p][j] * phi_s[l];
              }
            }
          }
        }
      }

      flux = model_.numericalFlux(normal, u_extended, u_n, qpg) - model_.numericalFlux(normal, u_s, u_n, qpg);
      flux *= (1.0 - capacity) * vWeightFactor_;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); j++) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s_transformed[insideView.tree().child(i).localIndex(j)] * flux * factor;
          r_n[outsideView.tree().child(i).localIndex(j)] += -phi_n_transformed[outsideView.tree().child(i).localIndex(j)] * flux * factor;
        }
      }
    }
  }

  template<class CutCell, class LocalView>
  void dodVolumeInternal(const CutCell& cutCell, const LocalView& localView,
                         const LocalVector& x, const LocalVector& x_extended,
                         std::size_t inflowFaceIndex,
                         K capacity, K weightFactor,
                         LocalVector& r) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    const auto& geometry = cutCell.geometry();

    auto vf = inflowLFVectorfields_[inflowFaceIndex];

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = geometry.global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, geometry.jacobianInverseTransposed(qp.position()));

      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u(0.0);

      fillInFunctionValue(localView, x_extended, phi, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(localView, x, phi, u, std::integral_constant<int, Model::m>());

      const K factor = qp.weight();

      auto flux = model_.flux(qpg, u_extended) - model_.flux(qpg, u);
      flux *= (1.0 - capacity) * weightFactor;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < gradphi.size(); ++j) {
          r[localView.tree().child(i).localIndex(j)] -= flux[i] * gradphi[j][0] * factor;
          r[localView.tree().child(i).localIndex(j)] -= (1.0 - capacity) * (u_extended[i] - u[i]) * 0.5 * (vf * gradphi[j][0]) * factor;
        }
      }
    }
  }

  template<class CutCell, class LocalView>
  void dodResidualInternal(const CutCell& cutCell, const LocalView& localView,
                           const LocalVector& x, const LocalVector& x_extended,
                           std::size_t inflowFaceIndex,
                           K capacity, K weightFactor,
                           LocalVector& r) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    const auto& geometry = cutCell.geometry();

    auto vf = inflowLFVectorfields_[inflowFaceIndex];

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = geometry.global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, geometry.jacobianInverseTransposed(qp.position()));

      Dune::FieldMatrix<K, Model::m, Model::dim> gradient(0.0);
      Dune::FieldMatrix<K, Model::m, Model::dim> extendedGradient(0.0);

      for (int k = 0; k < Model::m; ++k) {
        for (int l = 0; l < localView.tree().child(k).size(); ++l) {
          gradient[k] += x[localView.tree().child(k).localIndex(l)] * gradphi[l][0];
          extendedGradient[k] += x_extended[localView.tree().child(k).localIndex(l)] * gradphi[l][0];
        }
      }

      const K factor = qp.weight();

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < gradphi.size(); ++j) {
          LocalVector y(0.0);
          y[localView.tree().child(i).localIndex(j)] = 1.0;
          Dune::FieldVector<K, Model::m> v(0.0);
          fillInFunctionValue(localView, y, phi, v, std::integral_constant<int, Model::m>());
          auto flux = model_.flux(qpg, v);
          flux *= (1.0 - capacity) * weightFactor;

          for (std::size_t k = 0; k < Model::m; ++k) {
            r[localView.tree().child(i).localIndex(j)] -= flux[k] * extendedGradient[k] * factor;
          }

          r[localView.tree().child(i).localIndex(j)] -= (1.0 - capacity) * v[i] * 0.5 * (vf * extendedGradient[i]) * factor;
        }
      }
    }
  }

  template<class CutCell, class LocalView>
  void dodVolumeBoundary(const CutCell& cutCell, const LocalView& localView,
                         const LocalVector& x, const LocalVector& x_extended,
                         std::size_t inflowFaceIndex, std::size_t ghostCellFaceIndex,
                         K capacity, K weightFactor,
                         LocalVector& r) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    const auto& geometry = cutCell.geometry();

    auto vf = inflowLFVectorfields_[inflowFaceIndex];
    auto ghostVf = inflowReflectedLFVectorfields_[ghostCellFaceIndex];
    auto nu = -(ghostVf * reflectionVector_);

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order() + intorderadd_);

    for (auto&& qp : quadratureRule_) {
      auto qpg = geometry.global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, geometry.jacobianInverseTransposed(qp.position()));

      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u(0.0);

      fillInFunctionValue(localView, x_extended, phi, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(localView, x, phi, u, std::integral_constant<int, Model::m>());

      const K factor = qp.weight();

      auto flux = model_.flux(qpg, u_extended) - model_.flux(qpg, u);
      flux *= (1.0 - capacity) * weightFactor;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < gradphi.size(); ++j) {
          r[localView.tree().child(i).localIndex(j)] -= flux[i] * gradphi[j][0] * factor;
          r[localView.tree().child(i).localIndex(j)] -= (1.0 - capacity) * (u_extended[i] - u[i]) * 0.5 * (vf * gradphi[j][0]) * nu * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceCentralPart(const CutIntersection& cutIntersection,
                          const LocalView& insideView,
                          const LocalVector& x_s, const LocalVector& x_n,
                          const LocalVector& x_extended,
                          std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                          K capacity, K weightFactor,
                          LocalVector& r_s, LocalVector& r_n,
                          bool smallCutCellIsOutside) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qppos_s = geometryInInside.local(qp.position());
      const auto qpg = qp.position();

      insideFiniteElement.localBasis().evaluateFunction(qppos_s, phi_s);

      Dune::FieldVector<K, Model::m> u_s(0.0);
      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u_n(0.0);

      fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended, phi_s, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_n, phi_s, u_n, std::integral_constant<int, Model::m>());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      const auto factor = qp.weight();

      auto centralFlux = model_.centralFlux(normal, u_extended, u_n, qpg);
      centralFlux *= (1.0 - capacity) * weightFactor;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s[j] * centralFlux[i] * factor;
          r_n[insideView.tree().child(i).localIndex(j)] += -phi_s[j] * centralFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceDiffusivePart(const CutIntersection& cutIntersection,
                            const LocalView& insideView,
                            const LocalVector& x_s, const LocalVector& x_n,
                            const LocalVector& x_extended,
                            std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                            K capacity,
                            LocalVector& r_s, LocalVector& r_n,
                            bool smallCutCellIsOutside) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    auto vf = inflowLFVectorfields_[inflowFaceIndex];

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (auto&& qp : quadratureRule_) {
      const auto qpGlobal = qp.position();
      const auto qpInside = geometryInInside.local(qp.position());

      insideFiniteElement.localBasis().evaluateFunction(qpInside, phi_s);

      Dune::FieldVector<K, Model::m> u_s(0.0);
      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u_n(0.0);

      fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended, phi_s, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_n, phi_s, u_n, std::integral_constant<int, Model::m>());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      const auto factor = qp.weight();

      auto diffusiveFlux = (u_extended - u_n);
      diffusiveFlux *= 0.5 * (1.0 - capacity) * (vf * normal);

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s[j] * diffusiveFlux[i] * factor;
          r_n[insideView.tree().child(i).localIndex(j)] += -phi_s[j] * diffusiveFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceCentralPartHalfFlow(const CutIntersection& cutIntersection,
                                  const LocalView& insideView,
                                  const LocalVector& x,
                                  std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                                  K capacity, K weightFactor,
                                  LocalVector& r,
                                  bool smallCutCellIsOutside) const
  {
    const auto& finiteElement = insideView.tree().child(0).finiteElement();
    int intorder = 2 * finiteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(finiteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qppos_s = geometryInInside.local(qp.position());
      const auto qpg = qp.position();

      finiteElement.localBasis().evaluateFunction(qppos_s, phi);

      Dune::FieldVector<K, Model::m> u(0.0);
      fillInFunctionValue(insideView, x, phi, u, std::integral_constant<int, Model::m>());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      const auto factor = qp.weight();

      auto centralFlux = model_.centralFluxHalf(normal, u, qpg);
      centralFlux *= (1.0 - capacity) * weightFactor;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r[insideView.tree().child(i).localIndex(j)] += phi[j] * centralFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceCentralPartReflected(const CutIntersection& cutIntersection,
                                   const LocalView& insideView,
                                   const LocalVector& x_s, const LocalVector& x_n,
                                   const LocalVector& x_extended,
                                   std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                                   K capacity,
                                   K weightFactor,
                                   LocalVector& r_s, LocalVector& r_n,
                                   bool smallCutCellIsOutside) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qpGlobal = reflectPoint(qp.position());
      const auto qpInside = geometryInInside.local(qpGlobal);

      auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      normal = reflectVector(normal);

      insideFiniteElement.localBasis().evaluateFunction(qpInside, phi_s);

      Dune::FieldVector<K, Model::m> u_s(0.0);
      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u_n(0.0);

      fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended, phi_s, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_n, phi_s, u_n, std::integral_constant<int, Model::m>());

      const auto factor = qp.weight();

      // auto flux = model_.centralFlux(normal, u_extended, u_n, qpGlobal) - model_.centralFlux(normal, u_s, u_n, qpGlobal);
      auto flux = model_.centralFlux(normal, u_extended, u_n, qpGlobal);
      flux *= (1.0 - capacity) * weightFactor;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s[j] * flux[i] * factor;
          r_n[insideView.tree().child(i).localIndex(j)] += -phi_s[j] * flux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceCentralPartHalfFlowReflected(const CutIntersection& cutIntersection,
                                           const LocalView& insideView,
                                           const LocalVector& x,
                                           std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                                           K capacity, K weightFactor,
                                           LocalVector& r,
                                           bool smallCutCellIsOutside) const
  {
    const auto& finiteElement = insideView.tree().child(0).finiteElement();
    int intorder = 2 * finiteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(finiteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& cutCellGeometry = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qpGlobal = reflectPoint(qp.position());
      const auto qpInside = cutCellGeometry.local(qpGlobal);

      auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      normal = reflectVector(normal);

      finiteElement.localBasis().evaluateFunction(qpInside, phi);

      Dune::FieldVector<K, Model::m> u(0.0);
      fillInFunctionValue(insideView, x, phi, u, std::integral_constant<int, Model::m>());

      const auto factor = qp.weight();

      auto centralFlux = model_.centralFluxHalf(normal, u, qpGlobal);
      centralFlux *= (1.0 - capacity) * weightFactor;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r[insideView.tree().child(i).localIndex(j)] -= phi[j] * centralFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceDiffusivePartHalfFlow(const CutIntersection& cutIntersection,
                                    const LocalView& insideView,
                                    const LocalVector& x,
                                    std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                                    K capacity,
                                    LocalVector& r,
                                    bool smallCutCellIsOutside) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(insideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    auto vf = inflowLFVectorfields_[inflowFaceIndex];

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (auto&& qp : quadratureRule_) {
      const auto qpGlobal = qp.position();
      const auto qpInside = geometryInInside.local(qp.position());

      insideFiniteElement.localBasis().evaluateFunction(qpInside, phi);

      Dune::FieldVector<K, Model::m> u(0.0);
      fillInFunctionValue(insideView, x, phi, u, std::integral_constant<int, Model::m>());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());

      const auto factor = qp.weight();

      auto diffusiveFlux = u;
      diffusiveFlux *= 0.5 * (1.0 - capacity) * (vf * normal);

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r[insideView.tree().child(i).localIndex(j)] += phi[j] * diffusiveFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceDiffusivePartHalfFlowReflected(const CutIntersection& cutIntersection,
                                             const LocalView& insideView,
                                             const LocalVector& x,
                                             std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                                             K capacity,
                                             LocalVector& r,
                                             bool smallCutCellIsOutside) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi(insideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& cutCellGeometry = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    auto vf = inflowReflectedLFVectorfields_[inflowFaceIndex];

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;
    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (auto&& qp : quadratureRule_) {
      const auto qpGlobal = qp.position();
      const auto qpInside = cutCellGeometry.local(qpGlobal);

      auto normal = sign * cutIntersection.unitOuterNormal(qp.position());

      insideFiniteElement.localBasis().evaluateFunction(qpInside, phi);

      Dune::FieldVector<K, Model::m> u(0.0);
      fillInFunctionValue(insideView, x, phi, u, std::integral_constant<int, Model::m>());

      const auto factor = qp.weight();

      auto diffusiveFlux = u;
      diffusiveFlux *= 0.5 * (1.0 - capacity) * (vf * normal);

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r[insideView.tree().child(i).localIndex(j)] += phi[j] * diffusiveFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodFaceBoundaryDiffusivePart(const CutIntersection& cutIntersection,
                                    const LocalView& insideView,
                                    const LocalVector& x_s, const LocalVector& x_n,
                                    const LocalVector& x_extended,
                                    std::size_t outflowFaceIndex, std::size_t inflowFaceIndex,
                                    std::size_t ghostFaceIndex, K capacity,
                                    LocalVector& r_s, LocalVector& r_n,
                                    bool smallCutCellIsOutside) const
  {
    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    auto vf = inflowLFVectorfields_[inflowFaceIndex];
    auto ghostVf = inflowReflectedLFVectorfields_[ghostFaceIndex];

    auto nu = (ghostVf * reflectionVector_) / (vf * reflectionVector_);

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (auto&& qp : quadratureRule_) {
      const auto qpGlobal = qp.position();
      const auto qpInside = geometryInInside.local(qp.position());

      insideFiniteElement.localBasis().evaluateFunction(qpInside, phi_s);

      Dune::FieldVector<K, Model::m> u_s(0.0);
      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u_n(0.0);

      fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended, phi_s, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_n, phi_s, u_n, std::integral_constant<int, Model::m>());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      const auto factor = qp.weight();
      auto diffusiveFlux = (u_extended - u_n);
      diffusiveFlux *= 0.5 * (1.0 - capacity) * (vf * normal) * nu;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r_s[insideView.tree().child(i).localIndex(j)] += phi_s[j] * diffusiveFlux[i] * factor;
          r_n[insideView.tree().child(i).localIndex(j)] += -phi_s[j] * diffusiveFlux[i] * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView>
  void dodBoundary(const CutIntersection& cutIntersection, const LocalView& insideView,
                   const LocalVector& x, const LocalVector& x_extended,
                   const LocalVector& bc_coefficients, const LocalVector& bc_extended_coefficients,
                   std::size_t outflowFaceIndex, std::size_t inflowFaceIndex, K capacity,
                   LocalVector& r)
  {
    LocalVector x_transformed(0.0);
    LocalVector x_extended_transformed(0.0);
    LocalVector bc_transformed_coefficients(0.0);

    auto weight = weights_[inflowFaceIndex];

    weight.mv(x, x_transformed);
    weight.mv(x_extended, x_extended_transformed);
    weight.mv(bc_coefficients, bc_transformed_coefficients);

    const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

    int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

    using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi(insideFiniteElement.localBasis().size());
    std::vector<Dune::FieldVector<K, Model::m>> phi_transformed(insideView.tree().size());

    const auto& geometry = cutIntersection.geometry();
    const auto& geometryInInside = cutIntersection.geometryInInside();

    cutIntersection.quadratureRule(quadratureRule_, intorder);

    for (const auto& qp : quadratureRule_) {
      const auto qppos_s = geometryInInside.local(qp.position());
      const auto qpg = qp.position();

      insideFiniteElement.localBasis().evaluateFunction(qppos_s, phi);

      Dune::FieldVector<K, Model::m> u(0.0);
      Dune::FieldVector<K, Model::m> u_transformed(0.0);
      Dune::FieldVector<K, Model::m> u_extended(0.0);
      Dune::FieldVector<K, Model::m> u_extended_transformed(0.0);
      Dune::FieldVector<K, Model::m> bc_transformed(0.0);

      fillInFunctionValue(insideView, x, phi, u, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_transformed, phi, u_transformed, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended, phi, u_extended, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, x_extended_transformed, phi, u_extended_transformed, std::integral_constant<int, Model::m>());
      fillInFunctionValue(insideView, bc_transformed_coefficients, phi, bc_transformed, std::integral_constant<int, Model::m>());

      const auto normal = cutIntersection.unitOuterNormal(qp.position());
      const auto factor = qp.weight();

      auto flux = model_.numericalFlux(normal, u_extended_transformed, bc_transformed, qpg) - model_.numericalFlux(normal, u_transformed, bc_transformed, qpg);
      flux *= (1.0 - capacity) * uWeightFactor_;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
          r[insideView.tree().child(i).localIndex(j)] += phi[j] * flux[i] * factor;
        }
      }

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); j++) {
          phi_transformed[insideView.tree().child(i).localIndex(j)] = Dune::FieldVector<K, Model::m>(0.0);

          for (std::size_t k = 0; k < Model::m; ++k) {
            for (std::size_t l = 0; l < insideView.tree().child(k).size(); ++l) {
              phi_transformed[insideView.tree().child(i).localIndex(j)][k] += weight[insideView.tree().child(k).localIndex(l)][insideView.tree().child(i).localIndex(j)] * phi[l];
            }
          }
        }
      }

      const auto bc = boundaryCondition_(u, qpg, normal, t_);

      flux = model_.numericalFlux(normal, u_extended, bc, qpg) - model_.numericalFlux(normal, u, bc, qpg);
      flux *= (1.0 - capacity) * vWeightFactor_;

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < insideView.tree().child(i).size(); j++) {
          r[insideView.tree().child(i).localIndex(j)] += phi_transformed[insideView.tree().child(i).localIndex(j)] * flux * factor;
        }
      }
    }
  }

  template<class CutIntersection, class LocalView, class Flux>
  void computeFlowOperator(const CutIntersection& cutIntersection, LocalMatrix& out, const LocalView& localView, const Flux& flux, bool smallCutCellIsOutside) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());

    const auto cutCellGeometry = smallCutCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();
    cutIntersection.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order());
    const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    for (auto&& qp : quadratureRule_) {
      const auto qpGlobal = qp.position();
      const auto qpInside = cutCellGeometry.local(qp.position());

      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());

      finiteElement.localBasis().evaluateFunction(qpInside, phi);
      const auto factor = qp.weight();

      for (int j = 0; j < finiteElement.localBasis().size(); ++j) {
        for (int i = 0; i < Model::m; ++i) {
          Dune::FieldVector<K, Model::m> u(0.0);
          u[i] = phi[j];
          auto v = flux(normal, u, qpGlobal);

          for (int k = 0; k < finiteElement.localBasis().size(); ++k) {
            for (int l = 0; l < Model::m; ++l) {
              out[localView.tree().child(l).localIndex(k)][localView.tree().child(i).localIndex(j)] += factor * v[l] * phi[k];
            }
          }
        }
      }
    }
  }

  template<class Ip, class LocalView>
  void computeOutflowOperator(const Ip& ip, LocalMatrix& out, const LocalView& localView, bool smallCutCellIsOutside) const
  {
    computeFlowOperator(ip, out, localView, [this] (const auto& n, const auto& u, const auto& x) { return model_.outflow(n, u, x); }, smallCutCellIsOutside);
  }

  template<class Ip, class LocalView>
  void computeInflowOperator(const Ip& ip, LocalMatrix& out, const LocalView& localView, bool smallCutCellIsOutside) const
  {
    computeFlowOperator(ip, out, localView, [this] (const auto& n, const auto& u, const auto& x) { return model_.inflow(n, u, x); }, smallCutCellIsOutside);
  }

  template<class CutCell, class LocalView>
  void computeVolumeOperator(const CutCell& cutCell, const LocalView& localView, LocalMatrix& out) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    cutCell.quadratureRule(quadratureRule_, 2*finiteElement.localBasis().order());

    for (auto&& qp : quadratureRule_) {
      auto qpg = cutCell.geometry().global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, cutCell.geometry().jacobianInverseTransposed(qp.position()));

      const K factor = qp.weight();

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < finiteElement.localBasis().size(); ++j) {
          Dune::FieldVector<K, Model::m> u(0.0);
          u[i] = phi[j];
          const auto flux = model_.flux(qpg, u);

          for (std::size_t k = 0; k < Model::m; ++k) {
            for (std::size_t l = 0; l < finiteElement.localBasis().size(); ++l) {
              out[localView.tree().child(i).localIndex(j)][localView.tree().child(k).localIndex(l)] -= flux[k] * gradphi[l][0] * factor;
            }
          }
        }
      }
    }
  }

  template<class CutCell, class LocalView, class P0Weight, class Out>
  void computeResidualOperator(const CutCell& cutCell, const LocalView& localView, const P0Weight& uWeight, const P0Weight& vWeight, Out& out) const
  {
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType> phi(finiteElement.localBasis().size());
    std::vector<typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::JacobianType> gradphi(finiteElement.localBasis().size());

    cutCell.quadratureRule(quadratureRule_, 2 * finiteElement.localBasis().order());

    for (auto&& qp : quadratureRule_) {
      auto qpg = cutCell.geometry().global(qp.position());

      finiteElement.localBasis().evaluateFunction(qp.position(), phi);
      finiteElement.localBasis().evaluateJacobian(qp.position(), gradphi);
      transformGradients(gradphi, cutCell.geometry().jacobianInverseTransposed(qp.position()));

      const auto factor = qp.weight();

      for (std::size_t i = 0; i < Model::m; ++i) {
        for (std::size_t j = 0; j < finiteElement.localBasis().size(); ++j) {
          Dune::FieldVector<K, Model::m> u(0.0);
          u[i] = phi[j];
          Dune::FieldVector<K, Model::m> uTransformed(0.0);
          uWeight.mv(u, uTransformed);
          const auto flux = model_.flux(qpg, uTransformed);

          for (std::size_t k = 0; k < Model::m; ++k) {
            for (std::size_t l = 0; l < finiteElement.localBasis().size(); ++l) {
              for (int p = 0; p < Model::m; ++p) {
                out[localView.tree().child(k).localIndex(l)][localView.tree().child(i).localIndex(j)] -= flux[p] * gradphi[l][0] * factor * vWeight[p][k];
              }
            }
          }
        }
      }
    }
  }

  template<class CutIntersection, class FlowMatrix>
  void computeInflowP0Operator(const CutIntersection& cutIntersection, FlowMatrix& out, bool smallCutCellIsOutside) const
  {
    const int intorder = 5;
    cutIntersection.quadratureRule(quadratureRule_, intorder);
    K sign = smallCutCellIsOutside ? -1.0 : 1.0;

    for (auto&& qp : quadratureRule_) {
      const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
      const auto qpGlobal = qp.position();
      const auto inflow = model_.inflowMatrix(normal, qpGlobal);
      const auto factor = qp.weight();
      out += inflow * factor;
    }
  }

  template<class CutIntersectionIt, class CutCellInfo>
  K capacity(CutIntersectionIt it, CutIntersectionIt end, const CutCellInfo& info, int domainIndex) const
  {
    const auto cutCellVolume = info.volume;
    const K exp = 1. / Model::dim;
    const auto h = std::pow(info.boundingBox.entity().geometry().volume(), exp);

    std::vector<K> flow;

    for (; it != end; ++it) {
      flow.resize(flow.size() + Model::m);
      std::fill(flow.end() - Model::m, flow.end(), 0.0);

      if ((it->insideDomainIndex() == domainIndex) || (it->outsideDomainIndex() == domainIndex)) {
        const auto& geo = it->geometry();

        it->quadratureRule(quadratureRule_, 2 + intorderadd_);

        for (auto&& qp : quadratureRule_) {
          auto normal = it->unitOuterNormal(qp.position());

          if (it->neighbor() && it->domainIsOutside(domainIndex)) {
            normal *= -1.0;
          }

          auto eigenvalues = model_.eigenvalues(typename Model::State(0.0), normal, qp.position());
          const auto factor = qp.weight();

          for (std::size_t i = 0; i < Model::m; ++i) {
            *(flow.end() - (Model::m - i)) += std::abs(eigenvalues[i]) * factor;
          }
        }
      }
    }

    auto maxInflow = std::max_element(flow.begin(), flow.end());
    const K one = 1.0;
    // return std::min((cutCellVolume / (*maxInflow * dt_)) * capacityFactor_, one);
    K alpha = cutCellVolume / info.boundingBox.entity().geometry().volume();
    return std::min(alpha / tau_, one);
  }

  K capacity(K maxFlow, K cutCellVolume, K h) const
  {
    return std::min((cutCellVolume / (maxFlow * cfl_ * h)) * capacityFactor_, 1.0);
  }

  template<class CutIntersection, class LocalView>
  K maxInflow(CutIntersection cutIntersection, const LocalView& localView, const LocalVector& x, int domainIndex) const
  {
    bool smallCellIsOutside = cutIntersection.neighbor() && cutIntersection.domainIsOutside(domainIndex);

    const auto& geo = cutIntersection.geometry();
    const auto& geometryInNeighbor = smallCellIsOutside ? cutIntersection.geometryInOutside() : cutIntersection.geometryInInside();

    cutIntersection.quadratureRule(quadratureRule_, 2 + intorderadd_);

    const auto& finiteElement = localView.tree().child(0).finiteElement();
    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;
    std::vector<RangeType> phi_n(finiteElement.localBasis().size());

    Dune::FieldVector<K, Model::m> flow(0.0);

    for (auto&& qp : quadratureRule_) {
      auto normal = cutIntersection.unitOuterNormal(qp.position());

      if (smallCellIsOutside) {
        normal *= -1.0;
      }

      Dune::FieldVector<K, Model::m> u(0.0);
      finiteElement.localBasis().evaluateFunction(geometryInNeighbor.local(qp.position()), phi_n);
      fillInFunctionValue(localView, x, phi_n, u, std::integral_constant<int, Model::m>());

      auto eigenvalues = model_.eigenvalues(u, normal, qp.position());
      const auto factor = qp.weight();

      for (std::size_t i = 0; i < Model::m; ++i) {
        flow[i] += eigenvalues[i] * factor;
      }
    }

    return flow.infinity_norm();
  }

  template<class CutIntersection, class LocalView>
  K maxInflowOnBoundary(CutIntersection cutIntersection, const LocalView& localView, const LocalVector& x, int domainIndex) const
  {
    cutIntersection.quadratureRule(quadratureRule_, 2 + intorderadd_);

    const auto& geo = cutIntersection.geometry();
    const auto& geometryInInside = cutIntersection.geometryInInside();
    const auto& finiteElement = localView.tree().child(0).finiteElement();

    using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;

    std::vector<RangeType> phi_s(finiteElement.localBasis().size());

    Dune::FieldVector<K, Model::m> flow(0.0);

    for (auto&& qp : quadratureRule_) {
      auto normal = cutIntersection.unitOuterNormal(qp.position());

      Dune::FieldVector<K, Model::m> u(0.0);
      finiteElement.localBasis().evaluateFunction(geometryInInside.local(qp.position()), phi_s);
      fillInFunctionValue(localView, x, phi_s, u, std::integral_constant<int, Model::m>());
      const auto bc = boundaryCondition_(u, qp.position(), normal, t_);

      auto eigenvalues = model_.eigenvalues(bc, normal, qp.position());
      const auto factor = qp.weight();

      for (std::size_t i = 0; i < Model::m; ++i) {
        flow[i] += eigenvalues[i] * factor;
      }
    }

    return flow.infinity_norm();
  }

  void computeLFStabilityWeights(const std::vector<Dune::FieldVector<K, 2>>& normals, const std::vector<Dune::FieldVector<K, 2>>& tangents, const std::vector<K>& waveSpeeds, Dune::FieldVector<K, 3>& LFStabilitySolution) const
  {
    Dune::FieldMatrix<K, 3, 3> LFStabilitySystem(0.0);
    Dune::FieldVector<K, 3> LFStabilityRHS(0.0);

    LFStabilitySystem[0][1] = normals[0] * tangents[1];
    LFStabilitySystem[0][2] = normals[0] * tangents[2];
    LFStabilitySystem[1][0] = normals[1] * tangents[0];
    LFStabilitySystem[1][2] = normals[1] * tangents[2];
    LFStabilitySystem[2][0] = normals[2] * tangents[0];
    LFStabilitySystem[2][1] = normals[2] * tangents[1];

    auto rhsVF = waveSpeeds[0] * normals[0] + waveSpeeds[1] * normals[1] + waveSpeeds[2] * normals[2];

    LFStabilityRHS[0] = rhsVF * normals[0];
    LFStabilityRHS[1] = rhsVF * normals[1];
    LFStabilityRHS[2] = rhsVF * normals[2];

    auto tmp = LFStabilitySystem;

    pseudoInverse(LFStabilitySystem, LFStabilitySystem);

    LFStabilitySystem.mv(LFStabilityRHS, LFStabilitySolution);

    Dune::FieldVector<K, 3> residual(0.0);
    tmp.mv(LFStabilitySolution, residual);
    residual -= LFStabilityRHS;

    if (residual.infinity_norm() > 1e-14) {
      std::cout << "RHS not in image space" << std::endl;
    }
  }

  auto boundaryValue() const
  {
    // pass a zero normal vector for now, we figure out how to do it properly later
    return [=, this](const auto& u, const auto& x) { return boundaryCondition_(u, x, Dune::FieldVector<K, Model::dim>(0.0), t_); };
  }

  void setTime(K t)
  {
    t_ = t;
  }

  typename Model::Coordinate reflectPoint(const Model::Coordinate& x) const
  {
    return x + 2.0 * ((reflectionCenter_ - x) * reflectionVector_) * reflectionVector_;
  }

  typename Model::Coordinate reflectVector(const Model::Coordinate& v) const
  {
    return v - 2.0 * (v * reflectionVector_) * reflectionVector_;
  }

  typename Model::Coordinate reflectionVector() const
  {
    return reflectionVector_;
  }

  bool useMatrixStabilization() const
  {
    return useMatrixStabilization_;
  }

  void setTimestepSize(K dt)
  {
    dt_ = dt;
  }

private:
  Model model_;
  BoundaryCondition boundaryCondition_;

  mutable std::vector<QuadraturePoint> quadratureRule_;
  mutable std::vector<LocalMatrix> weights_;
  mutable std::vector<typename Model::Coordinate> inflowLFVectorfields_;
  mutable std::vector<typename Model::Coordinate> inflowReflectedLFVectorfields_;
  mutable typename Model::Coordinate reflectionVector_;
  mutable typename Model::Coordinate reflectionCenter_;

  K uWeightFactor_;
  K vWeightFactor_;
  K capacityFactor_;
  K dt_;
  K cfl_;
  K tau_;
  K t_;

  int intorderadd_;

  bool useMatrixStabilization_;
};

}

#endif
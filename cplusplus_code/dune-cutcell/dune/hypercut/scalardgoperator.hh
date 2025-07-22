#ifndef DUNE_HYPERCUT_SCALARDGOPERATOR_HH
#define DUNE_HYPERCUT_SCALARDGOPERATOR_HH

#include <cassert>
#include <type_traits>
#include <vector>

#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/hypercut/helper.hh>
#include <dune/hypercut/interpolation.hh>

#include <dune/subtriangulation/simpletpmctriangulation/localsubtriangulation.hh>

namespace Dune::Hypercut
{

template<class LST, class GlobalBasis, class VectorType, class LocalOp, class SmallCutCellClassifier>
class ScalarDGOperator
{
public:
  using LocalSubTriangulation = LST;
  using LocalOperator = LocalOp;
  using K = typename LST::ctype;

  using LocalVector = typename LocalOperator::LocalVector;
  using LocalMatrix = typename LocalOperator::LocalMatrix;

  using Domain = VectorType;

  static const int ComponentSpaceDimension = LocalOperator::ComponentSpaceDimension;

  ScalarDGOperator(LocalSubTriangulation& localSubTriangulation,
                 const GlobalBasis& globalBasis,
                 const LocalOperator& localOperator,
                 SmallCutCellClassifier& classifier,
                 K l2StabilityFactor = 5.0,
                 bool useBoundaryHack = false,
                 bool ignoreSmallCells = false) :
    localSubTriangulation_(localSubTriangulation),
    globalBasis_(globalBasis),
    localOperator_(localOperator),
    classifier_(classifier),
    l2StabilityFactor_(l2StabilityFactor),
    useBoundaryHack_(useBoundaryHack),
    ignoreSmallCells_(false)
  {
    for (int domainIndex = 0; domainIndex < localSubTriangulation_.domainConfiguration().numberOfDomains(); ++domainIndex) {
      insideViews_.push_back(globalBasis_.localView(domainIndex));
      outsideViews_.push_back(globalBasis_.localView(domainIndex));
      extendedViews_.push_back(globalBasis_.localView(domainIndex));
    }
  }

  template<class Vector>
  void operator()(K t, const Vector& x, Vector& y)
  {
    auto yBackend = Dune::Functions::istlVectorBackend(y);
    assemble(t, Dune::Functions::istlVectorBackend(x), yBackend);
  }

  template<class X, class Y>
  void assemble(K t, const X& x, const Y &y)
  {
    assembleSpaceDiscretization(t, x, y);
    applyInvertedMassMatrix(y);
  }

  template<class X, class Y>
  void assembleSpaceDiscretization(K t, const X& x, const Y &y)
  {
    localOperator_.setTime(t);

    classifier_.classify(localSubTriangulation_);

    const auto& indexSet = localSubTriangulation_.gridView().indexSet();

    // assemble the space discretization form
    for (const auto& element : elements(localSubTriangulation_.gridView())) {
      localSubTriangulation_.bind(element);
      localSubTriangulation_.createCutCells();
      localSubTriangulation_.createCutIntersections();

      std::vector<typename LocalOperator::LocalMatrix> corrections;

      typename LocalSubTriangulation::BoundingBox smallCellBoundingBox;
      Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> smallCellMassMatrix(0.0);

      // TODO: Find a proper place for this together with a proper handling
      K capacity = 1.0;

      bool smallCellIsBoundaryCell = false;

      for (auto cutCellIt = localSubTriangulation_.cutCellsBegin();
           cutCellIt != localSubTriangulation_.cutCellsEnd(); ++cutCellIt) {
        int currentDomainIndex = cutCellIt->domainIndex();

        auto& insideView = insideViews_[currentDomainIndex];
        insideView.bind(element);

        LocalVector insideLocalValue(0.0);
        LocalVector insideLocalResult(0.0);

        readLocalVector(x, insideView, insideLocalValue);

        bool isSmallCell = classifier_.isSmallCell(indexSet.index(element), currentDomainIndex);

        if (!ignoreSmallCells_ || !isSmallCell) {
          localOperator_.dgVolume(*cutCellIt, insideView, insideLocalValue, insideLocalResult);
        }

        if (isSmallCell) {
          // TODO: Same applies to the capacity
          assert(capacity == 1.0);

          capacity = localOperator_.capacity(localSubTriangulation_.cutIntersectionBegin(), localSubTriangulation_.cutIntersectionEnd(), localSubTriangulation_.cutCellInformation().information(element, currentDomainIndex), currentDomainIndex);

          localOperator_.computeWeights(*cutCellIt, localSubTriangulation_.cutIntersectionBegin(), localSubTriangulation_.cutIntersectionEnd(), insideView, currentDomainIndex, corrections);

          smallCellBoundingBox = localSubTriangulation_.cutCellInformation().information(element, currentDomainIndex).boundingBox;

          computeLocalMassMatrix(smallCellBoundingBox, insideView.tree().child(0).finiteElement().localBasis(), smallCellMassMatrix);
          smallCellMassMatrix.invert();

          assembleVolumeStabilization(cutCellIt, x, y, insideView, insideLocalValue, smallCellBoundingBox, smallCellMassMatrix, currentDomainIndex, capacity, insideLocalResult);
        }

        addLocalVector(y, insideView, insideLocalResult);
      }

      auto correctionIt = corrections.begin();
      std::size_t smallCellOutflowFaceIndex = 0;

      for (auto cutIntersectionIt = localSubTriangulation_.cutIntersectionBegin();
        cutIntersectionIt != localSubTriangulation_.cutIntersectionEnd(); ++cutIntersectionIt) {
        LocalVector insideLocalValue(0.0);
        LocalVector insideLocalResult(0.0);

        int insideDomainIndex = cutIntersectionIt->insideDomainIndex();

        auto& insideView = insideViews_[insideDomainIndex];
        insideView.bind(cutIntersectionIt->inside());

        readLocalVector(x, insideView, insideLocalValue);

        if (cutIntersectionIt->neighbor()) {
          int index = indexSet.index(cutIntersectionIt->outside());

          bool smallCutCellInside = classifier_.isSmallCell(indexSet.index(element), cutIntersectionIt->insideDomainIndex());
          bool smallCutCellOutside = cutIntersectionIt->insideDomainIndex() != cutIntersectionIt->outsideDomainIndex() && classifier_.isSmallCell(indexSet.index(element), cutIntersectionIt->outsideDomainIndex());
          bool neighborSmall = classifier_.isSmallCell(indexSet.index(cutIntersectionIt->outside()), cutIntersectionIt->outsideDomainIndex());

          if (visitFace(indexSet.index(cutIntersectionIt->inside()), index, cutIntersectionIt->insideDomainIndex(), cutIntersectionIt->outsideDomainIndex())) {
            auto& outsideView = outsideViews_[cutIntersectionIt->outsideDomainIndex()];
            outsideView.bind(cutIntersectionIt->outside());

            LocalVector outsideLocalValue(0.0);
            LocalVector outsideLocalResult(0.0);

            readLocalVector(x, outsideView, outsideLocalValue);

            if (!ignoreSmallCells_ || (!smallCutCellInside && !smallCutCellOutside && !neighborSmall)) {
              localOperator_.dgFace(*cutIntersectionIt, insideView, outsideView, insideLocalValue,
                                    outsideLocalValue, insideLocalResult, outsideLocalResult);
            }

            addLocalVector(y, outsideView, outsideLocalResult);
          }

          if (smallCutCellInside || smallCutCellOutside) {
            assembleFaceStabilization(cutIntersectionIt, smallCellOutflowFaceIndex, smallCellBoundingBox, smallCellMassMatrix, smallCutCellInside, x, y, insideView, insideLocalValue, correctionIt, capacity, insideLocalResult);

            ++smallCellOutflowFaceIndex;
          }
        } else {
          bool isSmallCell = classifier_.isSmallCell(indexSet.index(element), cutIntersectionIt->insideDomainIndex());

          if (!ignoreSmallCells_ || !isSmallCell) {
            localOperator_.dgBoundary(*cutIntersectionIt, insideView, insideLocalValue, insideLocalResult);
          }

          if (isSmallCell) {
            assembleBoundaryStabilization(cutIntersectionIt, smallCellOutflowFaceIndex, smallCellBoundingBox, smallCellMassMatrix, x, y, insideView, insideLocalValue, capacity, insideLocalResult);
            ++smallCellOutflowFaceIndex;
          }
        }

        addLocalVector(y, insideView, insideLocalResult);
      }
    }
  }

  template<class CutCellIt, class X, class Y, class LocalView, class BoundingBox>
  void assembleVolumeStabilization(CutCellIt cutCellIt, const X& x, Y& y, const LocalView& insideView, const typename LocalOperator::LocalVector& insideLocalValue, const BoundingBox& boundingBox, const Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension>& massMatrix, int currentDomainIndex, K capacity, typename LocalOperator::LocalVector& insideLocalResult)
  {
    std::size_t inflowFaceIndex = 0;

    for (auto smallCellIntersectionIt = localSubTriangulation_.cutIntersectionBegin();
      smallCellIntersectionIt != localSubTriangulation_.cutIntersectionEnd(); ++smallCellIntersectionIt) {
      if (smallCellIntersectionIt->insideDomainIndex() == currentDomainIndex || smallCellIntersectionIt->outsideDomainIndex() == currentDomainIndex) {
        typename LocalOperator::LocalVector extendedLocalValue(0.0);
        typename LocalOperator::LocalVector extendedLocalResult(0.0);

        Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> extendedBaseChange(0.0);

        int inflowNeighborDomainIndex = -1;

        std::vector<typename LocalOperator::LocalVector> extendedBasisCoefficients;

        // TODO: Fix this nested structure, we need a proper ordering of the faces
        if (smallCellIntersectionIt->neighbor()) {
          loadExtendedData(*smallCellIntersectionIt, boundingBox, currentDomainIndex, insideView, massMatrix, inflowNeighborDomainIndex, x, extendedLocalValue, extendedBaseChange);

          extendedBasisCoefficients.resize(extendedViews_[inflowNeighborDomainIndex].tree().size());

          for (int i = 0; i < LocalOperator::Model::m; ++i) {
            for (int j = 0; j < extendedViews_[inflowNeighborDomainIndex].tree().child(i).size(); ++j) {
              typename LocalOperator::LocalVector basisVector(0.0);
              basisVector[extendedViews_[inflowNeighborDomainIndex].tree().child(i).localIndex(j)] = 1.0;
              extendedBasisCoefficients[extendedViews_[inflowNeighborDomainIndex].tree().child(i).localIndex(j)] = transformBasisComponents(basisVector, extendedBaseChange, insideView, extendedViews_[inflowNeighborDomainIndex], std::integral_constant<int, LocalOperator::Model::m>());
            }
          }
        } else {
          assert(false);
          extendedBasisCoefficients.resize(0);
          interpolateFunctionLocal(localOperator_.boundaryValue(), insideLocalValue, boundingBox, insideView, massMatrix, extendedLocalValue, 2, std::integral_constant<int, LocalOperator::Model::m>());
        }

        typename LocalOperator::LocalVector insideLocalResultTmp(0.0);

        localOperator_.dodVolume(*cutCellIt, insideView, extendedViews_[inflowNeighborDomainIndex], insideLocalValue, extendedLocalValue, inflowFaceIndex, extendedBasisCoefficients, capacity, insideLocalResultTmp, extendedLocalResult);

        insideLocalResult += insideLocalResultTmp;

        if (smallCellIntersectionIt->neighbor()) {
          addLocalVector(y, extendedViews_[inflowNeighborDomainIndex], extendedLocalResult);
        }

        ++inflowFaceIndex;
      }
    }
  }

  template<class CutIntersectionIt, class X, class Y, class LocalView, class CorrectionIt>
  void assembleFaceStabilization(CutIntersectionIt cutIntersectionIt, std::size_t outflowFaceIndex, const typename LocalSubTriangulation::BoundingBox& smallCellBoundingBox, const Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension>& massMatrix, bool smallCutCellInside, const X& x, Y& y, const LocalView& insideView, const typename LocalOperator::LocalVector& insideLocalValue, CorrectionIt& correctionIt, K capacity, typename LocalOperator::LocalVector& insideLocalResult) {
    int smallCutCellDomainIndex = smallCutCellInside ? cutIntersectionIt->insideDomainIndex() : cutIntersectionIt->outsideDomainIndex();

    auto& outsideView = outsideViews_[cutIntersectionIt->outsideDomainIndex()];
    outsideView.bind(cutIntersectionIt->outside());

    Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> neighborBaseChange(0.0);

    if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
      const auto& neighborBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->outside(), cutIntersectionIt->outsideDomainIndex()).boundingBox;
      const auto& cutCellBasis =  insideView.tree().child(0).finiteElement().localBasis();
      interpolateBasis(smallCellBoundingBox, neighborBoundingBox, cutCellBasis, outsideViews_[cutIntersectionIt->outsideDomainIndex()].tree().child(0).finiteElement().localBasis(), massMatrix, neighborBaseChange);
    } else {
      const auto& neighborBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->inside(), cutIntersectionIt->insideDomainIndex()).boundingBox;
      const auto& cutCellBasis =  outsideView.tree().child(0).finiteElement().localBasis();
      interpolateBasis(smallCellBoundingBox, neighborBoundingBox, cutCellBasis, insideViews_[cutIntersectionIt->insideDomainIndex()].tree().child(0).finiteElement().localBasis(), massMatrix, neighborBaseChange);
    }

    std::size_t inflowFaceIndex = 0;

    for (auto smallCellIntersectionIt = localSubTriangulation_.cutIntersectionBegin();
      smallCellIntersectionIt != localSubTriangulation_.cutIntersectionEnd(); ++smallCellIntersectionIt) {
      if (smallCellIntersectionIt->insideDomainIndex() == smallCutCellDomainIndex || smallCellIntersectionIt->outsideDomainIndex() == smallCutCellDomainIndex) {
        LocalVector outsideLocalValue(0.0);
        readLocalVector(x, outsideView, outsideLocalValue);

        typename LocalOperator::LocalVector extentedLocalValue(0.0);

        int inflowNeighborDomainIndex = -1;

        Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> extendedBaseChange(0.0);

        if (smallCellIntersectionIt->neighbor()) {
          if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
            loadExtendedData(*smallCellIntersectionIt, smallCellBoundingBox, cutIntersectionIt->insideDomainIndex(), insideView, massMatrix, inflowNeighborDomainIndex, x, extentedLocalValue, extendedBaseChange);
          } else {
            loadExtendedData(*smallCellIntersectionIt, smallCellBoundingBox, cutIntersectionIt->outsideDomainIndex(), outsideView, massMatrix, inflowNeighborDomainIndex, x, extentedLocalValue, extendedBaseChange);
          }
        } else {
          if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
            interpolateFunctionLocal(localOperator_.boundaryValue(), insideLocalValue, smallCellBoundingBox, insideView, massMatrix, extentedLocalValue, 4, std::integral_constant<int, LocalOperator::Model::m>());
          } else {
            interpolateFunctionLocal(localOperator_.boundaryValue(), outsideLocalValue, smallCellBoundingBox, outsideView, massMatrix, extentedLocalValue, 4, std::integral_constant<int, LocalOperator::Model::m>());
          }
        }

        typename LocalOperator::LocalVector outsideLocalResult(0.0);
        typename LocalOperator::LocalVector insideLocalResultTmp(0.0);

        if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
          localOperator_.dodFace(*cutIntersectionIt, insideView, outsideView, insideLocalValue, outsideLocalValue, extentedLocalValue, outflowFaceIndex, inflowFaceIndex, neighborBaseChange, capacity, insideLocalResultTmp, outsideLocalResult, false);
          applyL2CorrectionFlux(y, *correctionIt, insideView, outsideView, inflowNeighborDomainIndex, outsideLocalValue, extentedLocalValue, capacity, neighborBaseChange, extendedBaseChange, smallCellIntersectionIt->neighbor());
        } else {
          localOperator_.dodFace(*cutIntersectionIt, outsideView, insideView, outsideLocalValue, insideLocalValue, extentedLocalValue, outflowFaceIndex, inflowFaceIndex, neighborBaseChange, capacity, outsideLocalResult, insideLocalResultTmp, true);
          applyL2CorrectionFlux(y, *correctionIt, outsideView, insideView, inflowNeighborDomainIndex, insideLocalValue, extentedLocalValue, capacity, neighborBaseChange, extendedBaseChange, smallCellIntersectionIt->neighbor());
        }

        insideLocalResult += insideLocalResultTmp;

        addLocalVector(y, outsideView, outsideLocalResult);

        ++inflowFaceIndex;
        ++correctionIt;
      }
    }
  }

  template<class CutIntersectionIt, class X, class Y, class LocalView>
  void assembleBoundaryStabilization(CutIntersectionIt cutIntersectionIt, std::size_t outflowFaceIndex, const typename LocalSubTriangulation::BoundingBox& smallCellBoundingBox, const Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension>& massMatrix, const X& x, Y& y, const LocalView& insideView, const typename LocalOperator::LocalVector& insideLocalValue, K capacity, typename LocalOperator::LocalVector& insideLocalResult) {
    std::size_t inflowFaceIndex = 0;

    for (auto smallCellIntersectionIt = localSubTriangulation_.cutIntersectionBegin();
      smallCellIntersectionIt != localSubTriangulation_.cutIntersectionEnd(); ++smallCellIntersectionIt) {
      if (smallCellIntersectionIt->insideDomainIndex() == cutIntersectionIt->insideDomainIndex() || smallCellIntersectionIt->outsideDomainIndex() == cutIntersectionIt->insideDomainIndex()) {
        if (inflowFaceIndex != outflowFaceIndex) {
          typename LocalOperator::LocalVector extentedLocalValue(0.0);
          typename LocalOperator::LocalVector extendedBoundaryData(0.0);

          int inflowNeighborDomainIndex = -1;

          Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> extendedBaseChange(0.0);

          typename LocalOperator::LocalVector boundaryData(0.0);

          interpolateFunctionLocalReflected(insideLocalValue, smallCellBoundingBox, insideView, massMatrix, [this] (const auto& x) { return localOperator_.reflectPoint(x); }, [this] (const auto& x) { return localOperator_.reflectionVector(); }, boundaryData, 2, std::integral_constant<int, LocalOperator::Model::m>());

          if (smallCellIntersectionIt->neighbor()) {
            loadExtendedData(*smallCellIntersectionIt, smallCellBoundingBox, cutIntersectionIt->insideDomainIndex(), insideView, massMatrix, inflowNeighborDomainIndex, x, extentedLocalValue, extendedBaseChange);

            if (localOperator_.useMatrixStabilization()) {
              interpolateFunctionLocalReflected(extentedLocalValue, smallCellBoundingBox, insideView, massMatrix, [this] (const auto& x) { return localOperator_.reflectPoint(x); }, [this] (const auto& x) { return localOperator_.reflectionVector(); }, extendedBoundaryData, 2, std::integral_constant<int, LocalOperator::Model::m>());
            }
          } else {
            assert(false);
            interpolateFunctionLocal(localOperator_.boundaryValue(), insideLocalValue, smallCellBoundingBox, insideView, massMatrix, extentedLocalValue, 2, std::integral_constant<int, LocalOperator::Model::m>());
          }

          localOperator_.dodBoundary(*cutIntersectionIt, insideView, insideLocalValue, extentedLocalValue, boundaryData, extendedBoundaryData, outflowFaceIndex, inflowFaceIndex, capacity, insideLocalResult);
        }

        ++inflowFaceIndex;
      }
    }
  }

  template<class Y>
  void applyInvertedMassMatrix(Y& y)
  {
    // TODO: Move this somehow into the local operator
    std::vector<typename LocalOperator::QuadraturePoint> quadratureRule;

    // apply the inverted mass matrix
    for (const auto& element : elements(localSubTriangulation_.gridView())) {
      localSubTriangulation_.bindOnVolume(element);
      localSubTriangulation_.createCutCells();

      for (auto cutCellIt = localSubTriangulation_.cutCellsBegin(); cutCellIt != localSubTriangulation_.cutCellsEnd(); ++cutCellIt) {
        int currentDomainIndex = cutCellIt->domainIndex();
        auto& insideView = insideViews_[currentDomainIndex];
        insideView.bind(element);

        LocalVector insideLocalValue(0.0);
        readLocalVector(y, insideView, insideLocalValue);

        Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> massMatrix(0.0);
        cutCellIt->quadratureRule(quadratureRule, 2 * insideView.tree().child(0).finiteElement().localBasis().order());
        computeMassMatrix(quadratureRule, insideView.tree().child(0).finiteElement().localBasis(), massMatrix);
        massMatrix.invert();

        // TODO: I don't like the fact here that we have to access the component number
        auto insideLocalResult = transformBasisComponents(insideLocalValue, massMatrix, insideView, insideView, std::integral_constant<int, LocalOperator::Model::m>());
        insideLocalResult *= -1.0;

        writeLocalVector(y, insideView, insideLocalResult);
      }
    }
  }

  bool visitFace(int insideIndex, int outsideIndex, int insideDomainIndex, int outsideDomainIndex) const
  {
    return insideIndex < outsideIndex || (insideIndex == outsideIndex && insideDomainIndex != outsideDomainIndex);
  }

  template<class CutIntersection, class BoundingBox, class LocalView, class Vector, class LocalVector, class MassMatrix, class BaseChange>
  void loadExtendedData(const CutIntersection& inflowFace, const BoundingBox& cutCellBoundingBox, int smallCutCellDomainIndex, const LocalView& smallCutCellLocalView, const MassMatrix& massMatrix, int& inflowNeighborDomainIndex, const Vector& x, LocalVector& extendedLocalValue, BaseChange& extendedBaseChange)
  {
    const auto& cutCellBasis = smallCutCellLocalView.tree().child(0).finiteElement().localBasis();

    if (inflowFace.insideDomainIndex() == inflowFace.outsideDomainIndex() || inflowFace.insideDomainIndex() == smallCutCellDomainIndex) {
      const auto& inflowElement = inflowFace.outside();
      inflowNeighborDomainIndex = inflowFace.outsideDomainIndex();
      extendedViews_[inflowNeighborDomainIndex].bind(inflowElement);

      const auto& inflowNeighborBoundingBox = localSubTriangulation_.cutCellInformation().information(inflowElement, inflowNeighborDomainIndex).boundingBox;
      interpolateBasis(cutCellBoundingBox, inflowNeighborBoundingBox, cutCellBasis, extendedViews_[inflowNeighborDomainIndex].tree().child(0).finiteElement().localBasis(), massMatrix, extendedBaseChange);
    } else {
      const auto& inflowElement = inflowFace.inside();
      inflowNeighborDomainIndex = inflowFace.insideDomainIndex();
      extendedViews_[inflowNeighborDomainIndex].bind(inflowElement);

      const auto& inflowNeighborBoundingBox = localSubTriangulation_.cutCellInformation().information(inflowElement, inflowNeighborDomainIndex).boundingBox;
      interpolateBasis(cutCellBoundingBox, inflowNeighborBoundingBox, cutCellBasis, extendedViews_[inflowNeighborDomainIndex].tree().child(0).finiteElement().localBasis(), massMatrix, extendedBaseChange);
    }

    readLocalVector(x, extendedViews_[inflowNeighborDomainIndex], extendedLocalValue);
    auto tmp = transformBasisComponents(extendedLocalValue, extendedBaseChange, smallCutCellLocalView, extendedViews_[inflowNeighborDomainIndex], std::integral_constant<int, LocalOperator::Model::m>());
    extendedLocalValue = tmp;
  }

  template<class Vector, class LocalView, class BaseChange>
  void applyL2CorrectionFlux(Vector& y, const LocalMatrix& correctionOperator, const LocalView& insideView, const LocalView& outsideView, int inflowNeighborDomainIndex, const LocalVector& outsideLocalValue, const LocalVector& extendedLocalValue, K capacity, const BaseChange& neighborBaseChange, const BaseChange& extendedBaseChange, bool hasNeighbor) const
  {
    LocalVector tmp = outsideLocalValue;
    tmp = transformBasisComponents(tmp, neighborBaseChange, insideView, outsideView, std::integral_constant<int, LocalOperator::Model::m>());

    LocalVector correction(0.0);
    correctionOperator.mv(extendedLocalValue - tmp, correction);

    LocalVector neighborCorrection(0.0);

    projectOntoNeighborBasis(insideView, outsideView, neighborBaseChange, correction, neighborCorrection);

    addLocalVector(y, outsideView, l2StabilityFactor_ * (1.0 - capacity) * neighborCorrection);

    if (hasNeighbor) {
      LocalVector inflowNeighborCorrection(0.0);

      projectOntoNeighborBasis(insideView, extendedViews_[inflowNeighborDomainIndex], extendedBaseChange, correction, inflowNeighborCorrection);

      addLocalVector(y, extendedViews_[inflowNeighborDomainIndex], -l2StabilityFactor_ * (1.0 - capacity) * inflowNeighborCorrection);
    }
  }

  template<class LocalView>
  void projectOntoNeighborBasis(const LocalView& smallCellView, const LocalView& neighborView, const Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension>& baseChange, const LocalVector& v, LocalVector& res) const
  {
    for (int i = 0; i < LocalOperator::Model::m; ++i) {
      for (int j = 0; j < neighborView.tree().child(i).size(); ++j) {
        LocalVector basisVector(0.0);
        basisVector[neighborView.tree().child(i).localIndex(j)] = 1.0;
        basisVector = transformBasisComponents(basisVector, baseChange, smallCellView, neighborView, std::integral_constant<int, LocalOperator::Model::m>());
        res[neighborView.tree().child(i).localIndex(j)] = v * basisVector;
      }
    }
  }


  Domain zero() const
  {
    VectorType vector;
    auto vectorBackend = Dune::Functions::istlVectorBackend(vector);
    vectorBackend.resize(globalBasis_);
    vector = 0.0;
    return vector;
  }

private:
  LocalSubTriangulation& localSubTriangulation_;
  const GlobalBasis& globalBasis_;
  std::vector<typename GlobalBasis::LocalView> insideViews_;
  std::vector<typename GlobalBasis::LocalView> outsideViews_;
  std::vector<typename GlobalBasis::LocalView> extendedViews_;
  LocalOperator localOperator_;
  SmallCutCellClassifier& classifier_;
  K l2StabilityFactor_;
  bool useBoundaryHack_;
  bool ignoreSmallCells_;
};

} // end namespace Dune::Hypercut

#endif
#ifndef DUNE_HYPERCUT_BETASEMINORM_HH
#define DUNE_HYPERCUT_BETASEMINORM_HH

#include <cassert>
#include <type_traits>
#include <vector>

#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/hypercut/helper.hh>
#include <dune/hypercut/interpolation.hh>

#include <dune/subtriangulation/simpletpmctriangulation/localsubtriangulation.hh>

namespace Dune::Hypercut
{
    template<class Model_, int ComponentSpaceDimension_, class K_>
    class LocalSeminorm
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

    LocalSeminorm(const Model& model, int intorderadd = 1) :
        model_(model),
        intorderadd_(intorderadd)
    {

    }

    template<class CutIntersection, class LocalView>
    K dgFace(const CutIntersection& cutIntersection,
             const LocalView& insideView, const LocalView& outsideView,
             const LocalVector& x_s, const LocalVector& x_n) const
    {
        const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();
        const auto& outsideFiniteElement = outsideView.tree().child(0).finiteElement();

        int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

        using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;

        std::vector<RangeType> phi_s(insideFiniteElement.localBasis().size());
        std::vector<RangeType> phi_n(outsideFiniteElement.localBasis().size());

        const auto& geometry = cutIntersection.geometry();
        const auto& geometryInInside = cutIntersection.geometryInInside();
        const auto& geometryInOutside = cutIntersection.geometryInOutside();

        cutIntersection.quadratureRule(quadratureRule_, intorder);

        K res = 0.0;

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
            const auto flux = model_.seminorm(normal, u_s - u_n, qpg);

            for (std::size_t i = 0; i < Model::m; ++i) {
                for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
                    res += 0.5 * flux * (u_s - u_n) * factor;
                }
            }
        }

        return res;
    }

    template<class CutIntersection, class LocalView, class Solution>
    K dgBoundary(const CutIntersection& cutIntersection, const LocalView& localView,
                 const LocalVector& x, const Solution& solution) const
    {
        const auto& finiteElement = localView.tree().child(0).finiteElement();

        int intorder = 2 * finiteElement.localBasis().order() + intorderadd_;
        using RangeType = typename std::decay_t<decltype(finiteElement.localBasis())>::Traits::RangeType;

        std::vector<RangeType> phi(finiteElement.localBasis().size());

        const auto& geometry = cutIntersection.geometry();
        const auto& geometryInInside = cutIntersection.geometryInInside();

        cutIntersection.quadratureRule(quadratureRule_, intorder);

        K res = 0.0;

        for (const auto& qp : quadratureRule_) {
            const auto qppos_s = geometryInInside.local(qp.position());
            const auto qpg = qp.position();

            finiteElement.localBasis().evaluateFunction(qppos_s, phi);

            Dune::FieldVector<K, Model::m> u(0.0);

            fillInFunctionValue(localView, x, phi, u, std::integral_constant<int, Model::m>());

            auto normal = cutIntersection.unitOuterNormal(qp.position());
            const K factor = qp.weight();
            const auto flux = model_.outflow(normal, solution(qpg) - u, qpg);

            for (std::size_t i = 0; i < Model::m; ++i) {
                for (std::size_t j = 0; j < localView.tree().child(i).size(); ++j) {
                    res += 0.5 * flux * (solution(qpg) - u) * factor;
                }
            }
        }

        return res;
    }

    template<class CutIntersection, class LocalView, class Weight, class Solution>
    K dodFace(const CutIntersection& cutIntersection,
              const LocalView& insideView, const LocalView& outsideView,
              const LocalVector& x_s,
              const LocalVector& x_n,
              const LocalVector& x_extended,
              const Weight& weight, K capacity,
              bool smallCutCellIsOutside,
              const Solution& solution,
              K weightedBetaMean) const
    {
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

        K res = 0.0;

        for (const auto& qp : quadratureRule_) {
            const auto qppos_s = geometryInInside.local(qp.position());
            const auto qppos_n = geometryInOutside.local(qp.position());
            const auto qpg = qp.position();

            insideFiniteElement.localBasis().evaluateFunction(qppos_s, phi_s);
            outsideFiniteElement.localBasis().evaluateFunction(qppos_n, phi_n);

            Dune::FieldVector<K, Model::m> u_s(0.0);
            Dune::FieldVector<K, Model::m> u_extended(0.0);
            Dune::FieldVector<K, Model::m> u_n(0.0);

            fillInFunctionValue(insideView, x_s, phi_s, u_s, std::integral_constant<int, Model::m>());
            fillInFunctionValue(insideView, x_extended, phi_s, u_extended, std::integral_constant<int, Model::m>());
            fillInFunctionValue(outsideView, x_n, phi_n, u_n, std::integral_constant<int, Model::m>());

            const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
            const auto factor = qp.weight();

            auto flux = model_.seminorm(normal, u_s - u_n, qpg);
            flux *= (1.0 - capacity);

            auto fluxExt = model_.outflow(normal, (weightedBetaMean - u_extended) - (solution(qpg) - u_n), qpg);
            fluxExt *= (1.0 - capacity);

            for (std::size_t i = 0; i < Model::m; ++i) {
                for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
                    res += 0.5 * fluxExt * ((weightedBetaMean - u_extended) - (solution(qpg) - u_n)) * factor;
                    res += -0.5 * flux * (u_s - u_n) * factor;
                }
            }
        }

        return res;
    }

    template<class CutIntersection, class LocalView, class Weight>
    K dodBoundary(const CutIntersection& cutIntersection, const LocalView& insideView,
                    const LocalVector& x, const LocalVector& x_extended,
                    const LocalVector& x_transformed, const LocalVector& x_extended_transformed,
                    const LocalVector& bc_transformed_coefficients,
                    const Weight& weight, K capacity)
    {
        const auto& insideFiniteElement = insideView.tree().child(0).finiteElement();

        int intorder = 2 * insideFiniteElement.localBasis().order() + intorderadd_;

        using RangeType = typename std::decay_t<decltype(insideFiniteElement.localBasis())>::Traits::RangeType;

        std::vector<RangeType> phi(insideFiniteElement.localBasis().size());
        std::vector<Dune::FieldVector<K, Model::m>> phi_transformed(insideView.tree().size());

        const auto& geometry = cutIntersection.geometry();
        const auto& geometryInInside = cutIntersection.geometryInInside();

        cutIntersection.quadratureRule(quadratureRule_, intorder);

        K res = 0.0;

        for (const auto& qp : quadratureRule_) {
            const auto qppos_s = geometryInInside.local(qp.position());
            const auto qpg = qp.position();

            insideFiniteElement.localBasis().evaluateFunction(qppos_s, phi);

            Dune::FieldVector<K, Model::m> u(0.0);
            Dune::FieldVector<K, Model::m> u_extended(0.0);

            fillInFunctionValue(insideView, x, phi, u, std::integral_constant<int, Model::m>());
            fillInFunctionValue(insideView, x_extended, phi, u_extended, std::integral_constant<int, Model::m>());

            const auto normal = cutIntersection.unitOuterNormal(qp.position());
            const auto factor = qp.weight();

            auto flux = model_.outflow(u_extended - u);
            flux *= (1.0 - capacity);

            for (std::size_t i = 0; i < Model::m; ++i) {
                for (std::size_t j = 0; j < insideView.tree().child(i).size(); ++j) {
                    res += flux * (u_extended - u) * factor;
                }
            }
        }

        return res;
    }

    // We assume a P0 scheme for a scalar equation here
    template<class CutIntersection, class Solution>
    K weightedBetaMean(const CutIntersection& cutIntersection, const Solution& u, bool smallCutCellIsOutside) const
    {
        const K sign = smallCutCellIsOutside ? -1.0 : 1.0;

        cutIntersection.quadratureRule(quadratureRule_, 4);

        K mean = 0.0;
        K volume = 0.0;

        for (const auto& qp : quadratureRule_) {
            const auto qpg = qp.position();
            const auto normal = sign * cutIntersection.unitOuterNormal(qp.position());
            const auto factor = qp.weight();

            auto inflowspeed = -model_.inflowMatrix(normal, qpg);


            if (inflowspeed > Epsilon) {
                mean += inflowspeed * u(qpg) * factor;
                volume += inflowspeed * factor;
            }
        }

        if (volume > Epsilon) {
            return mean / volume;
        }

        return 0.0;
    }


    private:
        Model model_;

        mutable std::vector<QuadraturePoint> quadratureRule_;

        int intorderadd_;
    };

    template<class LST, class GlobalBasis, class VectorType, class LocalOp, class SmallCutCellClassifier>
    class DodSeminorm
    {
    public:
    using LocalSubTriangulation = LST;
    using LocalOperator = LocalOp;
    using K = typename LST::ctype;

    using LocalSeminorm = LocalSeminorm<typename LocalOperator::Model,
                                        LocalOperator::ComponentSpaceDimension, K>;

    using LocalVector = typename LocalOperator::LocalVector;
    using LocalMatrix = typename LocalOperator::LocalMatrix;

    using Domain = VectorType;

    static const int ComponentSpaceDimension = LocalOperator::ComponentSpaceDimension;

    DodSeminorm(LocalSubTriangulation& localSubTriangulation,
                const GlobalBasis& globalBasis,
                const LocalOperator& localOperator,
                const LocalSeminorm& localSeminorm,
                const SmallCutCellClassifier& classifier) :
        localSubTriangulation_(localSubTriangulation),
        globalBasis_(globalBasis),
        localOperator_(localOperator),
        localSeminorm_(localSeminorm),
        classifier_(classifier)
    {
        for (int domainIndex = 0; domainIndex < localSubTriangulation_.domainConfiguration().numberOfDomains(); ++domainIndex) {
            insideViews_.push_back(globalBasis_.localView(domainIndex));
            outsideViews_.push_back(globalBasis_.localView(domainIndex));
            extendedViews_.push_back(globalBasis_.localView(domainIndex));
        }
    }

    template<class Vector, class Solution>
    K operator()(const Vector& x, const Solution& solution)
    {
        return assemble(Dune::Functions::istlVectorBackend(x), solution);
    }

    template<class X, class Solution>
    K assemble(const X& x, const Solution& solution)
    {
        classifier_.classify(localSubTriangulation_);

        const auto& indexSet = localSubTriangulation_.gridView().indexSet();

        K res = 0.0;

        // assemble the space discretization form
        for (const auto& element : elements(localSubTriangulation_.gridView())) {
            K tmp = res;

            localSubTriangulation_.bind(element);
            localSubTriangulation_.createCutCells();
            localSubTriangulation_.createCutIntersections();

            std::vector<typename LocalOperator::LocalMatrix> weights;
            std::vector<typename LocalOperator::LocalMatrix> corrections;
            std::vector<typename LocalOperator::LocalMatrix> residualOperators;

            // TODO: Find a proper place for this together with a proper handling
            K capacity = 1.0;

            for (auto cutCellIt = localSubTriangulation_.cutCellsBegin();
                cutCellIt != localSubTriangulation_.cutCellsEnd(); ++cutCellIt) {
                int currentDomainIndex = cutCellIt->domainIndex();

                bool isSmallCell = classifier_.isSmallCell(indexSet.index(element), currentDomainIndex);

                if (isSmallCell) {
                    // TODO: Right now we are only dealing with isolated cut cells. Hence at this point
                    // this should be empty since we only have one small cut cell per element
                    // In the future we need a better way to iterate over the snippets of a cut cell
                    assert(weights.empty());
                    // TODO: Same applies to the capacity
                    assert(capacity == 1.0);

                    capacity = localOperator_.capacity(localSubTriangulation_.cutIntersectionBegin(), localSubTriangulation_.cutIntersectionEnd(), localSubTriangulation_.cutCellInformation().information(element, currentDomainIndex), currentDomainIndex);
                }
            }

            assert(3 * weights.size() == corrections.size());

            auto correctionIt = corrections.begin();

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

                        // if (!smallCutCellInside && !smallCutCellOutside && !neighborSmall)
                        res += localSeminorm_.dgFace(*cutIntersectionIt, insideView, outsideView, insideLocalValue, outsideLocalValue);
                    }

                    if (smallCutCellInside || smallCutCellOutside) {
                        res += assembleFaceStabilization(cutIntersectionIt, smallCutCellInside, x, insideView, insideLocalValue, weights, correctionIt, capacity, solution);
                    }
                } else {
                    bool isSmallCell = classifier_.isSmallCell(indexSet.index(element), cutIntersectionIt->insideDomainIndex());

                    res += localSeminorm_.dgBoundary(*cutIntersectionIt, insideView, insideLocalValue, solution);
                }
            }
        }

        return res;
    }

    template<class CutIntersectionIt, class X, class LocalView, class CorrectionIt, class Solution>
    K assembleFaceStabilization(CutIntersectionIt cutIntersectionIt, bool smallCutCellInside, const X& x, const LocalView& insideView, const typename LocalOperator::LocalVector& insideLocalValue, const std::vector<typename LocalOperator::LocalMatrix>& weights, CorrectionIt& correctionIt, K capacity, const Solution& solution) {
        int smallCutCellDomainIndex = smallCutCellInside ? cutIntersectionIt->insideDomainIndex() : cutIntersectionIt->outsideDomainIndex();

        auto& outsideView = outsideViews_[cutIntersectionIt->outsideDomainIndex()];
        outsideView.bind(cutIntersectionIt->outside());

        Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> massMatrix(0.0);
        Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> neighborBaseChange(0.0);

        if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
            const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->inside(), smallCutCellDomainIndex).boundingBox;
            const auto& neighborBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->outside(), cutIntersectionIt->outsideDomainIndex()).boundingBox;
            const auto& cutCellBasis =  insideView.tree().child(0).finiteElement().localBasis();
            computeLocalMassMatrix(cutCellBoundingBox, cutCellBasis, massMatrix);
            massMatrix.invert();
            interpolateBasis(cutCellBoundingBox, neighborBoundingBox, cutCellBasis, outsideViews_[cutIntersectionIt->outsideDomainIndex()].tree().child(0).finiteElement().localBasis(), massMatrix, neighborBaseChange);
        } else {
            const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->outside(), smallCutCellDomainIndex).boundingBox;
            const auto& neighborBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->inside(), cutIntersectionIt->insideDomainIndex()).boundingBox;
            const auto& cutCellBasis =  outsideView.tree().child(0).finiteElement().localBasis();
            computeLocalMassMatrix(cutCellBoundingBox, cutCellBasis, massMatrix);
            massMatrix.invert();
            interpolateBasis(cutCellBoundingBox, neighborBoundingBox, cutCellBasis, insideViews_[cutIntersectionIt->insideDomainIndex()].tree().child(0).finiteElement().localBasis(), massMatrix, neighborBaseChange);
        }

        auto weightIt = weights.begin();

        K res = 0.0;

        for (auto smallCellIntersectionIt = localSubTriangulation_.cutIntersectionBegin();
            smallCellIntersectionIt != localSubTriangulation_.cutIntersectionEnd(); ++smallCellIntersectionIt) {
            if (smallCellIntersectionIt->insideDomainIndex() == smallCutCellDomainIndex || smallCellIntersectionIt->outsideDomainIndex() == smallCutCellDomainIndex) {
                LocalVector outsideLocalValue(0.0);
                readLocalVector(x, outsideView, outsideLocalValue);

                typename LocalOperator::LocalVector extentedLocalValue(0.0);
                typename LocalOperator::LocalVector extendedLocalResult(0.0);

                int inflowNeighborDomainIndex = -1;

                Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> extendedBaseChange(0.0);

                K weightedBetaMean = 0.0;

                if (smallCellIntersectionIt->neighbor()) {
                    if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
                        const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->inside(), smallCutCellDomainIndex).boundingBox;
                        loadExtendedData(*smallCellIntersectionIt, cutCellBoundingBox, cutIntersectionIt->insideDomainIndex(), insideView, massMatrix, inflowNeighborDomainIndex, x, extentedLocalValue, extendedBaseChange);

                        weightedBetaMean = localSeminorm_.weightedBetaMean(*smallCellIntersectionIt, solution, false);
                    } else {
                        const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->outside(), smallCutCellDomainIndex).boundingBox;
                        loadExtendedData(*smallCellIntersectionIt, cutCellBoundingBox, cutIntersectionIt->outsideDomainIndex(), outsideView, massMatrix, inflowNeighborDomainIndex, x, extentedLocalValue, extendedBaseChange);

                        weightedBetaMean = localSeminorm_.weightedBetaMean(*smallCellIntersectionIt, solution, true);
                    }
                } else {
                    if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
                        const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->inside(), cutIntersectionIt->insideDomainIndex()).boundingBox;
                        interpolateFunctionLocal(localOperator_.boundaryValue(), insideLocalValue, cutCellBoundingBox, insideView, massMatrix, extentedLocalValue, 4, std::integral_constant<int, LocalOperator::Model::m>());
                    } else {
                        const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->outside(), cutIntersectionIt->outsideDomainIndex()).boundingBox;
                        interpolateFunctionLocal(localOperator_.boundaryValue(), outsideLocalValue, cutCellBoundingBox, outsideView, massMatrix, extentedLocalValue, 4, std::integral_constant<int, LocalOperator::Model::m>());
                    }
                }

                typename LocalOperator::LocalVector insideLocalValueTransformed(0.0);
                typename LocalOperator::LocalVector outsideLocalValueTransformed(0.0);
                typename LocalOperator::LocalVector extentedLocalValueTransformed(0.0);

                weightIt->mv(extentedLocalValue, extentedLocalValueTransformed);

                typename LocalOperator::LocalVector outsideLocalResult(0.0);
                typename LocalOperator::LocalVector insideLocalResultTmp(0.0);

                if (smallCutCellDomainIndex == cutIntersectionIt->insideDomainIndex()) {
                    LocalVector tmp = outsideLocalValue;
                    tmp = transformBasisComponents(tmp, neighborBaseChange, insideView, outsideView, std::integral_constant<int, LocalOperator::Model::m>());
                    weightIt->mv(tmp, outsideLocalValueTransformed);
                    weightIt->mv(insideLocalValue, insideLocalValueTransformed);

                    if (weightedBetaMean > 0.0) {
                        res += localSeminorm_.dodFace(*cutIntersectionIt, insideView, outsideView, insideLocalValue, outsideLocalValue, extentedLocalValue, *weightIt, capacity, false, solution, weightedBetaMean);
                    }
                } else {
                    LocalVector tmp = insideLocalValue;
                    tmp = transformBasisComponents(tmp, neighborBaseChange, outsideView, insideView, std::integral_constant<int, LocalOperator::Model::m>());
                    weightIt->mv(tmp, insideLocalValueTransformed);
                    weightIt->mv(outsideLocalValue, outsideLocalValueTransformed);

                    if (weightedBetaMean > 0.0) {
                        res += localSeminorm_.dodFace(*cutIntersectionIt, outsideView, insideView, outsideLocalValue, insideLocalValue, extentedLocalValue, *weightIt, capacity, true, solution, weightedBetaMean);
                    }
                }

                ++weightIt;
                ++correctionIt;
            }
        }

        return res;
    }

    template<class CutIntersectionIt, class X, class Y, class LocalView>
    void assembleBoundaryStabilization(CutIntersectionIt cutIntersectionIt, const X& x, Y& y, const LocalView& insideView, const typename LocalOperator::LocalVector& insideLocalValue, const std::vector<typename LocalOperator::LocalMatrix>& weights, const std::vector<typename LocalOperator::LocalMatrix>& residualOperators, K capacity, typename LocalOperator::LocalVector& insideLocalResult) {
        auto weightIt = weights.begin();
        auto residualOperatorIt = residualOperators.begin();

        for (auto smallCellIntersectionIt = localSubTriangulation_.cutIntersectionBegin();
        smallCellIntersectionIt != localSubTriangulation_.cutIntersectionEnd(); ++smallCellIntersectionIt) {
        if (smallCellIntersectionIt->insideDomainIndex() == cutIntersectionIt->insideDomainIndex() || smallCellIntersectionIt->outsideDomainIndex() == cutIntersectionIt->insideDomainIndex()) {
            typename LocalOperator::LocalVector extentedLocalValue(0.0);
            typename LocalOperator::LocalVector extendedLocalResult(0.0);

            int inflowNeighborDomainIndex = -1;

            Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> massMatrix(0.0);
            Dune::FieldMatrix<K, LocalOperator::ComponentSpaceDimension, LocalOperator::ComponentSpaceDimension> extendedBaseChange(0.0);

            typename LocalOperator::LocalVector boundaryData(0.0);

            const auto& cutCellBoundingBox = localSubTriangulation_.cutCellInformation().information(cutIntersectionIt->inside(), cutIntersectionIt->insideDomainIndex()).boundingBox;
            computeLocalMassMatrix(cutCellBoundingBox, insideView.tree().child(0).finiteElement().localBasis(), massMatrix);
            massMatrix.invert();

            interpolateFunctionLocal(localOperator_.boundaryValue(), insideLocalValue, cutCellBoundingBox, insideView, massMatrix, boundaryData, 2, std::integral_constant<int, LocalOperator::Model::m>());

            if (smallCellIntersectionIt->neighbor()) {
                loadExtendedData(*smallCellIntersectionIt, cutCellBoundingBox, cutIntersectionIt->insideDomainIndex(), insideView, massMatrix, inflowNeighborDomainIndex, x, extentedLocalValue, extendedBaseChange);
            } else {
                interpolateFunctionLocal(localOperator_.boundaryValue(), insideLocalValue, cutCellBoundingBox, insideView, massMatrix, extentedLocalValue, 2, std::integral_constant<int, LocalOperator::Model::m>());
            }

            typename LocalOperator::LocalVector insideLocalValueTransformed(0.0);
            typename LocalOperator::LocalVector extentedLocalValueTransformed(0.0);
            typename LocalOperator::LocalVector boundaryDataTransformed(0.0);

            weightIt->mv(insideLocalValue, insideLocalValueTransformed);
            weightIt->mv(extentedLocalValue, extentedLocalValueTransformed);
            weightIt->mv(boundaryData, boundaryDataTransformed);

            // localOperator_.dodBoundary(*cutIntersectionIt, insideView, insideLocalValue, extentedLocalValue, insideLocalValueTransformed, extentedLocalValueTransformed, boundaryDataTransformed, *weightIt, capacity, insideLocalResult);

            ++weightIt;
            ++residualOperatorIt;
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

    private:
        LocalSubTriangulation& localSubTriangulation_;
        const GlobalBasis& globalBasis_;
        LocalOperator localOperator_;
        LocalSeminorm localSeminorm_;
        std::vector<typename GlobalBasis::LocalView> insideViews_;
        std::vector<typename GlobalBasis::LocalView> outsideViews_;
        std::vector<typename GlobalBasis::LocalView> extendedViews_;
        SmallCutCellClassifier classifier_;
    };
}

#endif
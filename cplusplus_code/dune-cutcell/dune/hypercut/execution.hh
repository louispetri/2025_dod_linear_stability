#ifndef DUNE_HYPERCUT_EXECUTION_HH
#define DUNE_HYPERCUT_EXECUTION_HH

#include <type_traits>

#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/functionspacebases/monomialbasis.hh>
#include <dune/functions/functionspacebases/subdomainbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/multidomainbasis.hh>

#include <dune/subtriangulation/simpletpmctriangulation.hh>
#include <dune/subtriangulation/simpletpmctriangulation/exampledomains.hh>

#include <dune/hypercut/classification.hh>
#include <dune/hypercut/helper.hh>
#include <dune/hypercut/interpolation.hh>
#include <dune/hypercut/localoperator.hh>
#include <dune/hypercut/runtimeparameter.hh>
#include <dune/hypercut/timesteppingmethods.hh>
#include <dune/hypercut/vtkoutput.hh>

namespace Dune::Hypercut
{
    template<class K, int m>
    struct Result
    {
        Dune::FieldVector<K, m> infError_ = Dune::FieldVector<K, m>(0.0);
        Dune::FieldVector<K, m> l1Error_ = Dune::FieldVector<K, m>(0.0);
        Dune::FieldVector<K, m> l2Error_ = Dune::FieldVector<K, m>(0.0);
        K operatorNorm_;
        K minCellFraction_ = -1.0;
        K maxCellFraction_ = -1.0;
    };

    template<class GlobalOperator, class = void>
    struct HasLocalOperator
    {
        using Result = std::false_type;
    };

    template<class GlobalOperator>
    struct HasLocalOperator<GlobalOperator, std::void_t<typename GlobalOperator::LocalOperator>>
    {
        using Result = std::true_type;
    };

    template<class GlobalOperator, class T>
    struct ComponentNumber
    {
        static const int m = -1;
    };

    template<class GlobalOperator>
    struct ComponentNumber<GlobalOperator, std::false_type>
    {
        static const int m = GlobalOperator::Model::m;
    };

    template<class GlobalOperator>
    struct ComponentNumber<GlobalOperator, std::true_type>
    {
        static const int m = GlobalOperator::LocalOperator::Model::m;
    };

    template<class GlobalOperator, class GlobalBasis, class SubTriangulation, class Vector, class Configuration, class InitialData, class Solution, class K, int localSpaceDimension>
    K simulation(GlobalOperator& globalOperator, const GlobalBasis& basis, SubTriangulation& subTriangulation, Vector& vector, const Configuration& configuration, const InitialData& initialData, const Solution& solution, K T, K deltaT, std::integral_constant<int, localSpaceDimension>)
    {
        const int m = ComponentNumber<GlobalOperator, typename HasLocalOperator<GlobalOperator>::Result>::m;

        auto ts = timesteppingMethod(globalOperator, configuration.timesteppingMethod_);
        auto vectorBackend = Dune::Functions::istlVectorBackend(vector);

        Dune::Hypercut::interpolate(subTriangulation, basis, initialData,
                                    vectorBackend, std::integral_constant<int, m>(), std::integral_constant<int, localSpaceDimension>());

        if (configuration.writeVtk_) {
            Dune::Hypercut::writeMultiDomainFunction(subTriangulation, Dune::Functions::istlVectorBackend(vector), basis, "result_0",
                                                    std::integral_constant<int, m>(), std::integral_constant<int, localSpaceDimension>(),
                                                    solution, 0.0);
        }

        K t = 0.0l;
        int step = 1;

        Vector result;
        auto resultBackend = Dune::Functions::istlVectorBackend(result);
        resultBackend.resize(basis);

        while (t < T - deltaT) {
            result = 0.0;

            ts(vector, deltaT);
            t += deltaT;

            if (configuration.writeVtk_) {
                Dune::Hypercut::writeMultiDomainFunction(subTriangulation, Dune::Functions::istlVectorBackend(vector), basis, std::string("result_") + std::to_string(step),
                                                        std::integral_constant<int, m>(), std::integral_constant<int, localSpaceDimension>(),
                                                        solution, t);
            }

            ++step;
        }

        deltaT = T - t;
        result = 0.0;
        ts(vector, deltaT);
        t += deltaT;

        if (configuration.writeVtk_) {
            Dune::Hypercut::writeMultiDomainFunction(subTriangulation, Dune::Functions::istlVectorBackend(vector), basis, std::string("result_") + std::to_string(step),
                                                    std::integral_constant<int, m>(), std::integral_constant<int, localSpaceDimension>(),
                                                    solution, t);
        }

        return t;
    }

    template<class Basis, class VectorBackend, class LocalSubTriangulation, class LocalOperator, class Model, class Solution>
    void computeErrors(const Basis& basis, const VectorBackend& vectorBackend, LocalSubTriangulation& localSubTriangulation, const LocalOperator& localOperator, const Model& model, const Solution& solution, typename Model::K t, const Configuration<typename Model::K>& configuration, Result<typename Model::K, Model::m>& result)
    {
        using K = typename Model::K;
        using ErrorFunction = Dune::Hypercut::ComputeInfError<std::decay_t<decltype(basis)>, std::decay_t<decltype(vectorBackend)>, decltype(solution), K, Model, typename std::decay_t<decltype(localOperator)>::LocalVector>;
        using L1ErrorFunction = Dune::Hypercut::ComputeL1Error<std::decay_t<decltype(basis)>, std::decay_t<decltype(vectorBackend)>, decltype(solution), K, Model, typename std::decay_t<decltype(localOperator)>::LocalVector>;
        using L2ErrorFunction = Dune::Hypercut::ComputeSquaredL2Error<std::decay_t<decltype(basis)>, std::decay_t<decltype(vectorBackend)>, decltype(solution), K, Model, typename std::decay_t<decltype(localOperator)>::LocalVector>;

        if (configuration.computeLInfNorm_) {
            ErrorFunction errorFunction(basis, vectorBackend, solution, t, 6);
            result.infError_ = cutcellVolumesAccumulate(localSubTriangulation, errorFunction, Dune::FieldVector<K, Model::m>(0.0));

            for (std::size_t i = 0; i < Model::m; ++i) {
                result.infError_[i] = localSubTriangulation.gridView().grid().comm().max(result.infError_[i]);
            }
        }

        if (configuration.computeL1Norm_) {
            L1ErrorFunction l1ErrorFunction(basis, vectorBackend, solution, t, 4);
            result.l1Error_ = cutcellVolumesAccumulate(localSubTriangulation, l1ErrorFunction, Dune::FieldVector<K, Model::m>(0.0));

            for (std::size_t i = 0; i < Model::m; ++i) {
                result.l1Error_[i] = localSubTriangulation.gridView().grid().comm().sum(result.l1Error_[i]);
            }
        }

        if (configuration.computeL2Norm_) {
            L2ErrorFunction l2ErrorFunction(basis, vectorBackend, solution, t, 2);
            result.l2Error_ = cutcellVolumesAccumulate(localSubTriangulation, l2ErrorFunction, Dune::FieldVector<K, Model::m>(0.0));

            for (std::size_t i = 0; i < Model::m; ++i) {
                result.l2Error_[i] = localSubTriangulation.gridView().grid().comm().sum(result.l2Error_[i]);
            }

            for (int i = 0; i < Model::m; ++i) {
                result.l2Error_[i] = std::sqrt(result.l2Error_[i]);
            }
        }
    }

    template<class GlobalOperator>
    typename GlobalOperator::K computeOperatorNorm(GlobalOperator& globalOperator)
    {
        using K = typename GlobalOperator::K;

        auto vector = globalOperator.zero();
        auto result = globalOperator.zero();
        auto vectorBackend = Dune::Functions::istlVectorBackend(vector);
        auto resultBackend = Dune::Functions::istlVectorBackend(result);

        Dune::DynamicMatrix<K> spaceDiscretizationMatrix(vector.size(), vector.size(), 0.0);
        Dune::DynamicMatrix<K> spaceDiscretizationTransposedMatrix(vector.size(), vector.size(), 0.0);
        Dune::DynamicMatrix<K> invertedMassMatrix(vector.size(), vector.size(), 0.0);
        Dune::DynamicMatrix<K> massMatrix(vector.size(), vector.size(), 0.0);
        Dune::DynamicMatrix<K> normalMatrix(vector.size(), vector.size(), 0.0);

        for (std::size_t i = 0; i < vector.size(); ++i) {
            vector[i] = 1.0;
            result = 0.0;

            globalOperator.assembleSpaceDiscretization(0.0, vectorBackend, resultBackend);

            for (std::size_t j = 0; j < vector.size(); ++j) {
                spaceDiscretizationMatrix[j][i] = result[j];
            }

            vector[i] = 0.0;
        }

        for (std::size_t i = 0; i < vector.size(); ++i) {
            vector[i] = 1.0;

            globalOperator.applyInvertedMassMatrix(vectorBackend);

            for (std::size_t j = 0; j < vector.size(); ++j) {
                invertedMassMatrix[j][i] = vector[j];
            }

            vector = 0.0;
        }

        spaceDiscretizationTransposedMatrix = spaceDiscretizationMatrix.transposed();

        normalMatrix = invertedMassMatrix;
        normalMatrix.rightmultiply(spaceDiscretizationTransposedMatrix);
        normalMatrix.rightmultiply(invertedMassMatrix);
        normalMatrix.rightmultiply(spaceDiscretizationMatrix);

        Dune::DynamicVector<std::complex<K>> eigenvals (vector.size());
        Dune::DynamicMatrixHelp::eigenValuesNonSym(normalMatrix, eigenvals);

        return std::sqrt(eigenvals.infinity_norm());
    }

    template<class Geometry, class Model, class OperatorFactory, class BoundaryCondition, class InitialData, int order>
    Result<typename Model::K, Model::m> testCase(const Geometry& geometry, const Model& model, const OperatorFactory& operatorFactory, const BoundaryCondition& boundaryCondition, const InitialData& initialData, const Configuration<typename Model::K>& configuration, std::integral_constant<int, order>)
    {
        using K = typename Model::K;
        using GV = typename Geometry::GridView;
        using SubTriangulation = Dune::SubTriangulation::SimpleTpmcTriangulation<GV, GV>;
        using CutCellIndexMap = Dune::Functions::CutCellIndexMap<SubTriangulation>;
        using Classifier = Dune::Hypercut::TriangleVolumeFractionCutCellClassifier<Dune::SubTriangulation::LocalSubTriangulation<GV, GV>>;
        using ComponentSpace = Dune::Functions::MonomialPreBasis<GV, order, K>;

        const auto& gridView = geometry.levelGridView();
        SubTriangulation subTriangulation(gridView, geometry.levelGridView(), geometry.domainConfiguration(), false, configuration.levelSetValueTolerance_);
        auto& localSubTriangulation = subTriangulation.localSubTriangulation();
        CutCellIndexMap cutCellIndexMap(subTriangulation);
        auto cutCellIndexMapPointer = &cutCellIndexMap;
        auto basis = Dune::Functions::BasisFactory::makeMultiDomainBasis(gridView,
            [=] (const auto& gridView, int domainIndex) { return Dune::Functions::BasisFactory::power<Model::m>(
            Dune::Functions::BasisFactory::subdomain<ComponentSpace, CutCellIndexMap>(cutCellIndexMapPointer, domainIndex), Dune::Functions::BasisFactory::FlatLexicographic())(gridView); },
            geometry.domainConfiguration().numberOfDomains(), Dune::Functions::BasisFactory::FlatLexicographic());

        const int localSpaceDimension = ComponentSpace::maxNodeSize();
        using LocalOperator = Dune::Hypercut::LocalOperator<Model, BoundaryCondition, localSpaceDimension, K>;

        Classifier classifier(configuration.smallCellThreshold_);

        if (configuration.writeVtk_) {
            Dune::Hypercut::writeGrid(subTriangulation, localSubTriangulation, classifier, "setup");
        }

        const K cflFactor = configuration.cflSafetyFactor_;
        const K cfl = cflFactor * (1.0 / (2.0 * order + 1.0));
        const K capacityFactor = (1.0 / (2.0 * order + 1.0));
        const K tau = configuration.tau_;

        LocalOperator localOperator(model, boundaryCondition, cfl, capacityFactor, tau, 1.0);
        K deltaT = cfl * geometry.backgroundEdgeLength();
        localOperator.setTimestepSize(deltaT);
        auto globalOperator = operatorFactory(localSubTriangulation, basis, localOperator, classifier, cfl, capacityFactor);

        auto vector = globalOperator.zero();
        auto vectorBackend = Dune::Functions::istlVectorBackend(vector);

        K t = simulation(globalOperator, basis, subTriangulation, vector, configuration, [&] (const auto& x) { return initialData(x, 0.0); }, [&](const auto& x, K t) { return  initialData(x, t);}, configuration.T_, deltaT, std::integral_constant<int, localSpaceDimension>());

        auto solution = [=](const auto& x, K t) { return initialData(x, t);};
        Result<K, Model::m> result;
        computeErrors(basis, vectorBackend, localSubTriangulation, localOperator, model, solution, t, configuration, result);

        if (configuration.computeOperatorNorm_) {
            result.operatorNorm_ = computeOperatorNorm(globalOperator);
        }

        result.minCellFraction_ = classifier.minSmallCellFraction();
        result.maxCellFraction_ = classifier.maxSmallCellFraction();

        result.minCellFraction_ = geometry.gridView().grid().comm().min(classifier.minSmallCellFraction());
        result.maxCellFraction_ = geometry.gridView().grid().comm().min(classifier.maxSmallCellFraction());

        return result;
    }
} // end namespace Dune::Hypercut

#endif
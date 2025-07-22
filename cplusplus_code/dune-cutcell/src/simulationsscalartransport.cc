#include <cmath>
#include <iostream>
// #include <format>
#include <omp.h>

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynmatrixev.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bvector.hh>

#include <dune/hypercut/scalartransportmodel.hh>
#include <dune/hypercut/linearsystemmodel.hh>
#include <dune/hypercut/scalardgoperator.hh>
#include <dune/hypercut/timesteppingmethods.hh>
#include <dune/hypercut/betaseminorm.hh>
#include <dune/hypercut/execution.hh>

const int dim = 2;
const int order = 2;
const int m = 1;
using K = double;
using Coordinate = Dune::FieldVector<K, dim>;
using Configuration = Dune::Hypercut::Configuration<K>;
using Result = Dune::Hypercut::Result<K, m>;

Dune::FieldVector<K, m> rotatingGaussPeak(const Coordinate& x, Coordinate center, K radius, const Coordinate& rotationCenter, K angle)
{
    center -= rotationCenter;
    Coordinate rotatedCenter(0.0);
    rotatedCenter[0] = std::cos(angle) * center[0] - std::sin(angle) * center[1];
    rotatedCenter[1] = std::sin(angle) * center[0] + std::cos(angle) * center[1];
    rotatedCenter += rotationCenter;

    K d = (x-rotatedCenter).two_norm();
    K w = radius / 4.0;
    return std::exp(-d*d/w/w);
}

Dune::FieldVector<K, m> oneDimBump(Coordinate x, const Coordinate& vf, const Coordinate& center, K t, K bumpStart, K bumpEnd)
{
    x -= t * vf;
    K dot = (x - center) * (vf / vf.two_norm());

    if (dot > bumpStart && dot < bumpEnd) {
      return std::exp(1.0 / ((dot - bumpStart) * (dot - bumpEnd))) / std::exp (1.0 / (bumpStart * bumpEnd));
    }

    return 0.0;
}

Dune::FieldVector<K, m> gaussPeak(const Coordinate& x, const Coordinate& center, K radius)
{
    K d = (x-center).two_norm();
    K w = radius / 4.0;
    Dune::FieldVector<K, m> res(0.0);
    res[0] = std::exp(-d*d/w/w);
    return res;
}

template<class Geometry, class Vectorfield, class BoundaryCondition, class InitialData>
Result scalarTransportTestCase(const Geometry& geometry, const Vectorfield& vectorfield, const BoundaryCondition& boundaryCondition, const InitialData& initialData, const Configuration& configuration)
{
  using Model = Dune::Hypercut::ScalarTransportModel<Vectorfield>;
  using Vector = Dune::BlockVector<K>;

  Model model(vectorfield);

  auto operatorFactory = [] (auto& localSubTriangulation, const auto& basis, const auto& localOperator, auto& classifier, K cfl, K capacityFactor)
  {
    using GlobalOperator = Dune::Hypercut::ScalarDGOperator<std::decay_t<decltype(localSubTriangulation)>, std::decay_t<decltype(basis)>, Vector, std::decay_t<decltype(localOperator)>, std::decay_t<decltype(classifier)>>;
    return GlobalOperator(localSubTriangulation, basis, localOperator, classifier, 1.0);
  };

  return testCase(geometry, model, operatorFactory, boundaryCondition, initialData, configuration, std::integral_constant<int, order>());
}

void rampTest(const Dune::ParameterTree& parameterTree, const Configuration& configuration)
{
  using Geometry = Dune::SubTriangulation::Ramp<K>;
  using Vectorfield = Dune::Hypercut::VaryingRampVectorfield<K>;
  using BoundaryCondition = Dune::Hypercut::RampBoundary<Vectorfield>;

  K rampAngle = parameterTree.get<double>("rampTest.rampAngle");
  K offset = parameterTree.get<double>("rampTest.offset");
  Coordinate center(0.0);
  center[0] = offset;

  rampAngle = (rampAngle / 180.0l) * M_PI;
  Dune::FieldVector<K, dim> normal { std::cos(rampAngle + 0.5 * M_PI), std::sin(rampAngle + 0.5 * M_PI) };

  Geometry geometry(configuration.numberOfCellsPerAxis_, rampAngle + 0.5 * M_PI, offset);
  Vectorfield vectorfield(center, normal, rampAngle);
  BoundaryCondition boundaryCondition(vectorfield, center, rampAngle);

  auto initialData = [=] (const auto& x, K t) {
    return boundaryCondition(Dune::FieldVector<K, m>(0.0), x, Coordinate(0.0), t);
  };

  Result result = scalarTransportTestCase(geometry, vectorfield, boundaryCondition, initialData, configuration);

  std::ofstream file(std::string("error") + std::to_string(configuration.numberOfCellsPerAxis_)
                     + std::string("x") + std::to_string(configuration.numberOfCellsPerAxis_));
  file << "{\n";
  file << "  \"order\": " << order << ",\n";
  file << "  \"numberOfCellsPerAxis\": " << configuration.numberOfCellsPerAxis_ << ",\n";
  file << "  \"cflFactor\": " << configuration.cflSafetyFactor_ << ",\n";
  file << "  \"L-infError\": [" << result.infError_ << "],\n";
  file << "  \"L2Error\": [" << result.l2Error_ << "],\n";
  file << "}";
}

void circleTest(const Dune::ParameterTree& parameterTree, const Configuration& configuration)
{
  using Geometry = Dune::SubTriangulation::Disk<K>;
  using Vectorfield = Dune::Hypercut::RotatingVectorfield<K>;
  using BoundaryCondition = Dune::Hypercut::OutgoingBoundaryCondition<m, K>;

  K radius = parameterTree.get<double>("circleTest.radius");
  K offset = parameterTree.get<double>("circleTest.offset");

  Geometry geometry(configuration.numberOfCellsPerAxis_, radius, offset);
  Vectorfield vectorfield(Coordinate(geometry.center()));
  BoundaryCondition boundaryCondition;

  auto initialData = [&] (const auto& x, K t) {
    Coordinate center = geometry.center();
    center[0] -= 0.9 * radius;
    return rotatingGaussPeak(x, center, 0.5, Coordinate(geometry.center()), t * 2.0 * M_PI);
  };

  Result result = scalarTransportTestCase(geometry, vectorfield, boundaryCondition, initialData, configuration);

  std::ofstream file(std::string("error") + std::to_string(configuration.numberOfCellsPerAxis_)
                     + std::string("x") + std::to_string(configuration.numberOfCellsPerAxis_));
  file << "{\n";
  file << "  \"order\": " << order << ",\n";
  file << "  \"numberOfCellsPerAxis\": " << configuration.numberOfCellsPerAxis_ << ",\n";
  file << "  \"cflFactor\": " << configuration.cflSafetyFactor_ << ",\n";
  file << "  \"L-infError\": [" << result.infError_ << "],\n";
  file << "  \"L2Error\": [" << result.l2Error_ << "],\n";
  file << "}";

  if (configuration.computeOperatorNorm_) {
    std::ofstream file(std::string("operator-norm") + std::to_string(configuration.numberOfCellsPerAxis_)
                       + std::string("-") + std::to_string(offset));
    file << "{\n";
    file << "  \"order\": " << order << ",\n";
    file << "  \"numberOfCellsPerAxis\": " << configuration.numberOfCellsPerAxis_ << ",\n";
    file << "  \"operatorNorm\": " << result.operatorNorm_ << ",\n";
    file << "  \"minCellFraction\": " << result.minCellFraction_ << ",\n";
    file << "  \"maxCellFraction\": " << result.maxCellFraction_ << "\n";
    file << "}";
  }
}

void openChannelTest(const Dune::ParameterTree& parameterTree, const Configuration& configuration)
{
  using Geometry = Dune::SubTriangulation::RotatedOpenChannel<K>;
  using Vectorfield = Dune::Hypercut::ConstantVectorfield<K>;
  using BoundaryCondition = Dune::Hypercut::ZeroBoundary<m, K>;

  K offset = parameterTree.get<double>("channelTest.offset");
  K angle = parameterTree.get<double>("channelTest.channelAngle");
  K vfAngle = parameterTree.get<double>("channelTest.vfAngle");
  K length = parameterTree.get<double>("channelTest.length");

  vfAngle = (vfAngle / 180.0l) * M_PI;
  Dune::FieldVector<K, dim> direction;
  direction[0] = std::cos(vfAngle);
  direction[1] = std::sin(vfAngle);

  Geometry geometry(configuration.numberOfCellsPerAxis_, length, angle, offset);
  Vectorfield vectorfield(direction);
  BoundaryCondition boundaryCondition;

  auto initialData = [=] (const auto& x, K t)
  {
    Coordinate center1(0.0);
    Coordinate center2(0.0);

    center1[0] = offset;
    center2[1] = offset;

    center1 += 0.2 * direction;
    center2 += 0.2 * direction;

    return gaussPeak(x - t * direction, center1, 0.2) + gaussPeak(x - t * direction, center2, 0.2);
  };

  Result result = scalarTransportTestCase(geometry, vectorfield, boundaryCondition, initialData, configuration);

  std::ofstream file(std::string("error") + std::to_string(configuration.numberOfCellsPerAxis_)
                     + std::string("x") + std::to_string(configuration.numberOfCellsPerAxis_));
  file << "{\n";
  file << "  \"order\": " << order << ",\n";
  file << "  \"numberOfCellsPerAxis\": " << configuration.numberOfCellsPerAxis_ << ",\n";
  file << "  \"vfAngle\": " << vfAngle << ",\n";
  file << "  \"cutAngle\": " << angle << ",\n";
  file << "  \"offset\": " << offset << ",\n";
  file << "  \"tau\": " << configuration.tau_ << ",\n";
  file << "  \"cflFactor\": " << configuration.cflSafetyFactor_ << ",\n";
  file << "  \"L-infError\": [" << result.infError_ << "],\n";
  file << "  \"L2Error\": [" << result.l2Error_ << "],\n";
  file << "  \"minCellFraction\": " << result.minCellFraction_ << ",\n";
  file << "  \"maxCellFraction\": " << result.maxCellFraction_ << "\n";
  file << "}";

  if (configuration.computeOperatorNorm_) {
    std::ofstream file(std::string("operator-norm") + std::to_string(configuration.numberOfCellsPerAxis_)
                       + std::string("-") + std::to_string(offset));
    file << "{\n";
    file << "  \"order\": " << order << ",\n";
    file << "  \"numberOfCellsPerAxis\": " << configuration.numberOfCellsPerAxis_ << ",\n";
    file << "  \"offset\": " << offset << ",\n";
    file << "  \"tau\": " << configuration.tau_ << ",\n";
    file << "  \"operatorNorm\": " << result.operatorNorm_ << ",\n";
    file << "  \"minCellFraction\": " << result.minCellFraction_ << ",\n";
    file << "  \"maxCellFraction\": " << result.maxCellFraction_ << "\n";
    file << "}";
  }
}

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    Configuration configuration;
    Dune::ParameterTree parameterTree;
    if (!readCommonParameters(parameterTree, configuration, helper, argc, argv)) {
      return 1;
    }

    openChannelTest(parameterTree, configuration);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }

  return 1;
}

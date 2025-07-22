#ifndef DUNE_HYPERCUT_RUMETIME_PARAMETER_HH
#define DUNE_HYPERCUT_RUMETIME_PARAMETER_HH

#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

namespace Dune::Hypercut {

template <class K>
struct Configuration
{
  K l2StabilityFactor_;
  K cflSafetyFactor_;
  K tau_;
  K T_;
  K smallCellThreshold_;
  K levelSetValueTolerance_;
  std::string timesteppingMethod_;
  int numberOfCellsPerAxis_;
  int intorderadd_;
  bool writeVtk_;
  bool computeOperatorNorm_;
  bool computeLInfNorm_;
  bool computeL1Norm_;
  bool computeL2Norm_;
};

template <class K>
bool readCommonParameters(Dune::ParameterTree& parameterTree, Configuration<K>& configuration,
                          const Dune::MPIHelper& helper, int argc, char** argv)
{
  Dune::ParameterTreeParser::readOptions(argc, argv, parameterTree);

  if (!parameterTree.hasKey("config")) {
    std::cout << "No configuration file provided. Aborting.." << std::endl;
    return false;
  }

  std::string configFile = parameterTree.get("config", "");

  Dune::ParameterTreeParser::readINITree(configFile, parameterTree, false);

  configuration.numberOfCellsPerAxis_ = parameterTree.get<int>("xresolution");

  configuration.intorderadd_ = parameterTree.get<int>("discretization.intorderadd");

  configuration.l2StabilityFactor_ = parameterTree.get<double>("stabilization.l2StabilityFactor");

  configuration.cflSafetyFactor_    = parameterTree.get<double>("timestepping.cflSafetyFactor");
  configuration.T_                  = parameterTree.get<double>("timestepping.T");
  configuration.timesteppingMethod_ = parameterTree["timestepping.method"];
  configuration.tau_                = parameterTree.get<double>("timestepping.tau");

  configuration.levelSetValueTolerance_ = parameterTree.get<double>("levelSet.valueTolerance");

  configuration.writeVtk_            = parameterTree.get<bool>("output.vtk");
  configuration.computeOperatorNorm_ = parameterTree.get<bool>("output.operatorNorm");
  configuration.computeLInfNorm_     = parameterTree.get<bool>("output.linfNorm");
  configuration.computeL1Norm_       = parameterTree.get<bool>("output.l1Norm");
  configuration.computeL2Norm_       = parameterTree.get<bool>("output.l2Norm");

  configuration.smallCellThreshold_
      = parameterTree.get<double>("smallCellClassification.threshold");

  return true;
}

} // end namespace Dune::Hypercut

#endif
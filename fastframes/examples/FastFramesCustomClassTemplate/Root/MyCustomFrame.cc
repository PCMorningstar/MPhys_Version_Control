#include "MyCustomFrame/MyCustomFrame.h"

#include "FastFrames/DefineHelpers.h"
#include "FastFrames/UniqueSampleID.h"

ROOT::RDF::RNode MyCustomFrame::defineVariables(ROOT::RDF::RNode mainNode,
                                                const std::shared_ptr<Sample>& /*sample*/,
                                                const UniqueSampleID& /*id*/) {

  // You can also use the UniqueSampleID object to apply a custom defione
  // based on the sample and the subsample
  //   sample->name(): is the name of the sample defined in the config
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define

  return mainNode;
}

ROOT::RDF::RNode MyCustomFrame::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                      const std::shared_ptr<Sample>& /*sample*/,
                                                      const UniqueSampleID& /*id*/) {

  return mainNode;
}

ROOT::RDF::RNode MyCustomFrame::defineVariablesTruth(ROOT::RDF::RNode node,
                                                     const std::string& /*sample*/,
                                                     const std::shared_ptr<Sample>& /*sample*/,
                                                     const UniqueSampleID& /*sampleID*/) {
  return node;
}

ROOT::RDF::RNode MyCustomFrame::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}
  

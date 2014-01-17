//
// Open a parameter set file and print it.
//
// $Id: StandAloneParameterSet_module.cc,v 1.1 2014/01/17 19:41:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/01/17 19:41:51 $
//
// Contact person Rob Kutschke
//

#include "GeneralUtilities/inc/ParameterSetFromFile.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// C++ includes.
#include <iostream>

namespace mu2e {

  class StandAloneParameterSet : public art::EDAnalyzer {
  public:

    explicit StandAloneParameterSet(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event& e);

  };

  StandAloneParameterSet::StandAloneParameterSet(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset){

    std::string filename=pset.get<std::string>("psetFile");

    // Look for the specified file in FHICL_FILE_PATH.
    // If found, open it and turn it into a fhcil::ParameterSet.
    ParameterSetFromFile params(filename);

    std::cout <<"The file name is: " << filename << std::endl;
    std::cout << params.pSet().to_indented_string() << std::endl;

  }

  void StandAloneParameterSet::analyze( const art::Event& event ) {
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::StandAloneParameterSet);

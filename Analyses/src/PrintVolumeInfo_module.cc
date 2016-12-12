//
// Print the information from the VolumeInfo
//
// $Id: PrintVolumeInfo_module.cc,v 1.1 2013/12/20 20:05:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/12/20 20:05:12 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class PrintVolumeInfo : public art::EDAnalyzer {
  public:

    explicit PrintVolumeInfo(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e) override;

    void beginRun ( const art::Run& r) override;

  private:

    art::InputTag _g4Tag;

  };

  PrintVolumeInfo::PrintVolumeInfo(fhicl::ParameterSet const& pset):
    EDAnalyzer(pset),
    _g4Tag(pset.get<string>("g4Tag")){
  }

  void PrintVolumeInfo::analyze(const art::Event& ){}

  void PrintVolumeInfo::beginRun(const art::Run& run){

    art::Handle<PhysicalVolumeInfoCollection> volsHandle;
    run.getByLabel(_g4Tag,volsHandle);
    PhysicalVolumeInfoCollection vols(*volsHandle);

    for ( auto const& vol : vols ){
      cout << "Volume: " << vol << endl;
    }

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PrintVolumeInfo);

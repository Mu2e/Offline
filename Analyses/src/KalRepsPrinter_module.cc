//
// Print minimal information about KalReps.
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/PathnameWithNextVersion.hh"
#include "Mu2eUtilities/inc/decodeTrackPatRecType.hh"
#include "Mu2eUtilities/inc/TrackPatRecType.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "TrkDiag/inc/KalDiag.hh"
#include "TrkDiag/inc/TrkInfo.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

using namespace std;

namespace mu2e {

  class KalRepsPrinter : public art::EDAnalyzer {
  public:

    explicit KalRepsPrinter(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e) override;

  private:

    // Module label of the module that created a KalRepPtrCollection.
    art::InputTag tracksTag_;

    // Tool for creating analysis ntuples.
    KalDiag kaldiag_;

    // Output stream; either to cout or to a named file.
    std::string outputFilename_;
    std::unique_ptr<std::ofstream> outputStream_;
    std::ostream& out_;

  };
}  // end namespace mu2e

mu2e::KalRepsPrinter::KalRepsPrinter(fhicl::ParameterSet const& pset) :
  art::EDAnalyzer(pset),
  tracksTag_(pset.get<std::string>("tracksTag")),
  kaldiag_(pset.get<fhicl::ParameterSet>("kalDiag")),
  outputFilename_(pset.get<std::string>("outputFilename","")),
  outputStream_( outputFilename_.empty() ?
                 std::unique_ptr<std::ofstream>() :
                 std::make_unique<std::ofstream>(PathnameWithNextVersion(outputFilename_).pathname().c_str())),
  out_( outputStream_ ? *outputStream_ : std::cout ){
}

void mu2e::KalRepsPrinter::analyze(const art::Event& event) {

  auto ptrs = event.getValidHandle<KalRepPtrCollection>(tracksTag_);
  out_ << tracksTag_.encode()
       << "  Event: " << event.id()
       << "  Number of tracks: " << ptrs->size()
       << endl;

  for ( auto const& ptr : *ptrs ){
     TrackPatRecType type = decodeTrackPatRecType( ptr, event);

    if(!ptr->fitCurrent()) {
      throw cet::exception("BADINPUT")
        <<"KalRepsPrinter: do not know what to do with a fitCurrent==0 track\n";
    }

    TrkInfo track;
    kaldiag_.fillTrkInfo( ptr.get(), track);
    TrkFitInfo track_ent;
    kaldiag_.fillTrkFitInfo( ptr.get(), track_ent);

    out_ << "    : "
         << tracksTag_.label() << " "
         << ptr                << " "
         << type.name()
         << " Status: "        << track._status
         << " nHits:  "        << track._nhits
         << " nActive: "       << track._nactive
         << " ent.fitmom:  "   << track_ent._fitmom
         << " ent.fitmomerr: " << track_ent._fitmomerr
         << " ent.tandip: "    << track_ent._fitpar._td
         << " ent.trkqual: " << track._trkqual
         << endl;
  }

} // end analyze

DEFINE_ART_MODULE(mu2e::KalRepsPrinter);

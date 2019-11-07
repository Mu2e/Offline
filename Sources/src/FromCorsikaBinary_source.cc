// Transfer into framework CORSIKA binary files
//
// Original author: Stefano Roberto Soleti, 2019

#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>
#include <string>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Utilities/Globals.h" // FIXME-KJK: should not be necessary to use this
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/BranchType.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Provenance/canonicalProductName.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "Sources/inc/CosmicCORSIKA.hh"
#include "SeedService/inc/SeedService.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

using namespace std;

namespace mu2e {

    //================================================================
    class CorsikaBinaryDetail : private boost::noncopyable {
      std::string myModuleLabel_;
      art::SourceHelper const& pm_;
      unsigned runNumber_; // from ParSet
      art::SubRunID lastSubRunID_;
      std::set<art::SubRunID> seenSRIDs_;

      std::string currentFileName_;
      FILE *currentFile_ = nullptr;
      float garbage;

      unsigned currentSubRunNumber_; // from file
      // A helper function used to manage the principals.
      // This is boilerplate that does not change if you change the data products.
      void managePrincipals ( int runNumber,
                              int subRunNumber,
                              int eventNumber,
                              art::RunPrincipal*&    outR,
                              art::SubRunPrincipal*& outSR,
                              art::EventPrincipal*&  outE);
      unsigned getSubRunNumber(const std::string& filename) const;

      unsigned currentEventNumber_;

      art::ProductID particlesPID_;

    public:
      CorsikaBinaryDetail(const Parameters &conf,
                          art::ProductRegistryHelper &,
                          const art::SourceHelper &);

      void readFile(std::string const& filename, art::FileBlock*& fb);

      bool readNext(art::RunPrincipal* const& inR,
                    art::SubRunPrincipal* const& inSR,
                    art::RunPrincipal*& outR,
                    art::SubRunPrincipal*& outSR,
                    art::EventPrincipal*& outE);

      void closeCurrentFile();

      CosmicCORSIKA _corsikaGen;

    };

    //----------------------------------------------------------------
    CorsikaBinaryDetail::CorsikaBinaryDetail(const Parameters& conf,
                                             art::ProductRegistryHelper& rh,
                                             const art::SourceHelper& pm)
      : myModuleLabel_("FromCorsikaBinary")
      , pm_(pm)
      , runNumber_(conf().runNumber())
      , currentSubRunNumber_(-1U)
      , currentEventNumber_(0)
      , _corsikaGen(conf())
    {
      if(!art::RunID(runNumber_).isValid()) {
        throw cet::exception("BADCONFIG", " FromCorsikaBinary: ")
          << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
      }

      rh.reconstitutes<mu2e::GenParticleCollection,art::InEvent>(myModuleLabel_);

      // auto get_ProductID = [](auto const& typeLabel, auto const& processName) {
      //   if (!typeLabel.hasEmulatedModule()) {
      //     throw cet::exception("BADCONFIG", " FromCorsikaBinary: ")
      //     << " Must provided emulated module name for reconstituted product.\n";
      //   }
      //   auto const canonical_product_name = art::canonicalProductName(typeLabel.friendlyClassName(),
      //                                                                 typeLabel.emulatedModule(),
      //                                                                 typeLabel.productInstanceName(),
      //                                                                 processName);
      //   return art::ProductID{canonical_product_name};
      // };
      // particlesPID_ = get_ProductID(gpcTypeLabel, art::Globals::instance()->processName());

    }

    //----------------------------------------------------------------
    void CorsikaBinaryDetail::readFile(const std::string& filename, art::FileBlock*& fb) {

      currentFileName_ = filename;
      currentSubRunNumber_ = getSubRunNumber(filename);
      currentEventNumber_ = 0;

      currentFile_ = fopen(filename.c_str(), "r");
      _corsikaGen.openFile(currentFile_);
      fb = new art::FileBlock(art::FileFormatVersion(1, "CorsikaBinaryInput"), currentFileName_);
    }

    //----------------------------------------------------------------
    unsigned CorsikaBinaryDetail::getSubRunNumber(const std::string& filename) const {
      const std::string::size_type islash = filename.find("DAT") ;
      const std::string basename = (islash == std::string::npos) ? filename : filename.substr(islash + 3);
      unsigned sr(-1);
      std::istringstream is(basename);
      if(!(is>>sr)) {
        throw cet::exception("BADINPUTS")<<"Expect an unsigned integer at the beginning of input file name, got "<<basename<<"\n";
      }
      return sr;
    }

    //----------------------------------------------------------------
    void CorsikaBinaryDetail::closeCurrentFile() {
      currentFileName_ = "";
      fclose(currentFile_);
    }

    //----------------------------------------------------------------
    bool CorsikaBinaryDetail::readNext(art::RunPrincipal* const& inR,
                                       art::SubRunPrincipal* const& inSR,
                                       art::RunPrincipal*& outR,
                                       art::SubRunPrincipal*& outSR,
                                       art::EventPrincipal*& outE)
    {
      std::unique_ptr<GenParticleCollection> particles(new GenParticleCollection());
      bool still_data = (_corsikaGen.generate(*particles));
      if (!still_data) {
        return false;
      }

      managePrincipals(runNumber_, currentSubRunNumber_, ++currentEventNumber_, outR, outSR, outE);
      art::put_product_in_principal(std::move(particles), *outE, myModuleLabel_);

      return true;

    } // readNext()


  // Each time that we encounter a new run, a new subRun or a new event, we need to make a new principal
  // of the appropriate type.  This code does not need to change as the number and type of data products changes.
  void CorsikaBinaryDetail::managePrincipals ( int runNumber,
                                          int subRunNumber,
                                          int eventNumber,
                                          art::RunPrincipal*&    outR,
                                          art::SubRunPrincipal*& outSR,
                                          art::EventPrincipal*&  outE){

    art::Timestamp ts;

    art::SubRunID newID(runNumber, subRunNumber);
    // if(newID.runID() != lastSubRunID_.runID()) {
    //   // art takes ownership of the object pointed to by outR and will delete it at the appropriate time.
    //   outR = pm_.makeRunPrincipal(runNumber, ts);
    // }

    if(newID != lastSubRunID_) {
      // art takes ownership of the object pointed to by outSR and will delete it at the appropriate time.
      outR = pm_.makeRunPrincipal(runNumber, ts);
      outSR = pm_.makeSubRunPrincipal(runNumber,
                                      subRunNumber,
                                      ts);
    }
    lastSubRunID_ = newID;

    // art takes ownership of the object pointed to by outE and will delete it at the appropriate time.
    outE = pm_.makeEventPrincipal(runNumber, subRunNumber, eventNumber, ts, false);

  } // managePrincipals()
    //----------------------------------------------------------------

} // namespace mu2e

typedef art::Source<mu2e::CorsikaBinaryDetail> FromCorsikaBinary;
DEFINE_ART_INPUT_SOURCE(FromCorsikaBinary);

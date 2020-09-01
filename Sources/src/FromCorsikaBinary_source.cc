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
#include "CLHEP/Units/SystemOfUnits.h"

#include "fhiclcpp/types/Name.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
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
#include "MCDataProducts/inc/CosmicLivetime.hh"

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
      std::ifstream *currentFile_ = nullptr;
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

      unsigned currentEventNumber_;

      art::ProductID particlesPID_;

      float _area;  // m2
      float _lowE;  // GeV
      float _highE; // GeV
      float _fluxConstant;

      const float _mm22m2 = CLHEP::mm2 / CLHEP::m2;

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
      , _fluxConstant(conf().fluxConstant())
      , _corsikaGen(conf(), art::ServiceHandle<SeedService>{}->getInputSourceSeed())
    {
      if(!art::RunID(runNumber_).isValid()) {
        throw cet::exception("BADCONFIG", " FromCorsikaBinary: ")
          << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
      }

      rh.reconstitutes<mu2e::GenParticleCollection,art::InEvent>(myModuleLabel_);
      rh.reconstitutes<mu2e::CosmicLivetime,art::InEvent>(myModuleLabel_);
      _area = (conf().targetBoxXmax() + 2 * conf().showerAreaExtension() - conf().targetBoxXmin())
            * (conf().targetBoxZmax() + 2 * conf().showerAreaExtension() - conf().targetBoxZmin()) * _mm22m2; // m^2
    }

    //----------------------------------------------------------------
    void CorsikaBinaryDetail::readFile(const std::string& filename, art::FileBlock*& fb) {

      currentFileName_ = filename;
      currentEventNumber_ = 0;

      currentFile_ = new ifstream(currentFileName_);

      unsigned subrun = 0;
      float lowE, highE;
      _corsikaGen.openFile(currentFile_, subrun, lowE, highE);
      currentSubRunNumber_ = subrun;
      _lowE = lowE;
      _highE = highE;
      fb = new art::FileBlock(art::FileFormatVersion(1, "CorsikaBinaryInput"), currentFileName_);
    }

    //----------------------------------------------------------------
    void CorsikaBinaryDetail::closeCurrentFile() {
      currentFileName_ = "";
      currentFile_->close();
    }

    //----------------------------------------------------------------
    bool CorsikaBinaryDetail::readNext(art::RunPrincipal* const& inR,
                                       art::SubRunPrincipal* const& inSR,
                                       art::RunPrincipal*& outR,
                                       art::SubRunPrincipal*& outSR,
                                       art::EventPrincipal*& outE)
    {
      std::unique_ptr<GenParticleCollection> particles(new GenParticleCollection());
      unsigned int primaries;
      bool still_data = _corsikaGen.generate(*particles, primaries);

      if (!still_data) {
        return false;
      }

      managePrincipals(runNumber_, currentSubRunNumber_, ++currentEventNumber_, outR, outSR, outE);
      art::put_product_in_principal(std::move(particles), *outE, myModuleLabel_);
      std::unique_ptr<CosmicLivetime> livetime(new CosmicLivetime(primaries, _area, _lowE, _highE, _fluxConstant));
      art::put_product_in_principal(std::move(livetime), *outE, myModuleLabel_);
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

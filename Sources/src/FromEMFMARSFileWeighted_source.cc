// Transfer into framework complete MARS outputs, including
// per-particle weights.  While weighted particles can not be directly
// used to simulate detector hits, this information is useful when the
// same dataset needs to be re-used repeatedly with additional
// randomization.
//
// Original author: Andrei Gaponenko, 2012

#include "Sources/inc/ExtMonFNALMARSUtils.hh"

#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>

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

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFMARSWeightedDetail : private boost::noncopyable {
      std::string myModuleLabel_;
      art::SourceHelper const& pm_;
      unsigned runNumber_; // from ParSet
      art::SubRunID currentSubRunID_;
      std::set<art::SubRunID> seenSRIDs_;

      std::string currentFileName_;
      std::ifstream currentFile_;
      unsigned currentSubRunNumber_; // from file name

      // We may need to consume several lines of the input file to produce a single event.
      //
      // This data structure holds data from an input line that
      // was read from the file but not used to make a event yet.
      //
      // The var is pre-filled by readFile() then by readNext(), so
      // it is non-empty at the beginning of readNext() - unless we
      // ran out of data in the current file.
      MARSParticle currentLine_;

      MARSMu2eConverter cnv_;

      unsigned getSubRunNumber(const std::string& filename) const;

      unsigned currentEventNumber_;

      art::ProductID particlesPID_;
      art::ProductID marsInfoPID_;

    public:

      EMFMARSWeightedDetail(const fhicl::ParameterSet&,
                            art::ProductRegistryHelper&,
                            const art::SourceHelper&);

      void readFile(std::string const& filename, art::FileBlock*& fb);

      bool readNext(art::RunPrincipal* const& inR,
                    art::SubRunPrincipal* const& inSR,
                    art::RunPrincipal*& outR,
                    art::SubRunPrincipal*& outSR,
                    art::EventPrincipal*& outE);

      void closeCurrentFile();
    };

    //----------------------------------------------------------------
    EMFMARSWeightedDetail::EMFMARSWeightedDetail(const fhicl::ParameterSet& pset,
                                                 art::ProductRegistryHelper& rh,
                                                 const art::SourceHelper& pm)
      : myModuleLabel_("FromEMFMARSFileWeighted")
      , pm_(pm)
      , runNumber_(pset.get<unsigned>("runNumber"))
      , currentSubRunNumber_(-1U)
      , currentEventNumber_(0)
    {
      if(!art::RunID(runNumber_).isValid()) {
        throw cet::exception("BADCONFIG", " FromEMFMARSFileWeighted: ")
          << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
      }

      auto const& gpcTypeLabel = rh.reconstitutes<mu2e::GenParticleCollection,art::InEvent>(myModuleLabel_);
      auto const& marsTypeLabel = rh.reconstitutes<mu2e::MARSInfoCollection,art::InEvent>(myModuleLabel_);
      rh.reconstitutes<mu2e::GenParticleMARSAssns,art::InEvent>(myModuleLabel_);

      auto get_ProductID = [](auto const& typeLabel, auto const& processName) {
        if (!typeLabel.hasEmulatedModule()) {
          throw cet::exception("BADCONFIG", " FromEMFMARSFileWeighted: ")
          << " Must provided emulated module name for reconstituted product.\n";
        }
        auto const canonical_product_name = art::canonicalProductName(typeLabel.friendlyClassName(),
                                                                      typeLabel.emulatedModule(),
                                                                      typeLabel.productInstanceName(),
                                                                      processName);
        return art::ProductID{canonical_product_name};
      };
      particlesPID_ = get_ProductID(gpcTypeLabel, art::Globals::instance()->processName());
      marsInfoPID_ = get_ProductID(marsTypeLabel, art::Globals::instance()->processName());
    }

    //----------------------------------------------------------------
    void EMFMARSWeightedDetail::readFile(const std::string& filename, art::FileBlock*& fb) {
      if(currentLine_.protonNumber != -1U) {
        throw cet::exception("UNKNOWN")<<"FromEMFMARSFileWeighted: readFile() called while currentLine_ not empty. Something is wrong.\n";
      }

      currentFileName_ = filename;
      currentSubRunNumber_ = getSubRunNumber(filename);

      if(!seenSRIDs_.insert(art::SubRunID(runNumber_, currentSubRunNumber_)).second) {
        ++runNumber_;
        const bool inserted [[gnu::unused]] = seenSRIDs_.insert(art::SubRunID(runNumber_, currentSubRunNumber_)).second;
        assert(inserted);
      }

      currentEventNumber_ = 0;

      if(!art::SubRunID(runNumber_, currentSubRunNumber_).isValid()) {
        throw cet::exception("BADINPUTS")<<"Got invalid subrun number="<<currentSubRunNumber_<<" from input filename="<<filename<<"\n";
      }
      currentFile_.open(filename.c_str());
      fb = new art::FileBlock(art::FileFormatVersion(1, "EMFMARSWeightedinput"), currentFileName_);

      readMARSLine(currentFile_, currentLine_);
    }

    //----------------------------------------------------------------
    unsigned EMFMARSWeightedDetail::getSubRunNumber(const std::string& filename) const {
      const std::string::size_type islash = filename.rfind('/') ;
      const std::string basename = (islash == std::string::npos) ? filename : filename.substr(islash + 1);

      unsigned sr(-1);
      std::istringstream is(basename);
      if(!(is>>sr)) {
        throw cet::exception("BADINPUTS")<<"Expect an unsigned integer at the beginning of input file name, got "<<basename<<"\n";
      }
      return sr;
    }

    //----------------------------------------------------------------
    void EMFMARSWeightedDetail::closeCurrentFile() {
      currentFileName_ = "";
      currentFile_.close();
    }

    //----------------------------------------------------------------
    bool EMFMARSWeightedDetail::readNext(art::RunPrincipal* const& inR,
                                         art::SubRunPrincipal* const& inSR,
                                         art::RunPrincipal*& outR,
                                         art::SubRunPrincipal*& outSR,
                                         art::EventPrincipal*& outE)
    {
      // currentLine_ pre-filled by readFile() or the previous readNext() call.
      if(currentLine_.protonNumber == -1U) { // no more data in the current file
        return false;
      }
      else { // do have more data

        std::unique_ptr<GenParticleCollection> particles(new GenParticleCollection());
        std::unique_ptr<MARSInfoCollection> info(new MARSInfoCollection());
        std::unique_ptr<GenParticleMARSAssns> assns(new GenParticleMARSAssns());

        art::Timestamp ts;

        art::SubRunID newID(runNumber_, currentSubRunNumber_);
        if(newID.runID() != currentSubRunID_.runID()) {
          outR = pm_.makeRunPrincipal(runNumber_, ts);
        }

        if(newID != currentSubRunID_) {
          outSR = pm_.makeSubRunPrincipal(runNumber_,
                                          currentSubRunNumber_,
                                          ts);
          currentSubRunID_ = newID;
        }

        outE = pm_.makeEventPrincipal(runNumber_, currentSubRunNumber_, ++currentEventNumber_, ts, false);

        const art::EDProductGetter* particlesGetter = outE->productGetter(particlesPID_);
        const art::EDProductGetter* infoGetter = outE->productGetter(marsInfoPID_);

        const unsigned nj = currentLine_.protonNumber;
        particles->push_back(cnv_.marsToMu2eParticle(currentLine_));
        info->push_back(MARSInfo(currentLine_.weight, currentLine_.protonNumber, currentSubRunNumber_, runNumber_));
        assns->addSingle(art::Ptr<GenParticle>(particlesPID_, particles->size()-1, particlesGetter),
                         art::Ptr<MARSInfo>(marsInfoPID_, info->size()-1, infoGetter)
                         );

        while(readMARSLine(currentFile_, currentLine_)) {
          if(currentLine_.protonNumber == nj) {
            particles->push_back(cnv_.marsToMu2eParticle(currentLine_));
            info->push_back(MARSInfo(currentLine_.weight, currentLine_.protonNumber, currentSubRunNumber_, runNumber_));
            assns->addSingle(art::Ptr<GenParticle>(particlesPID_, particles->size()-1, particlesGetter),
                             art::Ptr<MARSInfo>(marsInfoPID_, info->size()-1, infoGetter)
                             );
          }
          else {
            break;
          }
        }

        art::put_product_in_principal(std::move(particles), *outE, myModuleLabel_);
        art::put_product_in_principal(std::move(info), *outE, myModuleLabel_);
        art::put_product_in_principal(std::move(assns), *outE, myModuleLabel_);

        return true;

      } // else-have data

    } // readNext()

    //----------------------------------------------------------------

  } // namespace ExtMonFNAL
} // namespace mu2e

typedef art::Source<mu2e::ExtMonFNAL::EMFMARSWeightedDetail> FromEMFMARSFileWeighted;
DEFINE_ART_INPUT_SOURCE(FromEMFMARSFileWeighted);

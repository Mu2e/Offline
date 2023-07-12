// Transfer MARS outputs into framework unweighting the particles.
//
// Original author: Andrei Gaponenko, 2012

#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>

#include "fhiclcpp/types/Name.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "Offline/MCDataProducts/inc/GenParticle.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/Sources/inc/ExtMonFNALMARSUtils.hh"

#include "Offline/SeedService/inc/SeedService.hh"
//#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<" in "<<__func__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {


  //================================================================
  struct HelperEventEntry {
    double remainingWeight;
    GenParticle particle;
    HelperEventEntry(double w, const GenParticle& p)
      : remainingWeight(w), particle(p)
    {}
  };

  struct HelperEvent : public std::vector<HelperEventEntry> {};


  //================================================================
  class ExtMonFNALMARSDetail : private boost::noncopyable {
    std::string myModuleLabel_;
    art::SourceHelper const& pm_;
    unsigned runNumber_; // from ParSet
    art::SubRunID currentSubRunID_;
    std::set<art::SubRunID> seenSRIDs_;

    GlobalConstantsHandle<ParticleDataList> pdt_;

    std::string currentFileName_;
    std::ifstream currentFile_;
    unsigned currentSubRunNumber_; // from file name

    // We may need to consume several lines of the input file to produce a single event.
    // On the other hand, the same set of lines may generate more than one event (case w>1).
    //
    // This data structure holds data from an input line that
    // was read from the file but not used to make a event yet.
    //
    // The var is pre-filled by readFile() then by readNext(), so
    // it is non-empty at the beginning of readNext() - unless we
    // ran out of data in the current file.
    ExtMonFNAL::MARSParticle currentLine_;

    // Non-empty at readNext() if the previous readNext() got weight>1
    HelperEvent he_;

    unsigned getSubRunNumber(const std::string& filename) const;

    unsigned currentEventNumber_;

    //CLHEP::RandFlat random_;
    CLHEP::MTwistEngine random_;

    HelperEventEntry marsToMu2eParticle(const ExtMonFNAL::MARSParticle& mp);

    // convert particle ID, energy, time units
    ExtMonFNAL::MARSMu2eConverter marsMu2eConverter_;

    //----------------------------------------------------------------

  public:

    ExtMonFNALMARSDetail(const fhicl::ParameterSet&,
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
  ExtMonFNALMARSDetail::ExtMonFNALMARSDetail(const fhicl::ParameterSet& pset,
                                             art::ProductRegistryHelper& rh,
                                             const art::SourceHelper& pm)
    : myModuleLabel_("ExtMonFNALMARS")
    , pm_(pm)
    , runNumber_(pset.get<unsigned>("runNumber"))
    , pdt_()
    , currentSubRunNumber_(-1U)
    , currentEventNumber_(0)
      // , random_(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
      // , random_(art::ServiceHandle<art::RandomNumberGenerator>()->createEngine(art::ServiceHandle<SeedService>()->getSeed()))
      //, random_(art::ServiceHandle<SeedService>()->getSeed())
    , random_(pset.get<int>("seed", 1))
  {
    rh.reconstitutes<mu2e::GenParticleCollection,art::InEvent>(myModuleLabel_);

    if(!art::RunID(runNumber_).isValid()) {
      throw cet::exception("BADCONFIG")<<" FromExtMonFNALMARSFile:  fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
    }
  }

  //----------------------------------------------------------------
  void ExtMonFNALMARSDetail::readFile(const std::string& filename, art::FileBlock*& fb) {
    AGDEBUG("filename="<<filename<<", currentLine_="<<currentLine_);
    if(currentLine_.protonNumber != -1U) {
      throw cet::exception("UNKNOWN")<<"FromExtMonFNALMARSFile: readFile() called while currentLine_ not empty. Something is wrong.\n";
    }

    currentFileName_ = filename;
    currentSubRunNumber_ = getSubRunNumber(filename);

    if(!seenSRIDs_.insert(art::SubRunID(runNumber_, currentSubRunNumber_)).second) {
      ++runNumber_;
      const bool inserted [[gnu::unused]]= seenSRIDs_.insert(art::SubRunID(runNumber_, currentSubRunNumber_)).second;
      assert(inserted);
    }

    currentEventNumber_ = 0;

    if(!art::SubRunID(runNumber_, currentSubRunNumber_).isValid()) {
      throw cet::exception("BADINPUTS")<<"Got invalid subrun number="<<currentSubRunNumber_<<" from input filename="<<filename<<"\n";
    }
    currentFile_.open(filename.c_str());
    fb = new art::FileBlock(art::FileFormatVersion(1, "ExtMonFNALMARSinput"), currentFileName_);

    ExtMonFNAL::readMARSLine(currentFile_, currentLine_);
  }

  //----------------------------------------------------------------
  unsigned ExtMonFNALMARSDetail::getSubRunNumber(const std::string& filename) const {
    const std::string::size_type islash = filename.rfind('/') ;
    const std::string basename = (islash == std::string::npos) ? filename : filename.substr(islash + 1);

    AGDEBUG("Got basename = "<<basename);

    unsigned sr(-1);
    std::istringstream is(basename);
    if(!(is>>sr)) {
      throw cet::exception("BADINPUTS")<<"Expect an unsigned integer at the beginning of input file name, got "<<basename<<"\n";
    }
    return sr;
  }

  //----------------------------------------------------------------
  void ExtMonFNALMARSDetail::closeCurrentFile() {
    AGDEBUG("ExtMonFNALMARSDetail: currentFileName_: "<<currentFileName_);
    currentFileName_ = "";
    currentFile_.close();
  }

  //----------------------------------------------------------------
  bool ExtMonFNALMARSDetail::readNext(art::RunPrincipal* const& inR,
                                      art::SubRunPrincipal* const& inSR,
                                      art::RunPrincipal*& outR,
                                      art::SubRunPrincipal*& outSR,
                                      art::EventPrincipal*& outE)
  {
    AGDEBUG("ExtMonFNALMARSDetail::readNext(): starting next event. he_.size()="<<he_.size()<<", currentLine_ = "<<currentLine_);

    if(he_.empty()) {

      // Need to process more lines
      if(currentLine_.protonNumber != -1U) { // do have more data

        const unsigned nj = currentLine_.protonNumber;
        he_.push_back(HelperEventEntry(marsToMu2eParticle(currentLine_)));
        AGDEBUG("Added "<<he_.back().particle);
        while(ExtMonFNAL::readMARSLine(currentFile_, currentLine_)) {
          if(currentLine_.protonNumber == nj) {
            he_.push_back(HelperEventEntry(marsToMu2eParticle(currentLine_)));
            AGDEBUG("Added "<<he_.back().particle);
          }
          else {
            break;
          }
        }
      }
      //----------------------------------------------------------------
    } // he_.emtpy()

    // Generate mu2e event using he_ infos, subtract from he_ weights.
    if(he_.empty()) {
      // at end of file
      return false;
    }
    else {
      art::Timestamp ts;

      art::SubRunID newID(runNumber_, currentSubRunNumber_);
      if(newID.runID() != currentSubRunID_.runID()) {
        AGDEBUG("Creating new run: r="<<runNumber_);
        outR = pm_.makeRunPrincipal(runNumber_, ts);
      }

      if(newID != currentSubRunID_) {
        outSR = pm_.makeSubRunPrincipal(runNumber_,
                                        currentSubRunNumber_,
                                        ts);
        currentSubRunID_ = newID;
      }

      outE = pm_.makeEventPrincipal(runNumber_, currentSubRunNumber_, ++currentEventNumber_, ts, false);

      AGDEBUG("After run/subrun/event principal creation: run="<<runNumber_<<", subrun="<<currentSubRunNumber_<<", event="<<currentEventNumber_);

      std::unique_ptr<GenParticleCollection> particles(new GenParticleCollection());

      HelperEvent leftovers;
      for(HelperEvent::const_iterator i=he_.begin(); i!=he_.end(); ++i) {
        // accept/reject
        AGDEBUG("remainingWeight = "<<i->remainingWeight);
        if(random_.flat() < i->remainingWeight) {
          AGDEBUG("Particle accepted");
          particles->push_back(i->particle);
          double remainingWeight = i->remainingWeight - 1;
          if(remainingWeight > 0) {
            leftovers.push_back(HelperEventEntry(remainingWeight, i->particle));
          }
        }
      }
      he_ = leftovers;

      art::put_product_in_principal(std::move(particles), *outE, myModuleLabel_);
      return true;

    } // else - got he_ data

  } // readNext()

  //----------------------------------------------------------------
  HelperEventEntry ExtMonFNALMARSDetail::marsToMu2eParticle(const ExtMonFNAL::MARSParticle& mp) {
    const int pdgId = marsMu2eConverter_.marsToMu2eParticleCode(mp.pid);

    const double mass = pdt_->particle(pdgId).mass();

    const double energy = mass + marsMu2eConverter_.marsToMu2eEnergy(mp.kineticEnergy);
    const double p3mag = sqrt((energy-mass)*(energy+mass));

    const CLHEP::HepLorentzVector p4(mp.dcx * p3mag,
                                     mp.dcy * p3mag,
                                     mp.dcz * p3mag,
                                     energy
                                     );

    HelperEventEntry res(mp.weight,
                         GenParticle(PDGCode::type(pdgId),
                                     GenId::MARS,
                                     marsMu2eConverter_.marsToMu2ePosition(mp.x, mp.y, mp.z),
                                     p4,
                                     marsMu2eConverter_.marsToMu2eTime(mp.tof)
                                     )
                         );

    return res;
  }

  //----------------------------------------------------------------

} // namespace mu2e

typedef art::Source<mu2e::ExtMonFNALMARSDetail> FromExtMonFNALMARSFile;
DEFINE_ART_INPUT_SOURCE(FromExtMonFNALMARSFile)

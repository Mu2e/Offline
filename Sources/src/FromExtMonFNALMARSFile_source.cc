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

#include "Offline/SeedService/inc/SeedService.hh"
//#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<" in "<<__func__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  // relevant information from a single line of the mars text file
  struct MARSParticle {
    unsigned protonNumber;
    int pid;
    double kineticEnergy;
    double weight;
    double x, y, z;
    double dcx, dcy, dcz;
    double tof;
    MARSParticle() : protonNumber(-1U), pid(), kineticEnergy(), weight(), x(), y(), z(), dcx(), dcy(), dcz(), tof() {}
  };

  std::ostream& operator<<(std::ostream& os, const MARSParticle& mp) {
    return os<<"MARSParticle("<<mp.protonNumber
             <<", pid="<<mp.pid
             <<", kE="<<mp.kineticEnergy
             <<", w="<<mp.weight
             <<", tof="<<mp.tof
             <<")";
  }

  bool readMARSLine(std::istream& file, MARSParticle& res) {
    std::string line;
    std::getline(file, line);

    AGDEBUG("readMARSLine(): got "<<line);

    if(line.empty()) {
      res.protonNumber = -1U;
      return false;
    }
    else {
      std::istringstream is(line);

      if(!(is>>res.protonNumber
           >>res.pid>>res.kineticEnergy>>res.weight
           >>res.x>>res.y>>res.z
           >>res.dcx>>res.dcy>>res.dcz
           >>res.tof)
         ) {

        throw cet::exception("BADINPUT")<<" parseMARSLine(): error processing line: "<<line<<"\n";
      }
      return true;
    }
  }

  //================================================================
  struct HelperEventEntry {
    double remainingWeight;
    GenParticle particle;
    HelperEventEntry(double w, const GenParticle& p)
      : remainingWeight(w), particle(p)
    {}
  };

  struct HelperEvent : public std::vector<HelperEventEntry> {};

  CLHEP::Hep3Vector marsToMu2ePosition(double x, double y, double z) {
    // Apply the offsets and convert cm to mm
    return CLHEP::Hep3Vector(10*(x+390.4), 10*(y+0), 10*(z-906.8));
  }

  double marsToMu2eTime(double t) {
    // seconds to ns
    return t*1.e9;
  }

  double marsToMu2eEnergy(double ke) {
    // GeV to MeV
    return ke * 1.e3;
  }

  // returns pdgId
  int marsToMu2eParticleCode(int marsPID) {
    static std::map<int,int> table;
    if(table.empty()) {
      table[1] = PDGCode::proton; // proton
      table[2] = PDGCode::n0; // neutron
      table[3] = PDGCode::pi_plus; // pi+
      table[4] = PDGCode::pi_minus; // pi-
      table[5] = PDGCode::K_plus; // K+
      table[6] = PDGCode::K_minus; // K-
      table[7] = PDGCode::mu_plus; // mu+
      table[8] = PDGCode::mu_minus; // mu-
      table[9] = PDGCode::gamma; // gamma
      table[10] = PDGCode::e_plus; // e-
      table[11] = PDGCode::e_minus; // e+
      table[12] = PDGCode::anti_proton; // p- (antiproton)
      table[13] = PDGCode::pi0; // pi0
      table[14] = PDGCode::deuteron; // deutron
      table[15] = PDGCode::tritium; // H3
      table[16] = PDGCode::He3; // He3
      table[17] = PDGCode::alpha; // He4
      table[18] = PDGCode::nu_mu; // numu
      table[19] = PDGCode::anti_nu_mu; // numubar
      table[20] = PDGCode::nu_e; // nue
      table[21] = PDGCode::anti_nu_e; // nuebar
      table[22] = PDGCode::K_L0; // K0L
      table[23] = PDGCode::K_S0; // K0S
      //table[24] = ; // K0
      //table[25] = ; // K0bar
      table[26] = PDGCode::Lambda0; // Lambda
      table[27] = PDGCode::anti_Lambda0; // Lambdabar
      table[28] = PDGCode::Sigma_plus; // Sigma+
      table[29] = PDGCode::Sigma0; // Sigma0
      table[30] = PDGCode::Sigma_minus; // Sigma-
      table[31] = PDGCode::anti_n0; // nbar
      table[32] = PDGCode::Xi0; // Xi0
      table[33] = PDGCode::Xi_minus; // Xi-
      table[34] = PDGCode::Omega_minus; // Omega-
      table[35] = PDGCode::anti_Sigma_minus; // anti-Sigma-
      table[36] = PDGCode::anti_Sigma0; // anti-Sigma0
      table[37] = PDGCode::anti_Sigma_plus; // anti-Sigma+
      table[38] = PDGCode::anti_Xi0; // anti-Xi0
      table[39] = PDGCode::anti_Xi_plus; // anti-Xi+
      table[40] = PDGCode::anti_Omega_plus; // anti-Omega+
    }
    std::map<int,int>::const_iterator ip = table.find(marsPID);
    if(ip == table.end()) {
      throw cet::exception("BADINPUT")<<" marsToMu2eParticleCode(): unknonw MARS particle code "<<marsPID<<"\n";
    }
    return ip->second;
  }

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
    MARSParticle currentLine_;

    // Non-empty at readNext() if the previous readNext() got weight>1
    HelperEvent he_;

    unsigned getSubRunNumber(const std::string& filename) const;

    unsigned currentEventNumber_;

    //CLHEP::RandFlat random_;
    CLHEP::MTwistEngine random_;

    HelperEventEntry marsToMu2eParticle(const MARSParticle& mp);

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

    readMARSLine(currentFile_, currentLine_);
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
        while(readMARSLine(currentFile_, currentLine_)) {
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
  HelperEventEntry ExtMonFNALMARSDetail::marsToMu2eParticle(const MARSParticle& mp) {
    const int pdgId = marsToMu2eParticleCode(mp.pid);

    const double mass = pdt_->particle(pdgId).mass();

    const double energy = mass + marsToMu2eEnergy(mp.kineticEnergy);
    const double p3mag = sqrt((energy-mass)*(energy+mass));

    const CLHEP::HepLorentzVector p4(mp.dcx * p3mag,
                                     mp.dcy * p3mag,
                                     mp.dcz * p3mag,
                                     energy
                                     );

    HelperEventEntry res(mp.weight,
                         GenParticle(PDGCode::type(pdgId),
                                     GenId::MARS,
                                     marsToMu2ePosition(mp.x, mp.y, mp.z),
                                     p4,
                                     marsToMu2eTime(mp.tof)
                                     )
                         );

    return res;
  }

  //----------------------------------------------------------------

} // namespace mu2e

typedef art::Source<mu2e::ExtMonFNALMARSDetail> FromExtMonFNALMARSFile;
DEFINE_ART_INPUT_SOURCE(FromExtMonFNALMARSFile);

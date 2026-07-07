//
// Plugin to read StepPoints in PS Vacuum and create ntuples
//
//
// Original author Robert Bernstein
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadAntiProtonSteps : public art::EDAnalyzer {
  public:

    typedef vector<int> Vint;
    typedef SimParticleCollection::key_type key_type;
    double pbarMass;
    double pbarMass2;

    explicit ReadAntiProtonSteps(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _psVacuumStepPoints(pset.get<string>("psVacuumStepPoints","AntiProtonSteps")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _ntAntiProtonSteps(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _diagLevel(pset.get<int>("diagLevel",0))
    {

      Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
      if( pdg_ids.size()>0 ) {
        cout << "ReadAntiProtonSteps: save following particle types in the ntuple: ";
        for( size_t i=0; i<pdg_ids.size(); ++i ) {
          pdg_save.insert(pdg_ids[i]);
          cout << pdg_ids[i] << ", ";
        }
        cout << endl;
      }

      nt = new float[1000];
      pbarMass = GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::anti_proton).mass();
      pbarMass2 = pbarMass*pbarMass;
    }

    virtual ~ReadAntiProtonSteps() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:

    // Name of the VD and TVD StepPoint collections
    std::string  _psVacuumStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    TNtuple* _ntAntiProtonSteps;

    float *nt; // Need this buffer to fill TTree ntpsVacuum

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    int _diagLevel;
  };

  void ReadAntiProtonSteps::beginJob(){

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntAntiProtonSteps = tfs->make<TNtuple>( "ntpbars", "AntiProtonSteps ntuple",
                                             "run:evt:trk:pdg:time:x:y:z:gtime:initialPbarCosTheta:momentum:initialProtonMomentum:currentKE:px:py:pz");
  }

  void ReadAntiProtonSteps::beginRun(art::Run const& run){

  }

  void ReadAntiProtonSteps::analyze(const art::Event& event) {

    ++_nAnalyzed;

    if (_diagLevel > 0)
      {
        std::cout << " \n \n \n hi, nAnalyzed = " << _nAnalyzed << std::endl;
      }


    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_psVacuumStepPoints,hits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());

    art::Handle<GenParticleCollection> genParticleHandle;
    event.getByLabel(_generatorModuleLabel, genParticleHandle);
    bool haveGenPart = genParticleHandle.isValid();
    GenParticleCollection const& genParticles(*genParticleHandle);
    if ( haveGenPart ) haveGenPart = !(genParticles.empty());
    CLHEP::HepLorentzVector initialProtonFourMomentum(0.,0.,0.,0.);
    //
    // look at the genParticles and see that the geantino with the initial proton is still there.
    if (_diagLevel > 1){
      std::cout << "found gen particle: " << haveGenPart << std::endl;
    }

    // in new style, we don't use the gen particle!
    if (haveGenPart)
      {
        for (auto iGen : genParticles )
          {
            if (_diagLevel > 0)
              {
                std::cout << " particle id " << iGen.pdgId() << " particle momentum " << iGen.momentum() << " position " << iGen.position() << std::endl;
              }
            if (iGen.pdgId() == PDGCode::proton)
              {
                //                initialProtonFourMomentum = iGen.momentum();
              }
          }
      }

    if (_diagLevel > 0){
      std::cout << "about to print size of hits" << std::endl;
      std::cout << " size of hits " << hits.isValid() << " " << hits->size()  << std::endl;
    }

    // Loop over all hits.
          if( hits.isValid() )
    //    if (hits.isValid() && hits->size() ==1)
      {
        for ( size_t i=0; i< hits->size(); ++i ){
          // Alias, used for readability.
          const StepPointMC& hit = (*hits)[i];

          // Get the hit information.
          const CLHEP::Hep3Vector& pos = hit.position();
          const CLHEP::Hep3Vector& mom = hit.momentum();

          // Get track info
          key_type trackId = hit.trackId();
          int pdgId = 0;
          bool goodPDG = false;
          CLHEP::HepLorentzVector startingFourMomentum(0.,0.,0.,0.);
          CLHEP::Hep3Vector startingPosition(0.,0.,0.);
          double currentKE(0.);
          if ( haveSimPart ){
            if( !simParticles->has(trackId) ) {
              pdgId = 0;
            } else {
              SimParticle const& sim = simParticles->at(trackId);
              pdgId = sim.pdgId();
              for (auto iPDG : pdg_save)
                {
                  if (pdgId == iPDG)
                    {
                      goodPDG = true;
                    }
                  if (goodPDG)
                    {
                      //
                      // what was the starting momentum and direction of the track?
                      // just use angle wrt z as a close-enough estimate
                      auto originalParticle = sim.originParticle();
                      startingFourMomentum = originalParticle.startMomentum();
                      startingPosition     = originalParticle.startPosition();
                      currentKE = sqrt( mom.mag()*mom.mag() + pbarMass2) - pbarMass;
                      initialProtonFourMomentum = (sim.parent())->endMomentum();
                      if (_diagLevel > 0)
                        {
                          std::cout << "\n inside hit number: " << i << std::endl;

                          std::cout << "starting momentum = " << startingFourMomentum << std::endl;
                          std::cout << "starting cos theta = " << startingFourMomentum.cosTheta() << std::endl;
                          std::cout << "starting position = " << startingPosition << std::endl;

                          std::cout << "_nAnalyzed, pdg, propertime, time, volume " << _nAnalyzed << " " << pdgId << " " << hit.properTime() << " " << hit.time() << " " << hit.volumeId() <<"\n"
                                    << std::setprecision(15) << " position" <<   " " << pos.x() << " " << pos.y() << " " << pos.z()
                                    << " \n current KE " << currentKE << std::endl;
                        }

                    }
                }
            }
          }

          if (goodPDG)//hack to get rid of genParticles
            {

              // Fill the ntuple.
              nt[0]  = event.id().run();
              nt[1]  = event.id().event();
              nt[2]  = trackId.asInt();
              nt[3]  = pdgId;
              nt[4]  = hit.time();
              nt[5]  = pos.x();
              nt[6]  = pos.y();
              nt[7]  = pos.z();
              nt[8] = hit.properTime();
              //              nt[9] = startingFourMomentum.cosTheta();
              //
              // get angle of initial proton to pbar; convert from HepDouble whatever that is
              nt[9] =    (initialProtonFourMomentum.vect().dot(startingFourMomentum.vect()))
                /(initialProtonFourMomentum.vect().mag()*startingFourMomentum.vect().mag());
              nt[10] = startingFourMomentum.vect().mag();
              nt[11] = initialProtonFourMomentum.vect().mag();
              nt[12] = currentKE;
              nt[13] = mom.x();
              nt[14] = mom.y();
              nt[15] = mom.z();

              if (_diagLevel > 0)
                {

                  std::cout << " \n filling ntuple" << std::endl;
                  std::cout << "_nAnalyzed, pdg, time, volume " << _nAnalyzed << " " << pdgId << " " << hit.time()  << " " << hit.volumeId() << "\n"
                            << std::setprecision(15) << " position" <<   " " << pos.x() << " " << pos.y() << " " << pos.z()
                            << " \n current KE " << currentKE << "\n  momentum vector = " << mom << std::endl;
                }
              _ntAntiProtonSteps->Fill(nt);
              if ( _nAnalyzed < _maxPrint){
                cout << "VD hit: "
                     << event.id().run()   << " | "
                     << event.id().event() << " | "
                     << hit.volumeId()     << " "
                     << pdgId              << " | "
                     << hit.time()         << " "
                     << mom.mag()
                     << endl;
              }
            }

        } // end loop over hits.
      }
  }

}  // end namespace mu2e

using mu2e::ReadAntiProtonSteps;
DEFINE_ART_MODULE(ReadAntiProtonSteps)

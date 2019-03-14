//
// Plugin to read sensitive detectors data and create ntuples
//
// 
#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>
 
using namespace std;
 
using CLHEP::Hep3Vector;
using CLHEP::keV;
 
namespace mu2e {
    const unsigned int ntsdet = 20; // maximum number of time VD
    class ReadEnergyDeposit : public art::EDAnalyzer {
 
        typedef vector<int> Vint;
        typedef SimParticleCollection::key_type key_type;
 
        // Name of the VD and TVD StepPoint collections
 
        art::InputTag _sdInputTag;
        art::InputTag _tsdInputTag;
        art::InputTag _simpInputTag;
        art::InputTag _generatorInputTag;
        art::InputTag _physInputTag;
 
        // Control printed output.
        int _nAnalyzed;
        int _maxPrint;
        int _debugout;
 
        TNtuple* _ntsd;
        TNtuple* _nttsd;
        TTree* _ntpart;
        TTree* _ntpart1;
        TNtuple* _ntsdext;
 
        float *nt; // Need this buffer to fill TTree ntsd
        float *ntext; // Need this buffer to fill TTree ntsdext
 
        bool write_ntsd;
        bool write_nttsd;
        bool write_ntpart;
        bool write_ntpart1;
        bool write_ntsdext;
 
        // Pointers to the physical volumes we are interested in
        // -- stopping target
        map<int,int> vid_stop;
 
        // List of particles of interest for the particles ntuple
        set<int> pdg_save;
 
        // List of particles to drop from time VD
        set<int> tsd_drop_pdg;
 
        // List of virtual detectors to be saved
        set<int> sd_save;
 
        // Module label of the g4 module that made the hits.
        std::string _g4ModuleLabel;
 
        // Virtual detector, which has to be crossed by particle before
        // it is saved in particles ntuple
        int _sd_required;
 
        // Save in the particles ntuple only those particles, which die 
        // after this time (in ns)
        double _timeCut;
 
        // Save in the particles ntuple only particle with momentum larger than this
        double _minMomentum;
 
        // Save only stopped particles in the particles ntuple
        bool _stopped_only;
 
        // Save all particles
        bool _save_all_pdg;
 
        // Should we add together proper time for the whole decay chain
        bool _add_proper_time;
 
        // If we are analyzing output of the staged simulation, look for 
        // real parent, navigating through the staged SimParticles
        bool _navigate_to_parent;
 
        public:
 
        explicit ReadEnergyDeposit(fhicl::ParameterSet const& pset);
        virtual ~ReadEnergyDeposit() { }
 
        virtual void beginJob();
        virtual void beginRun(art::Run const&);
 
        void analyze(const art::Event& e);
 
    };
 
    ReadEnergyDeposit::ReadEnergyDeposit(fhicl::ParameterSet const& pset) :
        art::EDAnalyzer(pset),
        _nAnalyzed(0),
        _maxPrint(pset.get<int>("maxPrint",0)),
        _debugout(pset.get<int>("debugOutput",0)),
        _ntsd(0), _nttsd(0), _ntpart(0), _ntpart1(0), _ntsdext(0),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")), // obsolete, left for backward compatibility
        _sd_required(pset.get<int>("requireVD",0)),
        _timeCut(pset.get<double>("timeCut",0.0)),
        _minMomentum(pset.get<double>("minMomentum",-1.0)),
        _stopped_only(pset.get<bool>("saveStopped",false)),
        _save_all_pdg(pset.get<bool>("saveAllPDG",false)),
        _add_proper_time(pset.get<bool>("addProperTime",false)),
        _navigate_to_parent(pset.get<bool>("navigateToParent",true))
    {
        _sdInputTag = pset.get<std::string>("sdStepPoints","g4run:virtualdetector");
        _tsdInputTag = pset.get<string>("tsdStepPoints","g4run:timeVD");
        _simpInputTag = pset.get<string>("simParticleColl","g4run");
        _generatorInputTag = pset.get<std::string>("generatorModuleLabel", "generate");
        _physInputTag = pset.get<std::string>("physicsVolumeColl", "g4run");
 
        write_ntsd    = pset.get<bool>("writeNTVD",true); 
        write_nttsd   = pset.get<bool>("writeNTTVD",true);
        write_ntpart  = pset.get<bool>("writeNTPART",true);
        write_ntpart1 = pset.get<bool>("writeNTPART1",true);
        write_ntsdext = pset.get<bool>("writeNTVDEXT",false); 
 
        if( _debugout>0 ) cout << "ReadEnergyDeposit: fill ntuples "
            << " NTVD=" << write_ntsd
                << " NTTVD=" << write_nttsd
                << " NTPART=" << write_ntpart
                << " NTPART1=" << write_ntpart1
                << " NTVDEXT=" << write_ntsdext
                << endl;
 
        Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
        if( pdg_ids.size()>0 ) {
//          cout << "ReadEnergyDeposit: save following particle types in the ntuple: ";
            for( size_t i=0; i<pdg_ids.size(); ++i ) {
                pdg_save.insert(pdg_ids[i]);
                //cout << pdg_ids[i] << ", ";
            }
            //cout << endl;
        }
 
        Vint const & tsd_drop_ids = pset.get<Vint>("tsdDropPDG", Vint());
        if( tsd_drop_ids.size()>0 ) {
            //cout << "ReadEnergyDeposit: drop following particle types from time VD ntuple: ";
            for( size_t i=0; i<tsd_drop_ids.size(); ++i ) {
                tsd_drop_pdg.insert(tsd_drop_ids[i]);
        //      cout << tsd_drop_ids[i] << ", ";
            }
            //cout << endl;
        }
 
        Vint const & sd_ids = pset.get<Vint>("saveVD", Vint());
        if( sd_ids.size()>0 ) {
        //  cout << "ReadEnergyDeposit: save data from the following virtual detectors: ";
            for( size_t i=0; i<sd_ids.size(); ++i ) {
                sd_save.insert(sd_ids[i]);
                //cout << sd_ids[i] << ", ";
            }
        //  cout << endl;
        }
 
        nt    = new float[1000];
        ntext = new float[1000];
 
    }
 
    void ReadEnergyDeposit::beginJob(){
 
        vid_stop.clear();
 
        // Get access to the TFile service.
 
        art::ServiceHandle<art::TFileService> tfs;
        if (write_ntsdext){
            _ntsdext = tfs->make<TNtuple>( "ntsdext", "Virtual Detectors ntuple (extended)",
                    "run:subrun:evt:trk:sdid:pdg:time:gtime:ke:"
                    "x:y:z:px:py:pz:"
                    "xl:yl:zl:pxl:pyl:pzl:"
                    "code:creation_code:parent_pdg:"
                    "originx:originy:originz:totE:ionE:noiE:e:p");
        }
        // Have to use TTree here, because one cannot use more than 100 variables in TNtuple
    }
 
    void ReadEnergyDeposit::beginRun(art::Run const& run){
 
        // Get pointers to the physical volumes we are interested
        art::Handle<PhysicalVolumeInfoCollection> physVolumes;
        run.getByLabel(_physInputTag, physVolumes);
        if( physVolumes.isValid() ) {
 
            for ( size_t i=0; i<physVolumes->size(); ++i ) {
                if( (*physVolumes)[i].name().compare(0,11,"TargetFoil_") == 0 ) {
 
                    vid_stop[i] = (*physVolumes)[i].copyNo();
                    //cout << "ReadEnergyDeposit: register stopping target volume " << i << " = "
                        //<< (*physVolumes)[i].name() << " " << (*physVolumes)[i].copyNo() << endl;
                }
            }
 
        }
 
    }
 
    void ReadEnergyDeposit::analyze(const art::Event& event) {
        ++_nAnalyzed;
        // Access virtual detectors geometry information
        // If not virtual detectors are defined, skip the rest
        GlobalConstantsHandle<ParticleDataTable> pdt;
 
        // Ask the event to give us a "handle" to the requested hits.
        art::Handle<StepPointMCCollection> hits;
        event.getByLabel(_sdInputTag,hits);
        if( _debugout>0 ) {
            cout << "ReadEnergyDeposit: hits collection is " << hits.isValid();
            if( hits.isValid() ) cout << " size=" << hits->size();
            cout << endl;
        }
        art::Handle<SimParticleCollection> simParticles;
        event.getByLabel(_simpInputTag, simParticles);
        bool haveSimPart = simParticles.isValid();
        if ( haveSimPart ) haveSimPart = !(simParticles->empty());
        if( _debugout>0 ) {
            cout << "ReadEnergyDeposit: simpart collection is " << simParticles.isValid();
            if( simParticles.isValid() ) cout << " size=" << simParticles->size();
            cout << endl;
        }
 
        art::Handle<G4BeamlineInfoCollection> g4beamlineData;
        event.getByLabel(_generatorInputTag, g4beamlineData);
        bool haveG4BL = g4beamlineData.isValid();
        if ( haveG4BL ) haveG4BL = (g4beamlineData->size()==1);
        if( _debugout>0 ) {
            cout << "ReadEnergyDeposit: beamline collection is " << g4beamlineData.isValid();
            if( g4beamlineData.isValid() ) cout << " size=" << g4beamlineData->size();
            cout << endl;
        }
 
        // Fill virtual detectors ntuple
 
        if( haveG4BL ) {
            G4BeamlineInfo const& extra = g4beamlineData->at(0);
            nt[18] = extra.weight();
            nt[19] = extra.time();
        } else {
            nt[18] = 0;
            nt[19] = 0;
        }
//cout<<"hits size000000000000000000000000000-----"<<hits->size()<<endl;
        // Loop over all VD hits.
        if( hits.isValid() ) for ( size_t i=0; i<hits->size(); ++i ){
            // Alias, used for readability.
            const StepPointMC& hit = (*hits)[i];
            // Get the hit information.
            // If virtual detector id is not in the list - skip it
            const CLHEP::Hep3Vector& pos = hit.position();
            const CLHEP::Hep3Vector& mom = hit.momentum();
            // Get track info
            key_type trackId = hit.trackId();
            int pdgId = 0;
            double mass(0.0);
            if ( haveSimPart ){
                if( !simParticles->has(trackId) ) {
                    pdgId = 0;
                } else {
                    SimParticle const& sim = simParticles->at(trackId);
                    pdgId = sim.pdgId();
                    mass = pdt->particle(pdgId).ref().mass();
                    // Fill the ntuple.
                    int parent_pdg = 0;
                    CLHEP::Hep3Vector origin(0.0,0.0,0.0);
                    int this_creation_code = sim.creationCode();
                    if (sim.isPrimary()){ //creation code 56=mu2ePrimary, means it was the first particle of this stage
                        //loop back through all parents until we find one that wasn't a mu2ePrimary
                        SimParticle const* sim_parent = &sim;
                        while( sim_parent && sim_parent->hasParent() ) {
                            sim_parent = simParticles->getOrNull(sim_parent->parentId());
                            if( sim_parent && sim_parent->endDefined() ) {
                                if (!sim_parent->isPrimary()){
                                    //when creation code is 56, the parent is actually the same particle.
                                    //so it's real parent would be the grandparent in this hierarchy
                                    SimParticle const* sim_grandparent = simParticles->getOrNull(sim_parent->parentId());
                                    if (sim_grandparent && sim_grandparent->endDefined()) parent_pdg = sim_grandparent->pdgId();
                                    //else throw cet::exception("ReadEnergyDeposit")<< " (Grand)parent not defined. \n";
 
                                    origin = sim_parent->startPosition();
                                    this_creation_code = sim_parent->creationCode();
                                    break;//we're just looking for the first one that's not a primary
                                }
                            }
                        }//end while
                    } else {
                        SimParticle const* sim_parent = simParticles->getOrNull(sim.parentId());
                        if( sim_parent && sim_parent->endDefined() ){
                            parent_pdg = sim_parent->pdgId();
                        }
                        origin = sim.startPosition();
                        this_creation_code = sim.creationCode();
                    }
                    // Fill the extended ntuple.
                    ntext[0]  = event.id().run(); 
                    ntext[1]  = event.id().subRun(); 
                    ntext[2]  = event.id().event(); 
                    ntext[3]  = trackId.asInt();
                    ntext[4]  = hit.volumeId(); 
                    ntext[5]  = pdgId;     
                    ntext[6]  = hit.time();
                    ntext[7]  = hit.properTime();
                    ntext[8]  = sqrt(mom.mag2()+mass*mass)-mass; // compute kinetic energy: this is what Geant cuts on
                    ntext[9]  = pos.x();
                    ntext[10] = pos.y();
                    ntext[11] = pos.z();
                    ntext[12] = mom.x();
                    ntext[13] = mom.y();
                    ntext[14] = mom.z(); 
                    ntext[15] = pos.x();
                    ntext[16] = pos.y();
                    ntext[17] = pos.z(); 
                    ntext[18] = mom.x();
                    ntext[19] = mom.y();
                    ntext[20] = mom.z();
                    ntext[21] = sim.creationCode();//geant4 creation code (might be mu2ePrimary i.e. not physics)
                    ntext[22] = this_creation_code;//This is supposed to be only physics creation codes, so trace back to parent when mu2ePrimary
                    ntext[23] = parent_pdg;
                    ntext[24] = origin.x();
                    ntext[25] = origin.y();
                    ntext[26] = origin.z();
                    ntext[27]=hit.totalEDep();
                    ntext[28]=hit.ionizingEdep();
                    ntext[29]=hit.nonIonizingEDep();
                    ntext[30]=sqrt(mom.mag2()+mass*mass); // compute kinetic energy: this is what Geant cuts on
                    ntext[31]=sqrt(mom.mag2()); // compute kinetic energy: this is what Geant cuts on
//                  cout<<hit.totalEDep()<<endl;
                    if (write_ntsdext){
                        _ntsdext->Fill(ntext);
                    }
                }
            }
 
            if ( _nAnalyzed < _maxPrint){
                cout << "VD hit: "
                    << event.id().run()   << " | "
                    << event.id().subRun()<< " | "
                    << event.id().event() << " | "
                    << hit.volumeId()     << " "
                    << pdgId              << " | "
                    << hit.time()         << " "
                    //<< sim.startGlobalTime()     << " "
                    << mom.mag()
                    << endl;
 
            }
        } // end loop over hits.
    }
}  // end namespace mu2e
 
using mu2e::ReadEnergyDeposit;
DEFINE_ART_MODULE(ReadEnergyDeposit);

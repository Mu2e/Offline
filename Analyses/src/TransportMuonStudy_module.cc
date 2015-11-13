//
// A module to study background rates in the detector subsystems.
//
// $Id: TransportMuonStudy_module.cc,v 1.1 2014/01/15 17:21:19 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/01/15 17:21:19 $
//
// Original author Giovanni Tassielli
//

#include "CLHEP/Units/PhysicalConstants.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StatusG4.hh"
#include "GeometryService/inc/VirtualDetector.hh"
//#include <cmath>
#include <iostream>
//#include <set>
#include <memory>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

namespace mu2e {

class TransportMuonStudy : public art::EDAnalyzer {
public:

        explicit TransportMuonStudy(fhicl::ParameterSet const& pset):
        art::EDAnalyzer(pset),
        _diagLevel(pset.get<int>("diagLevel",0)),
        //_trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
        //_swiresStepPoints(pset.get<string>("swiresStepPoints","trackerSWires")),
        //_fwiresStepPoints(pset.get<string>("fwiresStepPoints","itrackerFWires")),
        //_makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _primaryMu(pset.get<bool>("runFromMuonsFile", false)),
        //_extractElectronsData(pset.get<string>("elextractModuleLabel")),
        //_minimumEnergyTracker(pset.get<double>("minimumEnergyTracker",0.0001)), // MeV
        _vdStepPoints(pset.get<string>("vdStepPoints","virtualdetector")),
        _partDatafile(pset.get<std::string>("partDatafile", "muonData.txt")),
        _nAnalyzed(0),
        _tNtup(0),
        _nBadG4Status(0),
        _nOverflow(0),
        _nKilled(0),
        _totalcputime(0),
        _totalrealtime(0)
        {
                cout << "Module TransportMuonStudy is starting" << endl;
        }
        virtual ~TransportMuonStudy() {
        }
        virtual void beginJob();
        virtual void endJob();

        void analyze(art::Event const& e );

        fstream* wirestxt;

private:

        void doMuons(art::Event const& evt);

        // Diagnostic level
        int _diagLevel;

        // Label of the module that made the hits.
        //std::string _makerModuleLabel;

        // Label of the generator.
        std::string _generatorModuleLabel;

        // Label of the G4 module
        std::string _g4ModuleLabel;

        bool _primaryMu;

        std::string _vdStepPoints;

        std::string _partDatafile;

        //number of analyzed events
        int _nAnalyzed;

        TNtuple* _tNtup;

        int _nBadG4Status, _nOverflow, _nKilled;
        float _totalcputime, _totalrealtime;
        int _nVaribles;

        std::vector<VirtualDetectorId::enum_type> _storeAdditionalVDinfo;

};

void TransportMuonStudy::beginJob( ) {
}

void TransportMuonStudy::analyze(art::Event const& evt ) {

        ++_nAnalyzed;

        //*****test code******
        static int ncalls(0);
        ++ncalls;

        art::Handle<StatusG4> g4StatusHandle;
        evt.getByLabel( _g4ModuleLabel, g4StatusHandle);
        StatusG4 const& g4Status = *g4StatusHandle;

        if ( g4Status.status() > 1 ) {
                ++_nBadG4Status;
                mf::LogError("G4")
                << "Aborting TransportMuonStudy::analyze due to G4 status\n"
                << g4Status;
                return;
        }

        if (g4Status.overflowSimParticles()) {
                ++_nOverflow;
                mf::LogError("G4")
                << "Aborting TransportMuonStudy::analyze due to overflow of particles\n"
                << g4Status;
                return;
        }

        if (g4Status.nKilledStepLimit() > 0) {
                ++_nKilled;
                mf::LogError("G4")
                << "Aborting TransportMuonStudy::analyze due to nkilledStepLimit reached\n"
                << g4Status;
                return;
        }

        _totalcputime += g4Status.cpuTime();
        _totalrealtime += g4Status.realTime();

        art::ServiceHandle<GeometryService> geom;

        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;

                if (!_primaryMu) { wirestxt = new fstream(_partDatafile.c_str(),ios::out); }

                art::ServiceHandle<art::TFileService> tfs;

                std::string generalInfo ="run:subrun:evt:pdg";
                std::string startInfo =":stX:stY:stZ:stPx:stPy:stPz:stTime";
                std::string endInfo   =":enX:enY:enZ:enPx:enPy:enPz:enTime";
                std::string c5OutInfo   =":c5X:c5Y:c5Z:c5Px:c5Py:c5Pz:c5Time";
                std::string c5InInfo   =":c5InX:c5InY:c5InZ:c5InPx:c5InPy:c5InPz:c5InTime";
                std::string c3OutInfo   =":c3OutX:c3OutY:c3OutZ:c3OutPx:c3OutPy:c3OutPz:c3OutTime";
                std::string c3InInfo   =":c3InX:c3InY:c3InZ:c3InPx:c3InPy:c3InPz:c3InTime";
                std::string c1OutInfo   =":c1OutX:c1OutY:c1OutZ:c1OutPx:c1OutPy:c1OutPz:c1OutTime";
                std::string c1InInfo   =":c1InX:c1InY:c1InZ:c1InPx:c1InPy:c1InPz:c1InTime";
                _nVaribles = std::count(generalInfo.begin(), generalInfo.end(), ':')+1;
                _nVaribles += std::count(startInfo.begin(), startInfo.end(), ':');
                _nVaribles += std::count(endInfo.begin(), endInfo.end(), ':');
                _nVaribles += std::count(c5OutInfo.begin(), c5OutInfo.end(), ':');
                _nVaribles += std::count(c5InInfo.begin(), c5InInfo.end(), ':');
                _nVaribles += std::count(c3OutInfo.begin(), c3OutInfo.end(), ':');
                _nVaribles += std::count(c3InInfo.begin(), c3InInfo.end(), ':');
                _nVaribles += std::count(c1OutInfo.begin(), c1OutInfo.end(), ':');
                _nVaribles += std::count(c1InInfo.begin(), c1InInfo.end(), ':');
                _tNtup        = tfs->make<TNtuple>( "ParticleData", "Particle Data", (generalInfo+startInfo+endInfo+c5OutInfo+c5InInfo+c3OutInfo+c3InInfo+c1OutInfo+c1InInfo).c_str() );
                //list of VD for which we want to save the information after the c5Out (Coll5_Out).
                _storeAdditionalVDinfo.push_back(VirtualDetectorId::Coll5_In);
                _storeAdditionalVDinfo.push_back(VirtualDetectorId::Coll32_Out);
                _storeAdditionalVDinfo.push_back(VirtualDetectorId::Coll31_In);
                _storeAdditionalVDinfo.push_back(VirtualDetectorId::Coll1_Out);
                _storeAdditionalVDinfo.push_back(VirtualDetectorId::Coll1_In);

        }

        doMuons(evt);

} // end of analyze

void TransportMuonStudy::endJob() {

        if (!_primaryMu) {
                wirestxt->close();
                delete wirestxt;
        }

        cout << "TransportMuonStudy::endJob Number of events skipped "
             << "due to G4 completion status: "
             << _nBadG4Status
             << "\nTransportMuonStudy::endJob Number of overflow events "
             << "due to too many particles in G4: "
             << _nOverflow
             << "\nTransportMuonStudy::endJob Number of events with killed particles "
             << "due to too many steps in G4: "
             << _nKilled
             << "\nTransportMuonStudy::endJob total CpuTime "
             << _totalcputime
             << "\nTransportMuonStudy::endJob total RealTime "
             << _totalrealtime
             << endl;
}

void TransportMuonStudy::doMuons(art::Event const& evt) {

        // Get handles to the generated and simulated particles.
        art::Handle<GenParticleCollection> genParticles;
        evt.getByLabel(_generatorModuleLabel, genParticles);

        art::Handle<G4BeamlineInfoCollection> g4beamlineData;
        evt.getByLabel(_generatorModuleLabel, g4beamlineData);

        art::Handle<SimParticleCollection> simParticles;
        evt.getByLabel(_g4ModuleLabel, simParticles);

        // Handle to information about G4 physical volumes.
        art::Handle<PhysicalVolumeInfoCollection> volumes;
        evt.getRun().getByLabel(_g4ModuleLabel, volumes);

        //Handle to VD steps
        art::Handle<StepPointMCCollection> vdHits;
        evt.getByLabel(_g4ModuleLabel,_vdStepPoints,vdHits);


        // Some files might not have the SimParticle and volume information.
        bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

        // Other files might have empty collections.
        if ( haveSimPart ){
                haveSimPart = !(simParticles->empty() || volumes->empty());
        }

        float *tntpArray = new float [_nVaribles];
        for (int iv=0; iv<_nVaribles; ++iv) {
                tntpArray[iv]=0.0;
        }

        // cout << "event " << evt.id().event() << ": first fills" << endl;

        int idx(0);
        tntpArray[idx++] = evt.run();
        tntpArray[idx++] = evt.subRun();
        if (_primaryMu) {
                if (g4beamlineData.isValid()) {
                        G4BeamlineInfo const& extra = g4beamlineData->back();
                        tntpArray[idx++] = extra.eventId();
                } else {
                        tntpArray[idx++] = evt.id().event();
                }
        } else {
                tntpArray[idx++] = evt.id().event();
        }

        //Loop over vd hits
        // cout << "event " << evt.id().event() << ": looking at VD steps" << endl;
        bool checkMuon(false);
        if (vdHits.isValid()) {
                for (size_t i=0; i<vdHits->size(); ++i) {
                        const StepPointMC& hit = (*vdHits)[i];
                        SimParticle const* tmpPart = simParticles->getOrNull(hit.trackId());
                        if (tmpPart==NULL) { continue; }
                        checkMuon = false;
                        if (_primaryMu && tmpPart->isPrimary()) { checkMuon = true; }
                        if (!_primaryMu && tmpPart->isSecondary()) { checkMuon = true; }
                        if (checkMuon) {
                                int id = hit.volumeId();
                                if (id == VirtualDetectorId::ST_In) {
                                        if (tmpPart->pdgId()==PDGCode::mu_minus
                                                        && tmpPart->endGlobalTime()<1000
                                        ) {
                                                bool foundInC5(false);
                                                tntpArray[idx++] = tmpPart->pdgId();
                                                tntpArray[idx++] = tmpPart->startPosition().x()-3904.0;
                                                tntpArray[idx++] = tmpPart->startPosition().y();
                                                tntpArray[idx++] = tmpPart->startPosition().z()+7929.0;
                                                tntpArray[idx++] = tmpPart->startMomentum().x();
                                                tntpArray[idx++] = tmpPart->startMomentum().y();
                                                tntpArray[idx++] = tmpPart->startMomentum().z();
                                                tntpArray[idx++] = tmpPart->startGlobalTime();
                                                tntpArray[idx++] = hit.position().x()-3904.0; //tmpPart->endPosition().x()-3904.0;
                                                tntpArray[idx++] = hit.position().y();        //tmpPart->endPosition().y();
                                                tntpArray[idx++] = hit.position().z()+7929.0; //tmpPart->endPosition().z()+7929.0;
                                                tntpArray[idx++] = hit.momentum().x();        //tmpPart->endMomentum().x();
                                                tntpArray[idx++] = hit.momentum().y();        //tmpPart->endMomentum().y();
                                                tntpArray[idx++] = hit.momentum().z();        //tmpPart->endMomentum().z();
                                                tntpArray[idx++] = hit.time();                //tmpPart->endGlobalTime();

                                                for (size_t ic=0; ic<vdHits->size(); ++ic) {
                                                        const StepPointMC& hitc = (*vdHits)[ic];
                                                        int idc = hitc.volumeId();
                                                        if (idc == VirtualDetectorId::Coll5_Out && hitc.trackId()==hit.trackId()) {
                                                                SimParticle const* tmpPartc = simParticles->getOrNull(hitc.trackId());
                                                                if (tmpPartc==NULL) { continue; }
                                                                tntpArray[idx++] = hitc.position().x()-3904.0;
                                                                tntpArray[idx++] = hitc.position().y();
                                                                tntpArray[idx++] = hitc.position().z()+7929.0;
                                                                tntpArray[idx++] = hitc.momentum().x();
                                                                tntpArray[idx++] = hitc.momentum().y();
                                                                tntpArray[idx++] = hitc.momentum().z();
                                                                tntpArray[idx++] = hitc.time();
                                                                foundInC5=true;
                                                                break;
                                                        }
                                                }

                                                //store additional Coll_xx VDs informations
                                                for ( std::vector<VirtualDetectorId::enum_type>::iterator vd_it=_storeAdditionalVDinfo.begin(); vd_it!=_storeAdditionalVDinfo.end(); ++vd_it ) {
                                                        for (size_t ic=0; ic<vdHits->size(); ++ic) {
                                                                const StepPointMC& hitc = (*vdHits)[ic];
                                                                int idc = hitc.volumeId();
                                                                if (idc == *vd_it && hitc.trackId()==hit.trackId()) {
                                                                        SimParticle const* tmpPartc = simParticles->getOrNull(hitc.trackId());
                                                                        if (tmpPartc==NULL) { continue; }
                                                                        tntpArray[idx++] = hitc.position().x()-3904.0;
                                                                        tntpArray[idx++] = hitc.position().y();
                                                                        tntpArray[idx++] = hitc.position().z()+7929.0;
                                                                        tntpArray[idx++] = hitc.momentum().x();
                                                                        tntpArray[idx++] = hitc.momentum().y();
                                                                        tntpArray[idx++] = hitc.momentum().z();
                                                                        tntpArray[idx++] = hitc.time();
                                                                        break;
                                                                }
                                                        }
                                                }
                                                //end store additional Coll_xx VDs informations

                                                if (foundInC5) {
                                                        if (!_primaryMu) {
                                                                (*wirestxt)<<tmpPart->startPosition().x()-3904.0<<" "
                                                                                <<tmpPart->startPosition().y()<<" "
                                                                                <<tmpPart->startPosition().z()+7929.0<<" "
                                                                                <<tmpPart->startMomentum().x()<<" "
                                                                                <<tmpPart->startMomentum().y()<<" "
                                                                                <<tmpPart->startMomentum().z()<<" "
                                                                                <<tmpPart->startGlobalTime()<<" "
                                                                                <<PDGCode::mu_minus<<" "
                                                                                <<(evt.subRun()+1)*100000+evt.id().event()<<" "
                                                                                <<1<<" "
                                                                                //<<PDGCode::mu_minus
                                                                                //<<1<<std::endl;
                                                                                <<tmpPart->parent()->pdgId()<<" "
                                                                                <<tmpPart->weight()<<std::endl;
                                                        }
                                                        _tNtup->Fill(tntpArray);
                                                        break;
                                                }
                                        }
                                }
                        }
                }
        }

        //    _tNtup->Fill(tntpArray);

        delete [] tntpArray;

} // end of doMuons

}

using mu2e::TransportMuonStudy;
DEFINE_ART_MODULE(TransportMuonStudy);


//
// performance a remapping module of the StrawHit in a manner that they can be accessed by Z and Sector IDs
//
// $Id: CreateZRotStrawHitMap_module.cc,v 1.1 2011/10/28 00:19:14 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/10/28 00:19:14 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <map>


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//#include "CLHEP/Vector/TwoVector.h"
//#include "CLHEP/Geometry/Point3D.h"
//#include "CLHEP/Units/SystemOfUnits.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
//#include "ITrackerGeom/inc/Cell.hh"
//#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "ITrackerGeom/inc/ITracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
//#include "MCDataProducts/inc/GenId.hh"
//#include "MCDataProducts/inc/VisibleGenElTrack.hh"
//#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroup.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroupCollection.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMap.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMapCollection.hh"


using namespace std;

namespace mu2e {

  typedef art::Ptr<StrawHit> StrawHitPtr;

  class CreateZRotStrawHitMap : public art::EDProducer {
  public:

    explicit CreateZRotStrawHitMap(fhicl::ParameterSet const& pset);
    virtual ~CreateZRotStrawHitMap() {
//            if (_fakeCanvas)        delete _fakeCanvas;
    }

    virtual void beginJob();
    //void endJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:

    // Start: run time parameters

//    // The module label of this module.
//    std::string _moduleLabel;
//
//    // Label of the G4 module
//    std::string _g4ModuleLabel;
//
//    // Name of the tracker StepPoint collection
//    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

//    // Label of the generator.
//    std::string _generatorModuleLabel;

    // Label of the module that made the hits.
    std::string _timeRejecterModuleLabel;

    // Label of the module that made the hits.
    std::string _geomRejecterModuleLabel;

    // End: run time parameters

    unsigned int  runID, eventID, evNHit;

    inline unsigned int iRot(int device, int sector){
            return (unsigned int) (device%2)+2*sector;
    }

    inline unsigned int iZpos(int device, int sector){
            return (unsigned int) 2*device+(sector%2);
    }

  };

  CreateZRotStrawHitMap::CreateZRotStrawHitMap(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
    _timeRejecterModuleLabel(pset.get<string>("tRejecterModuleLabel")),
    _geomRejecterModuleLabel(pset.get<string>("gRejecterModuleLabel"))
//    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
//    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
//    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),

//    _fakeCanvas(0),
  {
          runID=eventID=evNHit=0;

          // Tell the framework what we make.
          produces<ZRotStrawHitMapCollection>();
}

  void CreateZRotStrawHitMap::beginJob(){

	  cout<<"Starting CreateZRotStrawHitMap jos!"<<endl;

  }

  void CreateZRotStrawHitMap::produce(art::Event & event ) {


    auto_ptr<ZRotStrawHitMapCollection> zsctmap(new ZRotStrawHitMapCollection);

    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();
    double zSttDist = ttr.getStation(1).midZ()-ttr.getStation(0).midZ();
    //cout<<"******************** zSttDist "<<zSttDist<<endl;

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    //TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();

    art::Handle<SctrSttnClusterGroupCollection> gclusgtHandle;
    event.getByLabel(_geomRejecterModuleLabel,gclusgtHandle);
    SctrSttnClusterGroupCollection const* gclustgs = gclusgtHandle.product();


    cout<<"----------------------------------------------------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;

    StrawId sid;
    int stn, layern, devicen, sectorn;
    unsigned int absSect, absZpos;

    //double ptMeV, rho;
    //double sel_rho;
    //double B=1.0;
    //CLHEP::Hep2Vector radDir;
    //HepGeom::Point3D<double> CirCenter;

    runID=event.run();
    eventID=event.event();
    evNHit=hits->size();

    //---------------- for time algorithm ----------------
    //size_t nTimeClusPerEvent = tclusts->size();
    //nPeaksFound=nTimeClusPerEvent;


//
//            for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {
//                    TrackerHitTimeCluster const&  tclust(tclusts->at(ipeak));
//                    cout<<"\t Hit in peak "<<tclust._selectedTrackerHits.size()<<endl;
//
//                    elHitFound=0;
//                    std::vector<StrawHitPtr> tmpElHits;
//
//                    for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin(); iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
//                            // Access data
//
//                            //for ( int iElHit=0; iElHit<signaLEltrk->getNumOfHit(); iElHit++) {
//                            //        GenElHitData& genEl = signaLEltrk->getHit(iElHit);
//                            for ( std::vector<StrawHitPtr>::iterator iGoodSelHit_it=signElGoodHit.begin(); iGoodSelHit_it!=signElGoodHit.end(); ++iGoodSelHit_it ) {
//                                    if ( iTCHit->key()==iGoodSelHit_it->key() ) {
//                                            elHitFound++;
//                                            tmpElHits.push_back(*iTCHit);
//                                            break;
//                                    }
//                            }
//
//                            //                //                cout<<"\t\t iHit in peak at "<<*iTCHit<<endl;
//                            //                //StrawHit        const&      hit(hits->at(*iTCHit));
//                            //                StrawHit        const&      hit=*(*iTCHit);
//                            //                StrawIndex si = hit.strawIndex();
//                            //                const Straw & str = tracker.getStraw(si);
//                            //
//                            //                // cout << "Getting informations about cells" << endl;
//                            //
//                            //                sid = str.id();
//                    }
//
//                    if ( elHitFound>0 ) {
//                            foundElHitInPeak.insert( pair<unsigned int, size_t> (elHitFound,ipeak) );
//                            signElGoodHitsInTmpks.insert( pair< size_t, std::vector<StrawHitPtr> > (ipeak,tmpElHits) );
//                    }
//            }
//
//            if ( !foundElHitInPeak.empty() ) {
//                    _hNBestTimePeakEv->Fill(sel_ptMeV);
//                    BestFoundElInPeak_it=foundElHitInPeak.rbegin();
//                    bestTPcElNHit=BestFoundElInPeak_it->first;
//                    bestTPNHit=tclusts->at(BestFoundElInPeak_it->second)._selectedTrackerHits.size();
//                    _hNhitBestTimePeakEv->Fill(sel_ptMeV,bestTPcElNHit);
//                    _hLostHitBstTmPkEv->Fill(sel_ptMeV, (signElGoodHit.size() - bestTPcElNHit) );
//                    _hNoiseHitBstTmPkEv->Fill(sel_ptMeV, (bestTPNHit - bestTPcElNHit) );
//            }

    //---------------- for geom algorithm ----------------

    cout<<"N of group of Ttracker cluster of Hit found that could be tracks: "<<gclustgs->size()<<endl;
    //nPotentTracks=gclustgs->size();

    SctrSttnClusterGroupCollection::const_iterator gclustgs_it;
    std::vector<SectorStationCluster>::const_iterator iclustg_it;
    int igroup=0, iclust;


    //int sClGrpTotalHit, iHitInCl, iCl, elHitFoundInCl;

    for ( gclustgs_it=gclustgs->begin(); gclustgs_it!=gclustgs->end(); ++gclustgs_it ) {
            cout<<"i-th group :"<<igroup<<" : "<<endl<<*gclustgs_it;
            iclust=0;
            zsctmap->push_back(ZRotStrawHitMap());
            ZRotStrawHitMap &tmpMap = zsctmap->back();
            tmpMap._min_Time = gclustgs_it->_relatedTimeCluster->_minHitTime;
            tmpMap._max_Time = gclustgs_it->_relatedTimeCluster->_maxHitTime;
            if (gclustgs_it->_coupling==SctrSttnClusterGroup::good) {
                    tmpMap._min_HStep = (gclustgs_it->_meanPitch - 3.0*gclustgs_it->_sigmaPitch)*zSttDist;
                    tmpMap._max_HStep = (gclustgs_it->_meanPitch + 3.0*gclustgs_it->_sigmaPitch)*zSttDist;

            }
            for ( iclustg_it=gclustgs_it->_selectedClusters.begin(); iclustg_it!=gclustgs_it->_selectedClusters.end(); ++iclustg_it ){
                    cout<<"i-th clust :"<<iclust<<" : "<<endl<<*iclustg_it;
                    for ( std::vector<StrawHitPtr>::const_iterator iSClGrpHit_it=iclustg_it->_selectedTrackerHits.begin(); iSClGrpHit_it!=iclustg_it->_selectedTrackerHits.end(); ++iSClGrpHit_it ) {
                            StrawHit        const&      hit=*(*iSClGrpHit_it);
                            StrawIndex si = hit.strawIndex();
                            const Straw & str = tracker.getStraw(si);

                            // cout << "Getting informations about cells" << endl;

                            sid     = str.id();
                            stn     = sid.getStraw();
                            layern  = sid.getLayer();
                            devicen = sid.getDevice();
                            sectorn = sid.getSector();

                            absSect = iRot(devicen,sectorn);
                            absZpos = iZpos(devicen,sectorn);

                            tmpMap.AddHit( absZpos, absSect, *iSClGrpHit_it );
                    }
                    iclust++;
            }

    }


    event.put(zsctmap);

//    cerr << "Double click in the canvas_Fake to continue:" ;
//    _fakeCanvas->cd();
//    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
//    printEvN->SetTextFont(62);
//    printEvN->SetTextSizePixels(180);
//    printEvN->Draw();
//    _fakeCanvas->WaitPrimitive();
//    cerr << endl;
//    delete printEvN;

    cout<<"--------------------------- End of Analysis --------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;


  } // end analyze

//  void CreateZRotStrawHitMap::endJob(){
//
//    // cd() to correct root directory. See note 3.
//    TDirectory* save = gDirectory;
//    _directory->cd();
//
//    // Write canvas.  See note 4.
////    _canvas->Write();
//
//    // cd() back to where we were.  See note 3.
//    save->cd();
//
//  }

}  // end namespace mu2e

using mu2e::CreateZRotStrawHitMap;
DEFINE_ART_MODULE(CreateZRotStrawHitMap);

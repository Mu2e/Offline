//
//
//  $Id: ITGasLayerSD.cc,v 1.9 2011/05/20 22:22:22 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/05/20 22:22:22 $
//
// 

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "ITrackerGeom/inc/ITracker.hh"
//#include "GeometryService/inc/GeometryService.hh"
//#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  G4ThreeVector ITGasLayerSD::_mu2eDetCenter;

  ITGasLayerSD::ITGasLayerSD(G4String name, const SimpleConfig& config) :
                  G4VSensitiveDetector(name),
                  _collection(0),
                  _debugList(0),
                  _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
                  _currentSize(0)


  {
          // Get list of events for which to make debug printout.
          string key("g4.itgaslayerSDEventList");
          if ( config.hasName(key) ){
                  vector<int> list;
                  config.getVectorInt(key,list);
                  _debugList.add(list);
          }

          art::ServiceHandle<GeometryService> geom;

          if ( !geom->hasElement<ITracker>() ) {
                  throw cet::exception("GEOM")
                  << "Expected I Trackers but found neither.\n";
          }


          //    _superlayer=atoi(name.substr(5,2).c_str());
          //    _ring=atoi(name.substr(8,2).c_str());

//          _ittype=itracker->geomType();

          //    _nwires=itracker->nSWire()+_superlayer*itracker->nSDeltaWire();
          //    try {
          //        itracker->getCellGeometryHandle()->SelectCell(_superlayer,_ring,0);
          //            _nwires=itracker->getCellGeometryHandle()->GetITLayer()->nCells();
          //            _Dphi=CLHEP::twopi/_nwires;
          //    }catch (cet::exception e) {
          //            cerr<<e;
          //            _nwires=0;
          //            _Dphi=0.0;
          //    }

          //    G4String HCname;
          //    HCname="StepPointG4Collection_";
          //    HCname+=name;
          //    collectionName.insert(HCname);
  }


  ITGasLayerSD::~ITGasLayerSD(){ }

  void ITGasLayerSD::Initialize(G4HCofThisEvent* HCE){

          _currentSize=0;


////          std::cout<<SensitiveDetectorName<<std::endl;
//          _collection = new StepPointG4Collection
//                          (SensitiveDetectorName,collectionName[0]);
//          HCE->AddHitsCollection( HCE->GetNumberOfCollections(), _collection );

  }


//  void ITGasLayerSD::EndOfEvent(G4HCofThisEvent*){
//
//    if (verboseLevel>0) {
//      G4int NbHits = _collection->entries();
//      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
//             << " hits in the Drift Chamber chambers: " << G4endl;
//      for (G4int i=0;i<NbHits;i++) (*_collection)[i]->Print();
//    }
//  }

  void ITGasLayerSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize
                            << " Drift Chamber hits were generated in the event."
                            << endl
                            << "Only " << _sizeLimit << " are saved in output collection."
                            << endl;
      cout << "Total of " << _currentSize
           << " Drift Chamber hits were generated in the event."
           << endl
           << "Only " << _sizeLimit << " are saved in output collection."
           << endl;
    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " hits in the Drift Chamber: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }
  }

//  G4ThreeVector ITGasLayerSD::GetTrackerOrigin(const G4TouchableHandle & touchableHandle) {
//
//    // how deep in the hierachy is the tracker a.k.a tracker depth
//    // (depends on the tracker version)
//    // for LTracker 1,2,3 it is 3
//    // the tracker version is set in the constructor
//
//    size_t td = 3;
//
////    art::ServiceHandle<GeometryService> geom;
////    if ( geom->hasElement<TTracker>() ) {
////      td =_TrackerVersion +1;
////    }
//
//    //    cout << "Debugging: tracker depth/version: " << td << "/" << _TrackerVersion << endl;
//
//    size_t hdepth = touchableHandle->GetHistoryDepth();
//
//    G4ThreeVector cdo;
//
//    for (size_t dd=0; dd!=hdepth; ++dd) {
//
//      if (dd>=td) cdo += touchableHandle->GetVolume(dd)->GetTranslation();
//
////       cout << "Debugging: det depth name copy#: " << dd << " " <<
////         touchableHandle->GetVolume(dd)->GetName() << " " <<
////         touchableHandle->GetVolume(dd)->GetCopyNo() << endl;
//
////       G4LogicalVolume* lvp = touchableHandle->GetVolume(dd)->GetLogicalVolume();
////       G4int nd = lvp->GetNoDaughters();
////       for (G4int d = 0;d!=nd; ++d) {
////         cout << "Debugging: daughter: " << lvp->GetDaughter(d)->GetName() << " " <<
////           lvp->GetDaughter(d)->GetCopyNo() << endl;
////       }
//
////       cout << "Debugging det origin: " << touchableHandle->GetVolume(dd)->GetTranslation() <<
////         " " << cdo << endl;
//
//    }
//
//    return cdo;
//
//  }

  void ITGasLayerSD::beforeG4Event(StepPointMCCollection& outputHits) {
    _collection = &outputHits;
    return;
  } // end of beforeG4Event

} //namespace mu2e

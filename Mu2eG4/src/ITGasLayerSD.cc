// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"

using namespace std;

namespace mu2e {

  ITGasLayerSD::ITGasLayerSD(G4String name) :G4VSensitiveDetector(name){

    _superlayer=atoi(name.substr(5,2).c_str());
    _ring=atoi(name.substr(8,2).c_str());
    GeomHandle<ITracker> itracker;
//    _nwires=itracker->nSWire()+_superlayer*itracker->nSDeltaWire();
    try {
        itracker->getCellGeometryHandle()->SelectCell(_superlayer,_ring,0);
    	_nwires=itracker->getCellGeometryHandle()->GetITLayer()->nCells();
    	_Dphi=M_2PI/_nwires;
    }catch (cms::Exception e) {
    	cerr<<e;
    	_nwires=0;
    	_Dphi=0.0;
    }
    G4String HCname;
    HCname="StepPointG4Collection_";
    HCname+=name;
    collectionName.insert(HCname);
  }


  ITGasLayerSD::~ITGasLayerSD(){ }

  void ITGasLayerSD::Initialize(G4HCofThisEvent* HCE){

//	  std::cout<<SensitiveDetectorName<<std::endl;
	  _collection = new StepPointG4Collection
			  (SensitiveDetectorName,collectionName[0]);
	  HCE->AddHitsCollection( HCE->GetNumberOfCollections(), _collection );

  }
  

  void ITGasLayerSD::EndOfEvent(G4HCofThisEvent*){

    if (verboseLevel>0) { 
      G4int NbHits = _collection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
	     << " hits in the Drift Chamber chambers: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i]->Print();
    } 
  }
  
} //namespace mu2e


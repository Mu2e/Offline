#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
//Geom:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/WorldG4Maker.hh"
#include "GeometryService/inc/TrackerMaker.hh"
#include "GeometryService/inc/Mu2eHallMaker.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
//Mu2e Tracker Geom:
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/Mu2eCoordTransform.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
//TEve
#include "TEveEventDisplay/src/dict_classes/Geom_Interface.h"

using namespace mu2e;
namespace mu2e{

	Geom_Interface::Geom_Interface(){}

  // Function to descend and remove nodes above the DS - run after HideBuilding
  void Geom_Interface::InsideDS( TGeoNode * node, bool inDSVac ){
    std::string _name = (node->GetVolume()->GetName());
    if ( node->GetMotherVolume() ) {
      std::string motherName(node->GetMotherVolume()->GetName());
      if ( motherName == "DS2Vacuum" || motherName == "DS3Vacuum" ){
        inDSVac = true;
      }
      }
      if ( inDSVac && _name.find("VirtualDetector_TT_Mid") != 0 ) {
        node->SetVisibility(kTRUE);
      } else{
        node->SetVisibility(kFALSE);
      }
      int ndau = node->GetNdaughters();
      for ( int i=0; i<ndau; ++i ){
        TGeoNode * dau = node->GetDaughter(i);
        InsideDS( dau, inDSVac );
      }
  }

  //Function allows user to specifically hide a node by its GDML mateiral name. See Fix.gdml for the name.
  void Geom_Interface::hideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff) {
    std::string material(node->GetVolume()->GetMaterial()->GetName());
    if ( material.find(mat) != std::string::npos ) node->SetVisibility(onOff);
    int ndau = node->GetNdaughters();
    for ( int i=0; i<ndau; ++i ){
      TGeoNode * dau = node->GetDaughter(i);
      hideNodesByMaterial( dau, mat, onOff);
    }
  }

  //Function to hide element by its name
  void Geom_Interface::hideNodesByName(TGeoNode* node, const std::string& str, bool onOff, int _diagLevel) {
    std::string name(node->GetName());
    if ( name.find(str) != std::string::npos ){
      node->SetVisibility(onOff);
      if(_diagLevel > 0) std::cout <<"hiding "<< name << std::endl;
    }
    int ndau = node->GetNdaughters();
    for ( int i=0; i<ndau; ++i ){
      TGeoNode * dau = node->GetDaughter(i);
      hideNodesByName( dau, str, onOff, _diagLevel);
    }
  }
  
  void Geom_Interface::showNodesByName(TGeoNode* node, const std::string& str, bool onOff){
    std::string name(node->GetName());
    if ( name.find(str) != std::string::npos ){
      node->SetVisibility(onOff);
    }
    int ndau = node->GetNdaughters();
    for (int i=0; i<ndau; ++i){
      TGeoNode * dau = node->GetDaughter(i);
      showNodesByName(dau, str, onOff);
    }
  }

  //Function to show CRV. Run after InsideDS
  void Geom_Interface::InsideCRV( TGeoNode * node, bool inCRVVac ){
    static std::vector <std::string> substrings {"CRSAluminum", "CRV", "CRS"};
    for(auto& i: substrings) showNodesByName(node, i, kTRUE);
  }

  //Function to hide all elements which are not PS,TS, DS:
  void Geom_Interface::SolenoidsOnly(TGeoNode* node) {
    static std::vector <std::string> substrings  { "Ceiling",
    "backfill", "dirt", "concrete", "VirtualDetector",
    "pipeType","ExtShield", "PSShield"};
    //,"CRSAluminium","CRV","CRS", 
    for(auto& i: substrings) hideNodesByName(node,i,kFALSE, 0);
    static std::vector <std::string> materials { "MBOverburden", "CONCRETE"};
    for(auto& i: materials) hideNodesByMaterial(node,i,kFALSE);
  }

  //Funciton to hide building top
  void Geom_Interface::hideTop(TGeoNode* node, int _diagLevel) {
    TString name = node->GetName();
    if(_diagLevel > 0 and name.Index("Shield")>0) {
    std::cout << name << " " <<  name.Index("mBox_") << std::endl;
    }
    bool test = false;

    if(name.Index("mBox_45_")>=0) test = true;
    if(name.Index("mBox_46_")>=0) test = true;
    if(name.Index("mBox_47_")>=0) test = true;
    if(name.Index("mBox_48_")>=0) test = true;
    if(name.Index("mBox_49_")>=0) test = true;
    if(name.Index("mBox_74_")>=0) test = true;

    if(test) {
      std::cout << "turning off " << name << std::endl;
      node->SetVisibility(false);
    }

    // Descend recursively into each daughter TGeoNode.
    int ndau = node->GetNdaughters();
    for ( int i=0; i<ndau; ++i ){
      TGeoNode * dau = node->GetDaughter(i);
      hideTop( dau, _diagLevel );
    }
  }

  void Geom_Interface::TrackerVolumeHeirarchy( TGeoNode * node, std::vector<CLHEP::Hep3Vector> &TransformList ){
    std::string _name = (node->GetVolume()->GetName());
    if( _name == "HallAir") {
      cout<<"HallAir Origin IS "<<node->GetMotherVolume()->GetName();
      TGeoVolume *vol = node->GetVolume();
      TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
      Double_t master[3];
      const Double_t *local = shape->GetOrigin();
      if(shape!=NULL){
        gGeoManager->LocalToMaster(local,master);
        CLHEP::Hep3Vector hallToworld(master[0], master[1], master[2]);
        TransformList.push_back(hallToworld);
      }
    }
    if( _name == "DS3Vacuum") {
      cout<<"DS3 Origin IS "<<node->GetMotherVolume()->GetName();
      TGeoVolume *vol = node->GetVolume();
      TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
      Double_t master[3];
      const Double_t *local = shape->GetOrigin();
      if(shape!=NULL){
        gGeoManager->LocalToMaster(local,master);
        CLHEP::Hep3Vector DSTohall(master[0], master[1], master[2]);
        TransformList.push_back(DSTohall);
      }
    }
    if( _name == "TrackerMother") {
      cout<<"Tracker Origin IS "<<node->GetMotherVolume()->GetName();
      TGeoVolume *vol = node->GetVolume();
      TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
      Double_t master[3];
      const Double_t *local = shape->GetOrigin();
      if(shape!=NULL){
        gGeoManager->LocalToMaster(local,master);
        CLHEP::Hep3Vector TrackerToDS(master[0], master[1], master[2]);
        TransformList.push_back(TrackerToDS);
      }
    }
    // Descend into each daughter TGeoNode.
    int ndau = node->GetNdaughters();
    for ( int i=0; i<ndau; ++i ){
      TGeoNode * dau = node->GetDaughter(i);
      TrackerVolumeHeirarchy(dau, TransformList);
    }
 }
}

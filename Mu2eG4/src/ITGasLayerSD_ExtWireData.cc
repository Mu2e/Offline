// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD_ExtWireData.hh"

#include "TMath.h"

using namespace std;

namespace mu2e {

  ITGasLayerSD_ExtWireData::ITGasLayerSD_ExtWireData(G4String name) : ITGasLayerSD(name) {
  }

  G4bool ITGasLayerSD_ExtWireData::ProcessHits(G4Step* aStep,G4TouchableHistory*){

	int ring = _ring;
	if(ring==0) return false;
    G4double edep = aStep->GetTotalEnergyDeposit();

    // Eventually we will want this but not now.
    //if(edep==0.) return false;

    // Origin of the ITracker.  Need to get this from G4.
    static G4ThreeVector detectorOrigin( -3904., -7350., 6200.);

    // Position at start of step point, in world system and in
    // a system in which the center of the tracking detector is the origin.
    G4ThreeVector prePosWorld   = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector prePosTracker = prePosWorld - detectorOrigin;

    G4ThreeVector preMomWorld = aStep->GetPreStepPoint()->GetMomentum();

    GeomHandle<ITracker> itracker;

    G4ThreeVector pos(prePosTracker);
    pos*=0.1;

    Double_t phihit=TMath::ATan2(pos[0],pos[1]);
    Double_t hit[4]={pos[0],pos[1],pos[2],aStep->GetPreStepPoint()->GetGlobalTime()};


    int ring1=ring-2;
    int ring2=ring-1;

    float xywire[3];
    double phiwire1=0,phiwire0=0;
    double phirel1=0,phirel2=0,radius1=0,radius2=0;
    int wire1=-1,wire2=-1;
    double invcosstereo;

    if(ring1>=0){
       itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,1);
       invcosstereo=1./TMath::Cos(itracker->getCellGeometryHandle()->GetWireEpsilon());
       itracker->getCellGeometryHandle()->WirePosAtZ(pos[2]*invcosstereo,xywire);
       phiwire1=TMath::ATan2(xywire[1],xywire[0]);
       itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,0);
       itracker->getCellGeometryHandle()->WirePosAtZ(pos[2]*invcosstereo,xywire);
       phiwire0=TMath::ATan2(xywire[1],xywire[0]);
       phirel1=phihit-phiwire0;
       if(phiwire1<phiwire0) phirel1*=-1.;

       if(phirel1<-0.5*_Dphi) phirel1+=M_2PI;
       if(phirel1>=M_2PI-0.5*_Dphi) phirel1-=M_2PI;

       wire1=int((phirel1+0.5*_Dphi)/_Dphi);
       if(wire1<0||wire1>=_nwires)
     	  /*throw cms::Exception("GEOM")*/cerr<<"Wrong wire number "<<wire1<< "\n";

       phirel1=phirel1-wire1*_Dphi;
       if(phiwire1<phiwire0) phirel1*=-1.;

       radius1=TMath::Hypot(xywire[0],xywire[1]);
     }

     if(ring2<itracker->nRing()-1){
       itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring2,1);
       invcosstereo=1./TMath::Cos(itracker->getCellGeometryHandle()->GetWireEpsilon());
       itracker->getCellGeometryHandle()->WirePosAtZ(pos[2]*invcosstereo,xywire);
       phiwire1=TMath::ATan2(xywire[1],xywire[0]);
       itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring2,0);
       itracker->getCellGeometryHandle()->WirePosAtZ(pos[2]*invcosstereo,xywire);
       phiwire0=TMath::ATan2(xywire[1],xywire[0]);

       phirel2=phihit-phiwire0;
       if(phiwire1<phiwire0) phirel2*=-1.;

       if(phirel2<-0.5*_Dphi) phirel2+=M_2PI;
       if(phirel2>=M_2PI-0.5*_Dphi) phirel2-=M_2PI;

       wire2=int((phirel2+0.5*_Dphi)/_Dphi);
       if(wire2<0||wire2>=_nwires)
     	  /*throw cms::Exception("GEOM")*/cerr<<"Wrong wire number "<<wire2<<"\n";

       phirel2=phirel2-wire2*_Dphi;
       if(phiwire1<phiwire0) phirel2*=-1.;

       radius2=TMath::Hypot(xywire[0],xywire[1]);
     }

     double dist1,dist2;
     if(wire1>=0){
       itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,wire1);
       dist1=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
     }else{
       radius1=radius2-radius2*TMath::Sin(_Dphi/3./2.)*2.*TMath::Sqrt(3.)/2.;
       double phiw=phiwire0+wire2*(phiwire1-phiwire0)+TMath::Sign(0.5*_Dphi,phirel2);
       xywire[0]=TMath::Cos(phiw)*radius1;
       xywire[1]=TMath::Sin(phiw)*radius1;
       dist1=TMath::Hypot(hit[0]-xywire[0],hit[1]-xywire[1]);
     }

     if(wire2>=0){
       itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring2,wire2);
       dist2=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
     }else{
       radius2=radius1+radius1*TMath::Sin(_Dphi/3./2.)*2.*TMath::Sqrt(3.)/2.;
       double phiw=phiwire0+wire1*(phiwire1-phiwire0)+TMath::Sign(0.5*_Dphi,phirel1);
       xywire[0]=TMath::Cos(phiw)*radius2;
       xywire[1]=TMath::Sin(phiw)*radius2;
       dist2=TMath::Hypot(hit[0]-xywire[0],hit[1]-xywire[1]);
     }

     int wire;
     if(dist1<=dist2 && (ring1>=0 &&  wire1>=0) ){
       ring=ring1;
       wire=wire1;
     }else{
       ring=ring2;
       wire=wire2;
     }


     unsigned long det=itracker->getCellGeometryHandle()->computeDet(_superlayer,ring,wire);

    StepPointG4* newHit =
      new StepPointG4( aStep->GetTrack()->GetTrackID(),
     		 det,
 		     edep,
 		     prePosTracker,
 		     preMomWorld,
 		     aStep->GetPreStepPoint()->GetGlobalTime()
 		     );

    // The collection takes ownership of the hit. 
    _collection->insert( newHit );

    newHit->Print();
    //newHit->Draw();
    
    return true;
  }
  
} //namespace mu2e


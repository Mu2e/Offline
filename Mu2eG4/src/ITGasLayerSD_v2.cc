// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD_v2.hh"

#include "TMath.h"
#include "TRandom.h"

using namespace std;

namespace mu2e {

ITGasLayerSD_v2::ITGasLayerSD_v2(G4String name) : ITGasLayerSD(name){
}

G4bool ITGasLayerSD_v2::ProcessHits(G4Step* aStep,G4TouchableHistory*){

        int ring = _ring;
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

        double phihit=TMath::ATan2(pos[1],pos[0]);
        double hit[4]={pos[0],pos[1],pos[2],aStep->GetPreStepPoint()->GetGlobalTime()};


        int ring1;

        float xywire[3];
        double phiwire1=0;
        int wire1=-1;
        double cosstereo;
        double invcosstereo;

        try {
                //std::cout<<"S "<<_superlayer<<" R "<<ring<<std::endl;
                itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,0);

                double Dphi=_Dphi/3.0;

                cosstereo=TMath::Cos(itracker->getCellGeometryHandle()->GetWireEpsilon());
                invcosstereo=1.0/cosstereo;
                itracker->getCellGeometryHandle()->WirePosAtZ(pos[2]/**invcosstereo*/,xywire);

//                std::cout<<"0 wire center "<<xywire[0]<<" "<<xywire[1]<<" "<<xywire[2]<<std::endl;
//                std::cout<<"hit pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
//                std::cout<<"hit 0 wire dist "<<sqrt(pow(pos[0]-xywire[0],2)+pow(pos[1]-xywire[1],2)+pow(pos[2]-xywire[2],2))<<std::endl;

//                //---------------- test --------------------
//                //TRandom *rn = new TRandom();
//                pos[0]=xywire[0]+gRandom->Uniform(-6,6);
//                pos[1]=xywire[1]+gRandom->Uniform(-6,6);
//                hit[0]=pos[0];
//                hit[1]=pos[1];
//                prePosTracker[0]=pos[0];
//                prePosTracker[1]=pos[1];
//                phihit=TMath::ATan2(pos[1],pos[0]);
//                //------------------------------------------

                phiwire1=TMath::ATan2(xywire[1],xywire[0]);

//                std::cout<<"phihit "<<phihit <<" phiwire1 "<<phiwire1 <<" _Dphi "<<_Dphi<<std::endl;

                phihit-=phiwire1;
                int wire;
//                wire = TMath::Nint(phihit/_Dphi);
                wire = TMath::Nint(phihit/Dphi);

                if (wire<0) wire+=3*_nwires;

                double tdist, mdist;
                mdist=1.0e+9;
                int maxW = _nwires-1;
                int maxR = itracker->nRing()-1;
                switch (wire%3) {
                        case 1:
                                wire/=3;
                                itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,wire);
                                mdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                                if(ring==0 || ring==maxR) {
                                        if (mdist>itracker->getCellGeometryHandle()->GetCellRad()) return false;
                                }
                                else {
                                        ring1=ring;
                                        ++ring1;
                                        wire1=wire;
                                        itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,wire1);
                                        tdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                                        if (mdist>tdist) {
                                                mdist=tdist;
                                                ring=ring1;
                                                break;
                                        }
                                        ++wire1;
                                        if (wire1>maxW) wire1-=_nwires;
                                        itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,wire1);
                                        tdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                                        if (mdist>tdist) {
                                                mdist=tdist;
                                                ring=ring1;
                                                wire=wire1;
                                                break;
                                        }
                                }
                                break;
                        case 2:
                                wire/=3;
                                ring1=ring;
                                ++ring1;
                                wire1=wire;
                                ++wire;
                                if (wire>maxW) wire-=_nwires;
                                itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,wire);
                                mdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                                if (ring==0 || ring==maxR) {
                                        if (mdist>itracker->getCellGeometryHandle()->GetCellRad()) return false;
                                }else {
                                        itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,wire1);
                                        tdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                                        if (mdist>tdist) {
                                                mdist=tdist;
                                                ring=ring1;
                                                wire=wire1;
                                                break;
                                        }
                                        itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring1,wire/*wire1*/);
                                        tdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                                        if (mdist>tdist) {
                                                mdist=tdist;
                                                ring=ring1;
                                                break;
                                        }
                                }
                                break;
                        default:
                                wire/=3;
                                break;
                }

                if(ring==0) ++ring;

//                std::cout<<"selected wire "<<ring<<" "<<wire<<std::endl;

                unsigned long det=itracker->getCellGeometryHandle()->computeDet(_superlayer,ring,wire);

                StepPointG4* newHit =
                                new StepPointG4( aStep->GetTrack()->GetTrackID(),
                                                det,
                                                edep,
                                                prePosTracker,
                                                preMomWorld,
                                                aStep->GetPreStepPoint()->GetGlobalTime(),
                                                aStep->GetPreStepPoint()->GetProperTime(),
                                                aStep->GetStepLength()
                                );

                // The collection takes ownership of the hit.
                _collection->insert( newHit );

                newHit->Print();
                //newHit->Draw();

                return true;
        }catch (cms::Exception e) {
                cerr<<e;
                return false;
        }

}

} //namespace mu2e

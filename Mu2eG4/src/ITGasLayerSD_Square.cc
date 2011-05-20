#include "Mu2eG4/inc/ITGasLayerSD_Square.hh"

#include "TMath.h"
#include "cetlib/pow.h"
#include <string>

using namespace std;

using cet::square;

namespace mu2e {

  ITGasLayerSD_Square::ITGasLayerSD_Square(G4String name, const SimpleConfig& config) : ITGasLayerSD( name, config) { }


  ITGasLayerSD_Square::~ITGasLayerSD_Square() { }

  G4bool ITGasLayerSD_Square::ProcessHits(G4Step* aStep, G4TouchableHistory*){

          _currentSize += 1;

          if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
            if( (_currentSize - _sizeLimit)==1 ) {
              mf::LogWarning("G4") << "Maximum number of particles reached in ItrackerSD: "
                                    << _currentSize << endl;
            }
            return false;
          }

          G4double edep  = aStep->GetTotalEnergyDeposit();
          G4double nidep = aStep->GetNonIonizingEnergyDeposit();
          G4double step  = aStep->GetStepLength();
          G4double idep  = edep-nidep;

          if ( _debugList.inList() )  cout<<"edep "<<edep<<" nidep "<<nidep<<" step "<<step<<endl;
          // I am not sure why we get these cases but we do.  Skip them.
          if ( (edep == 0. || idep == 0.)/*&& step == 0.*/ ) {
                  if ( _debugList.inList() )  cout<<"Skipped"<<endl;
                  return false;
          }

          string volName = aStep->GetTrack()->GetVolume()->GetName();
          if ( _debugList.inList() )  cout<<"Step vol name "<<aStep->GetTrack()->GetVolume()->GetName()<<endl;

          _superlayer=atoi(volName.substr(5,2).c_str());
          _ring=atoi(volName.substr(8,2).c_str());

          try {
              itracker->getCellGeometryHandle()->SelectCell(_superlayer,_ring,0);
                  _nwires=itracker->getCellGeometryHandle()->GetITLayer()->nCells();
                  _Dphi=CLHEP::twopi/_nwires;
          }catch (cet::exception e) {
                  cerr<<e;
                  _nwires=0;
                  _Dphi=0.0;
          }

          int ring = _ring;
          //        if(ring==0) ring++;/*return false;*/
          //G4double edep = aStep->GetTotalEnergyDeposit();

          //cout<<"Step vol name "<<aStep->GetTrack()->GetVolume()->GetName()<<endl;

          // Eventually we will want this but not now.
          //if(edep==0.) return false;


          // Origin of the ITracker.  Need to get this from G4.
          //static G4ThreeVector detectorOrigin( -3904., -7350., 6200.);


          // Position at start of step point, in world system and in
          // a system in which the center of the tracking detector is the origin.
          G4ThreeVector prePosWorld   = aStep->GetPreStepPoint()->GetPosition();
          if ( _debugList.inList() )  std::cout<<"G4 hit pos in World"<<prePosWorld[0]<<" "<<prePosWorld[1]<<" "<<prePosWorld[2]<<std::endl;
          G4ThreeVector prePosTracker = prePosWorld - _mu2eDetCenter;

          G4ThreeVector preMomWorld = aStep->GetPreStepPoint()->GetMomentum();

          //GeomHandle<ITracker> itracker;

          G4ThreeVector pos(prePosTracker);

          double phihit=TMath::ATan2(pos[1],pos[0]);
          double hit[4]={pos[0],pos[1],pos[2],aStep->GetPreStepPoint()->GetGlobalTime()};

          float xywire[3];
          double phiwire1=0;
          int wire1=-1;
          double cosstereo;
          double invcosstereo;

          try {
                  if ( _debugList.inList() )  std::cout<<"S "<<_superlayer<<" R "<<ring<<std::endl;

                  itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,0);

                  cosstereo=TMath::Cos(itracker->getCellGeometryHandle()->GetWireEpsilon());
                  invcosstereo=1.0/cosstereo;
                  if ( _debugList.inList() )  std::cout<<"stereo "<<itracker->getCellGeometryHandle()->GetWireEpsilon()<<" invcosstereo "<<invcosstereo<<std::endl;
                  itracker->getCellGeometryHandle()->WirePosAtZ(pos[2]/**invcosstereo*/,xywire);

                  if ( _debugList.inList() )  std::cout<<"0 wire center "<<xywire[0]<<" "<<xywire[1]<<" "<<xywire[2]<<std::endl;
                  if ( _debugList.inList() )  std::cout<<"hit pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
                  if ( _debugList.inList() )  std::cout<<"hit to wire 0 dist "<<sqrt(square(pos[0]-xywire[0])+square(pos[1]-xywire[1])+square(pos[2]-xywire[2]))<<std::endl;

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

                  if ( _debugList.inList() )  std::cout<<"phihit "<<phihit <<" phiwire1 "<<phiwire1 <<" _Dphi "<<_Dphi<<std::endl;

                  phihit-=phiwire1;
                  int wire;
                  wire = TMath::Nint(phihit/_Dphi);
                  //                wire = TMath::Nint(phihit/Dphi);

                  if (wire<0) wire+=_nwires;

                  double tdist, mdist;
                  mdist=1.0e+9;
                  int maxW = _nwires-1;

                  itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,wire);
                  mdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                  wire1=wire;
                  --wire1;
                  if (wire1<0) wire1+=_nwires;
                  itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,wire1);
                  tdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                  if (mdist>tdist) {
                          mdist=tdist;
                          wire=wire1;
                  }else {
                          ++wire1;
                          ++wire1;
                          if (wire1>maxW) wire1-=_nwires;
                          itracker->getCellGeometryHandle()->SelectCell(_superlayer,ring,wire1);
                          tdist=itracker->getCellGeometryHandle()->DistFromWireCenter(hit);
                          if (mdist>tdist) {
                                  mdist=tdist;
                                  wire=wire1;
                          }
                  }
                  //if (mdist>itracker->getCellGeometryHandle()->GetCellRad()) return false;

                  if ( _debugList.inList() )  std::cout<<"Cirumscribed Cell Radius "<<itracker->getCellGeometryHandle()->GetCellRad()<<std::endl;
                  if ( _debugList.inList() )  std::cout<<"selected wire "<<ring<<" "<<wire<<std::endl;

                  unsigned long det=itracker->getCellGeometryHandle()->computeDet(_superlayer,ring,wire);

   //                StepPointG4* newHit =
   //                                new StepPointG4( aStep->GetTrack()->GetTrackID(),
   //                                                det,
   //                                                edep,
   //                                                prePosTracker,
   //                                                preMomWorld,
   //                                                aStep->GetPreStepPoint()->GetGlobalTime(),
   //                                                aStep->GetPreStepPoint()->GetProperTime(),
   //                                                aStep->GetStepLength()
   //                                );
   //
   //                // The collection takes ownership of the hit.
   //                _collection->insert( newHit );
   //
   //                newHit->Print();
   //                //newHit->Draw();

                  // We add the hit object to the framework itrackerHit collection created in produce

                  _collection->push_back( StepPointMC(aStep->GetTrack()->GetTrackID(),
                                                      det,
                                                      idep,/*edep,*/
                                                      aStep->GetNonIonizingEnergyDeposit(),
                                                      aStep->GetPreStepPoint()->GetGlobalTime(),
                                                      aStep->GetPreStepPoint()->GetProperTime(),
                                                      prePosTracker,
                                                      preMomWorld,
                                                      step,
                                                      ProcessCode()
                                                      ));

                  // Some debugging tests.
                  if ( !_debugList.inList() ) return true;

                  return true;
          }catch (cet::exception e) {
                  cerr<<e;
                  return false;
          }

   }


} //namespace mu2e

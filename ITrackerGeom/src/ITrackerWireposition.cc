///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  DCH wireposition class                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "ITrackerGeom/inc/ITrackerWireposition.hh"
#include "TError.h"
#include "TGeoMatrix.h"
#include "TMath.h"
#include "TSystem.h"
#include "cetlib/pow.h"
#include <iostream>
#include <sys/stat.h>

using namespace std;

using cet::diff_of_squares;
using cet::sum_of_squares;

namespace mu2e {
//_____________________________________________________________________________
ITrackerWireposition::ITrackerWireposition(Bool_t newFile, const char *WireDataFile )
{

  //
  // ITrackerWireposition default constructor
  //

  fwiredata= new ITrackerWiredata();
  if (newFile) {
    wirefile = new TFile(WireDataFile,"RECREATE");
    trwdata = new TTree("WireData","WireData");
    trwdata->Branch("WireDataMatrix","ITrackerWiredata",&fwiredata);

  }
  else {

    struct stat stFileInfo;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(WireDataFile,&stFileInfo);
    if(intStat == 0) {
      // We were able to get the file attributes
      // so the file obviously exists.

      wirefile = new TFile(WireDataFile,"READ");
      trwdata = (TTree*) wirefile->Get("WireData");
      trwdata->SetBranchAddress("WireDataMatrix",&fwiredata);
      for(int i=0;i< trwdata->GetEntries();i++){
        fWireDataAll[i]=new ITrackerWiredata();
              trwdata->SetBranchAddress("WireDataMatrix",&fWireDataAll[i]);
              trwdata->GetEntry(i);
      }


    }
    else {
      // We were not able to get the file attributes.
      // This may mean that we don't have permission to
      // access the folder which contains this file. If you
      // need to do that level of checking, lookup the
      // return values of stat which will give you
      // more details on why stat failed.

      wirefile = new TFile(WireDataFile,"RECREATE");
      trwdata = new TTree("WireData","WireData");
      trwdata->Branch("WireDataMatrix","ITrackerWiredata",&fwiredata);

    }
  }

  fSuperLayer=-1;
  fLayer=-1;
  fWire=-1;
  SelectedSL=-1;
  SelectedL=-1;
  SelectedCell=-1;

  selectedMat = 0x0;
  selectedAlfa = 0.0;
  selectedEpsilon = 0.0;
  selectedRadius = 0.0;

  DCHEndcapZ = 150.;
  SenseWireRadius = 0.5 * .0015;

}

//_____________________________________________________________________________
void ITrackerWireposition::AddSuperLayer(Int_t nCelLayer, Int_t nwire){
  fwiredata->NcelLayer = nCelLayer;
  fwiredata->epsilon = new Float_t [nCelLayer];
  fwiredata->alfa = new Float_t [nCelLayer];
  fwiredata->radius_z0 = new Float_t [nCelLayer];
  fwiredata->PosMatrix->Clear();
  fwiredata->PosMatrix->Expand( nCelLayer );
  for (Int_t i=0; i<nCelLayer; i++){
    fwiredata->PosMatrix->AddAt( new TObjArray(nwire), i );
  }
}

void ITrackerWireposition::InsertAlfaEpsilon(Int_t iCel, Float_t alfa, Float_t eps){
  fwiredata->epsilon[iCel]=eps;
  fwiredata->alfa[iCel]=alfa;
}

void ITrackerWireposition::InsertRadius(Int_t iCel, Float_t rad){
  fwiredata->radius_z0[iCel]=rad;
}

void ITrackerWireposition::InsertWireMatrix(Int_t nCelL, Int_t nw, TGeoHMatrix *matrix){
  ((TObjArray*)fwiredata->PosMatrix->At(nCelL))->AddAt( matrix, nw );
}

void ITrackerWireposition::WriteData(){
  wirefile->cd();
  trwdata->Write();
}

void ITrackerWireposition::FillData(){
  trwdata->Fill();
}



 void ITrackerWireposition::SelectWire(Int_t SupLayer, Int_t CelLayer, Int_t Wire){

  if (SupLayer==SelectedSL && CelLayer==SelectedL && Wire==SelectedCell) return;
  if (SupLayer!=SelectedSL) fwiredata=fWireDataAll[SupLayer];


  selectedMat = (TGeoHMatrix *) ((TObjArray*)fwiredata->PosMatrix->At(CelLayer))->At(Wire);
  selectedAlfa = fwiredata->alfa[CelLayer];
  selectedEpsilon = fwiredata->epsilon[CelLayer];
  selectedRadius = fwiredata->radius_z0[CelLayer];
  SelectedSL=SupLayer;
  SelectedL=CelLayer;
  SelectedCell=Wire;

}

void ITrackerWireposition::SelectWireDet(ULong_t det){

  // Return the SuperLayer
  fSuperLayer=(Int_t)(det*0.00001);

  //Return the Layer
  fLayer=(Int_t)((det)-((fSuperLayer)*100000));

  fLayer/=1000;

  //Return the Wire
  fWire=(((det)-((fSuperLayer)*100000))-fLayer*1000);

  fSuperLayer--;

  //Call the upper method
  SelectWire(fSuperLayer,fLayer,fWire);

}

Double_t ITrackerWireposition::DistFromWireCenter(Double_t *global){

  selectedMat->MasterToLocal(global, templocal);

  return TMath::Sqrt( sum_of_squares(templocal[0], templocal[1]) );
}

Double_t ITrackerWireposition::DistFromWire(Double_t *global){

  selectedMat->MasterToLocal(global, templocal);

  return TMath::Sqrt( sum_of_squares(templocal[0], templocal[1]) )
       - SenseWireRadius;
}

void ITrackerWireposition::WirePosAtEndcap(Float_t *right, Float_t *left){

//  if (fParam->GetEndCapType()==1) {
//    DCHEndcapZ = TMath::Sqrt( diff_of_squares(fParam->GetMaxEndCapDim()-fParam->GetEndCapWallThickness(), selectedRadius) );
//  }

  templocal[0] = 0.;
  templocal[1] = 0.;
  templocal[2] = DCHEndcapZ/TMath::Cos(selectedEpsilon);

  selectedMat->LocalToMaster(templocal, tempglobal);
  right[0] = tempglobal[0];
  right[1] = tempglobal[1];
  right[2] = tempglobal[2];

  templocal[2] *= -1.;
  selectedMat->LocalToMaster(templocal, tempglobal);
  left[0] = tempglobal[0];
  left[1] = tempglobal[1];
  left[2] = tempglobal[2];

}

void ITrackerWireposition::WirePosAtZ(Float_t z, Float_t *pos){

//  if (fParam->GetEndCapType()==1) {
//    DCHEndcapZ = TMath::Sqrt( diff_of_squares(fParam->GetMaxEndCapDim()-fParam->GetEndCapWallThickness(), selectedRadius) );
//  }

  if (TMath::Abs(z)>DCHEndcapZ) {
    cout<<"z exceeds maximum length\n";
    return;
  }

  templocal[0] = 0.;
  templocal[1] = 0.;
  templocal[2] = z;

  selectedMat->LocalToMaster(templocal, tempglobal);
  pos[0] = tempglobal[0];
  pos[1] = tempglobal[1];
  pos[2] = tempglobal[2];

}

} // namespace mu2e


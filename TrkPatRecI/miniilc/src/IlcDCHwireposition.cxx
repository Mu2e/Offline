/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
 
/////////////////////////////////////////////////////////////////////////////// 
//                                                                           // 
//  DCH wireposition class                                                   // 
//                                                                           // 
/////////////////////////////////////////////////////////////////////////////// 

 
#include <sys/stat.h> 

#include <TError.h> 
#include <TGeoMatrix.h> 
#include <TMath.h> 
#include <TSystem.h> 
#include "IlcDCHwireposition.h" 
#include <iostream>
using namespace std; 
//#include "IlcLog.h"
 
//ClassImp(IlcDCHwireposition) 
 
//_____________________________________________________________________________ 
IlcDCHwireposition::IlcDCHwireposition(int newFile, const char *WireDataFile ) 
{ 
  



  // 
  // IlcDCHwireposition default constructor 
  // 
  wirefile=0;
  trwdata=0;
  fwiredata= new IlcDCHwiredata();
  //IlcDCHwiredata*  fWireDataAll[20];
  if (newFile>0) {
    wirefile = new TFile(WireDataFile,"RECREATE");
    trwdata = new TTree("WireData","WireData");
    trwdata->Branch("WireDataMatrix","IlcDCHwiredata",&fwiredata);
    
  }
  else if(newFile==0){
 
    struct stat stFileInfo;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(WireDataFile,&stFileInfo);
    if(intStat == 0) {
      // We were able to get the file attributes
      // so the file obviously exists.
      
      wirefile = new TFile(WireDataFile,"READ");
      trwdata = (TTree*) wirefile->Get("WireData");
      //ricommentare la stringa sottostante
      trwdata->SetBranchAddress("WireDataMatrix",&fwiredata);
      for(int i=0;i< trwdata->GetEntries();i++){	
        fWireDataAll[i]=new IlcDCHwiredata();
      	trwdata->SetBranchAddress("WireDataMatrix",&fWireDataAll[i]);
      	trwdata->GetEntry(i);
//	fWireDataAll[i] = (IlcDCHwiredata*)fwiredata->Clone();
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
      trwdata->Branch("WireDataMatrix","IlcDCHwiredata",&fwiredata);

      // aggiungere la chamata alla creazione della geometria
      // WriteData();
    }
  }else{

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

//  IlcInfo(Form("Warning default Param instanced\n"));
  cout<<"Warning default Param instanced\n";
  fParam = new IlcDCHParam();

  DCHEndcapZ = fParam->GetLength()-fParam->GetEndCapWallThickness();
  SenseWireRadius = 0.5 * fParam->GetSWireDiameter();

//   fMatrixArray           = 0; 
//   fMatrixCorrectionArray = 0; 
//   
//   Init(); 
  
 
} 
 
//_____________________________________________________________________________ 
//IlcDCHwireposition::~IlcDCHwireposition() 
//{ 
  // 
  // IlcDCHwireposition destructor 
  // 
 // gSystem->Exec("echo punto2 distruttore;ps xu|grep ilcroot");

  //wirefile->Close();

  //if (fwiredata) {
  //cout<<fwiredata<<" cancello\n"; 
  //fwiredata->Delete();
  //delete fwiredata;
  
  
 //}

//   delete fMatrixArray; 
//   delete fMatrixCorrectionArray; 
 
//} 
 
void IlcDCHwireposition::AddSuperLayer(Int_t nCelLayer, Int_t nwire){
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

void IlcDCHwireposition::InsertAlfaEpsilon(Int_t iCel, Float_t alfa, Float_t eps){
  fwiredata->epsilon[iCel]=eps;
  fwiredata->alfa[iCel]=alfa;
}

void IlcDCHwireposition::InsertRadius(Int_t iCel, Float_t rad){
  fwiredata->radius_z0[iCel]=rad;
}

void IlcDCHwireposition::InsertWireMatrix(Int_t nCelL, Int_t nw, TGeoHMatrix *matrix){
  ((TObjArray*)fwiredata->PosMatrix->At(nCelL))->AddAt( matrix, nw );
}

void IlcDCHwireposition::WriteData(){
//  trwdata->Fill();
  wirefile->cd();
  trwdata->Write();
}

void IlcDCHwireposition::FillData(int suplay){
  if(trwdata) trwdata->Fill();
  if(suplay>=0){
    fWireDataAll[suplay]=fwiredata;
    fwiredata=new IlcDCHwiredata();
  }
}



 void IlcDCHwireposition::SelectWire(Int_t SupLayer, Int_t CelLayer, Int_t Wire){

//gSystem->Exec("echo punto1 chiamata di SelectWire;ps xu|grep ilcroot");
  
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

void IlcDCHwireposition::SelectWireDet(ULong_t det){

  // Return the SuperLayer 
  fSuperLayer=(det*0.00001);


  //Return the Layer
  fLayer=((det)-((fSuperLayer)*100000));
  
  fLayer*=0.001;

  //Return the Wire
  fWire=(((det)-((fSuperLayer)*100000))-fLayer*1000);
  
  fSuperLayer--;

  //Call the upper method
  SelectWire(fSuperLayer,fLayer,fWire);

}

void IlcDCHwireposition::SetParam(IlcDCHParam *par) {
  if (fParam) delete fParam;
  fParam = par;
  
  cout <<"Param selected\n";

  DCHEndcapZ = fParam->GetLength()-fParam->GetEndCapWallThickness();
  SenseWireRadius = 0.5 * fParam->GetSWireDiameter();
  
}
Double_t IlcDCHwireposition::DistFromWireCenter(Double_t *global){
  
  selectedMat->MasterToLocal(global, templocal);
  
  return (TMath::Sqrt( pow(templocal[0],2) + pow(templocal[1],2) ));
}

Double_t IlcDCHwireposition::DistFromWire(Double_t *global){
  
  selectedMat->MasterToLocal(global, templocal);
  
  return (TMath::Sqrt( pow(templocal[0],2) + pow(templocal[1],2) ) - SenseWireRadius);
}

void IlcDCHwireposition::WirePosAtEndcap(Float_t *right, Float_t *left){

  if (fParam->GetEndCapType()==1) {
    DCHEndcapZ = TMath::Sqrt( pow(fParam->GetMaxEndCapDim()-fParam->GetEndCapWallThickness(),2) - pow(selectedRadius,2) );
  }

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

void IlcDCHwireposition::WirePosAtZ(Float_t z, Float_t *pos){

  if (fParam->GetEndCapType()==1) {
    DCHEndcapZ = TMath::Sqrt( pow(fParam->GetMaxEndCapDim()-fParam->GetEndCapWallThickness(),2) - pow(selectedRadius,2) );
  }

  if (TMath::Abs(z)>DCHEndcapZ) {
//    IlcInfo("z exceeds maximum length\n");
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

#ifndef ITRACKERWIREPOSITION_HH
#define ITRACKERWIREPOSITION_HH

#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TGeoMatrix.h"

#include "ITrackerGeom/inc/ITrackerWiredata.hh"

namespace mu2e {

class ITrackerWireposition { 
 
  public: 
 
    ITrackerWireposition(Bool_t newFile = kFALSE, const char *WireDataFile = "DCHWireData.root");
    virtual ~ITrackerWireposition()
    {
      if (fwiredata) {delete fwiredata;}
       wirefile->Close();
    }
 
    void FillData();
    void WriteData();
    void AddSuperLayer(Int_t nCel, Int_t nwire);
    void InsertAlfaEpsilon(Int_t iCel, Float_t alfa, Float_t eps);
    void InsertRadius(Int_t iCel, Float_t rad);
    void InsertWireMatrix(Int_t nCelL, Int_t nw, TGeoHMatrix *matrix);
    void SelectWire(Int_t SupLayer, Int_t CelLayer, Int_t Wire);
    void SelectWireDet(ULong_t det);// Det Method
    void SetSuperLayer(Int_t suplay) { fSuperLayer = suplay; }
    void SetCelRing  (Int_t lay)       { fLayer  = lay; }
    void SetWire (Int_t wi)    { fWire = wi; }
    void Global2Local(Double_t *global, Double_t *local) { selectedMat->MasterToLocal(global, local); }
    void Local2Global(Double_t *local, Double_t *global) { selectedMat->LocalToMaster(local, global); }
    void WirePosAtEndcap(Float_t *right, Float_t *left);
    void WirePosAtZ(Float_t z, Float_t *pos);
    
    Int_t  GetSuperLayer()                { return  fSuperLayer; }
    Int_t  GetCelRing()                   { return  fLayer; }
    Int_t  GetWire()                      { return  fWire; }
    Float_t GetWireAlfa()                 { return selectedAlfa; }
    Float_t GetWireEpsilon()              { return selectedEpsilon; }
    Float_t GetLayerRad()                 { return selectedRadius; }
    
    Double_t DistFromWire(Double_t *global);
    Double_t DistFromWireCenter(Double_t *global);

  TGeoHMatrix *GetGeoMatrix(){return selectedMat;};

  private:
    ITrackerWiredata*  fWireDataAll[24];
    Int_t     fSuperLayer;           
    Int_t     fLayer;         
    Int_t     fWire;
    ITrackerWiredata *fwiredata;
    TFile *wirefile;
    TTree *trwdata;
    TGeoHMatrix *selectedMat;
    Float_t  selectedAlfa;
    Float_t  selectedEpsilon;
    Float_t  selectedRadius;
    Double_t templocal[3];
    Double_t tempglobal[3];
    Float_t  DCHEndcapZ;
    Double_t SenseWireRadius;
    Int_t SelectedSL, SelectedL, SelectedCell;
 
}; 
 
}  //namespace mu2e

#endif /*ITRACKERWIREPOSITION_HH*/

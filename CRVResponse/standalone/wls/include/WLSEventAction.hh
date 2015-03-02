#ifndef WLSEventAction_h
#define WLSEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

class TH1D;
class TH2D;
class TH3D;
class TFile;

class WLSEventAction : public G4UserEventAction
{
  public:

    WLSEventAction(int mode, int id=0);  //id is only needed for mode -1
    ~WLSEventAction();

  public:

    void   BeginOfEventAction(const G4Event*);
    void     EndOfEventAction(const G4Event*);

    static WLSEventAction* Instance()                      {return _fgInstance;}
    void                   SetOptPhotonStart(G4ThreeVector start) {_start=start;}
    void                   SetGeneratedOptPhotons(int n)          {_generatedPhotons=n;}
    
    G4ThreeVector          GetHistBinCenter(int binx, int biny, int binz); 
    double                 GetHistBinWidthX(int binx); 
    double                 GetHistBinWidthY(int biny); 
    double                 GetHistBinWidthZ(int binz); 
    void                   doStoreConstants(bool storeConstants) {_storeConstants = storeConstants;}

  private:

    static WLSEventAction*  _fgInstance;  
    TH1D*                   _histP[2][4];
    TH1D*                   _histT[2][4];
    TH3D*                   _histSurvivalProb[4][4];
    TH3D*                   _histTimeDifference[4][4];
    TH3D*                   _histFiberEmissions[4][4];
    G4ThreeVector           _start;
    int                     _generatedPhotons;
    int                     _mode;
    TFile*                  _fileLookupTable;
    bool                    _storeConstants;

    void                    Draw(const G4Event* evt) const;

    TH2D                    *_photonsVsIntegral, *_photonsVsPulseHeight;
    TH2D                    *_PEsVsIntegral, *_PEsVsPulseHeight;
};

#endif

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

    G4int GetEventNo();
    static WLSEventAction* Instance()                      {return _fgInstance;}
    TH1D*                  GetHistPE(int mode, int SiPM)   {return _histPE[mode][SiPM];}
    TH1D*                  GetHistT(int mode, int SiPM)    {return _histT[mode][SiPM];}
    void                   SetOptPhotonStart(G4ThreeVector start) {_start=start;}
    void                   SetGeneratedOptPhotons(int n)          {_generatedPhotons=n;}
    
    G4ThreeVector          GetHistBinCenter(int binx, int biny, int binz); 
    double                 GetHistBinWidthX(int binx); 
    double                 GetHistBinWidthY(int biny); 
    double                 GetHistBinWidthZ(int binz); 
    void                   doStoreConstants(bool storeConstants) {_storeConstants = storeConstants;}

  private:

    static WLSEventAction*  _fgInstance;  
    TH1D***                 _histPE;
    TH1D***                 _histT;
    TH3D***                 _histSurvivalProb;
    TH3D***                 _histTimeDifference;
    TH3D***                 _histFiberEmissions;
    G4ThreeVector           _start;
    int                     _generatedPhotons;
    int                     _mode;
    TFile*                  _fileLookupTable;
    bool                    _storeConstants;

    void                    Draw();
    void                    Draw(const G4Event* evt) const;
};

#endif

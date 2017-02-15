#ifndef WLSEventAction_h
#define WLSEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

class TH1D;
class TH2D;
class TH3D;
class TFile;
class TNtuple;

class WLSEventAction : public G4UserEventAction
{
    WLSEventAction();

  public:

    WLSEventAction(int mode, const std::string &singlePEWaveformFilename, int numberOfPhotons=-1, int simType=-1, int minBin=-1, bool verbose=false);  
                                                                                                                 //numberOfPhotons, simType, minBin, verbose is only needed for mode -1
    ~WLSEventAction();

  public:

    void   BeginOfEventAction(const G4Event*);
    void     EndOfEventAction(const G4Event*);

    static WLSEventAction* Instance()                      {return _fgInstance;}
    void                   SetGeneratedOptPhotons(int n)   {_generatedPhotons=n;}
    void                   SetStartZ(double startZ)        {_startZ=startZ;}

  private:

    static WLSEventAction*  _fgInstance;  
    TH1D*                   _histP[2][4];
    TH1D*                   _histT[2][4];
    TH1D*                   _histPE[4];
    TH3D*                   _histSurvivalProb[4][4];
    TH3D*                   _histTimeDifference[4][4];
    TH3D*                   _histFiberEmissions[4][4];
    TNtuple*                _ntuple;
    int                     _generatedPhotons;
    int                     _mode, _numberOfPhotons, _simType, _minBin;
    std::string             _singlePEWaveformFilename;
    bool                    _verbose;
    bool                    _storeConstants;
    double                  _startZ;

    void                    Draw(const G4Event* evt);

    std::vector<double>     _PEs[4];
    std::vector<double>     _recoPEs[4];
};

#endif

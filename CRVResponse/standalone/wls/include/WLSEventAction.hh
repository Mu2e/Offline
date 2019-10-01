#ifndef WLSEventAction_h
#define WLSEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include <vector>

#include "WLSSteppingAction.hh"

class TH1D;
class TFile;
class TNtuple;

class WLSEventAction : public G4UserEventAction
{
    WLSEventAction();

  public:

    WLSEventAction(WLSSteppingAction::simulationMode mode, const std::string &singlePEWaveformFilename, const std::string &photonMapFilename,  
                   int numberOfPhotons=-1, int simType=-1, unsigned int minBin=0, bool verbose=false);  
                   //numberOfPhotons, simType, minBin, currentBin, verbose is only needed for simulationMode::CreateLookupTables
    ~WLSEventAction();

  public:

    void   BeginOfEventAction(const G4Event*);
    void     EndOfEventAction(const G4Event*);

    static WLSEventAction* Instance()                      {return _fgInstance;}
    void                   SetGeneratedOptPhotons(int n)   {_generatedPhotons=n;}
    void                   SetBinNumber(unsigned int n)    {_currentBin=n;}

  private:

    static WLSEventAction*  _fgInstance;  

    WLSSteppingAction::simulationMode _mode;

    TH1D*                   _histP[2][4];
    TH1D*                   _histT[2][4];
    TH1D*                   _histPE[4];
    TNtuple*                _ntuple;
    int                     _numberOfPhotons, _simType; 
    unsigned int            _minBin, _currentBin;
    std::string             _singlePEWaveformFilename;
    std::string             _photonMapFilename;
    bool                    _verbose;

    int                     _generatedPhotons;  //set by WLSPrimaryGeneratorAction
    double                  _startZ;            //set by WLSPrimaryGeneratorAction

    void                    Draw(const G4Event* evt);

    std::vector<double>     _PEs[4];
    std::vector<double>     _recoPEs[4];
    std::vector<double>     _pulseTimes[4];
    std::vector<double>     _LETimes[4];
    std::vector<double>     _pulseWidths[4];
};

#endif

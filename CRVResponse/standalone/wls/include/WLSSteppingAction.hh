#ifndef WLSSteppingAction_h
#define WLSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include <vector>
#include <map>
#include <set>

#include "CLHEP/Random/Randomize.h"

class TFile;
class TNtuple;

namespace mu2eCrv
{
  class MakeCrvPhotons;
}

class WLSSteppingAction : public G4UserSteppingAction
{
  public:

    enum simulationMode {CreateLookupTables, UseGeantOnly, UseGeantAndLookupTables, Undefined};

    WLSSteppingAction(simulationMode mode, const std::string &lookupFileName = "", const std::string &visibleEnergyAdjustmentFileName = "");
                                                                                 //lookupFileName and visibleEnergyAdjustmentFileName
                                                                                 //only used for simulationMode::UseGeantAndLookupTables
    ~WLSSteppingAction();

    void                      UserSteppingAction(const G4Step*);
    static WLSSteppingAction* Instance() {return _fgInstance;}
    void                      Reset();
    const std::vector<double> &GetArrivalTimes(int SiPM);
    const std::vector<double> &GetArrivalTimesFromLookupTables(int SiPM);
    const std::vector<int>    &GetFiberEmissions(int SiPM);

  private:

    std::unique_ptr<mu2eCrv::MakeCrvPhotons> _crvPhotons;
    static WLSSteppingAction *_fgInstance;  
    std::vector<double>       _arrivalTimes[4];
    std::vector<double>       _arrivalTimesFromLookupTables[4];
    std::vector<int>          _fiberEmissions[4];
    simulationMode            _mode;

    std::map<int,int>         _wlsTrackParents;

    CLHEP::HepJamesRandom     _engine;
    CLHEP::RandFlat           _randFlat;
    CLHEP::RandGaussQ         _randGaussQ;
    CLHEP::RandPoissonQ       _randPoissonQ;

    void                      ShowVisibleEnergyTable(const G4Step *theStep);

    TNtuple*                  _ntuple;  //WLS fiber test
};

#endif

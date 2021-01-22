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

    struct PhotonInfo
    {
      double _arrivalTime;
      int _fiberEmissions;
      PhotonInfo(double arrivalTime, int fiberEmissions) : _arrivalTime(arrivalTime), _fiberEmissions(fiberEmissions) {}
    };

    WLSSteppingAction(simulationMode mode, const std::string &lookupFileName = ""); //lookupFileName and visibleEnergyAdjustmentFileName
                                                                                    //only used for simulationMode::UseGeantAndLookupTables
    ~WLSSteppingAction();

    void                      UserSteppingAction(const G4Step*);
    static WLSSteppingAction* Instance() {return _fgInstance;}
    void                      Reset();
    const std::vector<PhotonInfo> &GetPhotonInfo(int SiPM);
    const std::vector<double> &GetArrivalTimesFromLookupTables(int SiPM);

    void                      PrintFiberStats();

  private:

    std::unique_ptr<mu2eCrv::MakeCrvPhotons> _crvPhotons;
    static WLSSteppingAction *_fgInstance;  
    std::vector<PhotonInfo>   _photonInfo[4];
    std::vector<double>       _arrivalTimesFromLookupTables[4];
    simulationMode            _mode;

    std::map<int,int>         _wlsTrackParents;
    std::set<int>             _tracksHittingFiber;
    std::set<int>             _tracksGettingAbsorbedInFiber;
    int                       _zeroFiberEmissions;

    CLHEP::HepJamesRandom     _engine;
    CLHEP::RandFlat           _randFlat;
    CLHEP::RandGaussQ         _randGaussQ;
    CLHEP::RandPoissonQ       _randPoissonQ;

    TNtuple*                  _ntuple;  //WLS fiber test
};

#endif

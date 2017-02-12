#ifndef WLSSteppingAction_h
#define WLSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include <vector>
#include <map>

#include "CLHEP/Random/Randomize.h"

class TFile;
class TH3D;
class TNtuple;

namespace mu2eCrv
{
  class MakeCrvPhotonArrivals;
}

class WLSSteppingAction : public G4UserSteppingAction
{
  public:

    WLSSteppingAction(int mode, const std::string &lookupFileName = "");
    ~WLSSteppingAction();

    void                      UserSteppingAction(const G4Step*);
    static WLSSteppingAction* Instance() {return _fgInstance;}
    void                      Reset();
    const std::vector<double> &GetArrivalTimes(int i, int SiPM);
    const std::vector<int>    &GetFiberEmissions(int SiPM);

  private:

    std::unique_ptr<mu2eCrv::MakeCrvPhotonArrivals> _crvPhotonArrivals;
    static WLSSteppingAction *_fgInstance;  
    std::vector<double>       _arrivalTimes[2][4];
    std::vector<int>          _fiberEmissions[4];
    int                       _mode;

    std::map<int,int>         _wlsTracks;

    CLHEP::HepJamesRandom     _engine;
    CLHEP::RandFlat           _randFlat;
    CLHEP::RandGaussQ         _randGaussQ;
    CLHEP::RandPoissonQ       _randPoissonQ;

    void                      Test(const G4Step *theStep, int PDGcode);

    TNtuple*                  _ntuple;  //WLS fiber test
};

#endif

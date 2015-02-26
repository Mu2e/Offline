#ifndef MakeCrvPhotonArrivals_h
#define MakeCrvPhotonArrivals_h

#include <vector>
#include <map>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/Randomize.h"


class TFile;
class TH3D;

class MakeCrvPhotonArrivals
{
  public:

    MakeCrvPhotonArrivals(CLHEP::RandFlat &randFlat) : _actualHalfLength(NAN), _fileLookupTable(NULL), _randFlat(randFlat) {}

    ~MakeCrvPhotonArrivals();

    void                      LoadLookupTable(std::string filename);
    void                      MakePhotons(const CLHEP::Hep3Vector &stepStart,   //they need to be points
                                      const CLHEP::Hep3Vector &stepEnd,         //local to the CRV bar
                                      double timeStart, double timeEnd,
                                      int PDGcode, double beta, double charge,
                                      double energyDepositedTotal,
                                      double energyDepositedNonIonizing);
    int                       GetNumberOfPhotons(int SiPM);
    const std::vector<double> &GetArrivalTimes(int SiPM);
    void                      SetScintillationYield(double scintillationYield) {_scintillationYield=scintillationYield;}
    void                      SetScintillatorDecayTimeFast(double decayTime) {_scintillatorDecayTimeFast=decayTime;}
    void                      SetScintillatorDecayTimeSlow(double decayTime) {_scintillatorDecayTimeSlow=decayTime;}
    void                      SetFiberDecayTime(double decayTime) {_fiberDecayTime=decayTime;}

    void                      SetActualHalfLength(double actualHalfLength) {_actualHalfLength=actualHalfLength;}

  private:

//    int                       _photonArrivals[4];
    std::vector<double>       _arrivalTimes[4];
    double                    _scintillationYield;
    double                    _scintillatorDecayTimeFast, _scintillatorDecayTimeSlow; 
    double                    _fiberDecayTime;

    double                    _actualHalfLength;

    TFile*                    _fileLookupTable;
    TH3D***                   _histSurvivalProb;
    TH3D***                   _histTimeDifference;
    TH3D***                   _histFiberEmissions;
    float                     _halfThickness;
    float                     _halfWidth;
    float                     _halfLength;
    float                     _speedOfLightFiber;
    float                     _cerenkovRindex;
    float                     _cerenkovEinterval;
    float                     _ratioFastSlow;
    float                     _scintillatorDensity;
    float                     _scintillatorBirksConstant;
    float                     _fiberSeparation;
    float                     _holeRadius;
    float                     _fiberRadius;

    CLHEP::RandFlat           &_randFlat;

    bool   IsInsideScintillator(const CLHEP::Hep3Vector &p);
    int    IsInsideFiber(const CLHEP::Hep3Vector &p);
    double GetRandomTime(TH3D *timeDifference, double y, double z);
    int    GetRandomFiberEmissions(TH3D *fiberEmissions, double y, double z);
    double GetAverageNumberOfCerenkovPhotons(double beta, double charge);
    void   AdjustPosition(CLHEP::Hep3Vector &p);
    double VisibleEnergyDeposition(int PDGcode, double stepLength,
                                   double energyDepositedTotal,
                                   double energyDepositedNonIonizing);
};

#endif

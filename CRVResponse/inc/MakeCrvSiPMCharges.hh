/*
Author: Ralf Ehrlich
Based on Paul Rubinov's C# code
*/

#ifndef MakeCrvSiPMCharges_hh
#define MakeCrvSiPMCharges_hh

#include <memory>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include "CLHEP/Random/Randomize.h"

#include <TFile.h>
#include <TH2F.h>

namespace mu2eCrv
{

  struct Pixel
  {
    bool   _discharged; 
    double _t;          //time of last discharge (if _discharged is true), or NAN (if _discharged is false)
    Pixel() : _discharged(false), _t(NAN) {}
  };

  struct SiPMresponse
  {
    double _time;
    double _charge;       //in C
    double _chargeInPEs;  //in PEs
    size_t _photonIndex;  //index in the original photon vector
    bool   _darkNoise;
    SiPMresponse(double time, double charge, double chargeInPEs, size_t photonIndex, bool darkNoise) : 
                  _time(time), _charge(charge), _chargeInPEs(chargeInPEs), _photonIndex(photonIndex), _darkNoise(darkNoise) {}
  };

  struct ScheduledCharge
  {
    std::pair<int,int>  _pixelId;
    double              _time;
    size_t              _photonIndex; //index in the original photon vector
    bool                _darkNoise;   //this charge is dark noise and was not created by an "outside photon"
    ScheduledCharge(const std::pair<int,int> &pixelId, double time, size_t photonIndex, bool darkNoise) : 
                  _pixelId(pixelId), _time(time), _photonIndex(photonIndex), _darkNoise(darkNoise) {}
    bool operator<(const ScheduledCharge &r) const
    {
      return _time < r._time;
    };
    private:
    ScheduledCharge();
  };
  
  class MakeCrvSiPMCharges
  {
    int    _nPixelsX;
    int    _nPixelsY;
    double _overvoltage;      //in V  (operating overvoltage = bias voltage - breakdown voltage)
    double _blindTime;        //in ns
    double _microBunchPeriod; //in ns
    double _timeConstant;     //in ns
    double _capacitance;      //in F

    public:
    struct ProbabilitiesStruct
    {
      double _avalancheProbParam1;
      double _avalancheProbParam2;
      double _trapType0Prob;
      double _trapType1Prob;
      double _trapType0Lifetime;
      double _trapType1Lifetime;
      double _thermalRate;  //in ns^-1
      double _crossTalkProb; 
    };

    private:
    ProbabilitiesStruct                _probabilities;
    std::vector<std::pair<int,int> >   _inactivePixels;

    std::map<std::pair<int,int>,Pixel> _pixels;
    std::multiset<ScheduledCharge>     _scheduledCharges;

    std::vector<std::pair<int,int> > FindCrossTalkPixelIds(const std::pair<int,int> &pixelId);
    std::pair<int,int>               FindThermalNoisePixelId();
    std::pair<int,int>               FindFiberPhotonsPixelId();
    bool                             IsInactivePixelId(const std::pair<int,int> &pixelId);

    double GetAvalancheProbability(double v);
    double GenerateAvalanche(Pixel &pixel, const std::pair<int,int> &pixelId, double time, size_t photonIndex, bool darkNoise);
    double GetVoltage(const Pixel &pixel, double time);
    void   FillPhotonQueue(const std::vector<std::pair<double,size_t> > &photons);

    CLHEP::RandFlat     &_randFlat;
    CLHEP::RandPoissonQ &_randPoissonQ;
    double               _avalancheProbFullyChargedPixel;

    TFile *_photonMapFile;
    TH2F  *_photonMap;

    public:

    MakeCrvSiPMCharges(CLHEP::RandFlat &randFlat, CLHEP::RandPoissonQ &randPoissonQ, const std::string &photonMapFileName);
    ~MakeCrvSiPMCharges() {_photonMapFile->Close();}

    void SetSiPMConstants(int nPixelsX, int nPixelsY, double overvoltage, 
                          double blindTime, double microBunchPeriod, double timeConstant, 
                          double capacitance, ProbabilitiesStruct probabilities,
                          const std::vector<std::pair<int,int> > &inactivePixels);
    void Simulate(const std::vector<std::pair<double,size_t> > &photons, 
                  std::vector<SiPMresponse> &SiPMresponseVector);
  };

}

#endif

/*
Author: Ralf Ehrlich
Based on Paul Rubinov's C# code
*/

#ifndef MakeCrvSiPMResponses_hh
#define MakeCrvSiPMResponses_hh

#include <memory>
#include <set>
#include <map>
#include <vector>
#include "CLHEP/Random/Randomize.h"

namespace mu2eCrv
{

  struct Pixel
  {
    double _v;
    double _t;
    Pixel(double bias, double time) : _v(bias), _t(time) {}

    private:
    Pixel();  //disables the default constructor
  };

  struct SiPMresponse
  {
    double _time;
    double _charge;
    SiPMresponse(double time, double charge) : _time(time), _charge(charge) {}
  };

  struct ScheduledCharge
  {
    int    _cellid;
    double _time;
    ScheduledCharge(int cellid, double time) : _cellid(cellid), _time(time) {}
    bool operator<(const ScheduledCharge &r) const
    {
      return _time < r._time;
    };
    private:
    ScheduledCharge();
  };
  
  class MakeCrvSiPMResponses
  {
    int    _numberPixels;
    int    _numberPixelsAtFiber;
    double _bias;             //in V above breakdown
    double _blindTime;        //in ns
    double _microBunchPeriod; //in ns
    double _timeConstant;     //in ns

    double _time;             //in ns

    public:
    struct ProbabilitiesStruct
    {
      double _constGeigerProbCoef;
      double _constGeigerProbVoltScale;
      double _constTrapType0Prob;  //per unit voltage
      double _constTrapType1Prob;  //per unit voltage
      double _constTrapType0Lifetime;
      double _constTrapType1Lifetime;
      double _constThermalProb;         //per unit time
      double _constPhotonProduction;    //in 1/fC
    };

    private:
    ProbabilitiesStruct            _probabilities;

    std::map<int,Pixel>            _pixels;
    std::multiset<ScheduledCharge> _scheduledCharges;

    double GenerateAvalanche(Pixel &pixel, int cellid);
    void   RechargeCell(Pixel &pixel);
    void   FillPhotonQueue(const std::vector<double> &photons);

    CLHEP::RandFlat     &_randFlat;
    CLHEP::RandPoissonQ &_randPoissonQ;

    public:

    MakeCrvSiPMResponses(CLHEP::RandFlat &randFlat, CLHEP::RandPoissonQ &randPoissonQ) :
                         _randFlat(randFlat), _randPoissonQ(randPoissonQ) {}

    void SetSiPMConstants(int numberPixels, int numberPixelsAtFiber, double bias, 
                          double blindTime, double microBunchPeriod, double timeConstant, 
                          ProbabilitiesStruct probabilities);
    void Simulate(const std::vector<double> &photons, 
                  std::vector<SiPMresponse> &SiPMresponseVector);
  };

}

#endif

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

  struct Pixel
  {
    double _v;
    double _t;
    Pixel(double bias) : _v(bias), _t(0) {}

    private:
    Pixel();  //disables the default constructor
  };

  struct SiPMresponse
  {
    double _time;
    double _charge;
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
    double _bias;           //in V above breakdown
    double _timeStart;      //in ns
    double _timeEnd;        //in ns
    double _scaleFactor;    //based on a time step of 1.0ns

    double _time;           //in ns

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

    public:
    void SetSiPMConstants(double numberPixels, double bias, 
                          double timeStart, double timeEnd, double scaleFactor, 
                          ProbabilitiesStruct probabilities);
    void Simulate(const std::vector<double> &photons, 
                  std::vector<SiPMresponse> &SiPMresponseVector);
  };

#endif

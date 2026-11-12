/*
Author: Ralf Ehrlich
Based on Paul Rubinov's C# code
*/

#include "Offline/CRVResponse/inc/MakeCrvSiPMCharges.hh"

#include <cmath>
#include <iostream>
#include <algorithm>

//photon map gets created from the CRVPhoton.root file, which can be generated with the standalone program (in WLSSteppingAction)
//in ROOT: CRVPhotons->Draw("x/0.05+20:(fabs(y)-13)/0.05+20>>photonMap(40,0,40,40,0,40)","","COLZ")

//to get standalone version: compile with
//g++ MakeCrvSiPMCharges.cc -std=c++11 -I../../ -I$CLHEP_INCLUDE_DIR -L$CLHEP_LIB_DIR -lCLHEP -DSiPMChargesStandalone `root-config --cflags --glibs`

namespace mu2eCrv
{

double MakeCrvSiPMCharges::GetAvalancheProbability(double v)
{
  double avalancheProbability = _probabilities._avalancheProbParam1*(1 - exp(-v/_probabilities._avalancheProbParam2));
  return avalancheProbability;
}

std::pair<int,int> MakeCrvSiPMCharges::FindThermalNoisePixelId()
{
  int x=_randFlat.fire(_nPixelsX);
  int y=_randFlat.fire(_nPixelsY);
  return std::pair<int,int>(x,y);
}

std::pair<int,int> MakeCrvSiPMCharges::FindFiberPhotonsPixelId()
{
  double x,y;
  _photonMap->GetRandom2(x,y);
  return std::pair<int,int>(static_cast<int>(x),static_cast<int>(y));
}

double MakeCrvSiPMCharges::GenerateAvalanche(Pixel &pixel, const std::pair<int,int> &pixelId, double time, size_t photonIndex, bool darkNoise)
{
  double v = GetVoltage(pixel,time);

  if(_randFlat.fire() < GetAvalancheProbability(v))
  {
    //after pulses
    //for simplicity, it is assumed that all pixels are fully charged
    //the _trapType0Prob and _trapType1Prob are measured probabilities (from the Hamamatsu specs), however the actually production probabilities are higher, but are reduced by the avalanche probability
    //(measured probability = production probability * avalanche probability)
    //the production probabilities are needed here
    if(_randFlat.fire() < _probabilities._trapType0Prob/_avalancheProbFullyChargedPixel)
    {
      //create new Type0 trap (fast)
      double traptime = -_probabilities._trapType0Lifetime * log(_randFlat.fire());

      auto it = std::upper_bound(_scheduledCharges.begin()+_currentChargeIndex, _scheduledCharges.end(), time+traptime,
                                 [](double a, const ScheduledCharge& b) {return a < b._time;});
      _scheduledCharges.emplace(it, pixelId,time + traptime,photonIndex,darkNoise);
    }

    if(_randFlat.fire() < _probabilities._trapType1Prob/_avalancheProbFullyChargedPixel)
    {
      //create new Type1 trap (slow)
      double traptime = -_probabilities._trapType1Lifetime * log(_randFlat.fire());
      auto it = std::upper_bound(_scheduledCharges.begin()+_currentChargeIndex, _scheduledCharges.end(), time+traptime,
                                 [](double a, const ScheduledCharge& b) {return a < b._time;});
      _scheduledCharges.emplace(it, pixelId,time + traptime,photonIndex,darkNoise);
    }

    //cross talk can happen in all 4 neighboring pixels (distribute photons there for possible avalanches)
    //cross talk happens at the same time. so it can be inserted at the current index at the scheduledCharge vector
    if(pixelId.first>0)
    {
      if(_randFlat.fire() < _crossTalkProbabilitySinglePixel)
           _scheduledCharges.emplace(_scheduledCharges.begin()+_currentChargeIndex,std::pair<int,int>(pixelId.first-1,pixelId.second),time,photonIndex,darkNoise);
    }
    if(pixelId.first+1<_nPixelsX)
    {
      if(_randFlat.fire() < _crossTalkProbabilitySinglePixel)
           _scheduledCharges.emplace(_scheduledCharges.begin()+_currentChargeIndex,std::pair<int,int>(pixelId.first+1,pixelId.second),time,photonIndex,darkNoise);
    }
    if(pixelId.second>0)
    {
      if(_randFlat.fire() < _crossTalkProbabilitySinglePixel)
           _scheduledCharges.emplace(_scheduledCharges.begin()+_currentChargeIndex,std::pair<int,int>(pixelId.first,pixelId.second-1),time,photonIndex,darkNoise);
    }
    if(pixelId.second+1<_nPixelsY)
    {
      if(_randFlat.fire() < _crossTalkProbabilitySinglePixel)
           _scheduledCharges.emplace(_scheduledCharges.begin()+_currentChargeIndex,std::pair<int,int>(pixelId.first,pixelId.second+1),time,photonIndex,darkNoise);
    }

    //the pixel's overvoltage becomes 0, i.e. the pixel's voltage gets reduced to the breakdown voltage
    //the time when this happens gets recorded in the pixel variable
    pixel._discharged=true;
    pixel._t=time;
    pixel._counter=_counter;  //explanation for _counter is in header file

    double outputCharge = _capacitance*v;   //output charge = capacitance (of one pixel) * overvoltage
                                            //gain = outputCharge / elementary charge
    return outputCharge;
  }
  else return 0;  //no avalanche means no output charge
}

double MakeCrvSiPMCharges::GetVoltage(const Pixel &pixel, double time)
{
  if(!pixel._discharged || pixel._counter!=_counter) return _overvoltage;  //an old counter value indicates that it was discharged during the last call of Simulate()
                                                                           //and can now be considered to be fully charged

  double deltaT = time - pixel._t;   //time since last discharge
  double v = _overvoltage * (1.0-exp(-deltaT/_timeConstant));
  return v;
}

void MakeCrvSiPMCharges::SetSiPMConstants(int nPixelsX, int nPixelsY, double overvoltage, double timeConstant,
                                            double capacitance, ProbabilitiesStruct probabilities,
                                            const std::vector<std::pair<int,int> > &inactivePixels)
{
  _nPixelsX = nPixelsX;
  _nPixelsY = nPixelsY;
  _overvoltage = overvoltage;   //operating overvoltage = bias voltage - breakdown voltage
  _timeConstant = timeConstant;
  _capacitance = capacitance;  //capacitance per pixel
  _probabilities = probabilities;

  _avalancheProbFullyChargedPixel = GetAvalancheProbability(overvoltage);

  _pixels.resize(_nPixelsX*_nPixelsY);
  for(int pixelIndex=0; pixelIndex<_nPixelsX*_nPixelsY; ++pixelIndex)
  {
    for(size_t i=0; i<inactivePixels.size(); i++)
    {
      if(pixelIndex==inactivePixels[i].first*_nPixelsY+inactivePixels[i].second) _pixels[pixelIndex]._inactive=true;
    }
  }

  double probabilityNoCrossTalk = 1.0-_probabilities._crossTalkProb;              //prob that cross talk does not occur = 1 - prob that cross talk occurs
  double probabilityNoCrossTalkSinglePixel = pow(probabilityNoCrossTalk,1.0/4.0); //prob that cross talk does not occur at any of the 4 neighboring pixels
                                                                                  //=pow(prob that cross talk does not occur at a pixel,4)
  _crossTalkProbabilitySinglePixel = 1.0-probabilityNoCrossTalkSinglePixel;
  //the crossTalkProbabilitySinglePixel is the measured probability (based on the _crossTalkProb from the Hamamatsu specs),
  //however the actually production probability is higher, but is reduced by the avalanche probability
  //(measured probability = production probability * avalanche probability)
  //the production probability is needed here
  _crossTalkProbabilitySinglePixel /= _avalancheProbFullyChargedPixel;
}

void MakeCrvSiPMCharges::FillQueue(const std::vector<std::pair<double,size_t> > &photons, double startTime, double endTime)
{
//schedule charges caused by the CRV counter photons
  for(size_t i=0; i<photons.size(); i++)
  {
    std::pair<int,int> pixelId = FindFiberPhotonsPixelId();  //only pixels at fiber
    if(pixelId.first<0 || pixelId.second<0 || pixelId.first>=_nPixelsX || pixelId.second>=_nPixelsY) throw std::logic_error("Used pixel outside of SiPM.");
    _scheduledCharges.emplace_back(pixelId, photons[i].first, photons[i].second, false);
  }

//schedule random thermal charges
  double timeWindow = endTime-startTime;

  //for the dark noise simulation, it is assumed that all pixels are fully charged (for simplicity)

  //the actual thermal production rate gets scaled down by the avalanche probability,
  //i.e. only a fraction of the thermaly created charges lead to a dark noise pulse
  //thermal production rate * avalanche probability = thermal rate
  double thermalProductionRate = _probabilities._thermalRate/_avalancheProbFullyChargedPixel;

  //average number of thermaly created charges is thermalProductionRate*timeWindow
  int numberThermalCharges = _randPoissonQ.fire(thermalProductionRate * timeWindow);
  for(int i=0; i<numberThermalCharges; i++)
  {
    std::pair<int,int> pixelId = FindThermalNoisePixelId();  //all pixels
    double time = startTime + timeWindow * _randFlat.fire();
    _scheduledCharges.emplace_back(pixelId, time, 0, true);
  }
}

void MakeCrvSiPMCharges::Simulate(const std::vector<std::pair<double,size_t> > &photons,   //pair of photon time and index in the original photon vector
                                   std::vector<SiPMresponse> &SiPMresponseVector, double startTime, double endTime)
{
  ++_counter; //see header file for explanation

  _scheduledCharges.clear();
  FillQueue(photons, startTime, endTime);

  //sort the charges by time, all charges that get inserted later (cross-talk, after pulses) need to be inserted into the correct location
  std::stable_sort(_scheduledCharges.begin(), _scheduledCharges.end(),
  [](const ScheduledCharge& a, const ScheduledCharge& b) {return a._time < b._time;});

  _currentChargeIndex=0;

  while(1)
  {
    if(_currentChargeIndex>=_scheduledCharges.size()) break;  //no more scheduled charges

    //find the first charge
    const ScheduledCharge &currentCharge = _scheduledCharges[_currentChargeIndex];

    std::pair<int,int> pixelId = currentCharge._pixelId;
    int pixelIndex = pixelId.first*_nPixelsY+pixelId.second;

    double time = currentCharge._time;
    size_t photonIndex = currentCharge._photonIndex;
    bool darkNoise = currentCharge._darkNoise;
    ++_currentChargeIndex;

    if(_pixels.at(pixelIndex)._inactive) continue;

    if(time>endTime) continue; //this is relevant for afterpulses

    Pixel &pixel = _pixels.at(pixelIndex);
    double outputCharge = GenerateAvalanche(pixel, pixelId, time, photonIndex, darkNoise);   //the output charge (in Coulomb) of the pixel due to the avalanche
    double outputChargeInPEs = (outputCharge/_capacitance)/_overvoltage;                     //the output charge in units of single PEs of a fully charges pixel

    if(outputCharge>0) SiPMresponseVector.emplace_back(time, outputCharge, outputChargeInPEs, photonIndex, darkNoise);
  } //while(1)
}

MakeCrvSiPMCharges::MakeCrvSiPMCharges(CLHEP::RandFlat &randFlat, CLHEP::RandPoissonQ &randPoissonQ, const std::string &photonMapFileName) :
                                       _randFlat(randFlat), _randPoissonQ(randPoissonQ), _avalancheProbFullyChargedPixel(0)
{
  _photonMapFile = new TFile(photonMapFileName.c_str());
  if(!_photonMapFile || _photonMapFile->IsZombie()) throw std::logic_error("Could not open photon map file.");
  _photonMap = (TH2F*)_photonMapFile->FindObjectAny("photonMap");
  if(!_photonMap) throw std::logic_error("Could not find photon map.");
}

}

//sample program

#ifdef SiPMChargesStandalone
int main()
{
  std::vector<std::pair<double,size_t> > photonTimes;
  photonTimes.emplace_back(650,0);
  photonTimes.emplace_back(620,1);
  for(int i=0; i<100; i++) photonTimes.emplace_back(630,i+2);
  for(int i=0; i<100; i++) photonTimes.emplace_back(623,i+102);
  photonTimes.emplace_back(656,202);
  photonTimes.emplace_back(612,203);
  std::vector<mu2eCrv::SiPMresponse> SiPMresponseVector;

  mu2eCrv::MakeCrvSiPMCharges::ProbabilitiesStruct probabilities;
  probabilities._avalancheProbParam1 = 0.65;
  probabilities._avalancheProbParam2 = 2.7;
  probabilities._trapType0Prob = 0.0;
  probabilities._trapType1Prob = 0.0;
  probabilities._trapType0Lifetime = 5;
  probabilities._trapType1Lifetime = 50;
  probabilities._thermalRate = 1.0e-4;
  probabilities._crossTalkProb = 0.05;

  std::vector<std::pair<int,int> > inactivePixels = { {18,18}, {18,19}, {18,20}, {18,21},
                                                      {19,18}, {19,19}, {19,20}, {19,21},
                                                      {20,18}, {20,19}, {20,20}, {20,21},
                                                      {21,18}, {21,19}, {21,20}, {21,21} };

  CLHEP::HepJamesRandom engine(1);
  CLHEP::RandFlat randFlat(engine);
  CLHEP::RandPoissonQ randPoissonQ(engine);
  mu2eCrv::MakeCrvSiPMCharges sim(randFlat,randPoissonQ,"/cvmfs/mu2e.opensciencegrid.org/DataFiles/CRVConditions/v6_0/photonMap.root");
  sim.SetSiPMConstants(40, 40, 3.0, 13.3, 8.84e-14, probabilities, inactivePixels);

  sim.Simulate(photonTimes, SiPMresponseVector, 500, 1695);

  for(unsigned int i=0; i<SiPMresponseVector.size(); i++)
  {
    std::cout<<i<<"   "<<SiPMresponseVector[i]._time<<"   "<<SiPMresponseVector[i]._charge<<"   "<<SiPMresponseVector[i]._chargeInPEs;
    std::cout<<"     "<<(SiPMresponseVector[i]._darkNoise?"dark Noise":std::to_string(SiPMresponseVector[i]._photonIndex).c_str())<<std::endl;
  }

  return 0;
}

#endif

/*
Author: Ralf Ehrlich
Based on Paul Rubinov's C# code
*/

#include "MakeCrvSiPMResponses.hh"

#include <math.h>
#include <iostream>
#include <algorithm>


//to get standalone version: compile with
//g++ MakeCrvSiPMResponses.cc -std=c++11 -I../inc -I$CLHEP_INCLUDE_DIR -L$CLHEP_LIB_DIR -lCLHEP -DSiPMResponseStandalone

    // the idea here is that a pixel goes through the phases,
    // these next steps happen inside the cell itself:
    //   1) estimate if there is a thermaly generated free charge
    //   2) check if a trap releases a free charge
    //   3) for each free charge, check for avalanche by comparing to the Geiger prob (GP)
    //   4) check for avalanche:
    //   if there is an avelanche:
    //     5) compute size of avalanche based on OV
    //     6) release appropriate number of photons
    //     7) add appropriate traps
    //     8) reduce the overvoltage (OV) to zero
    //   else
    //     9) compute the recharge current
    //     10) increase the OV by the step charge Q

double MakeCrvSiPMResponses::GenerateAvalanche(Pixel &pixel, int cellid)
{
  double v = pixel._v;

  double GeigerProb = 1 - exp(-pow(v, _probabilities._constGeigerProbCoef)/_probabilities._constGeigerProbVoltScale);

  if(_randFlat.fire() < GeigerProb)
  {
    int trapsType0 = _randPoissonQ.fire(_probabilities._constTrapType0Prob/_bias * v);
    int trapsType1 = _randPoissonQ.fire(_probabilities._constTrapType1Prob/_bias * v);
    int photons    = _randPoissonQ.fire(_probabilities._constPhotonProduction/_bias * v);
  
    double trap0Lifetime = _probabilities._constTrapType0Lifetime;
    double trap1Lifetime = _probabilities._constTrapType1Lifetime;
    for(int i=0; i<trapsType0; i++)
    {
      //create new Type0 trap (fast)
      double life = -trap0Lifetime * log(_randFlat.fire());
      _scheduledCharges.emplace(cellid,_time + life); //constructs ScheduledCharge(cellid,_time+life)
    }
    for(int i=0; i<trapsType1; i++)
    {
      //create new Type1 trap (slow)
      double life = -trap1Lifetime * log(_randFlat.fire());
      _scheduledCharges.emplace(cellid,_time + life); //constructs ScheduledCharge(cellid,_time+life)
    }

    //cross talk (distribute photons over all pixels leading to release of free charges there)
    for(int i=0; i<photons; i++)
    {
      int cellidCrossTalk = _randFlat.fireInt(_numberPixels);
      _scheduledCharges.emplace_hint(_scheduledCharges.begin(),cellidCrossTalk,_time); 
                                                      //constructs ScheduledCharge(cellidCrossTalk,time)
    }

    //zero the _v
    pixel._v=0;
    pixel._t=_time;

    return v;  //avalanche voltage
  }
  else return 0;
}

void MakeCrvSiPMResponses::RechargeCell(Pixel &pixel)
{
  double v = pixel._v;

  if(v < _bias)
  {
    double timeStep = _time - pixel._t;
    v = pow(1-_scaleFactor, timeStep) * (v-_bias) + _bias;
// formula is based on the following iterative sequence 
// which gets called at every _timeStep = 1.0ns = const
//
// double rechargeCurrent = (SiPM._bias - _v) * SiPM._scaleFactor;
// _v += SiPM._timeStep * rechargeCurrent;

    if(v > _bias) { v = _bias; }

    pixel._v=v;
    pixel._t=_time;
  }
}

void MakeCrvSiPMResponses::SetSiPMConstants(double numberPixels, double bias,  
                                            double blindTime, double microBunchPeriod, double scaleFactor, 
                                            ProbabilitiesStruct probabilities)
{
  _numberPixels = numberPixels;
  _bias = bias;
  _blindTime = blindTime;
  _microBunchPeriod = microBunchPeriod;
  _scaleFactor = scaleFactor;
  _probabilities = probabilities;
}

    // for each time step, the flow is like this:
    // 0) add free charges from external photons, if any

    //these next steps happen inside the cell itself:
    //   1) estimate if there is a thermaly generated free charge
    //   2) check if a trap releases a free charge
    //   3) for each free charge, check for avalanche by comparing to the Gieger prob (GP)
    //   4) check for avalanche:
    //   if there is an avelanche:
    //     5) compute size of avalanche based on OV
    //     6) release appropriate number of photons
    //     7) add appropriate traps
    //     8) reduce the overvoltage (OV) to zero
    //   else
    //     9) compute the recharge current
    //     10) increase the OV by the step charge Q

void MakeCrvSiPMResponses::FillPhotonQueue(const std::vector<double> &photons)
{
//schedule charges caused by the CRV counter photons
//no check whether time>=_blindTime && time<_mircoBunchPeriod, since this should be done in the calling method
  std::vector<double>::const_iterator iter;
  for(iter=photons.begin(); iter!=photons.end(); iter++)
  {
    int cellid = _randFlat.fireInt(_numberPixels);
    _scheduledCharges.emplace(cellid, *iter);  //constructs ScheduledCharge(cellid, *iter)
  }

//schedule random thermal charges
  double timeWindow = _microBunchPeriod - _blindTime;
  for(int cellid=0; cellid<_numberPixels; cellid++)
  {
    int numberThermalCharges = _randPoissonQ.fire(_probabilities._constThermalProb * timeWindow);
    for(int i=0; i<numberThermalCharges; i++)
    {
      double time = _blindTime + timeWindow * _randFlat.fire();
      _scheduledCharges.emplace(cellid, time); //constructs ScheduledCharge(cellid, time)
    }
  }
}

void MakeCrvSiPMResponses::Simulate(const std::vector<double> &photons, 
                                   std::vector<SiPMresponse> &SiPMresponseVector)
{
  _pixels.clear();
  _scheduledCharges.clear();
  FillPhotonQueue(photons);

  while(1)
  {
    std::multiset<ScheduledCharge>::iterator currentCharge = _scheduledCharges.begin();
    if(currentCharge==_scheduledCharges.end()) break;  //no more scheduled charges

    int cellid = currentCharge->_cellid;
    _time = currentCharge->_time;
    _scheduledCharges.erase(currentCharge);

    double wrappedTime=fmod(_time,_microBunchPeriod);
    if(wrappedTime<_blindTime) continue;  //Current time is inside the blind time.
                                          //This may happen for charges which were added later (after pulses, cross talk of after pulses).
                                          //TODO: Do the voltages of all pixels need to be reset after a blind Time?

    //find pixel with cellid, create if it doesn't exist
    std::map<int,Pixel>::iterator p = _pixels.find(cellid);
    if(p==_pixels.end()) p=_pixels.emplace(cellid, Pixel(_bias, _time)).first; 
    // .first returns the iterator to the new pixel

    Pixel &pixel = p->second;
    RechargeCell(pixel);
    double output = GenerateAvalanche(pixel, cellid); //in units of PEs*biasVoltage

    if(output>0) SiPMresponseVector.emplace_back(wrappedTime, output/_bias); //output/bias is the charge is in units of PEs
                                                                             //output time is between blindTime and microBunchPeriod
  } //while(1)
}


//sample program

#ifdef SiPMResponseStandalone
int main()
{
  std::vector<double> photonTimes;
  photonTimes.push_back(50);
  photonTimes.push_back(20);
  photonTimes.push_back(30);
  photonTimes.push_back(30);
  photonTimes.push_back(30);
  photonTimes.push_back(23);
  photonTimes.push_back(56);
  photonTimes.push_back(12);
  std::vector<SiPMresponse> SiPMresponseVector;

  MakeCrvSiPMResponses::ProbabilitiesStruct probabilities;
  probabilities._constGeigerProbCoef = 1.0;
  probabilities._constGeigerProbVoltScale = 5.5;
  probabilities._constTrapType0Prob = 0.14;   //trap_prob*trap_type0_prob=0.2*0.7
  probabilities._constTrapType1Prob = 0.06;   //trap_prob*trap_type1_pron=0.2*0.3
  probabilities._constTrapType0Lifetime = 5;
  probabilities._constTrapType1Lifetime = 50;
  probabilities._constThermalProb = 6.31e-7; //ns^1     1MHz at SiPM --> 1MHz/#pixel=631Hz at Pixel --> 631s^-1 = 6.31-7ns^-1 
  probabilities._constPhotonProduction = 0.4;

  CLHEP::HepJamesRandom engine(1);
  CLHEP::RandFlat randFlat(engine);
  CLHEP::RandPoissonQ randPoissonQ(engine);
  MakeCrvSiPMResponses sim(randFlat,randPoissonQ);
  sim.SetSiPMConstants(1584, 2.4, 0.0, 1695, 0.08, probabilities);
  sim.Simulate(photonTimes, SiPMresponseVector);

  for(unsigned int i=0; i<SiPMresponseVector.size(); i++)
  std::cout<<i<<"   "<<SiPMresponseVector[i]._time<<"   "<<SiPMresponseVector[i]._charge<<std::endl;

  return 0;
}
#endif

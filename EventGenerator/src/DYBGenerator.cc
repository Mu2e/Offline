#include "EventGenerator/inc/DYBGenerator.hh"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <utility>

namespace mu2e
{
DYBGenerator::DYBGenerator(Direction direction, double minTheta, double maxTheta, double minEnergy, double maxEnergy, 
                           double minPhi, double maxPhi, int nBinsTheta, int nBinsEnergy) : 
                           _direction(direction), 
                           _minTheta(minTheta), _maxTheta(maxTheta), 
                           _minEnergy(minEnergy), _maxEnergy(maxEnergy), 
                           _minPhi(minPhi), _maxPhi(maxPhi),
                           _nBinsTheta(nBinsTheta), _nBinsEnergy(nBinsEnergy) 
{
  if(minTheta<0.0 || maxTheta>M_PI/2.0) throw std::logic_error("Theta needs to be within 0 and Pi/2");
  if(minEnergy<0.1 || maxEnergy>10000.0) throw std::logic_error("Energy needs to be within 0.1 GeV and 10,000 GeV");
  if(minTheta>=maxTheta) throw std::logic_error("minTheta must be less than maxTheta");
  if(minEnergy>=maxEnergy) throw std::logic_error("minEnergy must be less than maxEnergy");
  if(minPhi>=maxPhi) throw std::logic_error("minPhi must be less than maxPhi");

  switch(direction)
  {
    case NEGATIVE_Y: if(minPhi<0.0 || maxPhi>2.0*M_PI) throw std::logic_error("Phi needs to be within 0 and 2*Pi"); break;
    case POSITIVE_X: if(minPhi<-M_PI/2.0 || maxPhi>M_PI/2.0) throw std::logic_error("Phi needs to be within -1/2*Pi and 1/2*Pi"); break;
    case POSITIVE_Z: if(minPhi<0 || maxPhi>M_PI) throw std::logic_error("Phi needs to be within 0 and Pi"); break;
    case NEGATIVE_X: if(minPhi<M_PI/2.0 || maxPhi>3.0/2.0*M_PI) throw std::logic_error("Phi needs to be within 1/2*Pi and 3/2*Pi"); break;
    case NEGATIVE_Z: if(minPhi<M_PI || maxPhi>2.0*M_PI) throw std::logic_error("Phi needs to be within Pi and 2*Pi"); break;
            default: throw std::logic_error("Invalid direction");
  };
  _sinMinPhi=sin(minPhi);
  _sinMaxPhi=sin(maxPhi);
  _cosMinPhi=cos(minPhi);
  _cosMaxPhi=cos(maxPhi);

  _dTheta=(maxTheta-minTheta)/_nBinsTheta;
 
  double dLogEnergy=log(maxEnergy/minEnergy)/_nBinsEnergy;  //log(maxEnergy/minEnergy)=log(maxEnergy)-log(minEnergy)
  double logMinEnergy=log(minEnergy);
  _energyBinBorder.resize(_nBinsEnergy+1);
  _energyBinCenter.resize(_nBinsEnergy);
  _dEnergy.resize(_nBinsEnergy);
  for(int i=0; i<=_nBinsEnergy; i++)
  {
    _energyBinBorder[i]=exp(logMinEnergy+i*dLogEnergy);
    if(i>0)
    {
      _energyBinCenter[i-1]=(_energyBinBorder[i]+_energyBinBorder[i-1])/2.0;
      _dEnergy[i-1]=_energyBinBorder[i]-_energyBinBorder[i-1];
    }
  }

  PrepareMuonMap();
}
  
void DYBGenerator::PrepareMuonMap()
{
  double sumTheta=0;
  for(int iTheta=0; iTheta<_nBinsTheta; iTheta++)
  {
    double theta=_minTheta+_dTheta*(iTheta+0.5);
    double costh=cos(theta);

    double sumEnergy=0;
    for(int iEnergy=0; iEnergy<_nBinsEnergy; iEnergy++)
    {
      double energy=_energyBinCenter[iEnergy];
      sumEnergy+=f(costh, energy)*_dEnergy[iEnergy]*_dTheta;
      _FEnergy[std::pair<int,int>(iTheta,iEnergy)]=sumEnergy;
    }
    for(int iEnergy=0; iEnergy<_nBinsEnergy; iEnergy++) _FEnergy[std::pair<int,int>(iTheta,iEnergy)]/=sumEnergy;

    sumTheta+=sumEnergy;
    _FTheta[iTheta]=sumTheta;
  }
  for(int iTheta=0; iTheta<_nBinsTheta; iTheta++) _FTheta[iTheta]/=sumTheta;

  _rate = sumTheta;

  std::cout<<std::endl;
}

void DYBGenerator::GenerateMuon(double &theta, double &energy, double &phi, CLHEP::RandFlat &randFlat)
{
  double uTheta = randFlat.fire(); 
  int iTheta=0;
  for(; iTheta<_nBinsTheta; iTheta++)
  {
    if(_FTheta[iTheta]>=uTheta) break;
  }
  theta=_minTheta+_dTheta*(iTheta+randFlat.fire());

  double uEnergy = randFlat.fire(); 
  int iEnergy=0;
  for(; iEnergy<_nBinsEnergy; iEnergy++)
  {
    if(_FEnergy[std::pair<int,int>(iTheta,iEnergy)]>=uEnergy) break;
  }
  energy=_energyBinBorder[iEnergy]+randFlat.fire()*_dEnergy[iEnergy];

  double uPhi = randFlat.fire(); 
  switch(_direction)
  {
    case NEGATIVE_Y: phi=2.0*M_PI*uPhi; break;
    case POSITIVE_X: phi=asin(uPhi*(_sinMaxPhi-_sinMinPhi)+_sinMinPhi); break;
    case POSITIVE_Z: phi=acos(uPhi*(_cosMaxPhi-_cosMinPhi)+_cosMinPhi); break;
    case NEGATIVE_X: phi=    M_PI-asin(uPhi*(_sinMaxPhi-_sinMinPhi)+_sinMinPhi); break;
    case NEGATIVE_Z: phi=2.0*M_PI-acos(uPhi*(_cosMaxPhi-_cosMinPhi)+_cosMinPhi); break;
            default: throw std::logic_error("Invalid direction");
  }
}

double DYBGenerator::f(double costh, double energy)
{
  double p1 = 0.102573;
  double p2 = -0.068287;
  double p3 = 0.958633;
  double p4 = 0.0407253;
  double p5 = 0.817285;

  double c=sqrt((costh*costh+p1*p1+p2*pow(costh,p3)+p4*pow(costh,p5))/(1.0+p1*p1+p2+p4));
  double intensity=0.14*pow((energy+3.63698/pow(c,1.29685)),-2.7)*(1.0/(1.0+1.1*energy*c/115.0)+0.054/(1.0+1.1*energy*c/850.0));

  double probabilityDensity=0;
  double sinth=sqrt(1-costh*costh);
  switch(_direction)
  {
    case NEGATIVE_Y: probabilityDensity=intensity*costh*sinth*(_maxPhi-_minPhi); break;
    case POSITIVE_X: probabilityDensity=intensity*sinth*sinth*(_sinMaxPhi-_sinMinPhi); break;
    case POSITIVE_Z: probabilityDensity=intensity*sinth*sinth*(_cosMinPhi-_cosMaxPhi); break;
    case NEGATIVE_X: probabilityDensity=intensity*sinth*sinth*(_sinMinPhi-_sinMaxPhi); break;
    case NEGATIVE_Z: probabilityDensity=intensity*sinth*sinth*(_cosMaxPhi-_cosMinPhi); break;
            default: throw std::logic_error("Invalid direction");
  }

  return probabilityDensity;
}
}

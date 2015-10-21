#include "MakeCrvPhotonArrivals.hh"

#include <sstream>

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"

void MakeCrvPhotonArrivals::LoadLookupTable(const std::string &filename)
{
  std::ifstream lookupfile(filename,std::ios::binary);
  if(!lookupfile.good()) throw std::logic_error("Could not open lookup table file."+filename);

  _LC.Read(lookupfile);
  _LBD.Read(lookupfile);

  unsigned int nScintillatorBins = _LBD.getNScintillatorBins();
  unsigned int nFiberBins = _LBD.getNFiberBins();

  //0...scintillationInScintillator, 1...cerenkovInScintillator 2...cerenkovInFiber0, 3...cerenkovInFiber1
  _bins[0].resize(nScintillatorBins);
  _bins[1].resize(nScintillatorBins);
  _bins[2].resize(nFiberBins);
  _bins[3].resize(nFiberBins);

  std::cout<<"Reading CRV lookup tables "<<filename<<" ... "<<std::flush;
  for(unsigned int i=0; i<nScintillatorBins; i++) _bins[0][i].Read(lookupfile);
  for(unsigned int i=0; i<nScintillatorBins; i++) _bins[1][i].Read(lookupfile);
  for(unsigned int i=0; i<nFiberBins; i++) _bins[2][i].Read(lookupfile);
  for(unsigned int i=0; i<nFiberBins; i++) _bins[3][i].Read(lookupfile);
  std::cout<<"Done."<<std::endl;

  lookupfile.close();
}

MakeCrvPhotonArrivals::~MakeCrvPhotonArrivals()
{
}

void MakeCrvPhotonArrivals::MakePhotons(const CLHEP::Hep3Vector &stepStart,   //they need to be points
                          const CLHEP::Hep3Vector &stepEnd,     //local to the CRV bar
                          double timeStart, double timeEnd,
                          int PDGcode, double beta, double charge,
                          double energyDepositedTotal,
                          double energyDepositedNonIonizing)
{
  for(int SiPM=0; SiPM<4; SiPM++) _arrivalTimes[SiPM].clear();

  const CLHEP::Hep3Vector distanceVector = stepEnd-stepStart;
  double stepLength = distanceVector.mag();
  double theta = distanceVector.theta();  //0...+pi
  double phi = distanceVector.phi();      //-pi...+pi
  double r=0;  //distance from fiber center for fiber tables

  double energy = VisibleEnergyDeposition(PDGcode, stepLength, energyDepositedTotal, energyDepositedNonIonizing);

  double precision=0.1; //mm
  int    steps=static_cast<int>(stepLength/precision);
  double energyPortion=energy/steps;
  for(double distance=precision/2.0; distance<stepLength; distance+=precision)
  {
    CLHEP::Hep3Vector p = stepStart + distanceVector*distance/stepLength;  //start position of photons
    double t = timeStart + (timeEnd-timeStart)*distance/stepLength;

    bool isInScintillator = IsInsideScintillator(p);
    int fiber = IsInsideFiber(p,r);
    int binNumber;
    if(fiber==-1) binNumber=_LBD.findScintillatorBin(p.x(),p.y(),p.z());
    else binNumber=_LBD.findFiberBin(beta,theta,phi,r,p.z());
    if(binNumber==-1) break;  //if the particular bin wasn't found

    int nPhotonsScintillation=0;
    int nPhotonsCerenkov=0;
    if(isInScintillator)
    {
      nPhotonsScintillation = static_cast<int>(_scintillationYield*energyPortion+0.5);
      nPhotonsCerenkov = static_cast<int>(GetAverageNumberOfCerenkovPhotons(beta, charge, _LC.rindexScintillator, _LC.cerenkovEnergyIntervalScintillator)*precision+0.5);
    } 
    if(fiber>0) //this implies not in scintillator
    {
      nPhotonsCerenkov = static_cast<int>(GetAverageNumberOfCerenkovPhotons(beta, charge, _LC.rindexFiber, _LC.cerenkovEnergyIntervalFiber)*precision+0.5);
    } 
    int nPhotons = nPhotonsScintillation + nPhotonsCerenkov;

    //loop over all photons created at this point
    for(int i=0; i<nPhotons; i++)
    {
      int table = -1;
      if(isInScintillator)
      {
        if(i<nPhotonsScintillation) table=0;  //scintillation in scintillator
        else table=1;                         //cerenkov in scintillator
      }
      if(fiber!=-1) table = fiber+2;  //cerenkov in fiber
      if(table==-1) break;  //this can't actually happen

      //get the right bin
      const LookupBin &theBin = _bins[table][binNumber]; //TODO: can be optimized that one doesn't have to access the array for every photon

      //loop over all SiPMs
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        //photon arrival probability at SiPM
        double probability = theBin.arrivalProbability[SiPM];
        if(_randFlat.fire()<=probability)  //a photon arrives at the SiPM --> calculate arrival time
        {
          //start time of photons
          double arrivalTime = t;

          //add straigt line travel time
          double zSiPM=(SiPM%2==0?-_LC.halfLength:_LC.halfLength);
          double straightLineTravelTime=fabs(p.z()-zSiPM)/_LC.speedOfLightFiber;
          arrivalTime+=straightLineTravelTime;

          //add fiber decay times depending on the number of emissions
          int nEmissions = GetRandomFiberEmissions(theBin,SiPM);
          for(int iEmission=0; iEmission<nEmissions; iEmission++) arrivalTime+=-_fiberDecayTime*log(_randFlat.fire());

          //add scintillation decay time
          if(table==0)  //for scintillation in scintillator
          {
//            if(_randFlat.fire()<=_LC.ratioFastSlow)   //use user variable, instead of lookup value
            if(_randFlat.fire()<=_scintillatorRatioFastSlow)
              arrivalTime+=-_scintillatorDecayTimeFast*log(_randFlat.fire());
            else
              arrivalTime+=-_scintillatorDecayTimeSlow*log(_randFlat.fire());
          }

          //add additional time delay due to the photons bouncing around
          arrivalTime+=GetRandomTime(theBin,SiPM);

          _arrivalTimes[SiPM].push_back(arrivalTime);

        }// if a photon was created
      }//loop over all SiPMs
    }//loop over all photons at this point
  }//loop over all points along the track
}

bool MakeCrvPhotonArrivals::IsInsideScintillator(const CLHEP::Hep3Vector &p)
{
  if(fabs(p.x())>=_LC.halfThickness) return false;
  if(fabs(p.y())>=_LC.halfWidth) return false;
  if(fabs(p.z())>=_LC.halfLength) return false;

  CLHEP::Hep2Vector p2D(fabs(p.x()), fabs(p.y()));
  CLHEP::Hep2Vector fiberHole2D(0.0, _LC.fiberSeparation/2.0);
  double distance=(p2D-fiberHole2D).mag();
  if(distance<_LC.holeRadius) return false;

  return true;
}

int MakeCrvPhotonArrivals::IsInsideFiber(const CLHEP::Hep3Vector &p, double &r)
{
  CLHEP::Hep2Vector p2D(p.x(), p.y());
  CLHEP::Hep2Vector fiber0(0.0, -_LC.fiberSeparation/2.0);
  CLHEP::Hep2Vector fiber1(0.0, _LC.fiberSeparation/2.0);
  double distanceFiber0=(p2D-fiber0).mag();
  double distanceFiber1=(p2D-fiber1).mag();
  if(distanceFiber0<=_LC.fiberRadius) {r=distanceFiber0; return 0;}
  if(distanceFiber1<=_LC.fiberRadius) {r=distanceFiber1; return 1;}
  return -1;
}

double MakeCrvPhotonArrivals::GetRandomTime(const LookupBin &theBin, int SiPM)
{
  double rand=_randFlat.fire()*LookupBin::probabilityScale;
  double sumProb=0;
  int timeDelay=0;
  for(; timeDelay<LookupBin::nTimeDelays; timeDelay++)
  {
    sumProb+=theBin.timeDelays[SiPM][timeDelay];
    if(rand<=sumProb) break;
  }

  return timeDelay;
}

int MakeCrvPhotonArrivals::GetRandomFiberEmissions(const LookupBin &theBin, int SiPM)
{
  double rand=_randFlat.fire()*LookupBin::probabilityScale;
  double sumProb=0;
  int emissions=0;
  for(; emissions<LookupBin::nFiberEmissions; emissions++)
  {
    sumProb+=theBin.fiberEmissions[SiPM][emissions];
    if(rand<=sumProb) break;
  }

  return emissions;
}

//we have several CRV bars with different lengths
//in order to avoid different lookup tables, an adjusted position for side 1 is used for shorter CRV bars
//(moving it closer to the SiPM for side 1) 
//NOT USED ANYMORE. WE USE DIFFERENT LOOKUP TABLES FOR DIFFERENT LENGTHS
void MakeCrvPhotonArrivals::AdjustPosition(CLHEP::Hep3Vector &p, int SiPM) 
{
  if(isnan(_actualHalfLength)) return; //no adjustment
  if(_actualHalfLength>_LC.halfLength) 
    throw std::logic_error("Actual bar half length is larger than the half length of the lookup table.");

  double difference = _LC.halfLength - _actualHalfLength;
  if(SiPM%2==0) p.setZ(p.z()-difference);
  else p.setZ(p.z()+difference);
}

int MakeCrvPhotonArrivals::GetNumberOfPhotons(int SiPM)
{
  return _arrivalTimes[SiPM].size();
}

const std::vector<double> &MakeCrvPhotonArrivals::GetArrivalTimes(int SiPM)
{
  return _arrivalTimes[SiPM];
}

//average number of cerenkov photons per millimeter
double MakeCrvPhotonArrivals::GetAverageNumberOfCerenkovPhotons(double beta, double charge, double rindex, double cerenkovEnergyInterval) 
{ 
  const double Rfact = 369.81/(CLHEP::eV * CLHEP::cm); //from G4Cerenkov::GetAverageNumberOfPhotons() 

  if(beta<=1.0/rindex) return(0);  //particle too slow -> no Cerenkov radiation

  double n = 1.0 - 1.0/(rindex*rindex*beta*beta);
  n *= Rfact * charge/eplus * charge/eplus * cerenkovEnergyInterval;

  return n;		
}

//this mimics G4EmSaturation::VisibleEnergyDeposition
//but approximates the proton range as very large
//and uses a fit for the electron range which was obtained specifically for Polystyrene
//the error seems to be less than 1%
double MakeCrvPhotonArrivals::VisibleEnergyDeposition(int PDGcode, double stepLength,
                                            double energyDepositedTotal,
                                            double energyDepositedNonIonizing)
{
  if(energyDepositedTotal <= 0.0) { return 0.0; }

  double evis = energyDepositedTotal;

  if(PDGcode==22)
  {
    if(evis>0)
    {
      double eDepOverElectronRange = 27.0*exp(-0.247*pow(fabs(log(evis)+8.2),1.6))+0.177;
      evis /= (1.0 + _LC.scintillatorBirksConstant*eDepOverElectronRange);
    }
  }
  else 
  {
    // protections
    double nloss = energyDepositedNonIonizing;
    if(nloss < 0.0) nloss = 0.0;
    double eloss = energyDepositedTotal - nloss;

    // neutrons
    if(PDGcode==2112 || eloss < 0.0 || stepLength <= 0.0) 
    {
      nloss = energyDepositedTotal;
      eloss = 0.0;
    }

    // continues energy loss
    if(eloss > 0.0) { eloss /= (1.0 + _LC.scintillatorBirksConstant*eloss/stepLength); }
 
    evis = eloss + nloss;
  }

//  std::cout<<"Original/Visible Energy Deposition (manual): "<<energyDepositedTotal<<"/"<<evis<<"   PDGcode: "<<PDGcode<<std::endl;
  return evis;
}


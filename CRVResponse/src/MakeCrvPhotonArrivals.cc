#include "MakeCrvPhotonArrivals.hh"

#include <sstream>

#include <TFile.h>
#include <TH3D.h>
#include <TNtuple.h>

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"

void MakeCrvPhotonArrivals::LoadLookupTable(std::string filename)
{
  _fileLookupTable = new TFile(filename.c_str());
  if(_fileLookupTable==NULL) throw std::logic_error((filename+" not found.").c_str());

  _histSurvivalProb = new TH3D**[4];
  _histTimeDifference = new TH3D**[4];
  _histFiberEmissions = new TH3D**[4];
  for(int table=0; table<4; table++)  //0...scintillation, 1...cerenkov
  {                                   //2...cerenkovFiber0, 3...cerenkovFiber1
    _histSurvivalProb[table] = new TH3D*[4];
    _histTimeDifference[table] = new TH3D*[4];
    _histFiberEmissions[table] = new TH3D*[4];
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      std::string tableTag[4]={"scintillation","cerenkov","cerenkovFiber0","cerenkovFiber1"};
      std::stringstream s1, s2, s3;
      s1<<tableTag[table]<<"_SurvivalProb_SiPM_"<<SiPM;
      s2<<tableTag[table]<<"_TimeDifference_SiPM_"<<SiPM;
      s3<<tableTag[table]<<"_FiberEmissions_SiPM_"<<SiPM;
      _histSurvivalProb[table][SiPM] = dynamic_cast<TH3D*>(_fileLookupTable->Get(s1.str().c_str())); 
      _histTimeDifference[table][SiPM] = dynamic_cast<TH3D*>(_fileLookupTable->Get(s2.str().c_str())); 
      _histFiberEmissions[table][SiPM] = dynamic_cast<TH3D*>(_fileLookupTable->Get(s3.str().c_str())); 
       if(_histSurvivalProb[table][SiPM]==NULL) throw std::logic_error("didn't find the survival prob map");
       if(_histTimeDifference[table][SiPM]==NULL) throw std::logic_error("didn't find the time difference map");
       if(_histFiberEmissions[table][SiPM]==NULL) throw std::logic_error("didn't find the fiber emissions map");
    }
  }

  TNtuple *ntuple = dynamic_cast<TNtuple*>(_fileLookupTable->Get("CRVbarConstants"));
  if(ntuple==NULL) throw std::logic_error("didn't find CRV bar constants.");
  ntuple->SetBranchAddress("halfThickness",&_halfThickness);
  ntuple->SetBranchAddress("halfWidth",&_halfWidth);
  ntuple->SetBranchAddress("halfLength",&_halfLength);
  ntuple->SetBranchAddress("speedOfLightFiber",&_speedOfLightFiber);
  ntuple->SetBranchAddress("cerenkovRindex",&_cerenkovRindex);
  ntuple->SetBranchAddress("cerenkovEinterval",&_cerenkovEinterval);
  ntuple->SetBranchAddress("ratioFastSlow",&_ratioFastSlow);
  ntuple->SetBranchAddress("scintillatorDensity",&_scintillatorDensity);
  ntuple->SetBranchAddress("scintillatorBirksConstant",&_scintillatorBirksConstant);
  ntuple->SetBranchAddress("fiberSeparation",&_fiberSeparation);
  ntuple->SetBranchAddress("holeRadius",&_holeRadius);
  ntuple->SetBranchAddress("fiberRadius",&_fiberRadius);
  ntuple->GetEntry(0);
}

MakeCrvPhotonArrivals::~MakeCrvPhotonArrivals()
{
  if(_fileLookupTable) _fileLookupTable->Close();
  _fileLookupTable=NULL;

//TODO: delete ROOT dynamically allocated objects

}

void MakeCrvPhotonArrivals::MakePhotons(const CLHEP::Hep3Vector &stepStart,   //they need to be points
                          const CLHEP::Hep3Vector &stepEnd,     //local to the CRV bar
                          double timeStart, double timeEnd,
                          int PDGcode, double beta, double charge,
                          double energyDepositedTotal,
                          double energyDepositedNonIonizing)
{
  if(_fileLookupTable==NULL) throw std::logic_error("Lookup table needs to be loaded.");

  const CLHEP::Hep3Vector distanceVector = stepEnd-stepStart;
  double stepLength = distanceVector.mag();

  double energy = VisibleEnergyDeposition(PDGcode, stepLength, energyDepositedTotal, energyDepositedNonIonizing);

  double nPhotonsScintillator = _scintillationYield*energy;
  double nPhotonsCerenkov = GetAverageNumberOfCerenkovPhotons(beta, charge)*stepLength;
  double nPhotons = nPhotonsScintillator + nPhotonsCerenkov;
  double cerenkovFraction = nPhotonsCerenkov/nPhotons;

  for(int SiPM=0; SiPM<4; SiPM++)
  {
    _arrivalTimes[SiPM].clear();

    for(double i=0; i<nPhotons; i++)
    {
      CLHEP::Hep3Vector p = stepStart + distanceVector*i/nPhotons;
      if(SiPM%2==1) AdjustPosition(p); //side 1

      bool isScintillation = (_randFlat.fire()>=cerenkovFraction);
      bool isInScintillator = IsInsideScintillator(p);
      int table = -1;

      if(isScintillation)
      {
        if(isInScintillator) table = 0;
      }
      else //cerenkov
      {
        if(isInScintillator) table = 1;  //cerenkov in scintillator
        else
        {
          int fiber = IsInsideFiber(p);
          if(fiber!=-1)
          {
            table = fiber+2;  //cerenkov in fiber
            p.setX(0.0);
            p.setY(fiber==0?-_fiberSeparation/2.0:_fiberSeparation/2.0);
          }
        }
      }

      if(table==-1) continue;

      int bin = _histSurvivalProb[table][SiPM]->FindBin(p.x(),p.y(),p.z());
      double probability = _histSurvivalProb[table][SiPM]->GetBinContent(bin);

      if(_randFlat.fire()<=probability)  //a photon arrives at the SiPM
      {
        double arrivalTime = timeStart + (timeEnd-timeStart)*i/nPhotons;

        double zSiPM=(SiPM%2==0?-_halfLength:_halfLength);
        double straightLineTravelTime=abs(p.z()-zSiPM)/_speedOfLightFiber;
        arrivalTime+=straightLineTravelTime;

        int nEmissions = GetRandomFiberEmissions(_histFiberEmissions[table][SiPM], p.y(), p.z());
        for(int iEmission=0; iEmission<nEmissions; iEmission++) arrivalTime+=-_fiberDecayTime*log(_randFlat.fire());

        if(isScintillation)
        {
          if(_randFlat.fire()<=_ratioFastSlow)
            arrivalTime+=-_scintillatorDecayTimeFast*log(_randFlat.fire());
          else
            arrivalTime+=-_scintillatorDecayTimeSlow*log(_randFlat.fire());
        }

        arrivalTime+=GetRandomTime(_histTimeDifference[table][SiPM], p.y(), p.z());

        _arrivalTimes[SiPM].push_back(arrivalTime);
      }
    }
  }
}

bool MakeCrvPhotonArrivals::IsInsideScintillator(const CLHEP::Hep3Vector &p)
{
  if(abs(p.x())>=_halfThickness) return false;
  if(abs(p.y())>=_halfWidth) return false;
  if(abs(p.z())>=_halfLength) return false;

  CLHEP::Hep2Vector p2D(abs(p.x()), abs(p.y()));
  CLHEP::Hep2Vector fiberHole2D(0.0, _fiberSeparation/2.0);
  double distance=(p2D-fiberHole2D).mag();
  if(distance<_holeRadius) return false;

  return true;
}

int MakeCrvPhotonArrivals::IsInsideFiber(const CLHEP::Hep3Vector &p)
{
  CLHEP::Hep2Vector p2D(p.x(), p.y());
  CLHEP::Hep2Vector fiber0(0.0, -_fiberSeparation/2.0);
  CLHEP::Hep2Vector fiber1(0.0, _fiberSeparation/2.0);
  double distanceFiber0=(p2D-fiber0).mag();
  double distanceFiber1=(p2D-fiber1).mag();
  if(distanceFiber0<=_fiberRadius) return 0;
  if(distanceFiber1<=_fiberRadius) return 1;
  return -1;
}

double MakeCrvPhotonArrivals::GetRandomTime(TH3D *timeDifference, double y, double z)
{
  int binx = timeDifference->GetXaxis()->FindBin(y);
  int biny = timeDifference->GetYaxis()->FindBin(z);
  double integral=0;
  int binz=1;
  for(; binz<=timeDifference->GetNbinsZ(); binz++)
  {
    integral+=timeDifference->GetBinContent(binx,biny,binz);
  }

  double r=_randFlat.fire();
  double sum=0;
  binz=1;
  for(; binz<=timeDifference->GetNbinsZ(); binz++)
  {
    sum+=timeDifference->GetBinContent(binx,biny,binz);
    if(r<=sum/integral) break;
  }

  double time=timeDifference->GetZaxis()->GetBinCenter(binz);
  return time;
}

int MakeCrvPhotonArrivals::GetRandomFiberEmissions(TH3D *fiberEmissions, double y, double z)
{
  int binx = fiberEmissions->GetXaxis()->FindBin(y);
  int biny = fiberEmissions->GetYaxis()->FindBin(z);
  double integral=0;
  int binz=1;
  for(; binz<=fiberEmissions->GetNbinsZ(); binz++)
  {
    integral+=fiberEmissions->GetBinContent(binx,biny,binz);
  }

  double r=_randFlat.fire();
  double sum=0;
  binz=1;
  for(; binz<=fiberEmissions->GetNbinsZ(); binz++)
  {
    sum+=fiberEmissions->GetBinContent(binx,biny,binz);
    if(r<=sum/integral) break;
  }

  int nEmissions=static_cast<int>(fiberEmissions->GetZaxis()->GetBinCenter(binz));  //bin center is e.g. 2.5 --> n = 2
  return nEmissions;
}

//we have several CRV bars with different lengths
//in order to avoid different lookup tables, an adjusted position for side 1 is used for shorter CRV bars
//(moving it closer to the SiPM for side 1) 
void MakeCrvPhotonArrivals::AdjustPosition(CLHEP::Hep3Vector &p) 
{
  if(isnan(_actualHalfLength)) return; //no adjustment
  if(_actualHalfLength>_halfLength) 
    throw std::logic_error("Actual bar half length is larger than the half length of the lookup table.");

  double difference = _halfLength - _actualHalfLength;
  p.setZ(p.z()+difference);
}

int MakeCrvPhotonArrivals::GetNumberOfPhotons(int SiPM)
{
  return _arrivalTimes[SiPM].size();
}

const std::vector<double> &MakeCrvPhotonArrivals::GetArrivalTimes(int SiPM)
{
  return _arrivalTimes[SiPM];
}

double MakeCrvPhotonArrivals::GetAverageNumberOfCerenkovPhotons(double beta, double charge) 
{ 
  const double Rfact = 369.81/(CLHEP::eV * CLHEP::cm); //from G4Cerenkov::GetAverageNumberOfPhotons() 

  if(beta<=1.0/_cerenkovRindex) return(0);  //particle too slow -> no Cerenkov radiation

  double n = 1.0 - 1.0/(_cerenkovRindex*_cerenkovRindex*beta*beta);
  n *= Rfact * charge/eplus * charge/eplus * _cerenkovEinterval;

  return n;		
}

//this mimics G4EmSaturation::VisibleEnergyDeposition
//but approximates the proton range as large
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
      evis /= (1.0 + _scintillatorBirksConstant*eDepOverElectronRange);
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
    if(eloss > 0.0) { eloss /= (1.0 + _scintillatorBirksConstant*eloss/stepLength); }
 
    evis = eloss + nloss;
  }

//  std::cout<<"Visible Energy Deposition (manual): "<<evis<<"   PDGcode: "<<PDGcode<<std::endl;
  return evis;
}


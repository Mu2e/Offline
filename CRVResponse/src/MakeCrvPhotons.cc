#include "Offline/CRVResponse/inc/MakeCrvPhotons.hh"

#include <sstream>

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"

namespace mu2eCrv
{
void LookupConstants::Write(const std::string &filename)
{
  std::ofstream lookupfile(filename,std::ios::binary|std::ios::app);
  lookupfile.write(reinterpret_cast<char*>(this),sizeof(LookupConstants));
  lookupfile.close();
}
void LookupConstants::Read(std::ifstream &lookupfile)
{
  lookupfile.read(reinterpret_cast<char*>(this),sizeof(LookupConstants));
}

void LookupCerenkov::WriteMap(std::map<double,double> &m, std::ofstream &o)
{
  size_t n=m.size();
  o.write(reinterpret_cast<char*>(&n),sizeof(size_t));
  for(std::map<double,double>::const_iterator iter=m.begin(); iter!=m.end(); iter++)
  {
    o.write(reinterpret_cast<const char*>(&iter->first),sizeof(double));
    o.write(reinterpret_cast<const char*>(&iter->second),sizeof(double));
  }
}
void LookupCerenkov::ReadMap(std::map<double,double> &m, std::ifstream &i)
{
  size_t n;
  i.read(reinterpret_cast<char*>(&n),sizeof(size_t));
  double d1,d2;
  for(size_t j=0; j<n; j++)
  {
    i.read(reinterpret_cast<char*>(&d1),sizeof(double));
    i.read(reinterpret_cast<char*>(&d2),sizeof(double));
    m[d1]=d2;
  }
}
void LookupCerenkov::Write(const std::string &filename)
{
  std::ofstream lookupfile(filename,std::ios::binary|std::ios::app);
  WriteMap(photonsScintillator,lookupfile);
  WriteMap(photonsFiber,lookupfile);
  lookupfile.close();
}
void LookupCerenkov::Read(std::ifstream &lookupfile)
{
  ReadMap(photonsScintillator,lookupfile);
  ReadMap(photonsFiber,lookupfile);
}

void LookupBinDefinitions::WriteVector(std::vector<double> &v, std::ofstream &o)
{
  size_t n=v.size();
  o.write(reinterpret_cast<char*>(&n),sizeof(size_t));
  o.write(reinterpret_cast<char*>(v.data()),sizeof(double)*v.size());
}
void LookupBinDefinitions::ReadVector(std::vector<double> &v, std::ifstream &i)
{
  size_t n;
  i.read(reinterpret_cast<char*>(&n),sizeof(size_t));
  double *d=new double[n];
  i.read(reinterpret_cast<char*>(d),sizeof(double)*n);
  v.assign(d,d+n);
  delete [] d;
}
void LookupBinDefinitions::Write(const std::string &filename)
{
  std::ofstream lookupfile(filename,std::ios::binary|std::ios::app);
  WriteVector(xBins,lookupfile);
  WriteVector(yBins,lookupfile);
  WriteVector(zBins,lookupfile);
  WriteVector(betaBins,lookupfile);
  WriteVector(thetaBins,lookupfile);
  WriteVector(phiBins,lookupfile);
  WriteVector(rBins,lookupfile);
  lookupfile.close();
}
void LookupBinDefinitions::Read(std::ifstream &lookupfile)
{
  ReadVector(xBins,lookupfile);
  ReadVector(yBins,lookupfile);
  ReadVector(zBins,lookupfile);
  ReadVector(betaBins,lookupfile);
  ReadVector(thetaBins,lookupfile);
  ReadVector(phiBins,lookupfile);
  ReadVector(rBins,lookupfile);
}

unsigned int LookupBinDefinitions::getNScintillatorScintillationBins()
{
  unsigned int nXBins = xBins.size()-1;
  unsigned int nYBins = yBins.size()-1;
  unsigned int nZBins = zBins.size()-1;
  return nXBins*nYBins*nZBins;
}
unsigned int LookupBinDefinitions::getNScintillatorCerenkovBins()
{
  unsigned int nXBins = xBins.size()-1;
  unsigned int nYBins = yBins.size()-1;
  unsigned int nZBins = zBins.size()-1;
  unsigned int nBetaBins = betaBins.size()-1;
  return nXBins*nYBins*nZBins*nBetaBins;
}
unsigned int LookupBinDefinitions::getNFiberCerenkovBins()
{
  unsigned int nBetaBins = betaBins.size()-1;
  unsigned int nThetaBins = thetaBins.size()-1;
  unsigned int nPhiBins = phiBins.size()-1;
  unsigned int nRBins = rBins.size()-1;
  unsigned int nZBins = zBins.size()-1;
  return nBetaBins*nThetaBins*nPhiBins*nRBins*nZBins;
}

unsigned int LookupBinDefinitions::findBin(const std::vector<double> &v, const double &x, bool &notFound)
{
  if(v.size()<2 || !(v.front()<=x && x<=v.back()))
  {
    notFound=true;
    return(-1);
  }
  return(std::lower_bound(v.begin()+1,v.end(),x)-v.begin()-1);
}
int LookupBinDefinitions::findScintillatorScintillationBin(double x, double y, double z)
{
  bool notFound=false;
  unsigned int xBin=findBin(xBins,x,notFound);
  unsigned int yBin=findBin(yBins,y,notFound);
  unsigned int zBin=findBin(zBins,z,notFound);
  if(notFound) return(-1);

  unsigned int nYBins = yBins.size()-1;
  unsigned int nZBins = zBins.size()-1;
  return(zBin + yBin*nZBins + xBin*nYBins*nZBins);
}
int LookupBinDefinitions::findScintillatorCerenkovBin(double x, double y, double z, double beta)
{
  bool notFound=false;
  unsigned int xBin=findBin(xBins,x,notFound);
  unsigned int yBin=findBin(yBins,y,notFound);
  unsigned int zBin=findBin(zBins,z,notFound);
  unsigned int betaBin=findBin(betaBins,beta,notFound);
  if(notFound) return(-1);

  unsigned int nYBins = yBins.size()-1;
  unsigned int nZBins = zBins.size()-1;
  unsigned int nBetaBins = betaBins.size()-1;
  return(betaBin + zBin*nBetaBins + yBin*nZBins*nBetaBins + xBin*nYBins*nZBins*nBetaBins);
}
int LookupBinDefinitions::findFiberCerenkovBin(double beta, double theta, double phi, double r, double z)
{
  bool notFound=false;
  unsigned int betaBin=findBin(betaBins,beta,notFound);
  unsigned int thetaBin=findBin(thetaBins,theta,notFound);
  unsigned int phiBin=findBin(phiBins,phi,notFound);
  unsigned int rBin=findBin(rBins,r,notFound);
  unsigned int zBin=findBin(zBins,z,notFound);
  if(notFound) return(-1);

  unsigned int nThetaBins = thetaBins.size()-1;
  unsigned int nPhiBins = phiBins.size()-1;
  unsigned int nRBins = rBins.size()-1;
  unsigned int nZBins = zBins.size()-1;
  return(zBin + rBin*nZBins + phiBin*nRBins*nZBins + thetaBin*nPhiBins*nRBins*nZBins + betaBin*nThetaBins*nPhiBins*nRBins*nZBins);
}
bool LookupBinDefinitions::findScintillatorScintillationBinReverse(unsigned int bin, double &xbin, double &ybin, double &zbin)
{
  if(bin>=getNScintillatorScintillationBins()) return false;

  int nXBins=xBins.size()-1;   //e.g. 3 bins need 4 entries in the vector (for 3 bin boundaries)
  int nYBins=yBins.size()-1;
  int nZBins=zBins.size()-1;
  xbin = (bin / (nZBins*nYBins)) % nXBins;
  ybin = (bin / nZBins) % nYBins;
  zbin = bin % nZBins;
  return true;
}
bool LookupBinDefinitions::findScintillatorCerenkovBinReverse(unsigned int bin, double &xbin, double &ybin, double &zbin, double &betabin)
{
  if(bin>=getNScintillatorCerenkovBins()) return false;

  int nXBins=xBins.size()-1;   //e.g. 3 bins need 4 entries in the vector (for 3 bin boundaries)
  int nYBins=yBins.size()-1;
  int nZBins=zBins.size()-1;
  int nBetaBins=betaBins.size()-1;
  xbin = (bin / (nBetaBins*nZBins*nYBins)) % nXBins;
  ybin = (bin / (nBetaBins*nZBins)) % nYBins;
  zbin = (bin / nBetaBins) % nZBins;
  betabin = bin % nBetaBins;
  return true;
}
bool LookupBinDefinitions::findFiberCerenkovBinReverse(unsigned int bin, double &betabin, double &thetabin, double &phibin, double &rbin, double &zbin)
{
  if(bin>=getNFiberCerenkovBins()) return false;

  int nZBins=zBins.size()-1;
  int nBetaBins=betaBins.size()-1;
  int nThetaBins=thetaBins.size()-1;
  int nPhiBins=phiBins.size()-1;
  int nRBins=rBins.size()-1;
  betabin = (bin / (nZBins*nRBins*nPhiBins*nThetaBins)) % nBetaBins;
  thetabin = (bin / (nZBins*nRBins*nPhiBins)) % nThetaBins;
  phibin = (bin / (nZBins*nRBins)) % nPhiBins;
  rbin = (bin / nZBins) % nRBins;
  zbin = bin % nZBins;
  return true;
}

void LookupBin::WriteVector(std::vector<unsigned char> &v, std::ofstream &o)
{
  size_t n=v.size();
  o.write(reinterpret_cast<char*>(&n),sizeof(size_t));
  o.write(reinterpret_cast<char*>(v.data()),sizeof(unsigned char)*v.size());
}
void LookupBin::ReadVector(std::vector<unsigned char> &v, std::ifstream &i)
{
  size_t n;
  i.read(reinterpret_cast<char*>(&n),sizeof(size_t));
  unsigned char *d=new unsigned char[n];
  i.read(reinterpret_cast<char*>(d),sizeof(unsigned char)*n);
  v.assign(d,d+n);
  delete [] d;
}
void LookupBin::Write(const std::string &filename)
{
  //Lookup tables are created only for SiPM# 0 due to symmetry reasons
  std::ofstream lookupfile(filename,std::ios::binary|std::ios::app);
  lookupfile.write(reinterpret_cast<char*>(&binNumber),sizeof(unsigned int));
  lookupfile.write(reinterpret_cast<char*>(&arrivalProbability),sizeof(float));
  WriteVector(timeDelays, lookupfile);
  WriteVector(fiberEmissions, lookupfile);
  lookupfile.close();
}
void LookupBin::Read(std::ifstream &lookupfile, const unsigned int &i)
{
  //Lookup tables are created only for SiPM# 0 due to symmetry reasons
  lookupfile.read(reinterpret_cast<char*>(&binNumber),sizeof(unsigned int));
  lookupfile.read(reinterpret_cast<char*>(&arrivalProbability),sizeof(float));
  ReadVector(timeDelays, lookupfile);
  ReadVector(fiberEmissions, lookupfile);
  probabilityScaleTimeDelays=0;
  probabilityScaleFiberEmissions=0;
  for(size_t j=0; j<timeDelays.size(); ++j) probabilityScaleTimeDelays+=timeDelays[j];
  for(size_t j=0; j<fiberEmissions.size(); ++j) probabilityScaleFiberEmissions+=fiberEmissions[j];
  if(i!=binNumber) throw std::logic_error("Corrupt lookup table.");
}

  void MakeCrvPhotons::LoadLookupTable(const std::string &filename, int debug)
{
  _fileName = filename;
  std::ifstream lookupfile(filename,std::ios::binary);
  if(!lookupfile.good()) throw std::logic_error("Could not open lookup table file "+filename);

  _LC.Read(lookupfile);
  if(_LC.version1!=6) throw std::logic_error("This version of Offline expects a lookup table version 6.x.");
  if(_LC.reflector!=0 && _LC.reflector!=1 && _LC.reflector!=2) throw std::logic_error("Lookup tables can have either no reflector/absorber, or a reflector/absorber on the +z side.");

  _LCerenkov.Read(lookupfile);
  _LBD.Read(lookupfile);

  unsigned int nScintillatorScintillationBins = _LBD.getNScintillatorScintillationBins();
  unsigned int nScintillatorCerenkovBins      = _LBD.getNScintillatorCerenkovBins();
  unsigned int nFiberCerenkovBins             = _LBD.getNFiberCerenkovBins();

  //0...scintillationInScintillator, 1...cerenkovInScintillator 2...cerenkovInFiber
  _bins[0].resize(nScintillatorScintillationBins);
  _bins[1].resize(nScintillatorCerenkovBins);
  _bins[2].resize(nFiberCerenkovBins);

  if(debug>0) std::cout<<"Reading CRV lookup tables "<<filename<<" ... "<<std::flush;
  for(unsigned int i=0; i<nScintillatorScintillationBins; i++) _bins[0][i].Read(lookupfile,i);
  for(unsigned int i=0; i<nScintillatorCerenkovBins; i++)      _bins[1][i].Read(lookupfile,i);
  for(unsigned int i=0; i<nFiberCerenkovBins; i++)             _bins[2][i].Read(lookupfile,i);
  // a truncated table would otherwise leave the bin boundaries and yield curves holding
  // whatever was in the buffers, and the job would continue with silently wrong photon yields
  if(!lookupfile) throw std::logic_error("Error while reading CRV lookup table file "+filename+
                                         " (truncated or corrupt).");
  if(debug>0) std::cout<<"Done."<<std::endl;

  lookupfile.close();
}

MakeCrvPhotons::~MakeCrvPhotons()
{
}

void MakeCrvPhotons::MakePhotons(const CLHEP::Hep3Vector &stepStartTmp,   //they need to be points
                          const CLHEP::Hep3Vector &stepEndTmp,            //local to the CRV bar
                          double timeStart, double timeEnd,
                          double beta, double charge,
                          double visibleEnergyDeposited,
                          double trueTotalStepLength,   //may be longer than stepEnd-stepStart due to scattering
                                                        //is needed for the Cerenkov photons
                          int    reflector)
{
  if(_LC.reflector!=std::abs(reflector)) throw std::logic_error("Expected reflector/absorber doesn't match lookup table.");

  for(int SiPM=0; SiPM<_nSiPMs; SiPM++) _arrivalTimes[SiPM].clear();

  //coordinates are in local coordinates of the scintillator (x:thickness, y:width, z:length)
  CLHEP::Hep3Vector stepStart[_nSiPMs];
  CLHEP::Hep3Vector stepEnd[_nSiPMs];
  //local start/end positions for all for SiPMs derived from the symmetries of the counter,
  //because lookup tables are stored only for SiPM #0.
  for(int SiPM=0; SiPM<_nSiPMs; SiPM++)
  {
    stepStart[SiPM] = stepStartTmp;
    stepEnd[SiPM]   = stepEndTmp;
    if(reflector==-1 || reflector==-2)
    {
      stepStart[SiPM].setZ(-stepStart[SiPM].z());
      stepEnd[SiPM].setZ(-stepEnd[SiPM].z());
    }
  }
  stepStart[1].setZ(-stepStart[1].z());
  stepStart[2].setY(-stepStart[2].y());
  stepStart[3].setZ(-stepStart[3].z());
  stepStart[3].setY(-stepStart[3].y());
  stepEnd[1].setZ(-stepEnd[1].z());
  stepEnd[2].setY(-stepEnd[2].y());
  stepEnd[3].setZ(-stepEnd[3].z());
  stepEnd[3].setY(-stepEnd[3].y());

  for(int SiPM=0; SiPM<_nSiPMs; SiPM++)
  {
    //there are only lookup tables without reflector or with reflector on the +z side (i.e. at SiPMs #1 and #3)
    if(reflector!=0 && (SiPM==1 || SiPM==3)) continue;

    const CLHEP::Hep3Vector distanceVector = stepEnd[SiPM]-stepStart[SiPM];
    double totalStepLength = distanceVector.mag();
    double theta = distanceVector.theta();  //0...+pi

    double precision=0.1; //mm
    int    nSteps=std::max(static_cast<int>(totalStepLength/precision),1);

    double avgNPhotonsScintillation = _scintillationYield*visibleEnergyDeposited;
    double avgNPhotonsCerenkovInScintillator
           = GetAverageNumberOfCerenkovPhotons(beta, charge, _LCerenkov.photonsScintillator)*trueTotalStepLength;  //use the true path, since it  may be longer due to curved paths
    double avgNPhotonsCerenkovInFiber
           = GetAverageNumberOfCerenkovPhotons(beta, charge, _LCerenkov.photonsFiber)*trueTotalStepLength;  //use the true path, since it  may be longer due to curved paths

    int nPhotonsScintillationPerStep          = GetNumberOfPhotonsFromAverage(avgNPhotonsScintillation,nSteps);
    int nPhotonsCerenkovInScintillatorPerStep = GetNumberOfPhotonsFromAverage(avgNPhotonsCerenkovInScintillator,nSteps);
    int nPhotonsCerenkovInFiberPerStep        = GetNumberOfPhotonsFromAverage(avgNPhotonsCerenkovInFiber,nSteps);

    //first find all bin numbers, so that repeated bin numbers can be merged
    std::vector<int> binNumbers0, binNumbers1, binNumbers2;
    binNumbers0.reserve(nSteps);
    binNumbers1.reserve(nSteps);
    binNumbers2.reserve(nSteps);
    for(int step=0; step<nSteps; step++)
    {
      double stepFraction = (step+0.5)/nSteps;
      CLHEP::Hep3Vector p = stepStart[SiPM] + distanceVector*stepFraction;  //local start position of photons

      bool isInScintillator = IsInsideScintillator(p);
      double r=0;  //distance from fiber center for fiber tables
      double phi=0;  //angle w.r.t. the radius vector from the fiber center to the point in the 2D cross section plane
                     //0...+pi due to symmetry
      bool isInFiber = IsInsideFiber(p,distanceVector, r,phi);

      int binNumber0=-1;
      int binNumber1=-1;
      int binNumber2=-1;
      if(isInScintillator)
      {
        binNumber0=_LBD.findScintillatorScintillationBin(fabs(p.x()),p.y(),p.z());  //use only positive x values due to symmetry in x
        binNumber1=_LBD.findScintillatorCerenkovBin(fabs(p.x()),p.y(),p.z(),beta);  //use only positive x values due to symmetry in x
      }
      else if(isInFiber)
      {
        binNumber2=_LBD.findFiberCerenkovBin(beta,theta,phi,r,p.z());
      }
      binNumbers0.push_back(binNumber0);
      binNumbers1.push_back(binNumber1);
      binNumbers2.push_back(binNumber2);
    }

    for(int step=0; step<nSteps; step++)
    {
      int binNumber0=binNumbers0.at(step);
      int binNumber1=binNumbers1.at(step);
      int binNumber2=binNumbers2.at(step);

      //find out how many subsequent steps use the same bin numbers
      int repeatedBins=1; //merge steps if they use the same bin.
      for(int i=step+1; i<nSteps; ++i)
      {
        if(binNumbers0.at(i)==binNumber0 && binNumbers1.at(i)==binNumber1 && binNumbers2.at(i)==binNumber2) ++repeatedBins;
        else break;
      }

      double stepFraction = (step+0.5*(repeatedBins-1)+0.5)/nSteps;  //center of the sequence of multiple steps
      double photonStartTime = timeStart + (timeEnd-timeStart)*stepFraction;
      step+=repeatedBins-1;  //skip forward

      //lookup table numbers
      //- for scintillation in scintillator is 0
      //- for cerenkov in scintillator is 1
      //- for cerenkov in fiber is 2
      if(binNumber0>=0) ProcessPhotons(&_bins[0][binNumber0], nPhotonsScintillationPerStep*repeatedBins, photonStartTime, SiPM, reflector);
      if(binNumber1>=0) ProcessPhotons(&_bins[1][binNumber1], nPhotonsCerenkovInScintillatorPerStep*repeatedBins, photonStartTime, SiPM, reflector);
      if(binNumber2>=0) ProcessPhotons(&_bins[2][binNumber2], nPhotonsCerenkovInFiberPerStep*repeatedBins, photonStartTime, SiPM, reflector);
    }//loop over all steps along the track
  }//loop over all SiPMs

}

void MakeCrvPhotons::ProcessPhotons(const LookupBin *theBin, int nPhotons, double photonStartTime, int SiPM, int reflector)
{
  if(theBin==NULL || nPhotons<=0) return;

  //photon arrival probability at SiPM
  double probability = theBin->arrivalProbability;
  probability*=_photonYieldDeviation[SiPM];   //channel-specific deviation from nominal, e.g. due to scintillator variations or SiPM misalignments

  //number of photons arriving at the SiPM - calculated from a binomial distribution
  int nPhotonsArrived=0;
  //special CLHEP::RandBinomial cases: probability equal to 0 or 1 results in -1. (also when nPhotons equal to 0, but this case is checked above.)
  if(probability>=1.0) nPhotonsArrived=nPhotons;
  else if(probability>0.0)
  {
    //CLHEP::RandBinomial returns 0 for n*p<~0.01 and returns values that are too low up to n*p~0.03.
    //So it should be used only for large n*p, where it returns the correct results and is faster than the inversion.
    if(nPhotons*probability>15.0) nPhotonsArrived=static_cast<int>(_randBinomial.fire(nPhotons,probability));
    else //inverse transform sampling. takes more time for large n*p
    {
      double logq = std::log1p(-probability);
      double pk   = std::exp(nPhotons*logq);   //P(0)=q^n  //the n*p>15.0 cut prevents this from ever being zero (which could happen for very negative arguments)
      double cdf  = pk;                                    //if pk was 0, the while loop wouldn't be entered
      double r    = probability/(1.0-probability);         //(to prevent it from getting stuck, because pk would remain at 0 in the loop).
      double u    = _randFlat.fire();
      while(u>cdf && nPhotonsArrived<nPhotons && pk>0.0)
      {
        pk  *= r*(nPhotons-nPhotonsArrived)/(nPhotonsArrived+1.0);   //P(k+1)=P(k)*(n-k)/(k+1)*p/q
        cdf += pk;
        ++nPhotonsArrived;
      }
    }
  }

  //calculate arrival time for each photon
  for(long iPhoton=0; iPhoton<nPhotonsArrived; ++iPhoton)
  {
    //start time of photons
    double arrivalTime = photonStartTime;

    //add fiber decay times depending on the number of emissions
    int nEmissions = GetRandomFiberEmissions(theBin);
    for(int iEmission=0; iEmission<nEmissions; iEmission++) arrivalTime+=-_LC.WLSfiberDecayTime*log(_randFlat.fire());

    //add additional time delay due to the photons bouncing around
    arrivalTime+=GetRandomTime(theBin);

    if(reflector!=-1 && reflector!=-2) _arrivalTimes[SiPM].push_back(arrivalTime);
    else _arrivalTimes[SiPM+1].push_back(arrivalTime);

  }//photons arriving at the SiPM
}

bool MakeCrvPhotons::IsInsideScintillator(const CLHEP::Hep3Vector &p)
{
  if(fabs(p.x())>=_LC.halfThickness) return false;
  if(fabs(p.y())>=_LC.halfWidth) return false;
  if(fabs(p.z())>=_LC.halfLength) return false;

  CLHEP::Hep2Vector pos2D(fabs(p.x()), fabs(p.y()));
  CLHEP::Hep2Vector fiberPos2D(0.0, _LC.fiberSeparation/2.0);
  CLHEP::Hep2Vector fiberRadius2D(_LC.holeRadiusX, _LC.holeRadiusY);

  double ellipseSum=0;
  for(int i=0; i<2; i++)
  {
    double tmp=(pos2D[i]-fiberPos2D[i])/fiberRadius2D[i];
    ellipseSum+=tmp*tmp;
  }
  if(ellipseSum<1) return false; // <1 means inside ellipse, i.e. inside fiber hole. ellipse formula: (x-P_x)^2 / r_x^2 + (y-P_y)^2 / r_y^2 = 1

  CLHEP::Hep2Vector cornerPos2D(_LC.halfThickness-_LC.scintillatorCornerRadius,_LC.halfWidth-_LC.scintillatorCornerRadius);
  if(pos2D[0]>cornerPos2D[0] && pos2D[1]>cornerPos2D[1])
  {
    if((pos2D-cornerPos2D).mag()>_LC.scintillatorCornerRadius) return false;
  }

  return true;
}

bool MakeCrvPhotons::IsInsideFiber(const CLHEP::Hep3Vector &p, const CLHEP::Hep3Vector &dir, double &r, double &phi)
{
  if(fabs(p.z())>=_LC.halfLength) return false;

  CLHEP::Hep2Vector p2D(p.x(), p.y());
  CLHEP::Hep2Vector fiber0(0.0, -_LC.fiberSeparation/2.0);
  CLHEP::Hep2Vector fiber1(0.0, _LC.fiberSeparation/2.0);
  CLHEP::Hep2Vector radiusVector0=p2D-fiber0;  //radius vector from fiber0 center to p2D
  CLHEP::Hep2Vector radiusVector1=p2D-fiber1;  //radius vector from fiber1 center to p2D
  double distanceFiber0=radiusVector0.mag();
  double distanceFiber1=radiusVector1.mag();

  int fiberNumber=-1;
  if(distanceFiber0<=_LC.fiberRadius) {r=distanceFiber0; fiberNumber=0;}
  if(distanceFiber1<=_LC.fiberRadius) {r=distanceFiber1; fiberNumber=1;}
  if(fiberNumber==-1) return false;

  CLHEP::Hep2Vector dir2D(dir.x(),dir.y());
  if(fiberNumber==0) phi=dir2D.angle(radiusVector0);  //angle returns between 0 and pi
  if(fiberNumber==1) phi=dir2D.angle(radiusVector1);
  return true;
}

double MakeCrvPhotons::GetRandomTime(const LookupBin *theBin)
{
  //The lookup tables encodes probabilities as probability*mu2eCrv::LookupBin::probabilityScale(255),
  //so that the probabilities can be stored as integers. For example, the probability of 1 is stored as 255.
  //Due to rounding issues, the sum of all entries for this bin may not be 255.
  //This bin-specifc sum is the probabilityScaleTimeDelays.

  size_t timeDelay=0;
  double rand=_randFlat.fire()*theBin->probabilityScaleTimeDelays;
  double sumProb=0;
  size_t maxTimeDelay=theBin->timeDelays.size();
  for(timeDelay=0; timeDelay<maxTimeDelay; ++timeDelay)
  {
    sumProb+=theBin->timeDelays[timeDelay];
    if(rand<=sumProb) break;
  }

  return static_cast<double>(timeDelay);
}

int MakeCrvPhotons::GetRandomFiberEmissions(const LookupBin *theBin)
{
  //The lookup tables encodes probabilities as probability*mu2eCrv::LookupBin::probabilityScale(255),
  //so that the probabilities can be stored as integers. For example, the probability of 1 is stored as 255.
  //Due to rounding issues, the sum of all entries for this bin may not be 255.
  //This bin-specifc sum is the probabilityScaleFiberEmissions.

  size_t emissions=0;
  double rand=_randFlat.fire()*theBin->probabilityScaleFiberEmissions;
  double sumProb=0;
  size_t maxEmissions=theBin->fiberEmissions.size();
  for(emissions=0; emissions<maxEmissions; ++emissions)
  {
    sumProb+=theBin->fiberEmissions[emissions];
    if(rand<=sumProb) break;
  }

  return emissions;
}

int MakeCrvPhotons::GetNumberOfPhotonsFromAverage(double average, int nSteps)  //from G4Scintillation
{
  int nPhotons;
  if(average>10.0)
  {
    double sigma = std::sqrt(average);
    nPhotons = lrint(_randGaussQ.fire(average,sigma)/nSteps);
    if(nPhotons<0) nPhotons=0;
  }
  else
  {
    nPhotons = lrint(_randPoissonQ.fire(average)/nSteps);
  }
  return nPhotons;
}

int MakeCrvPhotons::GetNumberOfPhotons(int SiPM)
{
  return _arrivalTimes[SiPM].size();
}

const std::vector<double> &MakeCrvPhotons::GetArrivalTimes(int SiPM)
{
  return _arrivalTimes[SiPM];
}

//average number of cerenkov photons per millimeter
double MakeCrvPhotons::GetAverageNumberOfCerenkovPhotons(double beta, double charge, std::map<double,double> &photons)
{
  if(charge==0) return 0;

  double prevBeta=0;
  double prevNumberPhotons=0;
  std::map<double,double>::const_iterator i;
  for(i=photons.begin(); i!=photons.end(); i++)
  {
    if(beta<=i->first) //first beta in the Cerenkov map is always greater than 0.0 and has number of photons of 0.
    {
      if(i->first==0) return i->second;  //beta=0 in Cerenkov map shouldn't happen. if it did and was the first entry, i->first-prevBeta would result in a division by 0.
      double numberPhotons=prevNumberPhotons+(i->second-prevNumberPhotons)/(i->first-prevBeta)*(beta-prevBeta);
      numberPhotons*=fabs(charge/eplus);
      return numberPhotons;
    }
    prevBeta=i->first;
    prevNumberPhotons=i->second;
  }
  return photons.rbegin()->second*fabs(charge/eplus); //this shouldn't happen
}

} //namespace mu2e

#ifndef MakeCrvPhotonArrivals_h
#define MakeCrvPhotonArrivals_h

#include <vector>
#include <map>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/Randomize.h"

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>

struct LookupConstants
{
  double halfThickness, halfWidth, halfLength;
  double speedOfLightFiber;
  double rindexScintillator;
  double rindexFiber;
  double cerenkovEnergyIntervalScintillator;
  double cerenkovEnergyIntervalFiber;
  double ratioFastSlow;
  double scintillatorDensity;
  double scintillatorBirksConstant;
  double fiberSeparation, holeRadius, fiberRadius;
  void Write(const std::string &filename)
  {
    std::ofstream lookupfile(filename,std::ios::binary|std::ios::app);
    lookupfile.write(reinterpret_cast<char*>(this),sizeof(LookupConstants));
    lookupfile.close();
  }
  void Read(std::ifstream &lookupfile)
  {
    lookupfile.read(reinterpret_cast<char*>(this),sizeof(LookupConstants));
  }
};

struct LookupBinDefinitions
{
  std::vector<double> xBins;
  std::vector<double> yBins;
  std::vector<double> zBins;
  std::vector<double> betaBins;
  std::vector<double> thetaBins;
  std::vector<double> phiBins;
  std::vector<double> rBins;
  void WriteVector(std::vector<double> &v, std::ofstream &o)
  {
    size_t n=v.size();
    o.write(reinterpret_cast<char*>(&n),sizeof(size_t));
    o.write(reinterpret_cast<char*>(v.data()),sizeof(double)*v.size());
  }
  void ReadVector(std::vector<double> &v, std::ifstream &i)
  {
    size_t n;
    i.read(reinterpret_cast<char*>(&n),sizeof(size_t));
    double d[n];
    i.read(reinterpret_cast<char*>(d),sizeof(double)*n);
    v.assign(d,d+n);
  }
  void Write(const std::string &filename)
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
  void Read(std::ifstream &lookupfile)
  {
    ReadVector(xBins,lookupfile);
    ReadVector(yBins,lookupfile);
    ReadVector(zBins,lookupfile);
    ReadVector(betaBins,lookupfile);
    ReadVector(thetaBins,lookupfile);
    ReadVector(phiBins,lookupfile);
    ReadVector(rBins,lookupfile);
  }

  unsigned int getNScintillatorBins()
  {
    unsigned int nXBins = xBins.size()-1;
    unsigned int nYBins = yBins.size()-1;
    unsigned int nZBins = zBins.size()-1;
    return nXBins*nYBins*nZBins;
  }
  unsigned int getNFiberBins()
  {
    unsigned int nBetaBins = betaBins.size()-1;
    unsigned int nThetaBins = thetaBins.size()-1;
    unsigned int nPhiBins = phiBins.size()-1;
    unsigned int nRBins = rBins.size()-1;
    unsigned int nZBins = zBins.size()-1;
    return nBetaBins*nThetaBins*nPhiBins*nRBins*nZBins;
  }

  unsigned int findBin(const std::vector<double> &v, const double &x, bool &notFound)
  {
    for(unsigned int i=1; i<v.size(); i++)
    {
      if(v[i]>=x)
      {
        if(v[i-1]<=x) return(i-1);
      }
    }
    notFound=true;
    return(-1);
  }
  int findScintillatorBin(double x, double y, double z)
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
  int findFiberBin(double beta, double theta, double phi, double r, double z)
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
};

struct LookupBin
{
  static const int nTimeDelays=150;
  static const int nFiberEmissions=15;
  static const unsigned short probabilityScale=10000;  //still within unsigned short (2 bytes)

  float arrivalProbability[4];
  unsigned short timeDelays[4][nTimeDelays];
  unsigned short fiberEmissions[4][nFiberEmissions];
  void Write(const std::string &filename)
  {
    std::ofstream lookupfile(filename,std::ios::binary|std::ios::app);
    lookupfile.write(reinterpret_cast<char*>(arrivalProbability),sizeof(float)*4);
    lookupfile.write(reinterpret_cast<char*>(timeDelays),sizeof(unsigned short)*4*nTimeDelays);
    lookupfile.write(reinterpret_cast<char*>(fiberEmissions),sizeof(unsigned short)*4*nFiberEmissions);
    lookupfile.close();
  }
  void Read(std::ifstream &lookupfile)
  {
    lookupfile.read(reinterpret_cast<char*>(arrivalProbability),sizeof(float)*4);
    lookupfile.read(reinterpret_cast<char*>(timeDelays),sizeof(unsigned short)*4*nTimeDelays);
    lookupfile.read(reinterpret_cast<char*>(fiberEmissions),sizeof(unsigned short)*4*nFiberEmissions);
  }
};



class MakeCrvPhotonArrivals
{
  public:

    MakeCrvPhotonArrivals(CLHEP::RandFlat &randFlat) : _actualHalfLength(NAN), _randFlat(randFlat) {}

    ~MakeCrvPhotonArrivals();

    void                      LoadLookupTable(const std::string &filename);
    void                      MakePhotons(const CLHEP::Hep3Vector &stepStart,   //they need to be points
                                      const CLHEP::Hep3Vector &stepEnd,         //local to the CRV bar
                                      double timeStart, double timeEnd,
                                      int PDGcode, double beta, double charge,
                                      double energyDepositedTotal,
                                      double energyDepositedNonIonizing);
    int                       GetNumberOfPhotons(int SiPM);
    const std::vector<double> &GetArrivalTimes(int SiPM);
    void                      SetScintillationYield(double scintillationYield) {_scintillationYield=scintillationYield;}
    void                      SetScintillatorRatioFastSlow(double ratio) {_scintillatorRatioFastSlow=ratio;}
    void                      SetScintillatorDecayTimeFast(double decayTime) {_scintillatorDecayTimeFast=decayTime;}
    void                      SetScintillatorDecayTimeSlow(double decayTime) {_scintillatorDecayTimeSlow=decayTime;}
    void                      SetFiberDecayTime(double decayTime) {_fiberDecayTime=decayTime;}

    void                      SetActualHalfLength(double actualHalfLength) {_actualHalfLength=actualHalfLength;}

  private:

    std::vector<double>       _arrivalTimes[4];
    double                    _scintillationYield;
    double                    _scintillatorRatioFastSlow;
    double                    _scintillatorDecayTimeFast; 
    double                    _scintillatorDecayTimeSlow; 
    double                    _fiberDecayTime;

    double                    _actualHalfLength;

    LookupConstants           _LC;
    LookupBinDefinitions      _LBD;
    std::vector<LookupBin>    _bins[4];

    CLHEP::RandFlat           &_randFlat;

    bool   IsInsideScintillator(const CLHEP::Hep3Vector &p);
    int    IsInsideFiber(const CLHEP::Hep3Vector &p, double &r);
    double GetRandomTime(const LookupBin &theBin, int SiPM);
    int    GetRandomFiberEmissions(const LookupBin &theBin, int SiPM);
    double GetAverageNumberOfCerenkovPhotons(double beta, double charge, double rindex, double cerenkovEnergyInterval);
    void   AdjustPosition(CLHEP::Hep3Vector &p, int SiPM);
    double VisibleEnergyDeposition(int PDGcode, double stepLength,
                                   double energyDepositedTotal,
                                   double energyDepositedNonIonizing);

    public:
    void   DrawHistograms();
};

#endif

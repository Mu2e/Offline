#ifndef MakeCrvPhotons_h
#define MakeCrvPhotons_h

#include <vector>
#include <map>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/Randomize.h"

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>

namespace mu2eCrv
{

struct LookupConstants
{
  int    version1, version2;
  int    reflector;
  double halfThickness, halfWidth, halfLength;
  double fiberSeparation, holeRadiusX, holeRadiusY, fiberRadius, scintillatorCornerRadius;
  double scintillatorBirksConstant;
  double WLSfiberDecayTime;
  void Write(const std::string &filename);
  void Read(std::ifstream &lookupfile);
};

//number of photons per mm for a given beta
struct LookupCerenkov
{
  std::map<double,double> photonsScintillator;
  std::map<double,double> photonsFiber;
  void WriteMap(std::map<double,double> &m, std::ofstream &o);
  void ReadMap(std::map<double,double> &m, std::ifstream &i);
  void Write(const std::string &filename);
  void Read(std::ifstream &lookupfile);
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
  void WriteVector(std::vector<double> &v, std::ofstream &o);
  void ReadVector(std::vector<double> &v, std::ifstream &i);
  void Write(const std::string &filename);
  void Read(std::ifstream &lookupfile);

  unsigned int getNScintillatorScintillationBins();
  unsigned int getNScintillatorCerenkovBins();
  unsigned int getNFiberCerenkovBins();

  unsigned int findBin(const std::vector<double> &v, const double &x, bool &notFound);
  int findScintillatorScintillationBin(double x, double y, double z);
  int findScintillatorCerenkovBin(double x, double y, double z, double beta);
  int findFiberCerenkovBin(double beta, double theta, double phi, double r, double z);

  bool findScintillatorScintillationBinReverse(unsigned int bin, double &xbin, double &ybin, double &zbin);
  bool findScintillatorCerenkovBinReverse(unsigned int bin, double &xbin, double &ybin, double &zbin, double &betabin);
  bool findFiberCerenkovBinReverse(unsigned int bin, double &betabin, double &thetabin, double &phibin, double &rbin, double &zbin);
};

struct LookupBin
{
  static const int maxTimeDelays=150;
  static const int maxFiberEmissions=10;
  static const unsigned short probabilityScale=255;  //still within unsigned char (1 byte)
  //the lookup tables encodes probabilities as probability*probabilityScale(255), 
  //so that the probabilities can be stored as unsigned chars.
  //therefore, the probability of 1 is stored as 255.

  //Lookup tables are created only for SiPM# 0 due to symmetry reasons
  unsigned int binNumber;  //to check whether file was assembled correctly
  float arrivalProbability;
  std::vector<unsigned char> timeDelays;
  std::vector<unsigned char> fiberEmissions;
  unsigned int probabilityScaleTimeDelays;
  unsigned int probabilityScaleFiberEmissions;
  void WriteVector(std::vector<unsigned char> &v, std::ofstream &o);
  void ReadVector(std::vector<unsigned char> &v, std::ifstream &i);
  void Write(const std::string &filename);
  void Read(std::ifstream &lookupfile, const unsigned int &i);
};



class MakeCrvPhotons
{
  public:

    MakeCrvPhotons(CLHEP::RandFlat &randFlat, CLHEP::RandGaussQ &randGaussQ, CLHEP::RandPoissonQ &randPoissonQ) : 
                                                      _randFlat(randFlat), _randGaussQ(randGaussQ), _randPoissonQ(randPoissonQ) {}

    ~MakeCrvPhotons();

    const std::string         &GetFileName() const {return _fileName;}

    void                      LoadLookupTable(const std::string &filename);
    void                      MakePhotons(const CLHEP::Hep3Vector &stepStart,   //they need to be points
                                      const CLHEP::Hep3Vector &stepEnd,         //local to the CRV bar
                                      double timeStart, double timeEnd,
                                      double beta, double charge,
                                      double visibleEnergyDeposited,
                                      double trueStepLength,
                                      int reflector=0);
    int                       GetNumberOfPhotons(int SiPM);
    const std::vector<double> &GetArrivalTimes(int SiPM);
    void                      SetScintillationYield(double yield) {_scintillationYield=yield;}

  private:

    std::string               _fileName;
    int                       _reflector;

    std::vector<double>       _arrivalTimes[4];
    double                    _scintillationYield;

    LookupConstants           _LC;
    LookupCerenkov            _LCerenkov;
    LookupBinDefinitions      _LBD;
    std::vector<LookupBin>    _bins[3];   //scintillation in scintillator (0), Cerenkov in scintillator (1), Cerenkov in fiber (2)

    CLHEP::RandFlat           &_randFlat;
    CLHEP::RandGaussQ         &_randGaussQ;
    CLHEP::RandPoissonQ       &_randPoissonQ;

    bool   IsInsideScintillator(const CLHEP::Hep3Vector &p);
    bool   IsInsideFiber(const CLHEP::Hep3Vector &p, const CLHEP::Hep3Vector &dir, double &r, double &phi);
    double GetRandomTime(const LookupBin *theBin);
    int    GetRandomFiberEmissions(const LookupBin *theBin);
    double GetAverageNumberOfCerenkovPhotons(double beta, double charge, std::map<double,double> &photons);
    int    GetNumberOfPhotonsFromAverage(double average, int nSteps);

    public:
    void   DrawHistograms();
};

} //namespace mu2e

#endif

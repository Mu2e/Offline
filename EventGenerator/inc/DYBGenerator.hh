#include <vector>
#include <map>
#include "CLHEP/Random/RandFlat.h"

namespace mu2e
{ 
class DYBGenerator
{
  public: 
  enum Direction
  {
    UNDEFINED,
    NEGATIVE_Y,
    NEGATIVE_X,
    POSITIVE_X,
    NEGATIVE_Z,
    POSITIVE_Z
  };
  
  DYBGenerator(Direction direction, double minTheta, double maxTheta, double minEnergy, double maxEnergy, 
               double minPhi, double maxPhi, int nBinsTheta, int nBinsEnergy);
  ~DYBGenerator() {}

  void PrepareMuonMap();
  void GenerateMuon(double &theta, double &energy, double &phi, CLHEP::RandFlat &randFlat);

  Direction GetDirection() {return _direction;}
  double    GetRate()      {return _rate;}

  private:
  Direction _direction;

  double _minTheta, _maxTheta, _minEnergy, _maxEnergy, _minPhi, _maxPhi;
  int    _nBinsTheta, _nBinsEnergy;

  double _sinMinPhi, _sinMaxPhi, _cosMinPhi, _cosMaxPhi;
  double _dTheta;

  std::vector<double> _energyBinBorder;
  std::vector<double> _energyBinCenter;
  std::vector<double> _dEnergy;

  double _rate;
  std::map<int,double>                 _FTheta;
  std::map<std::pair<int,int>, double> _FEnergy;

  double f(double costh, double energy);
};

}

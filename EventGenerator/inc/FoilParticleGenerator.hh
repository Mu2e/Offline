#ifndef EventGenerator_FoilParticleGenerator_hh
#define EventGenerator_FoilParticleGenerator_hh

//
// Generate position (Mu2e coordinates) and time of a generic particle coming from target foils.
//
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Need to be improved at a later date.
//  - Limits on cos(theta) and phi but uniform within the range.


#include <memory>
#include <fstream>

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/RandomLimitedExpo.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"

//CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class FoilParticleGenerator {

  public:

    enum foilGen_enum {
      flatFoil, volWeightFoil, expoFoil, expoVolWeightFoil, muonFileInputFoil, lastFoil_enum
    };

    enum posGen_enum {
      flatPos, muonFileInputPos, lastPos_enum
    };

    enum timeGen_enum {
      flatTime, limitedExpoTime, negExp, lastTime_enum
    };
    
    
#define FOILALGO_NAMES							\
    "flatFoil", "volWeightFoil", "expoFoil", "expoVolWeightFoil" , "muonFileInputFoil" 
    
#define POSALGO_NAMES				\
    "flatPos", "muonFileInputPos"
    
#define TIMEALGO_NAMES				\
    "flatTime", "limitedExpoTime", "negExp"    
    
  public:
    
    FoilParticleGenerator( art::RandomNumberGenerator::base_engine_t& engine,
			   double tmin, double tmax, foilGen_enum foilAlgo,
			   posGen_enum  posAlgo, timeGen_enum  timeAlgo,
			   bool PTtoSTdelay, 
			   bool pPulseDelay, double pPulseShift,
			     std::string STinfilename, int linesToSkip = 0);
    
    ~FoilParticleGenerator();
    
    void generatePositionAndTime(CLHEP::Hep3Vector& pos, double& time, bool foldingTimeOption);
    
    int iFoil();
    
    double pulseDelay();
    double muDelay();

    static foilGen_enum findFoilGenByName (std::string const& name);
    static posGen_enum  findPosGenByName (std::string const& name);
    static timeGen_enum findTimeGenByName (std::string const& name);

    const static char *_foilName[];
    const static char* _posName[];
    const static char* _timeName[];

  private:

    // time generation range
    double _tmin, _tmax;

    double _maxtime;

    // foil, position and time random algorithm
    foilGen_enum  _foilAlgo;
    posGen_enum   _posAlgo;
    timeGen_enum  _timeAlgo;

    //number of foils
    int _nfoils;

    //extracted foil.
    int _ifoil;

    // Random numbers generators
    CLHEP::RandFlat     _randFlat;
    RandomLimitedExpo   _randTime;
    double _muTimeDecay;
    CLHEP::RandExponential     _randNegExpoTime;
    CLHEP::RandGeneral  _randFoils;
    CLHEP::RandGeneral  _randExpoFoils;
    CLHEP::RandGeneral  _delayTime;
    ProtonPulseRandPDF  _pulseTime;

    bool _DSFrame;

    //Include a delay in time due to the PT to ST path
    bool _PTtoSTdelay;
    bool _pPulseDelay;
    double _pPulseShift;
    std::string _STinfilename;

    //Number of lines to skip in the stopped muon input file
    int _ntoskip;

    double _muDelay, _pulseDelay;

    //Build a binned representation of foils volume
    std::vector<double> binnedFoilsVolume();
    std::vector<double> weightedBinnedFoilsVolume();
    std::vector<double> timePathDelay();

    // methods to extract foil, position and time w.r.t. the chosen algorithm
    int getFlatRndFoil() ;
    int getVolumeRndFoil() ;
    int getExpoRndFoil() ;
    int getVolumeAndExpoRndFoil() ;
    double includeTimeDelay();
    double includePulseDelay();
    CLHEP::Hep3Vector getFlatRndPos(TargetFoil const& theFoil) ;
    double getFlatRndTime() ;
    double getLimitedExpRndTime() ;
    void getInfoFromFile(CLHEP::Hep3Vector& pos, double& time);
    double getNegativeExpoRndTime();
    double getMuTimeDecay();
    void PointToBeginData();
    void SkipStartingLines();

    static std::fstream& STfile(std::string STinfilename);

  };

} // end namespace mu2e,

#endif /* EventGenerator_FoilParticleGenerator_hh */


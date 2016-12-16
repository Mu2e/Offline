#ifndef CRVAnalysisInfo_hh
#define CRVAnalysisInfo_hh

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

namespace mu2e
{

  struct CRVHitsReco   //information about a cluster of CRV coincidence triplets
  {
    int               _crvSectorType;   //CRV sector type
    CLHEP::Hep3Vector _pos;             //average position
    double            _timeWindowStart; //first hit time
    double            _timeWindowEnd;   //last hit time
    int               _PEs;                   //total number of PEs for this cluster
    int               _nCoincidenceHits;      //number of coincidence hits in this cluster
    CRVHitsReco(int crvSectorType, CLHEP::Hep3Vector pos, double timeWindowStart, double timeWindowEnd, int PEs, int nCoincidenceHits) :
                _crvSectorType(crvSectorType),
                _pos(pos),
                _timeWindowStart(timeWindowStart),
                _timeWindowEnd(timeWindowEnd),
                _PEs(PEs),
                _nCoincidenceHits(nCoincidenceHits)
                {}
  };

  struct CRVHitsMC  //information about the MC track which most likely caused the CRV coincidence triplets
  {
    bool              _valid;            //was an MC particle found that matches the coincidence triplets? 
    int               _pdgId;            //PDG ID of this MC particle
    int               _primaryPdgId;     //PDG ID of the primary particle of this MC particle (helps to determine whether it was a cosmic ray, etc.)
    std::string       _generator;        //generator of the primary particle
    CLHEP::Hep3Vector _pos;              //position of the MC particle when it "created" the first StepPointMC
    CLHEP::Hep3Vector _momentum;         //momentum of the MC particle when it "created" the first StepPointMC
    double            _time;             //time of the MC particle when it "created" the first StepPointMC
    double            _depositedEnergy;  //total energy deposited for this cluster (not just for this track)
    CRVHitsMC(bool valid, int pdgId, int primaryPdgId, const std::string &generator, 
              CLHEP::Hep3Vector pos, CLHEP::Hep3Vector momentum, double time, double depositedEnergy) :
              _valid(valid),
              _pdgId(pdgId),
              _primaryPdgId(primaryPdgId),
              _generator(generator),
              _pos(pos),
              _momentum(momentum),
              _time(time),
              _depositedEnergy(depositedEnergy)
              {}
  };

  typedef std::vector<CRVHitsReco> CRVHitsRecoCollection;  //this is the reco vector which will be stored in the main TTree 
  typedef std::vector<CRVHitsMC>   CRVHitsMCCollection;    //this is the MC vector which will be stored in the main TTree 

}
#endif



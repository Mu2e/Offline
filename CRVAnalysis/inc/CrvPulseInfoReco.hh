#ifndef CrvPulseInfoReco_hh
#define CrvPulseInfoReco_hh

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include <map>
#include "Rtypes.h"

namespace mu2e
{

  struct CrvPulseInfoReco   //information about CRV reco pulses
  {

    Float_t             _x, _y, _z;       //average position of counters
    Int_t               _barId;           //CRV counter ID
    Int_t               _sectorId;        //CRV sector ID
    Int_t               _SiPMId;          //SiPMId number
    Int_t               _PEs;             //PEs using pulse integral
    Int_t               _PEsPulseHeight;  //PEs using pulse height
    Float_t             _pulseHeight;     //Pulse height
    Float_t             _pulseBeta;       //Pulse beta
    Float_t             _pulseFitChi2;    //Pulse Fit chi2
    Float_t             _time;            //Time

    CrvPulseInfoReco(CLHEP::Hep3Vector pos, int barId, int sectorId, int SiPMId, int PEs, int PEsPulseHeight, float pulseHeight, float pulseBeta, float pulseFitChi2, float time) :
                _x(pos.x()), _y(pos.y()), _z(pos.z()),
                _barId(barId),
                _sectorId(sectorId),
                _SiPMId(SiPMId),
                _PEs(PEs),
                _PEsPulseHeight(PEsPulseHeight),
                _pulseHeight(pulseHeight),
                _pulseBeta(pulseBeta),
                _pulseFitChi2(pulseFitChi2),
                _time(time)
                {}
    CrvPulseInfoReco() :
                _x(0), _y(0), _z(0),
                _barId(-1),
                _sectorId(-1),
                _SiPMId(-1),
                _PEs(-1),
                _PEsPulseHeight(-1),
                _pulseHeight(-1),
                _pulseBeta(-1),
                _pulseFitChi2(-1),
                _time(-1)
                {}
  };

  typedef std::vector<CrvPulseInfoReco> CrvPulseInfoRecoCollection;  //this is the reco vector which will be stored in the main TTree

}
#endif

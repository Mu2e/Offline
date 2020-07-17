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


    struct myrecopulse
    {
      Int_t _PEs;
      Float_t _time;

      myrecopulse(int PEs, float time ) :
        _PEs(PEs), _time(time)
      {}
      myrecopulse() :
        _PEs(0), _time(0)
      {}
    };

    struct myrecodigis
    {
      Int_t _ADCs;
      Float_t _time;
      myrecodigis(int ADCs, float time ) :
        _ADCs(ADCs), _time(time)
      {}
      myrecodigis() :
        _ADCs(0), _time(0)
      {}
    };




    Float_t             _x, _y, _z;       //average position of counters
    Int_t               _barId;           //CRV counter ID
    Int_t               _sectorId;        //CRV sector ID
    Int_t               _SiPMId;          //SiPMId number
    std::vector<myrecopulse>  _recoPulses;
    std::vector<myrecodigis>  _recoDigis;

    CrvPulseInfoReco(CLHEP::Hep3Vector pos, int barId, int sectorId, int SiPMId) :
                _x(pos.x()), _y(pos.y()), _z(pos.z()),
                _barId(barId),
                _sectorId(sectorId),
                _SiPMId(SiPMId)
                {}
    CrvPulseInfoReco() :
                _x(0), _y(0), _z(0),
                _barId(-1),
                _sectorId(-1),
                _SiPMId(-1)
                {}
  };

  typedef std::map<int, CrvPulseInfoReco> CrvPulseInfoRecoCollection;  //this is the reco vector which will be stored in the main TTree
  //  typedef std::vector<std::pair<int, CrvPulseInfoReco>> CrvPulseInfoRecoCollection;  //this is the reco vector which will be stored in the main TTree

}
#endif

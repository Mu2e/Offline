//
// MC relationships of some objects.  Extracted from KalDiag
// Dave Brown, LBNL, 20 Jun 2016
//
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>
namespace mu2e
{
  using std::vector;

  MCRelationship::MCRelationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2) :
    MCRelationship(mcd1.earlyStrawGasStep()->simParticle(),mcd2.earlyStrawGasStep()->simParticle())
  {}

  MCRelationship::MCRelationship(StrawDigiMC const& mcd, SPPtr const& spp) :
    MCRelationship(mcd.earlyStrawGasStep()->simParticle(),spp)
  {}

  MCRelationship::MCRelationship(SPPtr const& sppi,SPPtr const& sppj) : _rel(none), _rem(-1) {
    if(sppi.isNonnull() && sppj.isNonnull()){
      if(sppi == sppj){
        _rel= same;
        _rem = 0;
      } else {
        SPPtr pi = sppi->originParticle().parent();
        SPPtr pj = sppj->originParticle().parent();
        if(pi.isNonnull() && pi == sppj){
          _rel = daughter;
          _rem = 1;
        } else if(pj.isNonnull() && pj == sppi) {
          _rel = mother;
          _rem = 1;
        } else if(pi.isNonnull() && pj.isNonnull() && pi == pj){
          _rel = sibling;
          _rem = 1;
        } else {
          vector<SPPtr > pvi, pvj;
          pvi.push_back(sppi);
          pvj.push_back(sppj);
          while(pi.isNonnull()){
            pvi.push_back(pi);
            pi = pi->originParticle().parent();
          }
          while(pj.isNonnull()){
            pvj.push_back(pj);
            pj = pj->originParticle().parent();
          }
          auto idau = find(pvi.begin(),pvi.end(),sppj);
          auto jdau = find(pvj.begin(),pvj.end(),sppi);
          if(idau != pvi.end()){
            _rel = udaughter;
            _rem = std::distance(pvi.begin(),idau);
          } else if(jdau != pvj.end()){
            _rel = umother;
            _rem = std::distance(pvj.begin(),jdau);
          } else {
            for(size_t jj=0;jj<pvj.size();++jj){
              auto icuz = find(pvi.begin(),pvi.end(),pvj[jj]);
              if(icuz != pvi.end()){
                _rel = usibling;
                _rem = jj + distance(pvi.begin(),icuz);
                break;
              }
            }
            for(size_t ii=0;ii<pvi.size();++ii){
              auto jcuz = find(pvj.begin(),pvj.end(),pvi[ii]);
              if(jcuz != pvj.end()){
                _rel = usibling;
                _rem = ii + distance(pvj.begin(),jcuz);
                break;
              }
            }
          }
        }
      }
    }
  }
}


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

  namespace {
    // A particle simulated in several stages has one SimParticle per stage; each later copy is a self-parent of
    // the copy before it. Relationships are between physical particles, so every Ptr is mapped to its first copy.
    MCRelationship::SPPtr firstCopy(MCRelationship::SPPtr sp) {
      while(sp->selfParent()) sp = sp->parent();
      return sp;
    }
    // the first copy of the parent of a first copy, or a null Ptr for a primary
    MCRelationship::SPPtr parentCopy(MCRelationship::SPPtr const& sp) {
      auto const& parent = sp->parent();
      return parent.isNonnull() ? firstCopy(parent) : parent;
    }
  }

  MCRelationship::MCRelationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2) :
    MCRelationship(mcd1.earlyStrawGasStep()->simParticle(),mcd2.earlyStrawGasStep()->simParticle())
  {}

  MCRelationship::MCRelationship(StrawDigiMC const& mcd, SPPtr const& spp) :
    MCRelationship(mcd.earlyStrawGasStep()->simParticle(),spp)
  {}

  MCRelationship::MCRelationship(SPPtr const& sppi_in,SPPtr const& sppj_in) : _rel(none), _rem(-1) {
    if(sppi_in.isNonnull() && sppj_in.isNonnull()){
      SPPtr const sppi = firstCopy(sppi_in);
      SPPtr const sppj = firstCopy(sppj_in);
      if(sppi == sppj){
        _rel= same;
        _rem = 0;
      } else {
        SPPtr pi = parentCopy(sppi);
        SPPtr pj = parentCopy(sppj);
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
            pi = parentCopy(pi);
          }
          while(pj.isNonnull()){
            pvj.push_back(pj);
            pj = parentCopy(pj);
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


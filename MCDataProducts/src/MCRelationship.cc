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
    // A particle simulated in several stages has one SimParticle per stage. Each later copy's parent is the copy
    // before it, and selfParent() is true for the later copy. Relationships are between physical particles, so
    // every Ptr is mapped to its origin particle, the earliest copy (the Ptr form of SimParticle::originParticle),
    // before comparing.
    // Example, a two-stage cosmic: P1 is the stage-1 primary, P2 its stage-2 copy (P2.parent() is P1), and E2 a
    // daughter made in stage 2 (E2.parent() is P2). MCRelationship(P2,P1) is 'same' and MCRelationship(E2,P1)
    // is 'daughter'. Comparing raw Ptrs gave 'none' for both: P2 != P1, and the walk up from P2 jumped straight
    // to P1's parent (null), so P1 itself was never compared.
    MCRelationship::SPPtr originPtr(MCRelationship::SPPtr sp) {
      while(sp->selfParent()) sp = sp->parent();
      return sp;
    }
    // the origin of the parent of an origin particle, or a null Ptr for a primary
    MCRelationship::SPPtr parentOriginPtr(MCRelationship::SPPtr const& sp) {
      auto const& parent = sp->parent();
      return parent.isNonnull() ? originPtr(parent) : parent;
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
      SPPtr const sppi = originPtr(sppi_in);
      SPPtr const sppj = originPtr(sppj_in);
      if(sppi == sppj){
        _rel= same;
        _rem = 0;
      } else {
        SPPtr pi = parentOriginPtr(sppi);
        SPPtr pj = parentOriginPtr(sppj);
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
            pi = parentOriginPtr(pi);
          }
          while(pj.isNonnull()){
            pvj.push_back(pj);
            pj = parentOriginPtr(pj);
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


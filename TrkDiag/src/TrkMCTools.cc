//
// Namespace for collecting tools used in MC truth evaluation
// Original author: Dave Brown (LBNL) 8/10/2016
//
#include "Offline/TrkDiag/inc/TrkMCTools.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"

#include <map>

namespace mu2e {
  namespace TrkMCTools {

    bool CEDigi(StrawDigiMC const& mcdigi) {
      bool conversion(false);
      auto const& spmcp = mcdigi.earlyStrawGasStep();
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      int gid(-1);
      if(spp->genParticle().isNonnull())
        gid = spp->genParticle()->generatorId().id();
      // a conversion electron is an electron from the CE generator.  The momentum requirement
      // removes cases where the CE loses a catastrophic amount of energy (ie albedo backsplash
      // from the calorimeter).
      conversion = (spp->pdgId() == 11 && gid == 2 && sqrt(spmcp->momentum().mag2())>90.0);
      return conversion;
    }

    unsigned primaryParticle(art::Ptr<SimParticle>& pspp, std::vector<StrawHitIndex> const& hits, const StrawDigiMCCollection* mcdigis) {
      unsigned retval(0);
      std::map<art::Ptr<SimParticle>, unsigned> spmap;

      for( auto hi : hits ) {
        StrawDigiMC const& mcdigi = mcdigis->at(hi);
        art::Ptr<SimParticle> spp = mcdigi.earlyStrawGasStep()->simParticle();
        auto ispp = spmap.find(spp);
        if(ispp != spmap.end())
          ++(ispp->second);
        else
          spmap[spp] = 1;
      }
      for( auto imap : spmap ) {
        if(imap.second > retval){
          retval = imap.second;
          pspp = imap.first;
        }
      }
      return retval;
    }

    unsigned countDigis(art::Ptr<SimParticle> const& pspp, const StrawDigiMCCollection* mcdigis) {
      unsigned retval(0);
      for( auto mcdigi : *mcdigis ) {
        art::Ptr<SimParticle> spp = mcdigi.earlyStrawGasStep()->simParticle();
        if(spp == pspp) {
          ++retval;
        }
      }
      return retval;
    }

    void findMCTrk(const KalSeed& kseed,art::Ptr<SimParticle>& spp, StrawDigiMCCollection const& mcdigis) {
      static art::Ptr<SimParticle> nullvec;
      spp = nullvec;
      std::vector<spcount> sct;
      findMCTrk(kseed,sct, mcdigis);
      if(sct.size()>0)
        spp = sct[0]._spp;
    }

    void findMCTrk(const KalSeed& kseed,std::vector<spcount>& sct, StrawDigiMCCollection const& mcdigis, bool saveall) {
      sct.clear();
      // find the SimParticles which contributed hits.
      // loop through the straw hits from the track
      static StrawHitFlag active(StrawHitFlag::active);
      for(const auto& tshs : kseed.hits()) {
        // loop over the hits and find the associated steppoints
        bool isactive = tshs.flag().hasAllProperties(active);
        StrawDigiMC const& mcdigi = mcdigis.at(tshs.index());
        // if mc info is not meaningful, then skip this digi.
        // this looks sketchy, but nowhere is an implicit association
        // between the sct and mcdigis collection actually assumed
        if (mcdigi.containsSimulation()){
          art::Ptr<SimParticle> spp = mcdigi.earlyStrawGasStep()->simParticle();
          // see if this particle has already been found; if so, increment, if not, add it
          bool found(false);
          for(auto& spc : sct ) {
            if(spc._spp == spp ){
              found = true;
              spc.append(spp,isactive);
              break;
            }
          }
          if(!found)sct.push_back(spcount(spp,isactive));
        }
      }
      // sort by # of contributions
      sort(sct.begin(),sct.end(),spcountcomp());
    }

    void findMCSteps(const StepPointMCCollection& mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids, std::vector<MCStepItr>& steps) {
      steps.clear();
      // Loop over the step points, and find the one corresponding to the given detector
      for( MCStepItr imcs =mcsteps.begin();imcs!= mcsteps.end();imcs++){
        if(vids.size() == 0 ||  (imcs->trackId() == trkid && find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
          steps.push_back(imcs);
        }
      }
      // sort these in time
      sort(steps.begin(),steps.end(),timecomp());
    }

    void primaryRelation(PrimaryParticle const& primary,
        StrawDigiMCCollection const& sdmccol, std::vector<StrawDigiIndex> const& indices,
        art::Ptr<SimParticle>& primarysim, unsigned& nprimary, MCRelationship& mcrel) {
      // reset
      primarysim = art::Ptr<SimParticle>();
      nprimary = 0;
      mcrel = MCRelationship();
      // loop over primary sim particles
      for( auto spp : primary.primarySimParticles()){
        unsigned count(0);
        art::Ptr<SimParticle> sp;
        for(auto sdi : indices) {
          // find relation of this digi to the primary
          MCRelationship prel(spp,sdmccol.at(sdi).earlyStrawGasStep()->simParticle());
          // count the highest relationship for these digis
          if(prel == mcrel && prel != MCRelationship::none)
            count++;
          else if(prel > mcrel){
            mcrel = prel;
            count = 1;
            sp = sdmccol.at(sdi).earlyStrawGasStep()->simParticle();
          }
        }
        if(count > nprimary){
          nprimary = count;
          primarysim = sp;
        }
      }
    }
  }
}

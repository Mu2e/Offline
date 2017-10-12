#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "TrackerMC/inc/StrawHitDiag.hh"

#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e {

   void StrawHitDiag::init()
   {
      art::ServiceHandle<art::TFileService> tfs;
      _shdiag =tfs->make<TTree>("shdiag","StrawHit diagnostics");
      _shdiag->Branch("shid",&_shid);          // strawhit ID
      _shdiag->Branch("peakfit",&_peakfit);    // ADC waveform fit information
      _shdiag->Branch("edep",&_edep,"edep/F"); // reconstructed deposited energy
      _shdiag->Branch("time",&_time,"time/F"); // reconstructed deposited energy
      _shdiag->Branch("dt",&_dt,"dt/F");       // reconstructed deposited energy
      _shdiag->Branch("mcinfo",&_shmcinfo);    // MC truth info
   }



   //------------------------------------------------------------------------------------------
   void StrawHitDiag::fill(const Straw& straw, const StrawHit &newhit,
                           const TrkChargeReco::PeakFitParams& params, 
                           const StrawDigiMCCollection* mcdigis, int isd)
   {       
       _shid    = SHID( straw.id() );
       _peakfit = params;
       _edep    = newhit.energyDep();
       _time    = newhit.time();
       _dt      = newhit.dt();

       if (mcdigis != 0) fillDiagMC( straw, (*mcdigis)[isd]);

       _shdiag->Fill();
   }

   void StrawHitDiag::fillDiagMC(const Straw& straw, const StrawDigiMC& mcdigi)
   {
      StrawDigi::TDCChannel itdc = StrawDigi::zero;
      if (!mcdigi.hasTDC(itdc)) itdc = StrawDigi::one;
      const art::Ptr<StepPointMC>& spmcp = mcdigi.stepPointMC(itdc);
      const art::Ptr<SimParticle>& spp   = spmcp->simParticle();

      // set of unique particles contributing to this hit 
      std::set<art::Ptr<SimParticle> > spptrs;   
      for(auto imcs = mcdigi.stepPointMCs().begin(); imcs!= mcdigi.stepPointMCs().end(); ++ imcs)
      {
         // if the simParticle for this step is the same as the one which fired the discrim, add the energy
         if( (*imcs)->simParticle() == spp ) _shmcinfo._threshenergy += (*imcs)->eDep();	
         spptrs.insert((*imcs)->simParticle()); 
      }

      CLHEP::Hep3Vector mcsep  = spmcp->position()-straw.getMidPoint();
      CLHEP::Hep3Vector dir    = spmcp->momentum().unit();
      CLHEP::Hep3Vector mcperp = (dir.cross(straw.getDirection())).unit();
      double dperp = mcperp.dot(mcsep);


      _shmcinfo._energy     = mcdigi.energySum();
      _shmcinfo._trigenergy = mcdigi.triggerEnergySum();     
      _shmcinfo._xtalk      = straw.index() != spmcp->strawIndex();
      _shmcinfo._threshenergy = 0.0;
      _shmcinfo._nmcpart = spptrs.size();
      _shmcinfo._pdg = spp->pdgId();
      _shmcinfo._proc = spp->creationCode();
      if(spp->genParticle().isNonnull())
        _shmcinfo._gen = spp->genParticle()->generatorId().id();
      else
        _shmcinfo._gen = -1;

      _shmcinfo._mom    = spmcp->momentum().mag();    
      _shmcinfo._dperp = fabs(dperp);
      _shmcinfo._ambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
      _shmcinfo._len = mcsep.dot(straw.getDirection());

   }

      


}


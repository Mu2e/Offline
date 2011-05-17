//
// Read a ToyGenParticleCollection that contains the position
// of pions that have ranged out in the stopping target foils.
// Generate new ToyGenParticleCollection that contains positrons
// from pi+ -> e+ nu decay that originate from the positions at
// which the pions stopped.
//
// $Id: EplusFromStoppedPion_module.cc,v 1.2 2011/05/17 22:22:46 wb Exp $
// $Author: wb $
// $Date: 2011/05/17 22:22:46 $
//
// Original author Rob Kutschke.
//
// Notes:
//
// 1) At the writing of this code, it was not possible to access the particle data table in
//    the constructor.  This is because the ConditionsService is not initialized until 
//    begin run time.  The service will be rewritten to initialize in its ctor those parts that can
//    be initialized without a run number.  At that time we can move the PDT related stuff to
//    the c'tor from the beginRun method.
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes.
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "GeneralUtilities/inc/TwoBodyKinematics.hh"
#include "GeneralUtilities/inc/pow.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"

// Root includes.
#include "TNtuple.h"
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class EplusFromStoppedPion : public art::EDProducer {
  public:
    explicit EplusFromStoppedPion(fhicl::ParameterSet const& pset);
    virtual ~EplusFromStoppedPion() { }

    void beginRun(art::Run& run);
    void produce( art::Event& event);

  private:

    // The label of the module that put the input collection into the event.
    std::string inputModuleLabel_;

    // Generator of unit vectors over a section of a unit sphere.
    RandomUnitSphere randomUnitSphere_;

    // Are histograms enabled.
    bool doHistograms_;

    // Limit on the number of entries in the ntuple (to manage file sizes).
    int maxNtup_;

    // The rest mass of an electron.
    double me_;
    
    // The momentum of the electron, in the pi rest frame, from pi -> e nu decay.
    double pe_;

    TNtuple* nt_;
    TH1F*    hzPos_;
    TH1F*    hcz_;
    TH1F*    hphi_;

    // Number of points in the scan along z.
    int nz_;

    // Limits and step size of z iteration, in mm.
    double z0_;
    double zend_;

  };

  EplusFromStoppedPion::EplusFromStoppedPion(fhicl::ParameterSet const& pset):

    // Run time arguments from the pset.
    inputModuleLabel_(pset.get<string>("inputModuleLabel")),

    randomUnitSphere_(createEngine( get_seed_value(pset)),
                      pset.get<double>("czmin"),
                      pset.get<double>("czmax")
                      ),

    doHistograms_(pset.get<bool>("doHistograms",true)),

    maxNtup_(pset.get<bool>("maxNtup",20000)),

    // Other data members used for generation.
    me_(),
    pe_(),

    // Data members used for histogramming.
    nt_(0),
    hzPos_(0),
    hcz_(0),
    hphi_(0),
    nz_(1443),
    z0_(2930.),
    zend_(17360.){

    // What does this module produce?
    produces<ToyGenParticleCollection>();

  }

  void EplusFromStoppedPion::beginRun(art::Run& run){

    // Get the positron and pi+ masses from the particle data table.  See note 1.
    ConditionsHandle<ParticleDataTable> pdt("ignored");
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_plus).ref();
    me_ = e_data.mass().value();

    const HepPDT::ParticleData& pi_data = pdt->particle(PDGCode::pi_plus).ref();
    double pimass = pi_data.mass().value();

    // Compute the momentum of the decay positron.
    double mneutrino(0.);
    TwoBodyKinematics decay( pimass, me_, mneutrino);
    pe_ = decay.p();

    // On the first run, book the histograms and ntuples.
    if ( !doHistograms_ ) return;
    if ( nt_ ) return;

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    hzPos_  =  tfs->make<TH1F>( "hzPos", "Z of Stopping Position;[mm]", 1800, 5400., 6400. );

    hcz_  = tfs->make<TH1F>( "hcz",  "Cos(theta) of generated momentum", 100, -1., 1.);
    hphi_ = tfs->make<TH1F>( "hphi", "Azimuth of generated momentum;(radians)", 100, -M_PI, M_PI );

    nt_  = tfs->make<TNtuple>( "nt", "Positions of the Stopped Pions",
                               "x:y:z:t:dz:r");

  }

  int findTarget( double z ){
    static double zmid=5471.;
    static double dz(50.);
    static int nTarget(17);
    static double z0=zmid-8.*dz;
    for ( int i=0; i<nTarget; ++i){
      double zt = z0 + i*dz;
      double delta = z-zt;
      if ( std::abs(delta) <= 0.1000001 ) return i;
    }
    return -1;
  }

  void
  EplusFromStoppedPion::produce(art::Event& event) {

    auto_ptr<ToyGenParticleCollection> output(new ToyGenParticleCollection);

    // Get handles to the generated and simulated particles.
    art::Handle<ToyGenParticleCollection> genHandle;
    event.getByLabel(inputModuleLabel_,genHandle);
    ToyGenParticleCollection const& gens(*genHandle);

    for ( ToyGenParticleCollection::const_iterator i=gens.begin();
          i!=gens.end(); ++i ){
      ToyGenParticle const& stoppedGen = *i;

      CLHEP::Hep3Vector p = randomUnitSphere_.fire(pe_);
      double ee = sqrt(pe_*pe_+ me_*me_);

      output->push_back( ToyGenParticle( PDGCode::e_plus, 
                                         GenId::fromG4BLFile,
                                         stoppedGen.position(),
                                         CLHEP::HepLorentzVector(p.x(), p.y(), p.z(), ee),
                                         stoppedGen.time() )
                         );

      if ( !doHistograms_ ) continue;

      ToyGenParticle const& gen = output->back();

      float buf[nt_->GetNvar()];

      int itarget = findTarget(gen.position().z());
      double dz = 5071. + itarget*50. - gen.position().z();

      // Local coordinates.
      double xl=gen.position().x()+3904.;
      double yl=gen.position().y();
      buf[0] = xl;
      buf[1] = yl;
      buf[2] = gen.position().z();
      buf[3] = gen.time();
      buf[4] = dz;
      buf[5] = sqrt(xl*xl + yl*yl);
      nt_->Fill(buf);

      hzPos_->Fill(gen.position().z());
      hcz_->Fill(gen.momentum().cosTheta());
      hphi_->Fill(gen.momentum().phi());


    }

    event.put(output);
 
  } // end of ::analyze.

}

using mu2e::EplusFromStoppedPion;
DEFINE_ART_MODULE(EplusFromStoppedPion);

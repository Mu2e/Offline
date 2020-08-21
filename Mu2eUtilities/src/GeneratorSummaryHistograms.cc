//
// Make histograms summarizing the information in the event generator.
//
//
// Contact person Rob Kutschke
//

#include <cmath>
#include <vector>

#include "CLHEP/Vector/LorentzVector.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "GeneralUtilities/inc/Binning.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"
#include "TH1.h"
#include "TH2.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
// Framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/Exception.h"

// Mu2e includes
#include "Mu2eUtilities/inc/GeneratorSummaryHistograms.hh"

namespace mu2e {

  GeneratorSummaryHistograms::GeneratorSummaryHistograms():
    detSys_(0),
    hMultiplicity_(0),
    hgenId_(0),
    hp_(0),
    hpt_(0),
    hke_(0),
    hcz_(0),
    hphi_(0),
    hradius_(0),
    hz_(0),
    htime_(0),
    hxy_(0),
    hrz_(0){
  }

  // Book histograms in the root directory for the current module.
  void  GeneratorSummaryHistograms::book( ){
    art::ServiceHandle<art::TFileService> tfs;
    book(*tfs);
  }

  // Book histograms in the subdirectory, given by the relativePath; that path is
  // relative to the root TFileDirectory for the current module.
  void GeneratorSummaryHistograms::book ( std::string const& relativePath ){

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( relativePath.c_str() );
    book (tfdir);
  }

  // Book the histograms.
  void GeneratorSummaryHistograms::book( art::TFileDirectory& tfdir ){

    // Fixme: both of the GeomHandles may someday change at round boundaries.  Update the code to deal with this.

    int nId(GenId::lastEnum);
    GeomHandle<StoppingTarget> target;
    Binning bins  = zBinningForFoils(*target,7);
    Binning bins2 = zBinningForFoils(*target,3);

    detSys_ = GeomHandle<DetectorSystem>().get();

    hMultiplicity_   = tfdir.make<TH1F>( "hMultiplicity",  "Generator Multiplicity",          20,    0.,    20. );
    hgenId_          = tfdir.make<TH1F>( "hgenId",         "GeneratorId",                    nId,    0.,   nId  );
    hp_              = tfdir.make<TH1F>( "hp",             "Momentum;(MeV)",                 150,    0.,   300. );
    hpt_             = tfdir.make<TH1F>( "hpt",            "Pt;(MeV)",                       150,    0.,   300. );
    hke_             = tfdir.make<TH1F>( "hke",            "Kinetic Energy;(MeV)",           150,    0.,   150. );
    hcz_             = tfdir.make<TH1F>( "hcz",            "Cos(theta)",                     100,   -1.,     1. );
    hphi_            = tfdir.make<TH1F>( "hphi",           "Azimuth;(radians)",              100, -M_PI,  M_PI  );
    hradius_         = tfdir.make<TH1F>( "hradius",        "Radius;(mm)",                    150,    0.,   150. );
    htime_           = tfdir.make<TH1F>( "htime",          "Time;(ns)",                      150,    0.,  1700. );

    hz_              = tfdir.make<TH1F>( "hz",             "z Tracker Coordinates;(mm)",
                                         bins.nbins(), bins.low(), bins.high() );

    hxy_             = tfdir.make<TH2F>( "hxyPos",         "Conversion Electron (x,y) at Production;(mm)",
                                         60,  -120., 120., 60, -120., 120. );
    hrz_             = tfdir.make<TH2F>( "hrzPos",         "Conversion Electron (z,r) at Production;(mm)",
                                         bins2.nbins(), bins2.low(), bins2.high(), 60, 0., 120. );

  } // end GeneratorSummaryHistograms::book

  void GeneratorSummaryHistograms::fill( GenParticleCollection const& gens ){

    hMultiplicity_->Fill( gens.size() );


    for ( GenParticleCollection::const_iterator i=gens.begin(), e=gens.end();
          i != e; ++i ){
      GenParticle const& gen(*i);

      // Position in the mu2e system.
      CLHEP::Hep3Vector const& mu2ePos(gen.position());

      // Position in the detector system.
      const CLHEP::Hep3Vector pos(detSys_->toDetector(mu2ePos));

      // Momentum
      CLHEP::Hep3Vector const&       p(gen.momentum().vect());
      CLHEP::HepLorentzVector const& p4(gen.momentum());


      double r(pos.perp());
      double z(pos.z());

      // Kinetic energy.
      double ke( p4.e() - p4.m() );

      // Direction cosine.
      double cz = p.cosTheta();

      hgenId_ ->Fill( gen.generatorId().id() );
      hp_     ->Fill( p.mag() );
      hpt_    ->Fill( p.perp() );
      hke_    ->Fill( ke );
      hcz_    ->Fill( cz );
      hphi_   ->Fill( p.phi()  );
      hradius_->Fill( r );
      hz_     ->Fill( z );
      htime_  ->Fill( gen.time() );
      hxy_    ->Fill( pos.x(), pos.y() );
      hrz_    ->Fill( z, r );

    } // end GeneratorSummaryHistograms::fill
  }

} // end namespace mu2e

// art tool to sample starting positions/times, according to some rule

#ifndef EventGenerator_PositionSamplerTool_hh
#define EventGenerator_PositionSamplerTool_hh

// stl
#include <utility>
#include <vector>

// art
#include <art/Framework/Principal/Event.h>

// canvas
#include <canvas/Persistency/Common/Ptr.h>

// clhep
#include <CLHEP/Vector/LorentzVector.h>

// mu2e
#include <Offline/MCDataProducts/inc/SimParticle.hh>

namespace mu2e{
  // in general, need to propagate SimParticles back to midstage generators
  // for provenance tracking
  typedef art::Ptr<SimParticle> SimParticlePtr ;
  typedef std::pair<SimParticlePtr,CLHEP::HepLorentzVector> ParticlePositionPair;
  typedef std::vector< SimParticlePtr > SimParticlePtrVector;

  class PositionSamplerTool{
    public:
      PositionSamplerTool(){ /**/ };
     ~PositionSamplerTool(){ /**/ };

      virtual ParticlePositionPair Sample(const SimParticlePtrVector&) = 0;
    protected:
      /**/
    private:
      /**/
  };
}; // namespace mu2e

#endif

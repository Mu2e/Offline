#ifndef BTrkHelper_ParticleInfo_hh
#define BTrkHelper_ParticleInfo_hh

//
// Supply information about particles to the btrk code.
//   - btrk code makes calls indexed by TrkParicle::type.
//   - PDT is indexed by PDGCode::type
// This code looks after the translation and caches results.
//

#include "Offline/GlobalConstantsService/inc/ParticleData.hh"

#include "BTrk/BaBar/ParticleInfoInterface.hh"

#include <string>
#include <map>

namespace mu2e {

  class ParticleDataList;

  class ParticleInfo : public ParticleInfoInterface {

  public:

    ParticleInfo();

    double mass  ( TrkParticle::type id ) const override {
      return getParticle(id)->mass();
    }

    double charge( TrkParticle::type id ) const override {
      return getParticle(id)->charge();
    }

    std::string name  ( TrkParticle::type id ) const override {
      return getParticle(id)->name();
    }

  private:

    //the following has to be mutable because BTrk holds this
    //with a const pointer

    // Guaranteed valid throughout the job.
    ParticleDataList const& pdt_;

    // Local cache of the information for particles that we care about;
    // indexed by TrkParticle::type, not by PDG::id.
    mutable std::map<TrkParticle::type,ParticleData const *> table_;

    // Find particle data in the local cache; fault to the full cache as needed.
    ParticleData const*  getParticle( TrkParticle::type ) const;

  };
}

#endif /* btrkHelper_ParticleInfo_hh */

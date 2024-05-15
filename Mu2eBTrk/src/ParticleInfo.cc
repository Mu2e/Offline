#include "Offline/Mu2eBTrk/inc/ParticleInfo.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

mu2e::ParticleInfo::ParticleInfo():pdt_(*GlobalConstantsHandle<ParticleDataList>()){
}

mu2e::ParticleData const*
mu2e::ParticleInfo::getParticle( TrkParticle::type id ) const{

  // First, check the local cache.
  auto q = table_.find(id);
  if ( q != table_.end() ) return q->second;

  // If not present in local cache, translate from TrkParticle::type
  // to PDGCode::type and fill the local cache.
  ParticleData const* p(nullptr);
  switch (id) {
    case TrkParticle::e_minus: {
      p = &pdt_.particle(PDGCode::e_minus);
      break;
    }

    case TrkParticle::e_plus: {
      p = &pdt_.particle(PDGCode::e_plus);
      break;
    }

    case TrkParticle::mu_minus: {
      p = &pdt_.particle(PDGCode::mu_minus);
      break;
    }

    case TrkParticle::mu_plus: {
      p = &pdt_.particle(PDGCode::mu_plus);
      break;
    }

    case TrkParticle::pi_minus: {
      p = &pdt_.particle(PDGCode::pi_minus);
      break;
    }

    case TrkParticle::pi_plus: {
      p = &pdt_.particle(PDGCode::pi_plus);
      break;
    }

    case TrkParticle::K_minus: {
      p = &pdt_.particle(PDGCode::K_minus);
      break;
    }

    case TrkParticle::K_plus: {
      p = &pdt_.particle(PDGCode::K_plus);
      break;
    }

    case TrkParticle::anti_p_minus: {
      p = &pdt_.particle(PDGCode::anti_proton);
      break;
    }

    case TrkParticle::p_plus: {
      p = &pdt_.particle(PDGCode::proton);
      break;
    }

    default: {
      throw cet::exception("RANGE")
        << "ParticleInfo::getParticle unrecognized TrkParticle type: "
        << id;
    }

  }

  table_[id] = p;

  return p;

}

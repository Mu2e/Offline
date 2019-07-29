#include "Mu2eBTrk/inc/ParticleInfo.hh"
#include "DataProducts/inc/PDGCode.hh"

mu2e::ParticleInfo::ParticleInfo():pdt_(){
}

HepPDT::ParticleData const*
mu2e::ParticleInfo::getParticle( TrkParticle::type id ) const{

  // First, check the local cache.
  auto q = table_.find(id);
  if ( q != table_.end() ) return q->second;

  // If not present in local cache, translate from TrkParticle::type
  // to PDGCode::type and fill the local cache.
  HepPDT::ParticleData const* p(nullptr);
  switch (id) {
    case TrkParticle::e_minus: {
      p = &pdt_->particle(PDGCode::e_minus).ref();
      break;
    }

    case TrkParticle::e_plus: {
      p = &pdt_->particle(PDGCode::e_plus).ref();
      break;
    }

    case TrkParticle::mu_minus: {
      p = &pdt_->particle(PDGCode::mu_minus).ref();
      break;
    }

    case TrkParticle::mu_plus: {
      p = &pdt_->particle(PDGCode::mu_plus).ref();
      break;
    }

    case TrkParticle::pi_minus: {
      p = &pdt_->particle(PDGCode::pi_minus).ref();
      break;
    }

    case TrkParticle::pi_plus: {
      p = &pdt_->particle(PDGCode::pi_plus).ref();
      break;
    }

    case TrkParticle::K_minus: {
      p = &pdt_->particle(PDGCode::K_minus).ref();
      break;
    }

    case TrkParticle::K_plus: {
      p = &pdt_->particle(PDGCode::K_plus).ref();
      break;
    }

    case TrkParticle::anti_p_minus: {
      p = &pdt_->particle(PDGCode::anti_p_minus).ref();
      break;
    }

    case TrkParticle::p_plus: {
      p = &pdt_->particle(PDGCode::p_plus).ref();
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

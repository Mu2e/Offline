// $Id: GenSimParticleLink.hh,v 1.1 2012/04/18 22:57:19 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/04/18 22:57:19 $
//
// Original author Gianni Onorato


#ifndef GenSimParticleLink_hh
#define GenSimParticleLink_hh

#include "canvas/Persistency/Common/Assns.h"

namespace mu2e {
  class GenParticle;
  class SimParticle;

  typedef art::Assns<GenParticle, SimParticle>  GenSimParticleLink;
}

#endif/*GenSimParticleLink_hh*/

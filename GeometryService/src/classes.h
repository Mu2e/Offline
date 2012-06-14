// $Id: classes.h,v 1.2 2012/06/14 20:33:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/06/14 20:33:09 $
//
// Original author Andrei Gaponenko
//

#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::WorldG4>;
template class art::Wrapper<mu2e::Mu2eEnvelope>;

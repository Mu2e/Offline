//
// Build a dictionary.
//
// $Id: classes.h,v 1.4 2014/04/08 04:25:46 murat Exp $
// $Author: murat $
// $Date: 2014/04/08 04:25:46 $
//
// Original author Rob Kutschke
//

#include <vector>
//#include <array>

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"

#include "CalPatRec/inc/AlgorithmID.hh"
#include "CalPatRec/inc/AlgorithmIDCollection.hh"

#include "CalPatRec/inc/CalTimePeak.hh"
#include "KalmanTests/inc/TrkDef.hh"

// Cannot use the typedefs in here - not sure why.

template class art::Ptr<mu2e::AlgorithmID>;
template class std::vector<art::Ptr<mu2e::AlgorithmID> >;
template class art::Wrapper<mu2e::AlgorithmIDCollection>;

template class art::Ptr<mu2e::CalTimePeak>;
template class std::vector<art::Ptr<mu2e::CalTimePeak> >;
template class art::Wrapper<mu2e::CalTimePeakCollection>;

template class art::Ptr<mu2e::hitIndex>;
template class std::vector<art::Ptr<mu2e::hitIndex> >;
template class art::Wrapper<mu2e::hitIndexCollection>;

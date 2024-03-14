///////////////////////////////////////////////////////////////////////////////
// PhiZSeedFinderAlg uses a face-based (instead of a panel-based) data organization
///////////////////////////////////////////////////////////////////////////////
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"

using CalPatRec::ChannelID;

namespace mu2e {

  using namespace PhiZSeedFinderTypes;

  PhiZSeedFinderAlg::PhiZSeedFinderAlg(const fhicl::Table<PhiZSeedFinderAlg::Config>& config, Data_t* Data) :
    _debugLevel            (config().debugLevel()       ),
    _diagLevel             (config().diagLevel()        ),
    _testOrder             (config().testOrder()        )
  {

    _data    = Data;
//-----------------------------------------------------------------------------
// as the algorithm is supposed to work only on time clusters, make sure that
// only one time bin is defined.
// could think of further simplification down the road
//-----------------------------------------------------------------------------
    _timeBin = 2000;

    printf("PhiZSeedFinderAlg created\n");
  }

//------------------------------------------------------------------------------
// all hits belong to the same time cluster, order them in Z and in time
//-----------------------------------------------------------------------------
  int PhiZSeedFinderAlg::orderHits(const TimeCluster* Tc) {
    ChannelID cx, co;
//-----------------------------------------------------------------------------
// vector of pointers to CH, ordered in time. Initial list is not touched
//-----------------------------------------------------------------------------
    _data->_nComboHits = Tc->nhits();

    _data->_v.resize(_data->_nComboHits);

    for (int i=0; i<_data->_nComboHits; i++) {
      const StrawHitIndex& ind = Tc->hits().at(i);
      _data->_v[i] = &(*_data->chcol)[ind];
    }

    std::sort(_data->_v.begin(), _data->_v.end(),
              [](const ComboHit*& a, const ComboHit*& b) { return a->time() < b->time(); });
//-----------------------------------------------------------------------------
// at this point hits in '_v' are already ordered in time
//-----------------------------------------------------------------------------
    for (int ih=0; ih<_data->_nComboHits; ih++) {
      const ComboHit* ch = _data->_v[ih];

      cx.Station                 = ch->strawId().station();
      cx.Plane                   = ch->strawId().plane() % 2;
      cx.Face                    = -1;
      cx.Panel                   = ch->strawId().panel();
//-----------------------------------------------------------------------------
// get Z-ordered location
//-----------------------------------------------------------------------------
      ChannelID::orderID(&cx, &co);

      int os       = co.Station;
      int of       = co.Face;
      int op       = co.Panel;

      if (_printErrors) {
        if ((os < 0) || (os >= kNStations     )) printf(" >>> ERROR: wrong station number: %i\n",os);
        if ((of < 0) || (of >= kNFaces        )) printf(" >>> ERROR: wrong face    number: %i\n",of);
        if ((op < 0) || (op >= kNPanelsPerFace)) printf(" >>> ERROR: wrong panel   number: %i\n",op);
      }
//-----------------------------------------------------------------------------
// prototype face-based hit storage
// hits are already time-ordered - that makes it easy to define fFirst
// for each face, define multiple time bins and indices of the first and the last
// hits in each bin
//-----------------------------------------------------------------------------
      FaceZ_t* fz  = &_data->fFaceData[os][of];
      int      loc = fz->fHitData.size();

      fz->fHitData.push_back(HitData_t(ch,of));
      int time_bin = int (ch->time()/_timeBin) ;

      if (fz->fFirst[time_bin] < 0) fz->fFirst[time_bin] = loc;
      fz->fLast[time_bin] = loc;
    }

    return 0;
  }


//-----------------------------------------------------------------------------
  void  PhiZSeedFinderAlg::run(const TimeCluster* Tc) {
    orderHits(Tc);
  }

}

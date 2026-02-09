#ifndef Data_Collections_h
#define Data_Collections_h
//Cosmics:
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
//Calo:
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
//MC Products:
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
//Kalman Tracks
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
//Tracker Hits:
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
//CRV:
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
//Art/FCL:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include <TObject.h>
#include <TROOT.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <iostream>
#include <vector>
#include <tuple>
using namespace CLHEP;

namespace mu2e{
        class Data_Collections
        {
    public:
      #ifndef __CINT__
      explicit Data_Collections(){};
      Data_Collections(const Data_Collections &){};
      Data_Collections& operator=(const Data_Collections &);

      //RecoDataProducts:
      const ComboHitCollection *chcol = 0;
      const TimeClusterCollection *tccol = 0;
      const CrvRecoPulseCollection* crvcoincol = 0;
      const CosmicTrackSeedCollection* cosmiccol = 0;
      const CaloClusterCollection* clustercol = 0;
      const CaloHitCollection* cryHitcol = 0;
      const HelixSeedCollection* hseedcol = 0;
      const KalSeedCollection* kalseedcol = 0;
      std::vector<const KalSeedCollection*> track_list;
      std::vector<std::string> track_labels;
      std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple;

      //MCDataProducts:
      const MCTrajectoryCollection *mctrajcol = 0;

      virtual ~Data_Collections(){};
      #endif
    ClassDef(Data_Collections,0);
        };

}

#endif

#ifndef CalPatRec_PhiZSeedFinder_types_hh
#define CalPatRec_PhiZSeedFinder_types_hh

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

#include "TObject.h"
#include "TClonesArray.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"

#include "Offline/CalPatRec/inc/CalPatRec_enums.hh"
#include "Offline/CalPatRec/inc/Pzz_t.hh"
#include "Offline/CalPatRec/inc/HitData_t.hh"

using CalPatRec::Pzz_t;
using CalPatRec::HitData_t;

namespace mu2e {
  class Panel;
  class SimParticle;
  class PhiZSeedFinderAlg;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace PhiZSeedFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type             {fhicl::Name("tool_type"             ), fhicl::Comment("tool type: PhiZSeedFinderDiag")  };
      fhicl::Atom<int>         mcTruth               {fhicl::Name("mcTruth"               ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>         diagLevel             {fhicl::Name("diagLevel"             ), fhicl::Comment("diagnostic level")               };
      fhicl::Atom<bool>        mcDiag                {fhicl::Name("mcDiag"                ), fhicl::Comment("MC diag")                        };
      fhicl::Atom<int>         printOTracker         {fhicl::Name("printOTracker"         ), fhicl::Comment("print ordered Tracker")          };
      fhicl::Atom<int>         printComboHits        {fhicl::Name("printComboHits"        ), fhicl::Comment("print combo hits")               };
      fhicl::Atom<int>         printGoodComboHits    {fhicl::Name("printGoodComboHits"    ), fhicl::Comment("print good combo hits")          };
      fhicl::Atom<int>         printShcol            {fhicl::Name("printShcol"            ), fhicl::Comment("if 1, print shColl"             )};

      fhicl::Table<McUtilsToolBase::Config> mcUtils  {fhicl::Name("mcUtils"               ), fhicl::Comment("MC Diag plugin"                 )};
    };

    struct FaceZ_t {
      int                     fID;         // 3*face+panel, for pre-calculating overlaps

      std::vector<HitData_t>  fHitData;
      int                     fFirst[100];   // ** FIXME - need larger dimension for off-spill cosmics...
      int                     fLast [100];

      Pzz_t                   fPanel[3];
      double                  z;           //

      Pzz_t*                  Panel(int I) { return &fPanel[I]; }

      int                     nHits        () { return fHitData.size(); }
      HitData_t*              hitData      (int I) { return &fHitData     [I]; }
    };

    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      const DiskCalorimeter*        calorimeter;

      art::InputTag                 chCollTag;
      art::InputTag                 sdmcCollTag;

      const ComboHitCollection*     chcol;
      const TimeClusterCollection*  tccol;

      PhiZSeedFinderAlg*            _finder;

      int                           debugLevel;              // printout level

      int                           _nTimeClusters;
      int                           _nComboHits;
      std::vector<const ComboHit*>  _v;                      // sorted
//-----------------------------------------------------------------------------
// try to avoid looping over panels
//-----------------------------------------------------------------------------
      FaceZ_t                       fFaceData   [kNStations][kNFaces];
      int                           stationUsed [kNStations];
//-----------------------------------------------------------------------------
// station #2 is the same as station #0 etc...
//-----------------------------------------------------------------------------
      int                           panelOverlap[2][12][12];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
      Data_t();
      ~Data_t();

      FaceZ_t* faceData   (int Station, int Face) { return &fFaceData[Station][Face]; }

      void     InitEvent(const art::Event* Evt, int DebugLevel);
      void     InitGeometry();

    };
  }
}
#endif

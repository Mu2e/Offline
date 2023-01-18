#ifndef __CalPatRec_DeltaFinder_types_hh__
#define __CalPatRec_DeltaFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

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

#include "Offline/CalPatRec/inc/DeltaFinder_structures.hh"
#include "Offline/CalPatRec/inc/DeltaSeed.hh"
#include "Offline/CalPatRec/inc/DeltaCandidate.hh"

namespace mu2e {
  class Panel;
  class SimParticle;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace DeltaFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type             {fhicl::Name("tool_type"             ), fhicl::Comment("tool type: DeltaFinderDiag")     };
      // fhicl::Atom<art::InputTag> spmcCollTag         {fhicl::Name("spmcCollTag"           ), fhicl::Comment("StepPointMC coll tag"      )     };
      fhicl::Atom<int>         mcTruth               {fhicl::Name("mcTruth"               ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>         diagLevel             {fhicl::Name("diagLevel"             ), fhicl::Comment("diagnostic level")               };
      fhicl::Atom<bool>        mcDiag                {fhicl::Name("mcDiag"                ), fhicl::Comment("MC diag")                        };
      fhicl::Atom<int>         printOTracker         {fhicl::Name("printOTracker"         ), fhicl::Comment("print ordered Tracker")          };
      fhicl::Atom<int>         printComboHits        {fhicl::Name("printComboHits"        ), fhicl::Comment("print combo hits")               };
      fhicl::Atom<int>         printElectrons        {fhicl::Name("printElectrons"        ), fhicl::Comment("print electrons")                };
      fhicl::Atom<int>         printElectronsHits    {fhicl::Name("printElectronsHits"    ), fhicl::Comment("print electron hits")            };
      fhicl::Atom<int>         printElectronsMinNHits{fhicl::Name("printElectronsMinNHits"), fhicl::Comment("minNhhits for printed electrons")};
      fhicl::Atom<float>       printElectronsMaxFReco{fhicl::Name("printElectronsMaxFReco"), fhicl::Comment("maxFReco for printed electrons" )};
      fhicl::Atom<float>       printElectronsMinMom  {fhicl::Name("printElectronsMinMom"  ), fhicl::Comment("min mom for printed electrons"  )};
      fhicl::Atom<float>       printElectronsMaxMom  {fhicl::Name("printElectronsMaxMom"  ), fhicl::Comment("max mom for printed electrons"  )};
      fhicl::Atom<int>         printDeltaSeeds       {fhicl::Name("printDeltaSeeds"       ), fhicl::Comment("if 1, print delta seeds"        )};
      fhicl::Atom<int>         printDeltaCandidates  {fhicl::Name("printDeltaCandidates"  ), fhicl::Comment("if 1, print delta candidates"   )};
      fhicl::Atom<int>         printShcol            {fhicl::Name("printShcol"            ), fhicl::Comment("if 1, print shColl"             )};

      fhicl::Table<McUtilsToolBase::Config> mcUtils{fhicl::Name("mcUtils"       ), fhicl::Comment("MC Diag plugin") };
    };

//-----------------------------------------------------------------------------
// data structure passed to the diagnostics plugin
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      const DiskCalorimeter*        calorimeter;

      art::InputTag                 chCollTag;
      art::InputTag                 chfCollTag;
      art::InputTag                 sdmcCollTag;

      const ComboHitCollection*     chcol;
      const StrawHitFlagCollection* chfColl;                 // input  combohit flags
      StrawHitFlagCollection*       outputChfColl;           // output combohit flags
      const TimeClusterCollection*  tpeakcol;

      int                           debugLevel;              // printout level

      TClonesArray*                 listOfSeeds[kNStations]; // seeds with the first station being this

      std::vector<DeltaCandidate>   listOfDeltaCandidates;

      PanelZ_t                      oTracker   [kNStations][kNFaces][kNPanelsPerFace];
      int                           stationUsed[kNStations];
//-----------------------------------------------------------------------------
// station #2 is the same as station #0 etc...
//-----------------------------------------------------------------------------
      int                           panelOverlap[2][12][12];

      float                         stationToCaloTOF[2][20];
      float                         faceTOF[80];
      float                         meanPitchAngle;

      int                           fNSeeds;
      int                           nseeds_per_station[kNStations];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
      Data_t();
      ~Data_t();

      DeltaCandidate* deltaCandidate(int I)              { return &listOfDeltaCandidates[I]; }
      DeltaSeed*      deltaSeed     (int Station, int I) {
        return (DeltaSeed*) listOfSeeds[Station]->UncheckedAt(I);
      }


      void InitEvent(const art::Event* Evt, int DebugLevel);
      void InitGeometry();

      int  NSeedsTot() { return fNSeeds; }

      int  NSeeds(int Station) {
        return listOfSeeds[Station]->GetEntriesFast();
      }

      DeltaSeed* NewDeltaSeed(int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1) {
        int loc = nseeds_per_station[Station];
        DeltaSeed* ds = new ((*listOfSeeds[Station])[loc]) DeltaSeed(loc,Station,Face0,Hd0,Face1,Hd1);
        nseeds_per_station[Station] += 1;
        fNSeeds++;
        return ds;
      }

      static void orderID           (ChannelID* X, ChannelID* Ordered);
      static void deOrderID         (ChannelID* X, ChannelID* Ordered);

      void       printHitData       (HitData_t*      HitData, const char* Option = "");
      void       printDeltaSeed     (DeltaSeed*      Seed   , const char* Option = "");
      void       printDeltaCandidate(DeltaCandidate* Delta  , const char* Option = "");

      void       testOrderID  ();
      void       testdeOrderID();
   };

//-----------------------------------------------------------------------------
// finally, utility functions still used by the diag tool
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result);
  }
}
#endif

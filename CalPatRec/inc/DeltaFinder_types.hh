#ifndef CalPatRec_DeltaFinder_types_hh
#define CalPatRec_DeltaFinder_types_hh

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

#include "Offline/CalPatRec/inc/DeltaFinder_structures.hh"
#include "Offline/CalPatRec/inc/DeltaSeed.hh"
#include "Offline/CalPatRec/inc/DeltaCandidate.hh"
#include "Offline/CalPatRec/inc/ProtonCandidate.hh"

namespace mu2e {
  class Panel;
  class SimParticle;
  class DeltaFinderAlg;
//-----------------------------------------------------------------------------
// assume that class T has, say, Init method, which reinitialize an object
// allocated before
// in particular, that assumes that T owns all its pointers (if any) and Init
// handles them correctly
// T is supposed to have a constructor T(int N) , wher N is the element index in the list
// light-weight and crippled version of TClonesArray
//-----------------------------------------------------------------------------
    template <class T> struct ManagedList {
      int              fN;
      int              fNAllocated;
      std::vector<T*>  fList;

      ManagedList() {
        fN          = 0;
        fNAllocated = 0;
      }

      ~ManagedList() {
        for (int i=0; i<fNAllocated; i++) delete fList[i];
      }

      T* at(int I) { return fList[I]; }

      void clear () { fN = 0; }

      void reserve (int N) { fList.reserve(N); }

      void reinitialize() {
        for (int i=0; i<fNAllocated; i++) delete fList[i];
        fList.clear();
        fN          = 0;
        fNAllocated = 0;
      }

      int N() { return fN; }

      T* New() {
        T* ds;
        if (fN < (int) fList.size()) {
//-----------------------------------------------------------------------------
// reuse already allocated slot
//-----------------------------------------------------------------------------
          ds = fList[fN];
        }
        else {
//-----------------------------------------------------------------------------
// allocate new slot
//-----------------------------------------------------------------------------
          ds = new T(fN);
          fList.push_back(ds);
          fNAllocated++;
        }
        fN++;
        return ds;
      }

    };
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace DeltaFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type             {fhicl::Name("tool_type"             ), fhicl::Comment("tool type: DeltaFinderDiag")     };
      fhicl::Atom<int>         mcTruth               {fhicl::Name("mcTruth"               ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>         diagLevel             {fhicl::Name("diagLevel"             ), fhicl::Comment("diagnostic level")               };
      fhicl::Atom<bool>        mcDiag                {fhicl::Name("mcDiag"                ), fhicl::Comment("MC diag")                        };
      fhicl::Atom<int>         printOTracker         {fhicl::Name("printOTracker"         ), fhicl::Comment("print ordered Tracker")          };
      fhicl::Atom<int>         printComboHits        {fhicl::Name("printComboHits"        ), fhicl::Comment("print combo hits")               };
      fhicl::Atom<int>         printGoodComboHits    {fhicl::Name("printGoodComboHits"    ), fhicl::Comment("print good combo hits")          };
      fhicl::Atom<int>         printElectrons        {fhicl::Name("printElectrons"        ), fhicl::Comment("print electrons")                };
      fhicl::Atom<int>         printElectronsHits    {fhicl::Name("printElectronsHits"    ), fhicl::Comment("print electron hits")            };
      fhicl::Atom<int>         printElectronsMinNHits{fhicl::Name("printElectronsMinNHits"), fhicl::Comment("minNhhits for printed electrons")};
      fhicl::Atom<float>       printElectronsMaxFReco{fhicl::Name("printElectronsMaxFReco"), fhicl::Comment("maxFReco for printed electrons" )};
      fhicl::Atom<float>       printElectronsMinMom  {fhicl::Name("printElectronsMinMom"  ), fhicl::Comment("min mom for printed electrons"  )};
      fhicl::Atom<float>       printElectronsMaxMom  {fhicl::Name("printElectronsMaxMom"  ), fhicl::Comment("max mom for printed electrons"  )};
      fhicl::Atom<int>         printDeltaSeeds       {fhicl::Name("printDeltaSeeds"       ), fhicl::Comment("if 1, print delta seeds"        )};
      fhicl::Atom<int>         printDeltaCandidates  {fhicl::Name("printDeltaCandidates"  ), fhicl::Comment("if 1, print delta candidates"   )};
      fhicl::Atom<int>         printShcol            {fhicl::Name("printShcol"            ), fhicl::Comment("if 1, print shColl"             )};
      fhicl::Atom<int>         printSeedNParents     {fhicl::Name("printSeedNParents"     ), fhicl::Comment("if>0, print seeds with N hits"  )};

      fhicl::Atom<int>         printMcProtons        {fhicl::Name("printMcProtons"        ), fhicl::Comment("if>0, print MC protons"         )};
      fhicl::Atom<int>         printProtonHits       {fhicl::Name("printProtonHits"       ), fhicl::Comment("if>0, print proton hits"        )};
      fhicl::Atom<int>         printProtonSeeds      {fhicl::Name("printProtonSeeds"      ), fhicl::Comment("if>0, print proton seeds"       )};
      fhicl::Atom<int>         printProtonCandidates {fhicl::Name("printProtonCandidates" ), fhicl::Comment("if>0, print proton candidates"  )};

      fhicl::Table<McUtilsToolBase::Config> mcUtils  {fhicl::Name("mcUtils"               ), fhicl::Comment("MC Diag plugin"                 )};
    };
//-----------------------------------------------------------------------------
// data structures passed to the diagnostics plugin
//-----------------------------------------------------------------------------
    struct Pzz_t {
      int                              fID;         // 3*face+panel, for pre-calculating overlaps
      double                           wx;          // direction cosines of the wire, all wires are assumed parallel
      double                           wy;
      double                           nx;          // direction cosines of the normal to the wires, pointing outwards
      double                           ny;
      float                            z;           // Z-coordinate of the face
    };

    struct FaceZ_t {
      int                     fID;         // 3*face+panel, for pre-calculating overlaps

      std::vector<HitData_t>  fHitData;
      int                     fFirst[100];   // ** FIXME - need larger dimension for off-spill cosmics...
      int                     fLast [100];

      std::vector<HitData_t*> fProtonHitData;
      int                     fPFirst[100];  // ** FIXME - need larger dimension for off-spill cosmics...
      int                     fPLast [100];

      Pzz_t                   fPanel[3];
      double                  z;           //

      Pzz_t*                  Panel(int I) { return &fPanel[I]; }

      int                     nHits        () { return fHitData.size(); }
      int                     nProtonHits  () { return fProtonHitData.size(); }

      HitData_t*              hitData      (int I) { return &fHitData     [I]; }
      HitData_t*              protonHitData(int I) { return fProtonHitData[I]; }
    };

    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      const DiskCalorimeter*        calorimeter;

      art::InputTag                 chCollTag;
      // art::InputTag                 chfCollTag;
      art::InputTag                 sdmcCollTag;

      const ComboHitCollection*     chcol;
      // const StrawHitFlagCollection* chfColl;                 // input  combohit flags
      ComboHitCollection*           outputChColl;

      // const StrawHitFlagCollection* chfColl;                 // input  combohit flags

      DeltaFinderAlg*               _finder;

      int                           debugLevel;              // printout level

      int                           _nComboHits;
      int                           _nStrawHits;
      std::vector<const ComboHit*>  _v;                      // sorted

      ManagedList<DeltaSeed>        fListOfSeeds       [kNStations];
      std::vector<DeltaSeed*>       fListOfProtonSeeds [kNStations];
      std::vector<DeltaSeed*>       fListOfComptonSeeds[kNStations];

      std::vector<DeltaCandidate>   fListOfDeltaCandidates;

      ManagedList<ProtonCandidate>  fListOfProtonCandidates;
//-----------------------------------------------------------------------------
// try to avoid looping over panels
//-----------------------------------------------------------------------------
      FaceZ_t                       fFaceData      [kNStations][kNFaces];
      int                           stationUsed    [kNStations];
//-----------------------------------------------------------------------------
// station #2 is the same as station #0 etc...
//-----------------------------------------------------------------------------
      int                           panelOverlap[2][12][12];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
      Data_t();
      ~Data_t();

      DeltaCandidate*  deltaCandidate (int I)             { return &fListOfDeltaCandidates [I]; }

      ProtonCandidate* protonCandidate(int I)             { return fListOfProtonCandidates.at(I); }

      DeltaSeed*       deltaSeed  (int Station, int    I) { return fListOfSeeds       [Station].at(I); }
      DeltaSeed*       ComptonSeed(int Station, int    I) { return fListOfComptonSeeds[Station][I]; }
      DeltaSeed*       protonSeed (int Station, int    I) { return fListOfProtonSeeds [Station][I]; }

      FaceZ_t*         faceData   (int Station, int Face) { return &fFaceData[Station][Face]; }


      void InitEvent(const art::Event* Evt, int DebugLevel);
      void InitGeometry();

      int  nSeedsTot();

      int  NSeeds       (int Station) { return fListOfSeeds[Station].N(); }
      int  NComptonSeeds(int Station) { return fListOfComptonSeeds[Station].size(); }
      int  nProtonSeeds (int Station) { return fListOfProtonSeeds [Station].size(); }

      int  nDeltaCandidates ()        { return fListOfDeltaCandidates.size (); }
      int  nProtonCandidates()        { return fListOfProtonCandidates.N(); }

      void AddProtonSeed (DeltaSeed* Seed, int Station) { fListOfProtonSeeds [Station].push_back(Seed); }
      void AddComptonSeed(DeltaSeed* Seed, int Station) { fListOfComptonSeeds[Station].push_back(Seed); }

      void addDeltaCandidate(DeltaCandidate* Delta) { fListOfDeltaCandidates.push_back(*Delta); }

      DeltaSeed*  newDeltaSeed(int Station) {
        DeltaSeed* ds = fListOfSeeds[Station].New();
        ds->SetStation(Station);
        return ds;
      }

      ProtonCandidate* newProtonCandidate() { return fListOfProtonCandidates.New(); }

      static void orderID           (ChannelID* X, ChannelID* Ordered);
      static void deOrderID         (ChannelID* X, ChannelID* Ordered);

      void        printHitData       (HitData_t*      HitData, const char* Option = "");
      void        printDeltaSeed     (DeltaSeed*      Seed   , const char* Option = "");
      void        printDeltaCandidate(DeltaCandidate* Delta  , const char* Option = "");

      void        testOrderID  ();
      void        testdeOrderID();
   };

//-----------------------------------------------------------------------------
// finally, utility functions still used by the diag tool
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result);
  }
}
#endif

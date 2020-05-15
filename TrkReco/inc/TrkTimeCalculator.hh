//
// Functor to calculate a track t0 from a TimeCluster.  This takes into account the specific offsets of each subsystem
//

#ifndef TrkReco_TrkTimeCalculator_HH
#define TrkReco_TrkTimeCalculator_HH

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"


namespace mu2e 
{
   class TrkTimeCalculator {

     public:

       struct Config
       {
          using Name    = fhicl::Name;
          using Comment = fhicl::Comment;
          fhicl::Atom<int>    debug{        Name("debugLevel"),          Comment("Debug flag"),0 }; 
          //fhicl::Atom<bool>   useFlag{      Name("useFlag"),      Comment("Which flags to use") }; 
          //fhicl::Atom<bool>   dontUseFlag{  Name("dontUseFlag"),  Comment("Which flags to skip") }; 
          fhicl::Atom<double> avgDriftTime{ Name("AverageDriftTime"),    Comment("Average time offset for straw hits"),22.5 }; 
          fhicl::Atom<bool>   useTOTdrift{  Name("UseTOTDrift"),         Comment("Switch to use the time drift"),true }; 
          fhicl::Atom<double> beta{         Name("StrawHitBeta"),        Comment("Beta of the particle-hypothesis used"),1.0 }; 
          fhicl::Atom<double> shErr{        Name("StrawHitTimeErr"),     Comment("Effective hit time res. without TOT"),9.7 }; 
          fhicl::Atom<double> caloZOffset{  Name("TrkToCaloTimeOffset"), Comment("Location of shower max WRT COG Z position"),-120.0 }; 
          fhicl::Atom<double> caloT0Offset{ Name("CaloTimeErr"),         Comment("Time offsets for downstream particls in the calorimeter"),-0.4 }; 
          fhicl::Atom<double> caloT0Err{    Name("CaloT0Err"),    Comment("Time resolution for calorimeter"),0.5 }; 
       };       

       // parameter set must be passed in on construction
       explicit TrkTimeCalculator(fhicl::ParameterSet const&);
       explicit TrkTimeCalculator(const Config& conf);

       virtual ~TrkTimeCalculator() {};

       // update the t0 value inside a TimeCluster      
       void updateT0(TimeCluster& tc, StrawHitCollection const& shcol);
       // same, taking a HelixSeed.  This uses the pitch to make a more sophisticated z correction
       void updateT0(HelixSeed& hs, StrawHitCollection const& shcol);

       double timeOfFlightTimeOffset(double hitz,double pitch) const; // z position in the tracker frame!!!
       double strawHitTimeErr() const { return _shErr; }
       double trkToCaloTimeOffset() const { return _caloT0Offset; }
       double caloClusterTimeErr() const { return _caloT0Err; }
       // same for a ComboHit
       double comboHitTime(ComboHit const& ch,double pitch);
       // calculate the t0 for a calo cluster.
       double caloClusterTime(CaloCluster const& cc,double pitch) const;


     private:

       int _debug;
       //StrawHitFlag _useflag, _dontuseflag;
       double _avgDriftTime;  
       bool   _useTOTdrift;
       double _beta;  
       double _shErr;
       double _caloZOffset; 
       double _caloT0Offset; 
       double _caloT0Err; 
   };
 }
#endif

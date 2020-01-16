//
// Two Niveau clustering algorithm fcl configuration
//
//  Bertrand Echenard (2017) CIT
//
#ifndef ClustererFclConfig_HH
#define ClustererFclConfig_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

     struct TNTClustererConfig
     {
         using Name = fhicl::Name;
         using Comment = fhicl::Comment;
         fhicl::Atom<float>            hitDistance {     Name("HitDistance"),      Comment("Minimum cluster hit distance")  };
         fhicl::Atom<float>            seedDistance{     Name("SeedDistance"),     Comment("Minimum distance for cluster seed")  };
         fhicl::Atom<float>            clusterDiameter{  Name("ClusterDiameter"),  Comment("Average cluster diameter")  };
         fhicl::Atom<float>            clusterTime{      Name("ClusterTime"),      Comment("Average cluster time spread")  };
         fhicl::Atom<float>            deltaTimeBinMin{  Name("DeltaTimeBinMin"),  Comment("Delta time for cluster lookup")  };
         fhicl::Atom<float>            maxHitTimeDiff{   Name("MaxHitTimeDiff"),   Comment("Maximum hit cluster tme difference")  };
         fhicl::Atom<float>            maxSumDistance{   Name("MaxSumDistance"),   Comment("Maximum sum pf hit-cluster distance for convergence")  };        
         fhicl::Atom<float>            minHitError{      Name("MinHitError"),      Comment("Min value of hit error")  };
         fhicl::Atom<float>            maxDistance{      Name("MaxDistance"),      Comment("Max hit-cluster distance")  };
         fhicl::Atom<float>            timeRMS{          Name("TimeRMS"),          Comment("Cluster time RMS")  };
         fhicl::Atom<unsigned>         maxCluIterations{ Name("MaxCluIterations"), Comment("Maximum number of cluster algo iterations") };
         fhicl::Atom<bool>             medianCentroid {  Name("MedianCentroid"),   Comment("Use median to calculate cluster centroid") };
         fhicl::Atom<bool>             comboInit{        Name("ComboInit"),        Comment("Start with combo hits") };
         fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
         fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
         fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
         fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };   
     };
}  

#endif

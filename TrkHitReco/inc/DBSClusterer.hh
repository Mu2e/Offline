//
// DBScan clustering
//
//  Bertrand Echenard (2025) CIT
//
#ifndef DBSClusterer_HH
#define DBSClusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/TrkHitReco/inc/BkgClusterer.hh"
#include "fhiclcpp/types/Sequence.h"


namespace mu2e {

  class DBSClusterer : public BkgClusterer
  {
    public:

      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<unsigned>         DBSminN{          Name("DBSminExpand"),     Comment("Min number neighbors for DBScan algo") };
        fhicl::Atom<float>            hitDeltaTime{     Name("DeltaTime"),        Comment("Max time difference between hits") };
        fhicl::Atom<float>            hitDeltaZ{        Name("DeltaZ"),           Comment("Max Z difference between hits") };
        fhicl::Atom<float>            hitDeltaXY{       Name("DeltaXY"),          Comment("Max XY difference between hits") };
        fhicl::Atom<unsigned>         minClusterHits{   Name("MinClusterHits"),   Comment("Min number hits in cluster") };
        fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
        fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
        fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
        fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };
      };


      DBSClusterer(const std::optional<Config> config);
      virtual ~DBSClusterer() {};

      void          init        ();
      virtual void  findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol, int iev);
      virtual float distance    (const BkgCluster& cluster,      const ComboHit& hit) const;


    private:
      unsigned findNeighbors    (unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors);
      void     calculateCluster (BkgCluster& cluster, const ComboHitCollection& chcol);
      void     dump             (const std::vector<BkgCluster>& clusters);

      unsigned                DBSminExpand_;
      float                   deltaTime_;
      float                   deltaZ_;
      float                   deltaXY2_;
      unsigned                minClusterHits_;
      StrawHitFlag            bkgmask_;
      StrawHitFlag            sigmask_;
      bool                    testflag_;
      int                     diag_;
  };
}

#endif

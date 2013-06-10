#ifndef CaloCluster_CaloClusterTools_hh
#define CaloCluster_CaloClusterTools_hh

//
// Tools for CaloCluster objects
//
// $Id: CaloClusterTools.hh,v 1.3 2013/06/10 04:11:36 murat Exp $
// $Author: murat $
// $Date: 2013/06/10 04:11:36 $
//
// Original author G. Pezzullo, A. Luca' & G. Tassielli
//

#include "RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {

  // Forward reference.
  class VaneCalorimeter;

  class CaloClusterTools{

  private:
    // The underlying persistent cluster.
    CaloCluster const& _cluster;

    // The geometry information.
    //    VaneCalorimeter const& _calorimeter;

  public:
    CaloClusterTools(CaloCluster const &clu);

    // The underlying persistent cluster.
    CaloCluster const& cluster() const { return _cluster; }

    // Properties computed from the persistent cluster information.
    double              timeErr() const; //ns
    double    timeFasterCrystal() const; //ns
    double timeFasterCrystalErr() const; //ns
    double         energyDepErr() const; //MeV
    double            showerDir() const;
    double         errShowerDir() const;
    int                   wSize() const;
    int                   vSize() const;
    int      cryEnergydepMaxRow() const;
    int   cryEnergydepMaxColumn() const;

    void print( std::ostream& ost, bool doEndl = true ) const;

  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const CaloClusterTools& c ){
    c.print(ost,false);
    return ost;
  }


} // namespace mu2e

#endif /* CaloCluster_CaloClusterTools_hh */

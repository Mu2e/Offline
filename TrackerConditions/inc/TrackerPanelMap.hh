// clang-format off
#ifndef TrackerConditions_TrackerPanelMap_hh
#define TrackerConditions_TrackerPanelMap_hh
//
// defines cross-mapping of different indices corresponding to a tracker panel
//

// Mu2e includes
#include "Offline/DbTables/inc/TrkPanelMap.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/type_traits.h"

// C++ includes
#include <map>
#include <memory>

namespace mu2e {

  class TrackerPanelMap : public ProditionsEntity {
  public:
    enum {
      kMaxPlanes = 100,
      kMaxPanels = 600
    };
      
    typedef std::shared_ptr<TrackerPanelMap> ptr_t;
    typedef std::shared_ptr<const TrackerPanelMap> cptr_t;
    constexpr static const char* cxname = {"TrackerPanelMap"};

    TrackerPanelMap();
    virtual ~TrackerPanelMap() {}

    const TrkPanelMap::Row* panel_map_by_mnid       (int MnID) const {
      return _tpm_by_mnid[MnID];
    };

    const TrkPanelMap::Row* panel_map_by_offline_ind(int UniquePlane, int Panel) const {
      return _tpm_by_offline[UniquePlane][Panel];
    }
    
    const TrkPanelMap::Row* panel_map_by_online_ind (int DtcID, int Link ) const {
      return _tpm_by_online[DtcID][Link];
    }

    void print(std::ostream&) const override;

    void add  (const TrkPanelMap::Row& Tpm);

  private:
                                        // panel map for a given run
    std::vector<TrkPanelMap::Row> _map;
//-----------------------------------------------------------------------------
// assume less than 100 planes, 600 panels
//-----------------------------------------------------------------------------
    const TrkPanelMap::Row*  _tpm_by_mnid   [kMaxPanels];    // indexed by 'minnesota ID'
    const TrkPanelMap::Row*  _tpm_by_offline[kMaxPlanes][6]; // indexed by the offline uniquePlane and panel 
    const TrkPanelMap::Row*  _tpm_by_online [kMaxPlanes][6]; // indexed by the DTC ID and link ID
  };

}

#endif

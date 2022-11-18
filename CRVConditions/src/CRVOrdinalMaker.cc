#include "Offline/CRVConditions/inc/CRVOrdinalMaker.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>

namespace mu2e {

CRVOrdinal::ptr_t CRVOrdinalMaker::fromFcl() {
  if (_config.verbose()) {
    std::cout << "CRVOrdinalMaker::fromFcl making nominal CRVOrdinal\n";
  }

  // both maps initialized to invalid
  CRVOrdinal::OfflineMap offMap;
  for (std::size_t i = 0; i < CRVId::nROC; i++) {
    for (std::size_t j = 0; j < CRVId::nFEBPerROC; j++) {
      for (std::size_t k = 0; k < CRVId::nChanPerFEB; k++) {
        offMap[i][j][k] = CRVId::nChannels;
      }
    }
  }
  CRVOrdinal::OnlineMap onMap(CRVId::nChannels,
                              CRVROC(0, 0, CRVId::nChanPerFEB));

  // This is very nominal numbering
  // there are some channel numbers which do not correspond to SiPMs
  // and these are not handled
  // TODO this will need to be updated
  for (std::size_t channel = 0; channel < CRVId::nChannels; channel++) {
    std::size_t subchannel = channel % CRVId::nChanPerFEB;
    std::size_t FEB = (channel / CRVId::nChanPerFEB) % CRVId::nFEBPerROC;
    std::size_t ROC = (channel / CRVId::nChanPerFEB / CRVId::nFEBPerROC);
    onMap.at(channel) = CRVROC(ROC, FEB, subchannel);
    offMap.at(ROC).at(FEB).at(subchannel) = channel;
  }

  auto ptr = std::make_shared<CRVOrdinal>(onMap, offMap);
  return ptr;

}  // end fromFcl

CRVOrdinal::ptr_t CRVOrdinalMaker::fromDb() {
  if (_config.verbose()) {
    std::cout << "CRVOrdinalMaker::fromDb making CRVOrdinal\n";
  }

  // no database dependence yet, so just return nominal
  return fromFcl();

}  // end fromDb

}  // namespace mu2e

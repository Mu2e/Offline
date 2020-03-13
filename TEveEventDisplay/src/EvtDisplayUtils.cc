#include  "TEveEventDisplay/src/dict_classes/NavState.h"
#include  "TEveEventDisplay/src/dict_classes/EvtDisplayUtils.h"
#include <string>

using namespace mu2e;

  EvtDisplayUtils::EvtDisplayUtils():fTbRun(0),fTbEvt(0){}

  void EvtDisplayUtils::PrevEvent(){
    NavState::Set(kPREV_EVENT);
  }
  void EvtDisplayUtils::NextEvent(){
    NavState::Set(kNEXT_EVENT);
  }
  void EvtDisplayUtils::GotoEvent(){
    int run = std::stoi(fTbRun->GetString());
    int event = std::stoi(fTbEvt->GetString());
    NavState::SetTarget(run, event);
    NavState::Set(kGOTO_EVENT);
  }


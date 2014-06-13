//-----------------------------------------------------------------------------
//  debug
//-----------------------------------------------------------------------------
#include "Stntuple/scripts/global_vars.h"

int debug(TDebugModule* Module) {
  Module = (TDebugModule*) g.x->AddModule("TDebugModule",0);
}


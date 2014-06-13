//-----------------------------------------------------------------------------
// trigger paths: way to select JET_20 of JET_50 ot JET_70 or JET_100 
//                or something else (ELECTRON_CENTRAL_18, for example)
//                or TAU_MET:111
//-----------------------------------------------------------------------------
#include "global_vars.h"

int setup_trigger_path(const char* TriggerPath) {
  TObjArray word;
  int       l1, l2, l3, nw;

  printf("[setup_trigger_path.C]: g.L3TrigPath=%s. g.DoMc=%i\n",TriggerPath,g.DoMc);

  if (g.L3TrigPath != "") {
    // emulate for MC
    //printf("g.L3TrigPath=%s, adding TriggerModule\n",TriggerPath);
    split(TriggerPath,':',&word);
    nw = word.GetEntries();
    printf("nw = %i\n",nw);

    if (nw == 1) {
				// only the trigger path specified, use default
				// mode
      if (g.DoMc == 1) then {
				// MC: emulation
	l1 = 2;
	l2 = 2;
	l3 = 2;
      }
      else {
				// data: check
	l1 = 1;
	l2 = 1;
	l3 = 1;
      }
    }
    else {
				// trigger emulation mode is specified
      char* bits = ((TObjString*) word[1])->GetString().Data();
      printf("bits=%s\n",bits);
      l1 = int(bits[0])-0x30;
      l2 = int(bits[1])-0x30;
      l3 = int(bits[2])-0x30;
    }

    char* path = ((TObjString*) word[0])->GetString().Data();
    printf("[setup_trigger_path]: path=%s,l1 l2 l3: %i %i %i\n",path,l1,l2,l3);
    m_trig = (TTriggerModule*) g.x->AddModule("TTriggerModule",1);
    m_trig->AddPath(path,l1,l2,l3);
  }
  return 0;
}

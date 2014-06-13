#ifndef Stntuple_obj_AbsEvent_hh
#define Stntuple_obj_AbsEvent_hh

#if defined (CDF1)
#  ifndef __CINT__
#    define AbsEvent THBookStnEvent
#  else
     class   THBookStnEvent;
     typedef THBookStnEvent AbsEvent;
#  endif
#elif defined (CDF)
#  include "Experiment/Experiment.hh"

#elif defined (MU2E)

#  ifndef __CINT__
#    include "art/Framework/Principal/Event.h"
     typedef art::Event AbsEvent;
#  else
     namespace art {
       class Event;
     }
     typedef art::Event   AbsEvent;

#  endif
#endif

#endif




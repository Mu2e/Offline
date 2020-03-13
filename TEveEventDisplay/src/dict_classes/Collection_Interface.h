#ifndef Collection_Interface_h
#define Collection_Interface_h

#include <TObject.h>
#include "TROOT.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

namespace mu2e{
	class Collection_Interface{
		   #ifndef __CINT__
		   explicit Collection_Interface();
		   public:
		  	 virtual ~Collection_Interface(){};
                   private:
		  	art::Event  *_event;
		  	art::Run    *_run;
		   #endif
		   ClassDef(Collection_Interface,0);
		   
	};

}//end namespace
#endif /*Collection_Interface_h*/

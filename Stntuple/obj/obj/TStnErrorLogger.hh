#ifndef TStnErrorLogger_hh
#define TStnErrorLogger_hh

#include "TQObject.h"
#include "RQ_OBJECT.h"

#include "TNamed.h"

class TStnErrorLogger : public TNamed {
  RQ_OBJECT("TStnErrorLogger")
public:
//------------------------------------------------------------------------------
//  function members
//------------------------------------------------------------------------------
  TStnErrorLogger(const char* Name  = "StntupleErrorLogger",
		  const char* Title = "Error Logger");
  virtual ~TStnErrorLogger();
					// ****** overloaded methods of 
					// TObject

  void  Report(Int_t ErrorCode, const char* Message);  //*SIGNAL*

  void  Clear(Option_t* opt="");

  void  Print(Option_t* opt="") const ;

  ClassDef(TStnErrorLogger,0)

};

#endif





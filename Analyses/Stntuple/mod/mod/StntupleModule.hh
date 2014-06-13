//--------------------------------------------------------------------------
// Stntuple/Stntuple/StntupleModule.hh
//
// StntupleModule: base class for STNTUPLE modules. See source file for 
// more details, CVS log for revision history
//
// Nov 23 2000 P.Murat
//------------------------------------------------------------------------
#ifndef Stntuple_mod_StntupleModule
#define Stntuple_mod_StntupleModule

#include "Stntuple/mod/THistModule.hh"

class TFolder;
class TStnEvent;
class TStnErrorLogger;
class TStnDataBlock;

class StntupleModule : public THistModule {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  static TStnEvent*       fgEvent;
  static TStnErrorLogger* fgErrorLogger;
  static TFolder*         fgStntupleFolder;
//-----------------------------------------------------------------------------
// function members
//-----------------------------------------------------------------------------
public:
					// constructors and destructor

  StntupleModule(fhicl::ParameterSet const& Pset, const char* Name);

  ~StntupleModule();
				        // ****** accessors

  TStnEvent*       Event        () { return fgEvent;       }
  TStnErrorLogger* ErrorLogger  () { return fgErrorLogger; }

					// ****** modifiers
					// ****** other methods

  static TStnDataBlock*  AddDataBlock(const char* branch_name, 
				      const char* class_name,
				      Int_t      (*f)(TStnDataBlock*,AbsEvent*,Int_t),
				      Int_t       buffer_size,
				      Int_t       split_level,
				      Int_t       compression);

  static Int_t SetResolveLinksMethod(const char* BlockName, 
				     Int_t      (*f)(TStnDataBlock*,AbsEvent*,Int_t));

  void      LogError(const char* Message);
  void      LogError(char* Message);

  //  ClassDef(StntupleModule,0)
};

#endif

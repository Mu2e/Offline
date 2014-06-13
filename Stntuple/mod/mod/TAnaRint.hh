#ifndef murat_inc_TAnaRint
#define murat_inc_TAnaRint

#include "TObject.h"
#include "TRint.h"

class TAnaRint : public TObject {
private:

  TAnaRint();
  ~TAnaRint();

  static TRint*     fgRint;
  static TAnaRint*  fgInstance;


  int   fInteractiveMode;

  class  Cleaner {
  public: 
    Cleaner();
    ~Cleaner();
  };
  friend class Cleaner;

public: 

  void   GetInteractiveMode(int& Mode) { Mode = fInteractiveMode; }

  static TAnaRint* Instance(int argc=0, char** argv=0);

  static void      Delete ();

  static TRint*    Rint() { return fgRint; }

  void   SetInteractiveMode(int Mode) { fInteractiveMode = Mode; }

  ClassDef(TAnaRint,0)

};

#endif

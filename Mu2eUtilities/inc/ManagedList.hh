//-----------------------------------------------------------------------------
// P.Murat : a templated, light-weight, and "crippled" version of TClonesArray
// assumptions:
// - class T has an initialization method, say, T::Init, which properly reinitializes
//   an object allocated before
// - an object of class T owns all its pointers (if any) and T::Init handles them correctly
// - class T has a constructor T(int N), where N is the element index in the list
//
// if those assumptions hold, a call to ManagedList::clear() 'empties' the list
// without reallocating the memory
//-----------------------------------------------------------------------------
#ifndef __Mu2eUtilities_ManagedList_hh
#define __Mu2eUtilities_ManagedList_hh

namespace mu2e {

  template <class T> struct ManagedList {
    int              fN;
    int              fNAllocated;
    std::vector<T*>  fList;

    ManagedList() {
      fN          = 0;
      fNAllocated = 0;
    }

    ~ManagedList() {
      for (int i=0; i<fNAllocated; i++) delete fList[i];
    }

    T* at(int I) { return fList[I]; }

    void clear () { fN = 0; }

    void reserve (int N) { fList.reserve(N); }

    void reinitialize() {
      for (int i=0; i<fNAllocated; i++) delete fList[i];
      fList.clear();
      fN          = 0;
      fNAllocated = 0;
    }

    int N() { return fN; }

    T* New() {
      T* ds;
      if (fN < (int) fList.size()) {
//-----------------------------------------------------------------------------
// reuse already allocated slot
//-----------------------------------------------------------------------------
        ds = fList[fN];
      }
      else {
//-----------------------------------------------------------------------------
// allocate new slot
//-----------------------------------------------------------------------------
        ds = new T(fN);
        fList.push_back(ds);
        fNAllocated++;
      }
      fN++;
      return ds;
    }
  };
}
#endif

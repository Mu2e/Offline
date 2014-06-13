#ifndef Stntuple_base_THexIndex_hh
#define Stntuple_base_THexIndex_hh

class THexIndex {
public:
  int  fL;
  int  fK;
  
  THexIndex() {}
  THexIndex(int L, int K) {
    fL = L;
    fK = K;
  }

  void Set(int L, int K) { fL = L; fK = K; }

  ~THexIndex() {}
					// hexagon vertices
  static THexIndex fgPos[6];
};

#endif

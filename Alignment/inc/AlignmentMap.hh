#ifndef Alignment_AlignmentMap_HH
#define Alignment_AlignmentMap_HH
// AlignmentMap.hh
// This is the definition file for the AlignmentMap class.
// This class holds a list of associations between GeomHandle objects and
// AlignmentObj objects for use in the Mu2e Alignment service.

#include <unordered_map>

class AlignmentMap {
public:
  AlignmentMap(){}

private:
  std::unordered_map<GeomHandle,AlignmentObj> _map;

};

#endif // Alignment_AlignmentMap_HH

///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/base/TVisNode.hh"

ClassImp(TVisNode)

//_____________________________________________________________________________
TVisNode::TVisNode(const char* name):
  fName(name)
{
  fClosestObject = NULL;
  fDebugLevel    = 0;
}


//_____________________________________________________________________________
TVisNode::~TVisNode() {
}


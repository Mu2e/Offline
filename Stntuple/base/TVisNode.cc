///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/base/TVisNode.hh"

ClassImp(TVisNode)

//_____________________________________________________________________________
TVisNode::TVisNode(const char* name):
  fName(name)
{
  fClosestObject = NULL;
}


//_____________________________________________________________________________
TVisNode::~TVisNode() {
}


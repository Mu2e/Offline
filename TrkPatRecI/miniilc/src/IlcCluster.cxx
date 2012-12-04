/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: IlcCluster.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------------------------
//               Implementation of the Cluster class
// that is the base for IlcTPCcluster, IlcVXDclusterV2 and IlcTRDcluster
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "IlcCluster.h"

//ClassImp(IlcCluster)
 
//_____________________________________________________________________________
IlcCluster::IlcCluster(): 
  fY(0), fZ(0), fSigmaY2(0), fSigmaZ2(0) 
{
  //
  //default constructor
  //
  fTracks[0]=fTracks[1]=fTracks[2]=-3141593; 
}

//_____________________________________________________________________________
IlcCluster::IlcCluster(Int_t *lab, Float_t *hit): 
  fY(hit[0]), fZ(hit[1]), fSigmaY2(hit[2]), fSigmaZ2(hit[3]) 
{
  //
  // Creates a simulated cluster
  //
  fTracks[0]  = lab[0];
  fTracks[1]  = lab[1];
  fTracks[2]  = lab[2];
}

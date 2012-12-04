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

/*
$Log: IlcComplexCluster.cxx,v $
Revision 1.1  2012/12/04 00:51:27  tassiell
merge between the itracker_dev_v00 Branch and the HEAD

Revision 1.1.2.1  2012/09/29 11:07:41  ignatov
add miniilc code and modules to use it
update KalFitI.cc to use in more consistent way

Revision 1.1  2008/05/07 14:03:19  cgatto
original version from TPC(only renamed)

Revision 1.1.1.1  2008/03/11 14:53:46  vitomeg
Initial commit

Revision 1.1.1.1  2006/11/16 16:27:50  garren
IlcRoot framework

Revision 1.1.1.1  2006/11/09 16:31:38  vitomeg
Lite Version of IlcRoot to convert in IlcRoot

Revision 1.1.1.1  2006/09/04 14:56:58  vitomeg
initial revision

Revision 1.8  2004/03/30 14:09:22  kowal2
Changes due to the coding conventions

Revision 1.7  2003/11/24 09:48:28  kowal2
Changes to obey the coding conventions

Revision 1.6  2003/11/24 09:43:03  kowal2
Obsolete - removed

Revision 1.5  2003/09/29 11:27:39  kowal2
new classes added

Revision 1.3  2002/11/15 14:27:45  hristov
First version of the parallel DCH tracking (M.Ivanov)

Revision 1.2  2001/02/05 14:43:13  hristov
Compare() declared const

Revision 1.1  2000/10/05 16:17:27  kowal2
New class replacing IlcCluster


*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber clusters objects                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//Begin_Html
/*
<img src="gif/IlcDCHCluster.gif">
*/
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
// **************
//aaaaaaaaa
//aaaaaaaaa
//aaaaaaaaa
//aaaaaaaa
//aaaaaaaaaa

#include "IlcComplexCluster.h"


//ClassImp(IlcComplexCluster)
//_____________________________________________________________________________
Int_t IlcComplexCluster::Compare(const TObject * o) const
{
  //
  // compare two clusters according y coordinata
  IlcComplexCluster *cl= (IlcComplexCluster *)o;
  if (fY<cl->fY) return -1;
  if (fY==cl->fY) return 0;
  return 1;  
}

Bool_t IlcComplexCluster::IsSortable() const
{
  //
  //make IlcComplexCluster sortabale
  return kTRUE; 
}

//ClassImp(IlcDCHExactPoint)
//ClassImp(IlcDCHClusterPoint)
//ClassImp(IlcDCHTrackerPoint)
//ClassImp(IlcDCHTrackPoint)
//ClassImp(IlcDCHTrackPoint2)
//ClassImp(IlcDCHTrackPointRef)



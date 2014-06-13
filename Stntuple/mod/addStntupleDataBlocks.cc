///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/base/TStnArrayI.hh"
#include "Stntuple/base/TStnArrayF.hh"

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/mod/addStntupleDataBlocks.hh"

extern void addStntupleDataBlocks() {
//-----------------------------------------------------------------------------
//  STNTUPLE data Blocks
//-----------------------------------------------------------------------------
  TStnArrayI                 arri;
  TStnArrayF                 arrf;
  TStnHeaderBlock            header_data;
  TStnClusterBlock           cluster_data;
  TStnTrackBlock             track_data;
}

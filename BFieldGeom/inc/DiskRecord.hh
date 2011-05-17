#ifndef BFieldGeom_DiskRecord_hh
#define BFieldGeom_DiskRecord_hh
//
// The MECO GMC format magnetic field files are in a binary format.
// This struct represents one such record, including the framing fields,
// as it appears on disk in the binary files.
//
// $Id: DiskRecord.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
//
// Definition of one record as it appears on the GMC input file.
//
// 1) Header:  exclusive byte count = 24
// 2) Data     24 bytes
// 3) Trailer: exclusive byte count = 24
//
// Total record length, including the framing fields is 32 bytes.

typedef struct {

  int   head;
  float x;
  float y;
  float z;
  float bx;
  float by;
  float bz;
  int   tail;

} DiskRecord;

#endif /* BFieldGeom_DiskRecord_hh */

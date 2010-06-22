#ifndef DISKRECORD_HH
#define DISKRECORD_HH
//
// The MECO GMC format magnetic field files are in a binary format.
// This struct represents one such record, including the framing fields,
// as it appears on disk in the binary files.
//
// $Id: DiskRecord.hh,v 1.1 2010/06/22 16:44:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/22 16:44:25 $
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

#endif

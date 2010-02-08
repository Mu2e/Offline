#ifndef DISKRECORD_HH
#define DISKRECORD_HH

//
// Definition of one record as it appears on the GMC input file.
//
// 1) Header byte count = 24
// 2) Data - 24 bytes
// 3) Trailer byte count = 24
//

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

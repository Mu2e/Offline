BEGIN { FS = " ";FILE=MOD ".dat"}
{# print $8 $9
  # print "testing" MOD
  if ( match($8,MOD) ) print $10 > FILE
}

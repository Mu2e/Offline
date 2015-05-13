#
# bash functions to provide easy reference to scripts in Offline
#

jsonMaker() {
  $MU2E_BASE_RELEASE/Tools/DH/jsonMaker.py "$@"
}
export jsonMaker

samRm() {
  $MU2E_BASE_RELEASE/Tools/DH/samRm.sh "$@"
}
export samRm

samGet() {
  $MU2E_BASE_RELEASE/Tools/DH/samGet.sh "$@"
}
export samGet

samPrestage() {
  $MU2E_BASE_RELEASE/Tools/DH/samPrestage.sh "$@"
}
export samPrestage

samDatasets() {
  $MU2E_BASE_RELEASE/Tools/DH/samDatasets.sh "$@"
}
export samDatasets

samToPnfs() {
  $MU2E_BASE_RELEASE/Tools/DH/samToPnfs.sh "$@"
}
export samToPnfs

samNoChildren() {
  $MU2E_BASE_RELEASE/Tools/DH/samNoChildren.sh "$@"
}
export samNoChildren

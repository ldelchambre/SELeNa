#!/bin/bash

################################################################################
echo "*** Delchambre (2019) ***"
eval "/home/ludovic/Bureau/Travail/GraL/SELeNa/simulation $1 $2 $3 $4 $5 $6 $7 $8 [-3:1000:3] [-3:1000:3]" | awk '
(FNR == 1) {
  for(i = 1; i <= NF; i++)
    col[$i] = i;
  printf("%13s %13s %13s %13s\n", $col["x"], $col["y"], $col["mag"], "mu");
}

(FNR > 1) {
  x = $col["x"];
  y = $col["y"];
  mag = $col["mag"];
  mu = 10^(-0.4*mag);
  printf("%+13.6e %+13.6e %+13.6e %+13.6e\n", x,y,mag,mu);
}
'
################################################################################

################################################################################
echo ""
echo "*** GravLens (Keeton) ***"
echo "" | awk -v b="$1" -v e="$2" -v te="$3" -v g="$4" -v tg="$5" -v s="$6" -v xs="$7" -v ys="$8" '
BEGIN {
  q = 1.0 - e;
  p = q * sqrt( 2.0 / ( 1.0 + q^2 ) );
  bp = b * p
  sp = s * p
  tep = 90.0 + te
  print "maingrid 1 -3 3 128 -3 3 128";
  print "setlens 1 1";
  print "alpha "bp"  0 0 "e" "tep" "g" "tg" "sp" 0 1";
  print "0 0 0 0 0 0 0 0 0 0";
  print "findimg "xs" "ys;
}
' | ./lensmodel | awk '
function abs(x) {
  return (x < 0) ? -x : x;
}

BEGIN {
  nimg = -1;
}

(nimg > 0 && NF == 4) {
  x=$1;
  y=$2;
  mu=abs($3);
  mag=-2.5*log(mu) / log(10);
  printf("%+13.6e %+13.6e %+13.6e %+13.6e\n", x,y,mag,mu);
}

($0 ~ /# [0-9]+ images:  x  y  mag  tdel/) {
  nimg = $2;
  printf("%13s %13s %13s %13s\n", "x", "y", "mag", "mu");
}
'
################################################################################

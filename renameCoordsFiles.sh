#!/bin/sh

dirPrefix=`pwd`;
dirStr=${dirPrefix}/$1;

cd ${dirStr};

for dir in ${dirStr}*;
do
  ddC=${dir}/Coordonnees_ADN_ch_0; 
  cd $ddC;
  echo "Renaming the ${ddC} directory.";
  rename nsweep\= nsweep\=0 nsweep\=?.txt;
  rename nsweep\= nsweep\=0 nsweep\=??.txt;
  rename nsweep\= nsweep\=0 nsweep\=???.txt;
  rename nsweep\= nsweep\=0 nsweep\=????.txt;
  rename nsweep\= nsweep\=0 nsweep\=?????.txt;
  rename nsweep\= nsweep\=0 nsweep\=??????.txt;
  rename nsweep\= nsweep\=0 nsweep\=???????.txt;
  rename nsweep\= nsweep\=0 nsweep\=????????.txt;
  rename nsweep\= nsweep\=0 nsweep\=?????????.txt;
  rename nsweep\= nsweep\=0 nsweep\=??????????.txt;
  cd ../../;
  for ddB in ${dir}/Binding_sites_*;
  do
    cd $ddB;
    echo "Renaming the ${ddB} directory.";
    rename nsweep\= nsweep\=0 nsweep\=?.txt;
    rename nsweep\= nsweep\=0 nsweep\=??.txt;
    rename nsweep\= nsweep\=0 nsweep\=???.txt;
    rename nsweep\= nsweep\=0 nsweep\=????.txt;
    rename nsweep\= nsweep\=0 nsweep\=?????.txt;
    rename nsweep\= nsweep\=0 nsweep\=??????.txt;
    rename nsweep\= nsweep\=0 nsweep\=???????.txt;
    rename nsweep\= nsweep\=0 nsweep\=????????.txt;
    rename nsweep\= nsweep\=0 nsweep\=?????????.txt;
    rename nsweep\= nsweep\=0 nsweep\=??????????.txt;
    cd ../../;
  done
done

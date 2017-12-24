#!/bin/bash

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Wrapper for bedtools unionbedg to create a multi-sample coverage
# matrix from the output of binCov.py

#Usage statement
usage(){
cat <<EOF

usage: makeMatrix.sh [-h] [-z] [-r RESTRICT] [-o OUTFILE] SAMPLES

Helper tool to automate creation of sorted coverage matrices from
binCov.py output bed files

Positional arguments:
  SAMPLES     List of samples and coverage files (tab-delimmed)

Optional arguments:
  -h  HELP        Show this help message and exit
  -z  BGZIP       Bgzip and tabix output file
  -r  RESTRICT    Restrict output to all bins in RESTRICT file
  -o  OUTFILE     Output file (default: stdout)

EOF
}

#Parse arguments
OUTFILE=/dev/stdout
RESTRICT=0
GZ=0
while getopts ":o:r:zh" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    z)
      GZ=1
      ;;
    r)
      RESTRICT=${OPTARG}
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
SAMPLES=$1

#Check for required input
if [ -z ${SAMPLES} ] || ! [ -e ${SAMPLES} ]; then
  usage
  exit 0
fi

#Check for gzip optioned only if output file specified
if [ ${OUTFILE} == "/dev/stdout" ] && [ ${GZ} == 1 ]; then
  echo -e "\nOUTFILE required for zip-compressed output"
  usage
  exit 0
fi

#Scrub ".gz" from output filename if provided by user
if [ ${GZ} == 1 ] && [ ${OUTFILE: -3} == ".gz" ]; then
  OUTFILE=$( echo "${OUTFILE}" | sed 's/\.gz//g' )
fi

#Use bedtools to create matrix
/usr/bin/env bedtools unionbedg \
  -header \
  -names $( while read ID cov; do echo "${ID}"; done < ${SAMPLES} | paste -s -d\  ) \
  -i $( while read ID cov; do echo "${cov}"; done < ${SAMPLES} | paste -s -d\  ) | \
  sed -e 's/chrom/Chr/g' -e 's/start/Start/g' -e 's/end/End/g' > ${OUTFILE}

#Restrict OUTFILE to RESTRICT bins (if optioned)
if [ ${RESTRICT} != "0" ]; then
  head -n1 ${OUTFILE} > ${OUTFILE}2
  fgrep -v "#" ${RESTRICT} | awk '{ print $1"_"$2"_"$3 }' | fgrep -wf - \
  <( paste <( sed '1d' ${OUTFILE} | awk '{ print $1"_"$2"_"$3 }' ) \
           <( sed '1d' ${OUTFILE} ) ) | cut -f2- >> ${OUTFILE}2
  mv ${OUTFILE}2 ${OUTFILE}
fi

#Add hash to header line of OUTFILE & sort
cat <( head -n1 ${OUTFILE} | awk '{ print "#"$0 }' ) \
    <( sed '1d' ${OUTFILE} | sort -Vk1,1 -k2,2n -k3,3n | uniq ) > ${OUTFILE}2
mv ${OUTFILE}2 ${OUTFILE}

#Bgzip & tabix OUTFILE, if optioned
if [ ${GZ}==1 ] && [ ${OUTFILE} != "/dev/stdout" ]; then
  bgzip -f ${OUTFILE}
  tabix -f ${OUTFILE}.gz
fi

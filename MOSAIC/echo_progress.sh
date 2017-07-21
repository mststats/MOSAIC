if [ "$#" -ne 5 ]; then
  echo "Usage: target firstind lastind chrs L"
  exit 1
fi
target=$1
firstind=$2
lastind=$3
chrs=$4
L=$5
grep round LOGS/"$target"_"$firstind"_"$L".out | tail -n 1
tail -n 2 LOGS/"$target"_"$firstind"_"$L".out; echo
List=`ls -t RESULTS/"$target"_"$L"way_"$firstind"-"$lastind"_"$chrs"_*.out`
set -- $List
tmp=$1
if [ -f $tmp ]
then
  tmp=`ls -1rt $tmp | tail -n 1`
  c=`tail -n 1 $tmp | wc -w`
  #let c=c-$L-1 # won't work on Ubuntu?
  c=`expr $c - $L - 1`
  tail -n 1 $tmp | cut -d " " -f $c
  #let l=L-1
  #awk '{print $(NF-$l)}' $tmp | tail -n 1
fi

#!/bin/bash
readonly program=$(basename $0)
function print_usage_and_exit() {
  echo >&2 "Usage: ./${program} CASE data/CASE"
  exit 1
}
if [ $# -ne 2 ]; then
  print_usage_and_exit
fi
#----------------------------------------------------------------

egrep ":  |2     ="  $1 | awk '{printf("%30s %30s\n", $1, $NF)}'> res1
egrep ":  |2     ="  $2 | awk '{print $NF}'> res2
paste res1 res2 |awk '{d=sqrt(($2-$3)^2);printf("%27s %22s %22s %e\n", $1, $2, $3, d)}' >res0


awk '{e=1;
    if($1=="rnorm^2"){
        if($2<1e-15){e=0}
    }else{
        if(NR>=26 && NR<=39){
                 if($4/$2<1e-6){e=0}
        }else{
                 if($4/$2<1e-15){e=0}
        }
    }
printf("%27s %22s %22s %9s %d\n",$1,$2,$3,$4,e) }' res0 > res
cat res
err_count=`cat res | awk '{a = a+ $NF}END{ print a}'`

rm res res0 res1 res2
exit $err_count

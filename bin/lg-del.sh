fn1='lm0.5705'
fn2='k3L52'
fn3='648'
fn4='684'

rm ../src/outputs/${fn1}${fn2}n${fn3}-${fn4}_output
for i in `seq 648 683`
do
  rm ../src/confs/nicolai_${fn1}*${fn2}n${i}.dat
  rm ../src/confs/phi_${fn1}*${fn2}n${i}*.dat
done

for i in `seq 71 107`
do
  ls -l ../src/confs/nicolai*n${i}.dat | grep "*アクセス*"
done

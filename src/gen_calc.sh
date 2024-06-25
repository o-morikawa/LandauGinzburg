#!/bin/bash
li=52
kval="3"
lval="0.5705"
numnic=36
#for i in `seq 0 71`
for i in `seq 26 35`
do
  num=`expr $i \* $numnic`
  num2=`expr $num + $numnic`
  j=`expr $i + 1`

###calc_str
  fn1="\"lm${lval}k${kval}L${li}n${num}-${num2}\""
  echo "#!/bin/bash" > calc_str${j}
  echo "export OMP_NUM_THREADS=18" >> calc_str${j}
  echo "fn1=${fn1}" >> calc_str${j}
  echo "" >> calc_str${j}
  echo "for fn in \$fn1" >> calc_str${j}
  echo "do" >> calc_str${j}
  echo "STARTTIME=\`date --iso-8601=seconds\`" >> calc_str${j}
  echo "(time ./\${fn}.out) &> tmp_out/\${fn}_output" >> calc_str${j}
  echo "ENDTIME=\`date --iso-8601=seconds\`" >> calc_str${j}
  echo "mv tmp_out/\${fn}_output outputs/\${fn}_output" >> calc_str${j}
  echo "" >> calc_str${j}
  echo "echo \"From: <morikawa@neumann.higgs.kyushu-u.ac.jp>\" > mail.txt" >> calc_str${j}
  echo "echo \"To: <o-morikawa@phys.kyushu-u.ac.jp>\" >> mail.txt" >> calc_str${j}
  echo "echo \"Subject: [LG-cutoff] Finished \${fn}\" >> mail.txt" >> calc_str${j}
  echo "echo \"\" >> mail.txt" >> calc_str${j}
  echo "echo \"O. M.\" >> mail.txt" >> calc_str${j}
  echo "echo \"\" >> mail.txt" >> calc_str${j}
  echo "echo \"\${fn}.out\" >> mail.txt" >> calc_str${j}
  echo "echo \"Status     : Simulation finished\" >> mail.txt" >> calc_str${j}
  echo "echo \"Start time : \$STARTTIME\" >> mail.txt" >> calc_str${j}
  echo "echo \"End time   : \$ENDTIME\" >> mail.txt" >> calc_str${j}
  echo "echo \"\" >> mail.txt" >> calc_str${j}
  echo "echo \"##\" >> mail.txt" >> calc_str${j}
  echo "cat outputs/\${fn}_output >> mail.txt" >> calc_str${j}
  echo "echo \"\" >> mail.txt" >> calc_str${j}
  echo "echo \"##\" >> mail.txt" >> calc_str${j}
  echo "echo \"\" >> mail.txt" >> calc_str${j}
  echo "echo \"Please do not reply to this message.\" >> mail.txt" >> calc_str${j}
  echo "cat mail.txt | sendmail -i -t" >> calc_str${j}
  echo "rm mail.txt" >> calc_str${j}
  echo "done" >> calc_str${j}
  chmod +x calc_str${j}

###lg.sh
#  echo "#!/bin/bash" > lg${j}.sh
#  echo "#PJM -L \"rscunit=ito-a\"" >> lg${j}.sh
#  echo "#PJM -L \"rscgrp=ito-ss\"" >> lg${j}.sh
#  echo "#PJM -L \"vnode=1\"" >> lg${j}.sh
#  echo "#PJM -L \"vnode-core=36\"" >> lg${j}.sh
#  echo "#PJM -L \"elapse=96:00:00\"" >> lg${j}.sh
#  echo "#PJM -j" >> lg${j}.sh
#  echo "#PJM -X" >> lg${j}.sh
#  echo "" >> lg${j}.sh
#  echo "export OMP_NUM_THREADS=18" >> lg${j}.sh
#  echo "./calc_str${j}" >> lg${j}.sh

#  pjsub lg${j}.sh
  ./calc_str${j}
done

fi="spt"

#pjstat -H day=1 -v
#ls -l ../src/outputs/${fi}* -tr
ls -l ../src/tmp_out/${fi}* -tr
#cat ../src/outputs/${fi}*output
cat ../src/tmp_out/${fi}*output
ls -l ../src/confs/nicolai_${fi}* -tr | tail
ls -l ../src/confs/nicolai_${fi}* | wc -l

#!/bin/bash
#PJM -L "rscunit=ito-a"        
#PJM -L "rscgrp=ito-s"
#PJM -L "vnode=2"
#PJM -L "vnode-core=36"
#PJM -L "elapse=:00:00"
#PJM -j
#PJM -X

export OMP_NUM_THREADS=3
./calc_str

# -- Add dNami python related files to path
export  INSTALLPATH=$PWD
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/src_py
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/pymod
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/generate/

# - post processing functions
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/../pst/utils

# -- 
#export  KMP_AFFINITY='verbose,granularity=fine,compact,1,0'
export  KMP_AFFINITY='verbose,granularity=fine,proclist=[10,11,12,13,14,15,16,17,18,19,20,0,1,2,3,4,5,6,7,8,9],explicit'

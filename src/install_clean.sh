rm .*.dblite
scons -c
export  INSTALLPATH=$PWD
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/src_py
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/pymod
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/generate/
export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/generate/archive/forConsLaw/	
python3 generate/genRhs.py
scons -j 4

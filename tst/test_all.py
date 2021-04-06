#  ===================================================
#
#                  Run test and compare to 
#                   reference results
#
# Details of each run are written out in 'test.log'
# 
# The following things are assumed:
# 
#    o each test case comes with a 'copy.sh' file 
#    o each compute specifies the number of precessors
#      as 'with_proc = [X,X,X]' 
#    o each test case prints an out.dat file and each
#      test case contains a reference.dat file 
#
#                             - 12/03/21  Stephen Winn 
# 
#  TO DO:
#    o add code generation tests
#    o refine pass/fail conditions for some tests
#    o check process.stdout for errors properly
#  ===================================================
# ---------------------------------------------------
import glob
import os 
import subprocess
import numpy as np

# ---------------------------------------------------

# -- Color class 
class bcolors:
    FAIL    = '\033[91m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ENDC    = '\033[0m'

# -- Get list of tests
test_dir  = os.getcwd()
test_log  = test_dir + '/test.log'
test_list = glob.glob('./*')
test_list.remove('./test_all.py')
try:
    test_list.remove('./test.log')
except:
    pass
test_stat = {} 

# -- Print the list
print(' ============================================' )
print(' List of tests:')
for test in test_list:
    print('    o ', test)
print(' ============================================' )
print(' ')
print(' ')
print(' ')

# -- Before first test, create wrk and remove files from previous test
try:
    os.mkdir('../wrk/')
    os.remove(test_log)
except:
    pass

# -- Go through tests
for test in test_list:

    print(f'{bcolors.WARNING}Running test ... {test}{bcolors.ENDC}')

    # -- Run copy script
    cmd = 'cd ' + test+ '; ./copy.sh >> ' + test_log 
    try:
        #process = subprocess.run(cmd,shell=True,capture_output=True,text=True,check=True)
        process = subprocess.run(cmd,shell=True)
    except Exception as e:
        print('Error at copy: ', e)
        test_stat[test] = f'{bcolors.FAIL}COPY FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Copied.')

    # -- Compile 
    cmd = 'cd ../src; ./install_clean.sh >> ' + test_log
    try:
        #process = subprocess.run(cmd,shell=True,capture_output=True,text=True,check=True)
        process = subprocess.run(cmd,shell=True)
    except Exception as e:
        print('Error at compilation: ', e)
        test_stat[test] = f'{bcolors.FAIL}COMPILATION FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Compiled.')

    # -- Run 
    # Get mpi proc number:
    lines = open('../wrk/compute.py').readlines()
    for line in lines:
        if 'with_proc' in line:
            break
    nproc = line.split('[')[1].split(']')[0]
    if ',' in nproc:
        nproc = nproc.split(',')
        nnproc = 1
        for dproc in nproc:
            nnproc *= int(dproc)
    else:
        nnproc = int(nproc)
    if nnproc == 1:
        print(' Single processor used.')
    else:
        print(f' {nnproc} processors used.')
    nnproc = str(nnproc)
    # Command 
    cmd = 'cd ../src; \
        export  INSTALLPATH=$PWD; \
        export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/src_py; \
        export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/pymod; \
        export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/generate/; \
        cd ../wrk; mpirun -n ' + nnproc + ' python3 compute.py >> ' + test_log
    try:
        #process = subprocess.run(cmd,shell=True,capture_output=True,text=True,check=True)
        process = subprocess.run(cmd,shell=True)
    except subprocess.CalledProcessError as e:
        print('Error at runtime: ', e)
        test_stat[test] = f'{bcolors.FAIL}RUN FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Ran.')

    # -- Compare and declare valid or not
    try:
        ref_dat = np.loadtxt(test+'/reference.dat')
    except Exception as e:
        print('Error:...', e)
    try:
        out_dat = np.loadtxt('../wrk/out.dat') 
    except Exception as e:
        print('Error:...', e)

    print(' Reference value: {}. Output data: {}'.format(ref_dat, out_dat))

    if np.abs(ref_dat-out_dat)/ref_dat < 1e-10:
        print(f'STATUS: {bcolors.OKGREEN}PASS{bcolors.ENDC}' )
        test_stat[test] = f'{bcolors.OKGREEN}PASS{bcolors.ENDC}'
    else:
        print(f'STATUS: {bcolors.FAIL}FAIL{bcolors.ENDC}' )
        test_stat[test] = f'{bcolors.FAIL}FAIL{bcolors.ENDC}'

    # -- Clean wrk for next test
    cmd = 'rm -r ../wrk/*'
    #subprocess.run(cmd,shell=True,capture_output=True,text=True,check=True)
    subprocess.run(cmd,shell=True)

    print('')

print('===================================================================')
print('All done. Synopsis:')

for key in test_stat.keys():
    print(' ' + key + ': ' + test_stat[key])


#  ===================================================
#
#                  Run test and compare to 
#                   reference results
#
# Details of each run are written out in 'test.log'
# 
# The following things are assumed:
# 
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
import shutil

# -- Color class 
class bcolors:
    FAIL    = '\033[91m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ENDC    = '\033[0m'

# --------------------------------------------------
# -- Check if the dependencies are met:

deps={ 'numpy'    : None,
       'scons'    : None,
       'mpi'      : None,
       'fcompiler': None, 
       'mpi4py'   : None
       }

# -- Numpy
try:
    import numpy as np
    deps['numpy'] = f'{bcolors.OKGREEN}OK{bcolors.ENDC}'
except:
    deps['numpy'] = f'{bcolors.FAIL}MISSING{bcolors.ENDC}'

# -- scons 
if shutil.which('scons'):
    deps['scons'] = f'{bcolors.OKGREEN}OK{bcolors.ENDC}'
else:
    deps['scons'] = f'{bcolors.FAIL}MISSING{bcolors.ENDC}'

# -- mpirun or mpiexec 
if shutil.which('mpirun') or shutil.which('mpiexec'):
    deps['mpi'] = f'{bcolors.OKGREEN}OK{bcolors.ENDC}'
else:
    deps['mpi'] = f'{bcolors.FAIL}MISSING{bcolors.ENDC}'

# -- mpirun or mpiexec 
if shutil.which('gfortran') or shutil.which('ifort'):
    deps['fcompiler'] = f'{bcolors.OKGREEN}OK{bcolors.ENDC}'
else:
    deps['fcompiler'] = f'{bcolors.FAIL}MISSING{bcolors.ENDC}'

# -- MPI4PY
try:
    import mpi4py 
    deps['mpi4py'] = f'{bcolors.OKGREEN}OK{bcolors.ENDC}'
except:
    deps['mpi4py'] = f'{bcolors.FAIL}MISSING{bcolors.ENDC}'


print(' --- Dependencies --- ')
for key in deps.keys():
    print('Item: {:10s}'.format(key) +  '  Status: {:10s}'.format(deps[key]))

print('')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print( f'{bcolors.WARNING}If all dependencies are not OK, the next steps will likely fail .... {bcolors.ENDC}') 
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('')
# ---------------------------------------------------

# -- Get list of tests
test_dir  = os.getcwd()
test_log  = test_dir + '/test.log'
test_list = sorted(glob.glob('./*'))
test_list.remove('./test_all.py')
try:
    test_list.remove('./test.log')
except:
    pass
test_stat = {} 

# -- Check if all tests have a reference.dat file:
for test in test_list:

    # Get folder contents:
    files = glob.glob(test + '/*')
    ref = [i for i in files if 'reference.dat' in i]
    if len(ref) == 0:
        test_list.remove(test)



# -- Print the list
print(' ============================================' )
print(' List of tests with references:')
for test in test_list:
    print('    o ', test)
print(' ============================================' )
print(' ')
print(' ')
print(' ')

# -- Before first test, create temporary work directory and remove files from previous test
tmp_dir = 'wrk_tmp'
try:
    os.mkdir('../'+tmp_dir+'/')
    os.remove(test_log)
except:
    pass

# -- Go through tests
for test in test_list:

    print(f'{bcolors.WARNING}Running test ... {test}{bcolors.ENDC}')

    # -- Run copy step 
    try:
        shutil.copyfile( test+'/rhs.py', '../src/generate/rhs.py' )
        shutil.copyfile( test+'/genRhs.py', '../src/generate/genRhs.py' )
    except Exception as e:
        print('Error at copy: ', e)
        test_stat[test] = f'{bcolors.FAIL}COPY FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Copied.')

    # -- Compile 
    cmd = './install_clean.sh >> ' + test_log
    try:
        process = subprocess.run(cmd, cwd='../src',shell=True,capture_output=True)
    except Exception as e:
        print('Error at compilation: ', e)
        test_stat[test] = f'{bcolors.FAIL}COMPILATION FAIL{bcolors.ENDC}'
        exit()
        continue
    else:
        print(' Compiled.')

    # -- Run 
    # Copy compute:
    try:
        shutil.copyfile( test+'/compute.py', '../'+tmp_dir+'/compute.py' )
    except Exception as e:
        print('Error at copy: ', e)
        test_stat[test] = f'{bcolors.FAIL}COPY FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Copied.')

    # Get mpi proc number:
    lines = open('../'+tmp_dir+'/compute.py').readlines()
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
    cwd = '/'.join(os.getcwd().split('/')[:-1]) + '/src/'
    cwd = cwd.replace(' ', '\ ')
    cmd = 'export  INSTALLPATH='+cwd+'; \
           export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/src_py; \
           export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/pymod;  \
           export  PYTHONPATH=$PYTHONPATH:$INSTALLPATH/generate; \
           mpirun -n ' + nnproc + ' python3 compute.py >> ' + test_log
    try:
        process = subprocess.run(cmd,shell=True,cwd='../'+tmp_dir)
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
        test_stat[test] = f'{bcolors.FAIL}COMPARE FAIL{bcolors.ENDC}'
        print('Error:...', e)
        continue

    try:
        out_dat = np.loadtxt('../'+tmp_dir+'/out.dat') 
    except Exception as e:
        test_stat[test] = f'{bcolors.FAIL}COMPARE FAIL{bcolors.ENDC}'
        print('Error:...', e)
        continue

    print(' Reference value: {}. Output data: {}'.format(ref_dat, out_dat))

    if ref_dat != 0.:
        delta = np.abs(ref_dat-out_dat)/ref_dat
    else:
        delta = np.abs(ref_dat-out_dat)

    if  delta < 1e-8:
        print(f'STATUS: {bcolors.OKGREEN}PASS{bcolors.ENDC}' )
        test_stat[test] = f'{bcolors.OKGREEN}PASS{bcolors.ENDC}'
    else:
        print(f'STATUS: {bcolors.FAIL}FAIL{bcolors.ENDC}' )
        test_stat[test] = f'{bcolors.FAIL}FAIL{bcolors.ENDC}'

    # -- Clean wrk for next test
    cmd = 'rm -r ../'+tmp_dir+'/*'
    subprocess.run(cmd,shell=True)

    print('')

print('===================================================================')

print('Cleaning up ...')

cmd = 'rm -r ../'+tmp_dir
subprocess.run(cmd,shell=True)

print('All done. Synopsis:')

for key in test_stat.keys():
    print('{:25} '.format(key)  +  ': ' + test_stat[key])


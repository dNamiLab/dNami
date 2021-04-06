#  ===================================================
#
#               Check if code compiles 
#
# Details of each run are written out in 'test.log'
# 
# The following things are assumed:
# 
#    o each test case comes with a 'copy.sh' file 
#
#                             - 15/03/21  Stephen Winn 
# 
#  ===================================================
# ---------------------------------------------------
import glob
import os 
import subprocess

# ---------------------------------------------------

# -- Color class 
class bcolors:
    FAIL    = '\033[91m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ENDC    = '\033[0m'

# -- Get list of exms
exm_dir  = os.getcwd()
exm_log  = exm_dir + '/exm.log'
exm_stat = {} 
exm_list = [
        '2d_vortex_advection',
        '3d_tgv',
        '3d_turbulent_counter_flow/uniform_grid'
        ]

# -- Print the list
print(' ============================================' )
print(' List of exms:')
for exm in exm_list:
    print('    o ', exm)
print(' ============================================' )
print(' ')
print(' ')
print(' ')

# -- Before first exm, create wrk and remove files from previous exm
try:
    os.mkdir('../wrk/')
    os.remove(exm_log)
except:
    pass

# -- Go through exms
for exm in exm_list:

    print(f'{bcolors.WARNING}Running exm ... {exm}{bcolors.ENDC}')

    # -- Run copy script
    cmd = 'cd ' + exm+ '; ./copy.sh >> ' + exm_log 
    try:
        #process = subprocess.run(cmd,shell=True,capture_output=True,text=True,check=True)
        process = subprocess.run(cmd,shell=True)
    except Exception as e:
        print('Error at copy: ', e)
        exm_stat[exm] = f'{bcolors.FAIL}COPY FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Copied.')

    # -- Compile 
    cmd = 'cd ../src; ./install_clean.sh >> ' + exm_log
    try:
        #process = subprocess.run(cmd,shell=True,capture_output=True,text=True,check=True)
        process = subprocess.run(cmd,shell=True)
    except Exception as e:
        print('Error at compilation: ', e)
        exm_stat[exm] = f'{bcolors.FAIL}COMPILATION FAIL{bcolors.ENDC}'
        continue
    else:
        print(' Compiled.')

    # -- Success if no errors
    exm_stat[exm] = f'{bcolors.OKGREEN}PASS{bcolors.ENDC}'



    print('')

print('All done. Synopsis:')

for key in exm_stat.keys():
    print(' ' + key + ': ' + exm_stat[key])


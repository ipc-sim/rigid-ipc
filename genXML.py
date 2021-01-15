import subprocess
import sys

if len(sys.argv) <= 1:
    print('need xml file name')
    exit()

mode = '0' # 0 for bullet, 1 for mujoco
if len(sys.argv) > 2:
    mode = sys.argv[2]

runCommand = 'build/FixingCollisions_sim 1 fixtures/3D/' + sys.argv[1] + '.json ' + mode
subprocess.call([runCommand], shell=True)

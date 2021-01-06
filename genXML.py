import subprocess
import sys

if len(sys.argv) <= 1:
    print('need xml file name')
    exit()

runCommand = 'build/FixingCollisions_sim 1 fixtures/3D/' + sys.argv[1] + '.json'
subprocess.call([runCommand], shell=True)

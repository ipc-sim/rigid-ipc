import subprocess
import sys

if len(sys.argv) <= 1:
    print('need xml file name')
    exit()

bulletOutputFilePath = '../bullet3/test.txt'
if len(sys.argv) > 2:
    bulletOutputFilePath = '../bullet3/' + sys.argv[2]

runCommand = 'build/FixingCollisions_sim 2 fixtures/3D/' + sys.argv[1] + '.json ' + bulletOutputFilePath
subprocess.call([runCommand], shell=True)

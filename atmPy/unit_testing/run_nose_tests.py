#!/opt/local/bin/python3

import subprocess
import os

def run():
    # folder = os.path.split(__file__)[0]
    folder = os.path.split(os.path.realpath(__file__))[0]
    os.chdir(folder)

    try:
        out = subprocess.check_output('nosetests-3.5 nose_tests.py', shell = True, stderr=subprocess.STDOUT)
        out = out.decode()
    except subprocess.CalledProcessError as e:
        out = e.stdout.decode()

    raus = open('testresults.log', 'w')
    raus.write(out)
    raus.close()

    out_lines = out.splitlines()

    title = 'atm-py unit test (%s)'%(out_lines[-1])

    print(out_lines)
    if out_lines[-1] == 'OK':
        message = out_lines[2]
    else:
        message = out_lines[-6]

    openthis = folder + '/testresults.log'
    os.system("terminal-notifier -title '%s' -message '%s' -execute 'open -a TextWrangler %s'"%(title,message,openthis))

if __name__ == "__main__":
    run()
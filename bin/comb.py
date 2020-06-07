import os

files = os.listdir('.')
files = [f for f in files
         if 'red' in f]
files.sort()

with open('full.ckin', 'a') as full:
    full.write('REACTIONS\n')
    for f in files:
        print(f)
        with open(f, 'r') as fobj:
            fstr = fobj.read()
            full.write(fstr)
    full.write('\nEND')

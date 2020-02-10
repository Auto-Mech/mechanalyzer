
with open('arrfit.out', 'r') as f:
    lines = f.readlines()

lines.reverse()
for line in lines:
    if line.startswith(' params'):
        bits1 = lines[lines.index(line)-1].split()
        bits2 = lines[lines.index(line)-2].split()
        break

print(bits1)
print(bits2)

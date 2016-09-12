import sys

ref = {}
f = open(sys.argv[1])
g = open(sys.argv[2])

key = ""
for line in g:
    if line[0] == '>':
        key = line[1:].strip()
        ref[key] = ""
    else:
        ref[key] += line.strip()
print ref['>CHROMOSOME_I']
g.close()

total = 0.0
count = 0
for line in f:
    if line[0] == '@':
        continue

    count += 1

f.close()

import sys


f = open(sys.argv[1])
string = ""
print f.readline().strip()


for line in f:
  if line[0] == '>':
    print len(string)
    string = ""
    continue
  else:
    string += line.strip()

if len(string) > 0:
  print len(string)

f.close()

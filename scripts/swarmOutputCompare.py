#!/usr/bin/python
import os, sys

def removeEmpty(s):
  while(s.count('') > 0):
    s.remove('')

def removeLine(arr, line):
  while (arr.count(line) > 0):
    arr.remove(line)

if (len(sys.argv) < 3):
  print "Usage: swarmOutputCompare.py output_file reference_file"
  print "Returns 0 if files are same other than floating point differeces."
  print "Returns 1 if files differ."
  sys.exit(0)

outfile = sys.argv[1]
reffile = sys.argv[2]

f = open(outfile, 'rb')
out = f.read().split('\n')
f.close() 
removeEmpty(out)
#In monte carlo test, "no error" appeared 3 times in output and 1 time in reference.
removeLine(out, "no error")

f = open(reffile, 'rb')
ref = f.read().split('\n')
f.close()
removeEmpty(ref)
removeLine(ref, "no error")

#If different # of lines, files differ.
if (len(out) != len(ref)):
  print "Files differ.  Output has "+str(len(out))+" lines while reference has "+str(len(ref))+" lines."
  sys.exit(1)

#Do line by line comparison
for j in range(len(out)):
  if (out[j] == ref[j]):
    #Lines match exactly, continue
    continue
  #If not exact match, split up into tokens
  outsplit = out[j].split()
  refsplit = ref[j].split()
  if (len(outsplit) != len(refsplit)):
    #If different number of tokens, lines differ.
    print "Files differ on line "+str(j+1)
    sys.exit(1)
  #Loop over and compare tokens
  for l in range(len(outsplit)):
    try:
      #If this fails, token is a string, handle, and compare directly
      outval = float(outsplit[l])
      refval = float(refsplit[l])
      #Compare float vals
      if (abs(outval-refval) > 1.e-5):
	print "Files differ on line "+str(j+1)
	sys.exit(1)
    except ValueError:
      #Tokens are strings, should match exactly
      if (outsplit[l] != refsplit[l]):
	print "Files differ on line "+str(j+1)
	sys.exit(1)
print "Files match."
sys.exit(0)

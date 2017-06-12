#!/usr/bin/python

import sys;


# --  take input of fileName
nArgs = len(sys.argv)
if ( nArgs != 2 ) :
  print "number of input arguments is not 2, nArguments: ", nArgs
  print "Format: python ",sys.argv[0]," <fileName> "
  quit()
print nArgs
inFileName = sys.argv[1];
# -- DEBUG
print "The file input name is : ",inFileName
# -- END Debug


# -- set outFile Names
dotIdx = inFileName.rindex('.');
filewoE = inFileName[:dotIdx];
print "FileName Basis: ", filewoE
outVertices  = filewoE + "Vertices.txt"
outEdges     = filewoE + "Edges.txt"
outLandmarks = filewoE + "Landmarks.txt"
print outVertices
print outEdges
print outLandmarks

# -- read input file and assign to output files
inFileStream = open(inFileName,'r')
outVFileStream = open(outVertices,'w')
outEFileStream = open(outEdges,'w')
outLFileStream = open(outLandmarks,'w')

for curLine in inFileStream:
  #curLine = filn.readline()
  # -- check length of line
  l = len(curLine)
  if (l < 10) :
    continue
  brkIdx = curLine.index(' ')
  lineDesp = curLine[:brkIdx]
  dataLine = curLine[brkIdx+1:]
  if (lineDesp == "VERTEX_SE2") :
    outVFileStream.write(dataLine)
  elif (lineDesp == "VERTEX_XY") :
    outLFileStream.write(dataLine)
  elif (lineDesp == "EDGE_SE2") :
    outEFileStream.write(dataLine)
  else :
    continue

# -- close the file descriptors
inFileStream.close()
outVFileStream.close()
outEFileStream.close()
outLFileStream.close()
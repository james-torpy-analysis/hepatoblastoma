#!/bin/bash

home_dir='/share/ScratchGeneral/jamtor'
in_dir='/share/ScratchGeneral/jamtor/projects/hepatoblastoma/andre/trimmed'

import glob
import pyfastx
import gzip
from itertools import islice
import regex
import os
import errno

read2_ptn='ATTGGAGTCCT'
read1_ptn='AGGACTCCAAT' # if insert is small this barcode can be sequenced into: AGGACTCCAATNNNNNNNNNNNN
read1_ptn_=regex.compile("("+read1_ptn+"){s<=1}")   # s<=1 allow one mismatch, e<=1 also allow indel
read2_ptn_len=len(read2_ptn)

inArr=glob.glob(in_dir + '/*_R1*.fastq.gz', recursive=False); force_ptn2=True; outFolder='no_barcode'
inArr=inArr[0:1]

def fastq(fname1, fname2):  
  fh1=gzip.open(fname1, 'rb')
  fh2=gzip.open(fname2, 'rb')
  while True:
    try:
      _name1, _seq1, _x1, _qual1 =[x.strip().decode('utf-8') for x in islice(fh1, 4)]
    except:
      _name1, _seq1, _x1, _qual1 = [None,None,None,None]
    try:
      _name2, _seq2, _x2, _qual2 =[x.strip().decode('utf-8') for x in islice(fh2, 4)]
    except:
      _name2, _seq2, _x2, _qual2 = [None,None,None,None]
    if _qual1 is None and _qual2 is None:
      break #end of both files
    elif _qual1 is None or _qual2 is None:
      raise ValueError(j,"Diff record count",_name1, _name2)
    _id1, _desc1=_name1[1:].split(' ')
    _id2, _desc2=_name2[1:].split(' ')
    if _id1 != _id2:
      raise ValueError(j, "id not matching", _id1," != ", id2)
    yield _id1, _desc1, _seq1, _qual1, _id2, _desc2, _seq2, _qual2

def get_match_count(_seq2):
  match_c=0
  for k in range(12,min(22,len(_seq2)) ):
    #print(' {}:{} {} {}'.format(k, k-12, _seq2[k], read2_ptn[k-12] ))
    if _seq2[k]==read2_ptn[k-12]: match_c+=1
  return match_c

for i, in1 in enumerate(inArr):
  print(' {}/{} {} '.format(i+1, len(inArr), in1 )) 
  in2=in1.replace('_R1','_R2')
  out1n=in1.replace('trimmed','trimmed/'+outFolder)
  out2n=in2.replace('trimmed','trimmed/'+outFolder)
  # make out directories:
  try:
    os.makedirs(os.path.dirname(out1n))
  except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass
  try:
    os.makedirs(os.path.dirname(out2n))
  except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass
  out1 = gzip.open(out1n, "wb")
  out2 = gzip.open(out2n, "wb")
  _skipped=0
  _all=0
  for _id1, _desc1, _seq1, _qual1, _id2, _desc2, _seq2, _qual2 in fastq(in1,in2):
    # search for read 2 pattern
    _all+=1
    _pos=12 if get_match_count(_seq2)>=5 else -1
    _tag=''
    if force_ptn2 and _pos<0: 
      _skipped+=1
      continue
      # remove ID and add into header
    if _pos>=0:
      _tag=_seq2[:_pos]
      _seq2=_seq2[_pos+read2_ptn_len:]
      _qual2=_qual2[_pos+read2_ptn_len:]
      if len(_seq2)<30 or len(_seq1)<30:
        _skipped+=1
        continue
    # search for read1 pattern
    r=read1_ptn_.search( _seq1)
    if r:
      _pos=r.start()
      _seq1=_seq1[:_pos]
      _qual1=_qual1[:_pos]
    else:
      # try to trim from end
      for j in reversed(range(1,len(read1_ptn)+1)):
        #print(' var:{} {} {}'.format( j, _seq1[-j:], read1_ptn[:j]))
        if _seq1[-j:]==read1_ptn[:j]:
          _seq1=_seq1[:-j]
          _qual1=_qual1[:-j]
          break
    if len(_seq1)<30: 
      _skipped+=1
      continue
    r=out1.write( '@{} {}\n{}\n+\n{}\n'.format( _id1+'_'+_tag,_desc1,_seq1,_qual1).encode('utf-8')  )
    r=out2.write( '@{} {}\n{}\n+\n{}\n'.format( _id2+'_'+_tag,_desc2,_seq2,_qual2).encode('utf-8')  )
    #print('@{} {}\n{}\n+\n{}\n'.format( _id1, _desc1, _seq1, _qual1 ))
    #print('@{} {}\n{}\n+\n{}\n'.format( _id2, _desc2, _seq2, _qual2 ))
    #break
  out1.close()
  out2.close()
  # determine percentageof reads skipped:
  if _skipped != 0:
    prop_skipped = _skipped/_all*100
  else:
    prop_skipped = 0
  print(' skipped:{} / {} ( {:.3f} %)'.format( _skipped, _all,  prop_skipped))

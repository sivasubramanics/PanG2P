import sys


### Give your input file path ###

file=open(sys.argv[1],'r')
fh=file.read()
file.close()




### Give your output file name ###

file1=open(sys.argv[2],'w')

#####################################################


fh=fh.strip()
id=fh.split('\n')
l=len(id)


x=0
i=0
j=1
p='p'

c=1

while x<(l-1):
 a=id[i]
 b=id[j]

 id1=a.split('\t')
 id2=b.split('\t')
 ac1=id1[3]
 ac2=id2[3]

 if (ac1==ac2) :
  p=id1[0]+"\t"+id1[1]+"\t"+id1[2]+"\t"+id1[3]+"_"+str(c)
  file1.write(p)
  file1.write('\n')
  c=c+1
  
 elif (ac1!=ac2):
  p=id1[0]+"\t"+id1[1]+"\t"+id1[2]+"\t"+id1[3]+"_"+str(c)
  file1.write(p)
  file1.write('\n')
  c=1


 i=i+1
 j=j+1
 x=x+1
 if j>(l-1):
  #print p
  m=id[(l-1)]
  id1=m.split('\t')
  p2=id1[0]+"\t"+id1[1]+"\t"+id1[2]+"\t"+id1[3]+"_L"
  file1.write(p2)
  file1.write('\n')

#print ('### your file is ready now...')

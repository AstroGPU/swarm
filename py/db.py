from bsddb.db import *
from struct import *

d = DB()
d.open("log.db")
c = d.cursor()

def process(s):
	(msgid, length) = unpack('ii',s[0:8]) # msgid, len
	if msgid == 1 :
		(time,sys,flags,nbod) = unpack('diii',s[8:28])
		bodies = {}
		for b in range(1,nbod+1):
			(x,y,z,vx,vy,vz,mass,body_id) = \
					unpack('dddddddi4x',s[(32+(b-1)*64):(32+(b)*64)])
			bodies[body_id] = ( [x,y,z] , [vx,vy,vz] , mass )
		return { 'time':time,'system':sys,'state':flags,'bodies':bodies }
	else:
		return None


for i in range(1,20):
	r = c.next()
	if(r == None):
		break;
	else:
		print process(r[1])




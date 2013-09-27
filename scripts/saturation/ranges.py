import math
initial=32
q=math.sqrt(2.0)
p=initial
x= []
for i in range(1,22):
	x.append(p)
	p = q * p
print(','.join(map(str,x)))

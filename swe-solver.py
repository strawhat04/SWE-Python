import numpy as np
from math import cos, pi, sin 
def F1(t,i):
	return q[t,i]

def F2(t,i):
	return h[t,i]*u[t,i]**2 + 0.5*g*h[t,i]**2

def BC(t):
	h[t,1]=htransient(t,0)
	h[t,0]=h[t,1]
	h[t,-1]=h[0,-2]

	u[t,0]=u[t,1]
	u[t,-1]=0*u[t,-1]

def profile():
	for i in range(1,len(h[0,:-1])):
		#if xgrid[i]>=2/5 and xgrid[i]<=3/5:
			#b[i]=( cos(10*pi*(xgrid[i]-0.5))+1 )/4
		b[i]=10.5+40*xgrid[i-1]/l - 10*sin(pi*(4*xgrid[i-1]/l + 0.5))

def h_init():
	for i in range(1,len(h[0,:-1])):
		h[0,i]=50.5-40*xgrid[i-1]/l + 10*sin(pi*(4*xgrid[i-1]/l + 0.5))
		# if xgrid[i]<=0.1:
		# 	h[0,i]=1-b[i]

		# if xgrid[i]>=0.1 and xgrid[i]<=0.2:
		# 	h[0,i]=1.2-b[i]

		# else:
		# 	h[0,i]=1-b[i]
			
def htransient(t,i):
	if t<=43200:
		return 64.5 + 4*sin(pi*(4*t/86400 - 0.5))

	else:
		return 60.5

#mesh code
l=600000
nx=501
dt=1
runtime=10800

xgrid=np.linspace(0,l,nx,dtype=float)

dx=l/(nx-1)

print("xgrid",xgrid)


g=9.8

lam=dt/dx



iniGRIDFunc=lambda m: [np.zeros(((int(runtime/dt))+1,xgrid.shape[0]+2)) for _ in range(m)]
h,q,u=iniGRIDFunc(3)
b=np.zeros((xgrid.shape[0]+2))
#input conditions
# h[0,:int(nx/2)]=1
# h[0,int(nx/2):]=0.5
#u[0,:]=0.2

profile()
h_init()


q[0,:]=np.multiply(h[0,:],u[0,:])

BC(0)
for t in range(1, len(h[:,0])):

	
	for i in range(1,len(h[0,:-1])):
		# if t==21:
		# 	print(t,h[t,i])

		h[t,i]=(h[t-1,i+1] + h[t-1,i-1])/2 - lam*(F1(t-1,i+1) - F1(t-1,i-1))/2 

		q[t,i]=(q[t-1,i+1] + q[t-1,i-1])/2 - lam*(F2(t-1,i+1) - F2(t-1,i-1))/2 - lam*g*h[t-1,i]*(b[i]-b[i-1]) 

	

	BC(t)
	for i in range(1,len(h[0,:-1])):
		if h[t,i]==0:
			u[t,i]=0
		else:
			u[t,i]=q[t,i]/h[t,i]

print(h)
print(u)

np.save('h',h)
np.save('u',u)
np.save('xgrid',xgrid)
np.save('b',b)
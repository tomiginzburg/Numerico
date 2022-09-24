import matplotlib.pyplot as plt
import numpy as np

def DifDiv(X,f):
	n = len(X)-1
	F=[f]
	for i in range(n):
		aux=[]
		for j in range(n-i,0,-1):
			aux.insert(0,(F[i][j]-F[i][j-1])/(X[j+i]-X[j-1]))
		F.append(aux)
	a=[F[k][0] for k in range(len(F))]
	
	f=0
	for i in range(len(a)):
		b=1
		for j in range(i):
			b*=(x-X[j])
		f+= a[i]*b

	return (a,f)

def Hermit(X,f,g):
	n = len(X)-1
	Z=[]
	F1=[]
	for h in range(n+1):
		Z.extend([X[h],X[h]])
		F1.extend([f[h],f[h]])
		F2=[]
	for w in range(n):
		F2.extend([g[w],(F1[2*w+2]-F1[2*w])/(Z[2*w+2]-Z[2*w])])
	F2.append(g[n])
	F=[F1,F2]
	for i in range(1,2*n+1):
		aux=[]
		for j in range(2*n-i+1,0,-1):
			aux.insert(0,(F[i][j]-F[i][j-1])/(Z[j+i]-Z[j-1]))
		F.append(aux)
	a=[F[k][0] for k in range(len(F))]
	f=0
	for i in range(len(a)):
		b=1
		for j in range(i):
			b*=(x-Z[j])
		f+= a[i]*b
	return (a,f)

X = [-1,0,1,2,3]
Y = [-2,4,3,0,-1.8]
Z = [10,2,-4,-2,0]

x = np.linspace(X[0],X[-1],100)

a,f = Hermit(X,Y,Z)
b,g = DifDiv(X,Y)

plt.plot(X,Y,'x')
plt.plot(x,f,'r')
plt.plot(x,g,'g')
plt.show()


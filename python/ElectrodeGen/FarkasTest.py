from ElectrodeGenerator import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

p1 = GetRandomPlatelet(1,[40,110,110], [4,2,0.5])
p2 = GetRandomPlatelet(2,[40,110,110], [4,2,0.5])

print(np.linalg.norm(p1.GetPos()-p2.GetPos()))

# Calculate the Matrix description of the platelet
p1.GetPlateletDef()
p2.GetPlateletDef()

#Build the matrizes that describe the linear programing problem
#A = np.vstack((np.vstack((p1.A, p1.A)),-np.identity(3)))
#Acsr = csr_matrix(A)
#b = np.append(np.vstack((p1.b, p1.b)), np.zeros(3))
#c = np.ones(3)
#dims = {'l': np.shape(A)[0], 'g': []}

# Build the matrizes that describe the linear programing problem
A = np.vstack((np.vstack((p1.A, p1.A)),-np.identity(3)))
Acsr = csr_matrix(np.vstack((np.vstack((p1.A[0], p1.A[2])),-np.identity(3))))
b = np.append(np.vstack((p1.b, p1.b)), np.zeros(3))
bcsr = np.append(np.vstack((p1.b[0], p1.b[2])), np.zeros(3))
c = np.ones(3)
dims = {'l': np.shape(Acsr)[0], 'g': []}

#A = csr_matrix(np.vstack((p1.A[0:1],-np.identity(3))))
#b = np.append(p1.b[0:1], np.zeros(3))
#c = np.ones(3)
#dims = {'l': np.shape(A)[0], 'g': []}

print('Shape A = ', np.shape(A))
print('Shape b = ', np.shape(b))
print('Shape c = ', np.shape(c))

#print('A = ', A)
print('b = ', b)
#print 'c = ', c

res = ecos.solve(c, Acsr, bcsr, dims)
print(res)



# plot the surface
fig = plt.figure()
ax = fig.gca(projection='3d')
xx, yy = np.meshgrid(range(0,110), range(0,110))

#ax.plot_surface(xx,yy, z1, alpha=0.5, color='red')
#ax.plot_surface(xx,yy, z3, alpha=0.5, color='yellow')
#ax.plot_surface(xx,yy, z4, alpha=0.5, color='cyan')
#ax.plot_surface(xx,yy, z5, alpha=0.5, color='black')

z0 = (A[0][0]*xx - A[0][1]*yy + b[0])/A[0][2]
z2 = (A[2][0]*xx - A[2][1]*yy + b[2])/A[2][2]
n0 = normvec(A[0])
n2 = normvec(A[2])
print('N0: ', n0,'  N2: ', n2)
print('Angle N0: ', angle(n0, [1,1,1])*57.295)
print('Angle N2: ', angle(n2, [1,1,1])*57.295)
#ax.quiver(xx,yy,z0, n0[0], n0[1], n0[2], length=0.5, color= 'red')
#ax.quiver(xx,yy,z2, n2[0], n2[1], n2[2], length=0.5, color = 'blue')
ax.plot_surface(xx,yy, z0, alpha=0.5, color='blue')
ax.plot_surface(xx,yy, z2, alpha=0.5, color='green')

#z1 = (A[1][0]*xx - A[1][1]*yy + b[1])*1./A[1][2]
#z3 = (A[3][0]*xx - A[3][1]*yy + b[3])*1./A[3][2]
#n1 = normvec(A[1])
#n3 = normvec(A[3])
#ax.quiver(xx,yy,z1, n1[0], n1[1], n1[2], length=0.5, color= 'red')
#ax.quiver(xx,yy,z3, n3[0], n3[1], n3[2], length=0.5, color = 'blue')
#ax.plot_surface(xx,yy, z1, alpha=0.5, color='blue')
#ax.plot_surface(xx,yy, z3, alpha=0.5, color='blue')

#z4 = (A[4][0]*xx - A[4][1]*yy + b[4])*1./A[4][2]
#z5 = (A[5][0]*xx - A[5][1]*yy + b[5])*1./A[5][2]
#n4 = normvec(A[4])
#n5 = normvec(A[5])
#ax.quiver(xx,yy,z4, n4[0], n4[1], n4[2], length=0.5, color= 'red')
#ax.quiver(xx,yy,z5, n5[0], n5[1], n5[2], length=0.5, color = 'blue')

opt = res['x']
ax.scatter(opt[0], opt[1], opt[2], color='red')

plt.show()

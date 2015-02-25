from ElectrodeGenerator import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

p1 = GetRandomPlatelet(1,[4,11,11], [4,2,0.5])
p2 = GetRandomPlatelet(2,[4,11,11], [4,2,0.5])

print np.linalg.norm(p1.GetPos()-p2.GetPos())

# Calculate the Matrix description of the platelet
p1.GetPlateletDef()
p2.GetPlateletDef()

## Build the matrizes that describe the linear programing problem
#A = csr_matrix(np.vstack((np.vstack((p1.A, p1.A)),-np.identity(3))))
#b = np.append(np.vstack((p1.b, p1.b)), np.zeros(3))
#c = np.ones(3)
#dims = {'l': np.shape(A)[0], 'g': []}

# Build the matrizes that describe the linear programing problem
A = csr_matrix(np.vstack((np.vstack((-p1.A[4], p1.A[5])),-np.identity(3))))
b = np.append(np.vstack((p1.b[0], p1.b[2])), np.zeros(3))
c = np.ones(3)
dims = {'l': np.shape(A)[0], 'g': []}

#A = csr_matrix(np.vstack((p1.A[0:2],-np.identity(3))))
#b = np.append(p1.b[0:2], np.zeros(3))
#c = np.ones(3)
#dims = {'l': np.shape(A)[0], 'g': []}
#print 'Position: ', p1.GetPos()

print 'Shape A = ', np.shape(A)
print 'Shape b = ', np.shape(b)
print 'Shape c = ', np.shape(c)

print 'A = ', A
#print 'b = ', b
#print 'c = ', c

res = ecos.solve(c, A, b, dims)
print res


xx, yy = np.meshgrid(range(30), range(30))


# calculate corresponding z
z1 = (-p1.A[0][0]*xx - p1.A[0][1]*yy - p1.b[0])*1./p1.A[0][2]
z2 = (-p1.A[1][0]*xx - p1.A[1][1]*yy - p1.b[1])*1./p1.A[1][2]
z3 = (-p1.A[2][0]*xx - p1.A[2][1]*yy - p1.b[2])*1./p1.A[2][2]
z4 = (-p1.A[3][0]*xx - p1.A[3][1]*yy - p1.b[3])*1./p1.A[3][2]
#z5 = (-p1.A[4][0]*xx - p1.A[4][1]*yy - p1.b[4])*1./p1.A[4][2]
#z6 = (-p1.A[5][0]*xx - p1.A[5][1]*yy - p1.b[5])*1./p1.A[5][2]
#
# plot the surface
plt3d = plt.figure()
ax = Axes3D(plt3d)
ax.plot_surface(xx,yy, z1, alpha=0.5, color='blue')
ax.plot_surface(xx,yy, z2, alpha=0.5, color='red')
ax.plot_surface(xx,yy, z3, alpha=0.5, color='green')
ax.plot_surface(xx,yy, z4, alpha=0.5, color='yellow')
#ax.plot_surface(xx,yy, z5, alpha=0.5, color='cyan')
#ax.plot_surface(xx,yy, z6, alpha=0.5, color='black')
plt.show()

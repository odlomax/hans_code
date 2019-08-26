import numpy as np
from matplotlib import pyplot as plt

# Fortran back-end interpolator
from vector_interpolate import vector_interp_3d

# Scipy interpolator
from scipy.interpolate import RegularGridInterpolator

# let's define a couple of 3d scalar fields
def exp_r(x,y,z):
    
    f=np.exp(-np.sqrt(x**2+y**2+z**2))
    
    return f

def exp2_r(x,y,z):
    
    f=np.exp(-(x**2+y**2+z**2))
    
    return f


# set some grid sizes
n_x=100
n_y=100
n_z=100

# make arrays of coordinate values
x_arr=np.linspace(0.,1.,n_x)
y_arr=np.linspace(0.,1.,n_y)
z_arr=np.linspace(0.,1.,n_z)

X,Y,Z=np.meshgrid(x_arr,y_arr,z_arr,indexing="ij")

# allocate memory for vector field (order="f" is important)
f_arr=np.zeros((2,n_x,n_y,n_z),dtype=np.float,order="f")

# get vector field
f_arr[0,...]=exp_r(X,Y,Z)
f_arr[1,...]=exp2_r(X,Y,Z)

# initialise Scipy interpolator
sp_interp_0=RegularGridInterpolator((x_arr,y_arr,z_arr),f_arr[0,...])
sp_interp_1=RegularGridInterpolator((x_arr,y_arr,z_arr),f_arr[1,...])

# initialise Fortran interpolator
f_interp=vector_interp_3d(x_arr,y_arr,z_arr,f_arr)


# do some imshows of the fields to make sure they're not garbage
plt.figure(0)
plt.imshow(f_interp.f_arr[0,...].mean(axis=2))
plt.figure(1)
plt.imshow(f_interp.f_arr[1,...].mean(axis=2))

# generate some random coordinates
n_point=10000
x=np.random.uniform(0.,1.,(n_point,3))

# allocate some arrays for field values
f=np.zeros((n_point,2))
g=np.zeros((n_point,2))



for i in range(n_point):
        
    # perform interpolations
    f[i,:]=f_interp(x[i,0],x[i,1],x[i,2])
    g[i,0]=sp_interp_0((x[i,0],x[i,1],x[i,2]))
    g[i,1]=sp_interp_1((x[i,0],x[i,1],x[i,2]))
    
    
# plot differences between interpolated values 
# should be of order 1e-16
plt.figure(2)
plt.plot(x[:,0],f[:,0]-g[:,0],".")
plt.plot(x[:,0],f[:,1]-g[:,1],".")
plt.show()

# plot field projections
plt.figure(3)
plt.plot(x[:,0],f[:,1],".")
plt.plot(x[:,0],f[:,0],".")
plt.show()

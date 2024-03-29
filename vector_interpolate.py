# class to perform efficient vector-field interpolation
import numpy as np
from f_module import m_interpolate

class vector_interp_3d:
    
    # 3d vector field interpolator
    
    def __init__(self,x_arr,y_arr,z_arr,f_list=None,f_ndarray=None):
        
        """
        
        Subroutine: initialise interpolator
        
        Arguments
        ---------
        
        x_arr[n_x]: float
            array of x values
        
        y_arr[n_y]: float
            array of y values
            
        z_arr[n_z]: float
            array of z values
            
        f_list[n_v]: float (optional)
            list of n_v scalar fields with shape [n_x,n_y,n_z]
        
        f_ndarray[n_v,n_x,n_y,n_z]: float (optional)
             vector field
        
        """
        
        # get sizes
        n_x=x_arr.size
        n_y=y_arr.size
        n_z=z_arr.size
        
        # associate array
        self.x_arr=x_arr
        self.y_arr=y_arr
        self.z_arr=z_arr
        if f_list is not None:
            
            # pack f_arr into column-major 4d array
            self.f_arr=np.zeros((len(f_list),n_x,n_y,n_z),dtype=np.float,order="f")
            for i in range(len(f_list)):
                self.f_arr[i,...]=f_list[i]
                
        elif f_ndarray is not None:
            
            # associate array
            if np.isfortran(f_ndarray):
                
                # already in column-major order
                self.f_arr=f_ndarray
            else:
                
                # convert to column-major order
                self.f_arr=np.asfortranarray(f_ndarray)
            
        else:
            print("Error, no f_arr given.")
            return
        
        
        # calc 1/dx, 1/dy, 1/dz
        self.inv_dx=(n_x-1)/(self.x_arr[-1]-self.x_arr[0])
        self.inv_dy=(n_y-1)/(self.y_arr[-1]-self.y_arr[0])
        self.inv_dz=(n_z-1)/(self.z_arr[-1]-self.z_arr[0])
        
        return
    
    def __call__(self,x,y,z):
        
        """
        
        Function: perform interpolation
        
        Arguments
        ---------
        
        x,y,z: float
            x, y and z positions (scalars)
        
        Result
        ------
        
        f[n_v]: float
            field value at (x,y,z)
        
        """
        
        # ENGAGE FORTRAN!
        f=m_interpolate.vector_interp_3d(x,y,z,self.inv_dx,self.inv_dy,self.inv_dz,
          self.x_arr,self.y_arr,self.z_arr,self.f_arr)
        
        return f
module m_interpolate

    ! module for linear interpolation functions
    
    implicit none
    
    ! set functions to private by default
    private
    public :: vector_interp_3d

    contains
    
    ! interpolate values from a 3d vector field grid
    pure function vector_interp_3d(x,y,z,inv_dx,inv_dy,inv_dz,x_arr,y_arr,z_arr,f_arr) result(f)
    
        ! argument declarations
        real(kind=8),intent(in) :: x                       ! x position
        real(kind=8),intent(in) :: y                       ! y position
        real(kind=8),intent(in) :: z                       ! z position
        real(kind=8),intent(in) :: inv_dx                  ! 1/dx
        real(kind=8),intent(in) :: inv_dy                  ! 1/dy
        real(kind=8),intent(in) :: inv_dz                  ! 1/dz
        real(kind=8),intent(in) :: x_arr(:)                ! array of x grid points (n_x_grid)
        real(kind=8),intent(in) :: y_arr(:)                ! array of y grid points (n_y_grid)
        real(kind=8),intent(in) :: z_arr(:)                ! array of z grid points (n_z_grid)
        real(kind=8),intent(in) :: f_arr(:,:,:,:)          ! array of f(x,y,z) grid points (n_v,n_x_grid,n_y_grid,n_z_grid)
        
        ! result declaration
        real(kind=8) :: f(size(f_arr,1))                   ! array of interpolated f(x,y,z) values (n_v,n_point)
        
        ! variable declarations
        integer :: i                                       ! x axis index
        integer :: j                                       ! y axis index
        integer :: k                                       ! z axis index
        real(kind=8) :: w_x                                ! x axis interpolation weight
        real(kind=8) :: w_y                                ! y axis interpolation weight
        real(kind=8) :: w_z                                ! z axis interpolation weight
        
        ! get index between 1 and n-1
        i=1+floor((x-x_arr(1))*inv_dx)
        i=min(max(i,1),size(x_arr)-1)
        ! set weight
        w_x=(x-x_arr(i))*inv_dx
        
        ! get index between 1 and n-1
        j=1+floor((y-y_arr(1))*inv_dy)
        j=min(max(j,1),size(y_arr)-1)
        ! set weight
        w_y=(y-y_arr(j))*inv_dy
        
        ! get index between 1 and n-1
        k=1+floor((z-z_arr(1))*inv_dz)
        k=min(max(k,1),size(z_arr)-1)
        ! get weight
        w_z=(z-z_arr(k))*inv_dz
        
        ! perform interpolation
        f=trilerp(w_x,w_y,w_z,&
            &f_arr(:,i,j,k),f_arr(:,i,j,k+1),f_arr(:,i,j+1,k),f_arr(:,i,j+1,k+1),&
            &f_arr(:,i+1,j,k),f_arr(:,i+1,j,k+1),f_arr(:,i+1,j+1,k),f_arr(:,i+1,j+1,k+1))
        
        
        return
    
    end function
    
    ! linear interpolation 
    elemental function lerp(w_x,f_0,f_1) result(f)
   
        ! argument declarations
        real(kind=8),intent(in) :: w_x                     ! (x-x_0)/(x_1-x_0)
        real(kind=8),intent(in) :: f_0                     ! f(x) interpolation points
        real(kind=8),intent(in) :: f_1
        
        ! result declaration
        real(kind=8) :: f                                  ! f(x)
        
        ! solve for f
        f=f_0+(f_1-f_0)*w_x
        
        return
   
    end function
    
    ! bilinear interpolation
    elemental function bilerp(w_x,w_y,f_00,f_01,f_10,f_11) result(f)
    
        ! argument declarations
        real(kind=8),intent(in) :: w_x                     ! (x-x_0)/(x_1-x_0)
        real(kind=8),intent(in) :: w_y                     ! (y-y_0)/(y_1-y_0) 
        real(kind=8),intent(in) :: f_00                    ! f(x,y) values
        real(kind=8),intent(in) :: f_01
        real(kind=8),intent(in) :: f_10
        real(kind=8),intent(in) :: f_11
        
        ! result declaration
        real(kind=8) :: f                                  ! f(x,y)
        
        ! solve for f
        f=lerp(w_y,lerp(w_x,f_00,f_10),lerp(w_x,f_01,f_11))
        
        return
        
    end function
    
    ! trilinear interpolation
    elemental function trilerp(w_x,w_y,w_z,&
        &f_000,f_001,f_010,f_011,f_100,f_101,f_110,f_111) result(f)
    
        ! argument declarations
        real(kind=8),intent(in) :: w_x                     ! (x-x_0)/(x_1-x_0)
        real(kind=8),intent(in) :: w_y                     ! (y-y_0)/(y_1-y_0)
        real(kind=8),intent(in) :: w_z                     ! (z-z_0)/(z_1-z_0)
        real(kind=8),intent(in) :: f_000                   ! f(x,y) values
        real(kind=8),intent(in) :: f_001
        real(kind=8),intent(in) :: f_010
        real(kind=8),intent(in) :: f_011
        real(kind=8),intent(in) :: f_100
        real(kind=8),intent(in) :: f_101
        real(kind=8),intent(in) :: f_110
        real(kind=8),intent(in) :: f_111
        
        ! result declaration
        real(kind=8) :: f                                  ! f(x,y)
        
        ! solve for f
        f=lerp(w_z,bilerp(w_x,w_y,f_000,f_010,f_100,f_110),&
            &bilerp(w_x,w_y,f_001,f_011,f_101,f_111))
        
        return
        
    end function
    
end module

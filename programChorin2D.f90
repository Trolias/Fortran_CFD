program Chorin2D
    implicit none

    ! Parameters definition
    integer :: nx, ny, it, ITERATIONS, POISSON_ITERATIONS, i
    real    :: L, h, w, Re, dt, rho, mu, nu, dx, dy, u_inlet, dp
    real,allocatable,dimension(:,:) :: u, v, P, u_new, v_new, P_new, gx, gy, u_temp, v_temp,b 
    real,allocatable,dimension(:) :: x, y, ua
    character(len=45) :: analytical, numerical
    
    analytical = 'Analytical_outlet.dat'
    numerical = 'Chorin_2ndOrder_Pressure_outlet.dat'

    Re = 100.0
    L = 1.0         ! Length of the channel
    h = 1.0         ! Height of the channel
    w = 1.0         ! Width of the channel

    ! Mesh construction and domain definition
    nx = 41
    ny = 41
    dx = L/(nx-1)
    dy = h/(ny-1)
    allocate(x(nx)); x = 0.0
    allocate(y(ny)); y = 0.0
    allocate(ua(ny)); ua = 0.0

    ! X-direction Points
    do i = 1, nx
        x(i) = (i-1)*dx
    end do

    ! Y-direction Points
    do i = 1, ny
        y(i) = (i-1)*dy
    end do

    ! Initialize Fields
    allocate(u(nx, ny)); u = 0.0
    allocate(u_new(nx, ny)); u_new = 0.0
    allocate(v(nx, ny)); v = 0.0
    allocate(v_new(nx, ny)); v_new = 0.0
    allocate(u_temp(nx,ny)); u_temp = 0.0
    allocate(v_temp(nx,ny)); v_temp = 0.0
    allocate(P(nx, ny)); P = 0.0
    allocate(P_new(nx, ny)); P_new = 0.0
    allocate(b(nx, ny)); b = 0.0
    allocate(gx(nx, ny)); gx = 0.0
    allocate(gy(nx, ny)); gy = 0.0

    ! Simulation Parameters
    rho = 1     ! Density of water [kg/m^3]
    mu = 0.1   ! Dynamic viscosity [Pa·s]
    nu = mu/rho     ! Kinematic viscosity
    u_inlet = Re*nu/(2.0*h)
    dp = (12.0*mu*L*u_inlet)/(h*h)

    ! TIME STEP 
    dt = 0.001

    ! Artificial Compressibility parameters
    POISSON_ITERATIONS = 300
    ITERATIONS = 1500

    ! Initial Condition
    ! u(2:nx-1, :) = u_inlet
    u(:,ny) = 1.0

! SOLUTION ALGORITHM
    do it=1,ITERATIONS
    ! SOLVE MOMENTUM EQUATIONS FOR ROTATIONAL PART 
    
        ! Momentum x-direction
        u_temp = u + dt * ( &
            - u*Upwind_x(u, nx, ny, dx) - v*Upwind_y(u, nx, ny, dy) &
            + nu*Diff_Central(u, nx, ny, dx, dy)                    &
        )
        ! Momentum y-direction
        v_temp = v + dt * ( &
            - u*Upwind_x(v, nx, ny, dx) - v*Upwind_y(v, nx, ny, dy) &
            + nu*Diff_Central(v, nx, ny, dx, dy)                    &
        )

    ! BOUNDARY CONSITIONS FOR TEMPORARY VELOCITIES
        ! ------------ X-DIRECTION ---------------
        u_temp(1, :) = 0.0
        u_temp(:, 1) = 0.0
        u_temp(:, ny) = 1.0
        u_temp(nx, :) = 0.0
        ! ------------ Y-DIRECTION ---------------
        v_temp(1, :) = 0.0
        v_temp(nx, :) = 0.0
        v_temp(:, 1) = 0.0
        v_temp(:, ny) = 0.0

    ! SOLVE PRESSURE POISSON
        ! Source term definition
        b = rho/dt * ( Upwind_x(u_temp, nx, ny, dx) + Upwind_y(v_temp, nx, ny, dy) )
        ! Pressure-Poisson
        call Pressure_Poisson_GaussSeidel(P,nx,ny,dx,dy,POISSON_ITERATIONS,b)
        
    ! CALCULATE DIVERGENCE FREE VELOCITY
        ! x-velocity
        u_new = u_temp - dt/rho * Grad_x(P, nx, ny, dx)
        v_new = v_temp - dt/rho * Grad_y(P, nx, ny, dy) 
    
    ! VELOCITY BOUNDARY CONDITIONS
    ! -------------- X-DIRECTION ---------------------
        u_new(1, :) = 0.0
        u_new(:, 1) = 0.0
        u_new(:, ny) = 1.0
        u_new(nx, :) = 0.0
    ! -------------- Y-DIRECTION ----------------------
        v_new(1, :) = 0.0
        v_new(nx, :) = 0.0
        v_new(:, 1) = 0.0
        v_new(:, ny) = 0.0

    ! UPDATE VELOCITIES
        u = u_new
        v = v_new

    end do 


    ! Analytical Solution
    do i = 1, ny
        ua(i) = (dp/(2.0*mu*L)) * y(i) * (h - y(i))
    end do
    ! Save data to files
    call Output
    call writeTofile(ua, ny, y, analytical)
    call writeTofile(u_new(nx, :), ny, y, numerical)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUBROUTINES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

    function Upwind_x(u, nx, ny, dx) result(res)
        integer, intent(in) :: nx, ny
        real, intent(in) :: u(nx, ny), dx
        real :: res(nx, ny)
        integer :: i, j
        res = 0.0
        do j = 2, ny-1
            do i = 2, nx-1
                ! res(i, j) = (u(i, j) - u(i-1, j)) / dx
                res(i, j) = (u(i+1, j) - u(i-1, j)) / (2.0*dx)
            end do
        end do
    end function Upwind_x

    function Upwind_y(u, nx, ny, dy) result(res)
        integer, intent(in) :: nx, ny
        real, intent(in) :: u(nx, ny), dy
        real :: res(nx, ny)
        integer :: i, j
        res = 0.0
        do j = 2, ny-1
            do i = 2, nx-1
                ! res(i, j) = (u(i, j) - u(i, j-1)) / dy
                res(i, j) = (u(i, j+1) - u(i, j-1)) / (2.0*dy)

            end do
        end do
    end function Upwind_y

    function Diff_Central(u, nx, ny, dx, dy) result(res)
        integer, intent(in) :: nx, ny
        real, intent(in) :: u(nx, ny), dx, dy
        real :: res(nx, ny)
        integer :: i, j
        res = 0.0
        do j = 2, ny-1
            do i = 2, nx-1
                res(i, j) = (u(i+1, j) - 2.0*u(i, j) + u(i-1, j))/(dx*dx) &
                          + (u(i, j+1) - 2.0*u(i, j) + u(i, j-1))/(dy*dy)
            end do
        end do
    end function Diff_Central

    function Grad_x(P, nx, ny, dx) result(res)
        integer, intent(in) :: nx, ny
        real, intent(in) :: P(nx, ny), dx
        real :: res(nx, ny)
        integer :: i, j
        res = 0.0
        do j = 2, ny-1
            do i = 2, nx-1
                ! res(i, j) = (P(i+1, j) - P(i, j)) / dx  
                res(i, j) = (P(i+1, j) - P(i-1, j)) / (2.0*dx)  

            end do
        end do
    end function Grad_x

    function Grad_y(P, nx, ny, dy) result(res)
        integer, intent(in) :: nx, ny
        real, intent(in) :: P(nx, ny), dy
        real :: res(nx, ny)
        integer :: i, j
        res = 0.0
        do j = 2, ny-1
            do i = 2, nx-1
                ! res(i, j) = (P(i, j+1) - P(i, j)) / dy 
                res(i, j) = (P(i, j+1) - P(i, j-1)) / (2.0*dy)  
            end do
        end do
    end function Grad_y

    subroutine Pressure_Poisson_GaussSeidel(P,nx,ny,dx,dy,iter,b)
        integer i,J,n
        integer nx,ny, iter
        real dx,dy
        real P(nx,ny), P_old(nx,ny), b(nx,ny)

        do n=1,iter
            P_old = P
            do i=2,nx-1
                do j=2,ny-1
                    P(i,j) = (  1.0 /( 2.0 * (dx*dx + dy*dy ) ) * &
                                ( &
                                 (dy*dy) * (p_old(i+1,j)+P(i-1,j)) &
                                +(dx*dx) * (p_old(i,j+1)+P(i,j-1)) & 
                                - b(i,j) * (dx*dx) * (dy*dy)           &
                                ) &
                             )
                end do
            end do
        ! BOUNDARY CONDITIONS
            P(1,1) = 0.0
            P(:, 1)  = P(:, 2)
            P(:, ny) = P(:,ny-1)
            P(1, :)  = P(2, :)
            P(nx,:)  = P(nx-1,:)
            
        ! 2nd-Order for Outletfor pipe problem
            ! do j = 2, ny-1
            !     ! P(nx, j) = 2.0*P(nx-1, j) - P(nx-2, j)
            ! end do

        ! Pressure Boundary Conditions at corners
            ! P(1, 1) = 0.5*(P(1, 2) + P(2, 1))
            ! P(nx, 1) = 0.5*(P(nx, 2) + P(nx-1, 1))
            ! P(1, ny) = 0.5*(P(2, ny) + P(1, ny-1))
            ! P(nx, ny) = 0.5*(P(nx-1, ny) + P(nx, ny-1))

            
        end do 


    end subroutine

! %%%%%%%%%%%%%%  WRITING DATA ROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine Output
        integer :: i, j
        open(1, file='Chorin_2ndOrder_Paraview.dat')
        write(1,*) 'TITLE="Solution"'
        write(1,*) 'VARIABLES="X" "Y" "u" "v" "P"'
        write(1,*) 'ZONE I=', nx, ', J=', ny, ', F=POINT'
        do j = 1, ny
            do i = 1, nx
                write(1, '(5(2X,E11.4))') x(i), y(j), u_new(i, j), v_new(i, j), P(i, j)
            end do
        end do
        close(1)
    end subroutine Output

    subroutine writeTofile(field, nx, vec_x, filename)
        integer, intent(in) :: nx
        real, intent(in) :: field(nx), vec_x(nx)
        character(len=*), intent(in) :: filename
        integer :: i
        open(15, file=filename)
        do i = 1, nx
            write(15, *) vec_x(i), field(i)
        end do
        close(15)
    end subroutine writeTofile

end program Chorin2D
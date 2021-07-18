module poisson_direct_dp

    use param    
    use iso_fortran_env, only : dp => real64
                                   
    implicit none

    private

    public  :: poisson_solver_2d_dp, &
               poisson_solver_3d_dp

    !====================================

    type poisson_solver_2d_dp

        private

#ifdef OPENMP
        integer                                :: nproc, np, resi
#endif
        integer                                :: order
        integer                                ::   nx,  ny
        real(dp), dimension(:),   allocatable  ::  evx, evy
        real(dp), dimension(:,:), allocatable  ::   ex,  ey

    contains

        procedure, public, pass(self) :: init => poisson_2d_init_dp
        procedure, public, pass(self) :: solve => poisson_2d_solve_dp
        procedure, public, pass(self) :: destroy => poisson_2d_destroy_dp

    end type poisson_solver_2d_dp

    !====================================

    type poisson_solver_3d_dp

        private

#ifdef OPENMP
        integer                                :: nproc, np, resi
#endif
        integer                                :: order
        integer                                ::   nx,  ny,  nz
        real(dp), dimension(:),   allocatable  ::  evx, evy, evz
        real(dp), dimension(:,:), allocatable  ::   ex,  ey,  ez

    contains

        procedure, public, pass(self) :: init => poisson_3d_init_dp
        procedure, public, pass(self) :: solve => poisson_3d_solve_dp
        procedure, public, pass(self) :: destroy => poisson_3d_destroy_dp

    end type poisson_solver_3d_dp


    real(dp), parameter               :: eps = epsilon(1.0_dp)

contains

    !2d section double precision

    subroutine poisson_2d_init_dp(self, np, step, bdy)
#ifdef OPENMP
        use omp_lib
#endif
        implicit none
        class(poisson_solver_2d_dp), intent(inout)  :: self
        integer,  dimension(2),      intent(in)     :: np
        real(dp), dimension(2),      intent(in)     :: step
        integer,  dimension(4),      intent(in)     :: bdy  !boundary condition

        !local
        real(dp), dimension(2)                      :: tmp

        self % order = maxval( nint( log10(step) ) )
        tmp          = step / 10 ** self % order

        !compute eigenvalues and eigenvectors
        self % nx = np(1)
        self % ny = np(2)
#ifdef OPENMP
        self % nproc = omp_get_max_threads()
        self % np    = np(2) / self % nproc
        self % resi  = self % nproc - mod(np(2), self % nproc)
#endif
        allocate( self % evx(np(1)),   &
                  self % evy(np(2)),   &
                  self % ex(np(1),np(1)), &
                  self % ey(np(2),np(2)) )

        call eigen_dp(np(1), tmp(1), self%evx, self%ex, bdy(1), bdy(2))
        call eigen_dp(np(2), tmp(2), self%evy, self%ey, bdy(3), bdy(4))
    end subroutine poisson_2d_init_dp

    subroutine poisson_2d_solve_dp(self, f, p)
#ifdef OPENMP
        use omp_lib
#endif
        !--------------------------------------------------------
        ! This subroutine solves the eq : del p = f
        !--------------------------------------------------------
        ! Parameters:
        ! f         forcing term of poisson eq
        ! p         solution of poisson eq
        !--------------------------------------------------------
        implicit none
        class(poisson_solver_2d_dp), intent(in)  :: self
        real(dp), dimension(:, :),   intent(in)  :: f
        real(dp), dimension(:, :),   intent(out) :: p       


        !local variables
        integer                                :: i, j
        integer                                :: nx, ny
        real(dp), dimension(:, :), allocatable :: f1, f2
#ifdef OPENMP
        integer                                :: nproc, np, resi, myid, jds, jde
        real(dp), dimension(:, :), allocatable :: ftmp

        nproc = self % nproc
        np    = self % np
        resi  = self % resi
#endif

        nx = self % nx
        ny = self % ny

        !check the dimensions are the same as we set in init
        if (nx /= size(f,1) .or. nx /= size(p,1)) stop "First  dimension is not the same"
        if (ny /= size(f,2) .or. ny /= size(p,2)) stop "Second dimension is not the same"

        !dimensions are Ok, then allocate work arrays
#ifdef OPENMP
        allocate(ftmp(nx, ny))

        !$omp parallel firstprivate(np) private(myid, i, j, jds, jde, f1, f2)
        myid = omp_get_thread_num()
        jds  = np *  myid + 1  + max(myid - resi,     0)
        jde  = np * (myid + 1) + max(myid - resi + 1, 0)
        np   = jde - jds + 1

        allocate(f1(nx, jds:jde), f2(nx, jds:jde))

        call dgemm('n','n',nx,np,ny,1.0_dp,      f,nx,self%ey(:,jds:jde),ny,0.0_dp,f1,nx)
        call dgemm('t','n',nx,np,nx,1.0_dp,self%ex,nx,                f1,nx,0.0_dp,f2,nx)

        do concurrent (i=1:nx,j=jds:jde, abs(self%evx(i) + self%evy(j)) > eps)
            f2(i,j) = f2(i,j) / (self%evx(i) + self%evy(j))
        end do

        call dgemm('n','n',nx,np,nx,1.0_dp,self%ex,nx,                f2,nx,0.0_dp,ftmp(:,jds:jde),nx)
        !$omp barrier
        call dgemm('n','t',nx,np,ny,1.0_dp,   ftmp,nx,self%ey(jds:jde,:),np,0.0_dp,   p(:,jds:jde),nx)

        deallocate(f1, f2)

        p(:,jds:jde) = p(:,jds:jde) * 10 ** (2 * self % order)

        !$omp end parallel

        deallocate(ftmp)
#else
        allocate(f1(nx, ny), f2(nx, ny))

        call dgemm('n','n',nx,ny,ny,1.0_dp,      f,nx,self%ey,ny,0.0_dp,f1,nx)
        call dgemm('t','n',nx,ny,nx,1.0_dp,self%ex,nx,     f1,nx,0.0_dp,f2,nx)

        do concurrent (j=1:ny,i=1:nx, abs(self%evx(i) + self%evy(j)) > eps)
            f2(i,j) = f2(i,j) / (self%evx(i) + self%evy(j))
        end do

        call dgemm('n','n',nx,ny,nx,1.0_dp,self%ex,nx,     f2,nx,0.0_dp,f1,nx)
        call dgemm('n','t',nx,ny,ny,1.0_dp,     f1,nx,self%ey,ny,0.0_dp, p,nx)

        deallocate(f1, f2)      

        p = p * 10 ** (2 * self % order)
#endif
    end subroutine poisson_2d_solve_dp

    subroutine poisson_2d_destroy_dp(self)
        implicit none
        class(poisson_solver_2d_dp),  intent(inout)  :: self

        deallocate( self % evx, &
                    self % evy, &
                    self %  ex, &
                    self %  ey )
    end subroutine poisson_2d_destroy_dp


    !3d section double precision

    subroutine poisson_3d_init_dp(self, np, step, bdy)
#ifdef OPENMP
        use omp_lib
#endif
        implicit none
        class(poisson_solver_3d_dp), intent(inout)  :: self
        integer,  dimension(3),      intent(in)     :: np
        real(dp), dimension(3),      intent(in)     :: step
        integer,  dimension(6),      intent(in)     :: bdy  !boundary condition

        !local
        real(dp), dimension(3)                      :: tmp

        self % order = maxval( nint( log10(step) ) )
        tmp          = step / 10 ** self % order

        !compute eigenvalues and eigenvectors
        self % nx = np(1)
        self % ny = np(2)
        self % nz = np(3)
#ifdef OPENMP
        self % nproc = omp_get_max_threads()
        self % np    = np(3) / self % nproc
        self % resi  = self % nproc - mod(np(3), self % nproc)
#endif
        allocate( self % evx(np(1)),   &
                  self % evy(np(2)),   &
                  self % evz(np(3)),   &
                  self % ex(np(1),np(1)), &
                  self % ey(np(2),np(2)), &
                  self % ez(np(3),np(3)) )

        call eigen_dp(np(1), tmp(1), self%evx, self%ex, bdy(1), bdy(2))
        call eigen_dp(np(2), tmp(2), self%evy, self%ey, bdy(3), bdy(4))
        call eigen_dp(np(3), tmp(3), self%evz, self%ez, bdy(5), bdy(6))
    end subroutine poisson_3d_init_dp

    subroutine poisson_3d_solve_dp(self, f, p)
#ifdef OPENMP
        use omp_lib
#endif
        !--------------------------------------------------------
        ! This subroutine solves the eq : del p = f
        !--------------------------------------------------------
        ! Parameters:
        ! f         forcing term of poisson eq
        ! p         solution of poisson eq
        !--------------------------------------------------------
        implicit none
        class(poisson_solver_3d_dp),  intent(in)  :: self
        real(dp), dimension(:, :, :), intent(in)  :: f
        real(dp), dimension(:, :, :), intent(out) :: p       


        !local variables
        integer                                   :: i, j, k
        integer                                   :: nx, ny, nz, nxny
        integer                                   :: nxnz, nynz
        real(dp), dimension(:, :, :), allocatable :: f1, f2, f3, f4
#ifdef OPENMP
        integer                                   :: nproc, np, resi, myid, kds, kde
        real(dp), dimension(:, :, :), allocatable :: ftmp
#endif


        nx = self % nx
        ny = self % ny
        nz = self % nz
        nxny = nx * ny
#ifdef OPENMP
        nproc = self % nproc
        np    = self % np
        resi  = self % resi
#else
        nxnz = nx * nz
        nynz = ny * nz
#endif

        !check the dimensions are the same as we set in init
        if (nx /= size(f,1) .or. nx /= size(p,1)) stop "First  dimension is not the same"
        if (ny /= size(f,2) .or. ny /= size(p,2)) stop "Second dimension is not the same"
        if (nz /= size(f,3) .or. nz /= size(p,3)) stop "Third  dimension is not the same"

#ifdef OPENMP
        allocate(ftmp(nx, ny, nz))

        !$omp parallel firstprivate(np) private(myid, i, j ,k, kds, kde, nxnz, nynz, f1, f2, f3, f4)
        myid = omp_get_thread_num()
        kds  = np *  myid + 1  + max(myid - resi,     0)
        kde  = np * (myid + 1) + max(myid - resi + 1, 0)
        np   = kde - kds + 1
        nxnz = nx * np
        nynz = ny * np

        !dimensions are Ok, then allocate work arrays
        allocate( f1(nx, ny, kds:kde), f2(nx, ny, kds:kde), &
                  f3(ny, nx, kds:kde), f4(ny, nx, kds:kde) )

        call dgemm('n','n',nxny,np,nz,1.0_dp,      f,nxny,self%ez(:,kds:kde),nz,0.0_dp,f1,nxny)
        call dgemm('t','n',nx,nynz,nx,1.0_dp,self%ex,  nx,                f1,nx,0.0_dp,f2,  nx)

        do concurrent (k=kds:kde)
            f3(:,:,k) = transpose(f2(:,:,k))
        end do

        call dgemm('t','n',ny,nxnz,ny,1.0_dp,self%ey,  ny,     f3,ny,0.0_dp,f4,  ny)

        do concurrent (j=1:ny,i=1:nx,k=kds:kde, abs(self%evx(i) + self%evy(j) + self%evz(k)) > eps)
            f4(j,i,k) = f4(j,i,k) / (self%evx(i) + self%evy(j) + self%evz(k))
        end do

        call dgemm('n','n',ny,nxnz,ny,1.0_dp,self%ey,  ny,     f4,ny,0.0_dp,f3,  ny)

        do concurrent (k=kds:kde)
            f2(:,:,k) = transpose(f3(:,:,k))
        end do

        call dgemm('n','n',nx,nynz,nx,1.0_dp,self%ex,  nx,                f2,nx,0.0_dp, ftmp(:,:,kds:kde),  nx)
        !$omp barrier
        call dgemm('n','t',nxny,np,nz,1.0_dp,   ftmp,nxny,self%ez(kds:kde,:),np,0.0_dp,    p(:,:,kds:kde),nxny)

        deallocate(f1, f2, f3, f4)

        p(:,:,kds:kde) = p(:,:,kds:kde) * 10 ** (2 * self % order)

        !$omp end parallel

        deallocate(ftmp)
#else
        !dimensions are Ok, then allocate work arrays
        allocate( f1(nx, ny, nz), f2(nx, ny, nz), &
                  f3(ny, nx, nz), f4(ny, nx, nz) )

        call dgemm('n','n',nxny,nz,nz,1.0_dp,      f,nxny,self%ez,nz,0.0_dp,f1,nxny)
        call dgemm('t','n',nx,nynz,nx,1.0_dp,self%ex,  nx,     f1,nx,0.0_dp,f2,  nx)

        do concurrent (k=1:nz)
            f3(:,:,k) = transpose(f2(:,:,k))
        end do

        call dgemm('t','n',ny,nxnz,ny,1.0_dp,self%ey,  ny,     f3,ny,0.0_dp,f4,  ny)

        do concurrent (j=1:ny,i=1:nx,k=1:nz, abs(self%evx(i) + self%evy(j) + self%evz(k)) > eps)
            f4(j,i,k) = f4(j,i,k) / (self%evx(i) + self%evy(j) + self%evz(k))
        end do

        call dgemm('n','n',ny,nxnz,ny,1.0_dp,self%ey,  ny,     f4,ny,0.0_dp,f3,  ny)

        do concurrent (k=1:nz)
            f2(:,:,k) = transpose(f3(:,:,k))
        end do

        call dgemm('n','n',nx,nynz,nx,1.0_dp,self%ex,  nx,     f2,nx,0.0_dp,f1,  nx)
        call dgemm('n','t',nxny,nz,nz,1.0_dp,     f1,nxny,self%ez,nz,0.0_dp, p,nxny)

        deallocate(f1, f2, f3, f4)

        p = p * 10 ** (2 * self % order)
#endif
    end subroutine poisson_3d_solve_dp

    subroutine poisson_3d_destroy_dp(self)
        implicit none
        class(poisson_solver_3d_dp),  intent(inout)  :: self
        deallocate( self % evx, &
                    self % evy, &
                    self % evz, &
                    self %  ex, &
                    self %  ey, &
                    self %  ez )
    end subroutine poisson_3d_destroy_dp

    subroutine eigen_dp(n, dx, ev, ee, first, last)
        implicit none
        integer,                   intent(in)  :: n, first, last
        real(dp),                  intent(in)  :: dx
        real(dp), dimension(n),    intent(out) :: ev
        real(dp), dimension(n, n), intent(out) :: ee

        !local variables
        integer                             :: info
        real(dp)                            :: ddxx
        real(dp), dimension(n)              :: d
        real(dp), dimension(n-1)            :: e
        real(dp), dimension(:), allocatable :: work
        integer,  dimension(:), allocatable :: iwork
        integer,  dimension(2*n)            :: isuppz
        integer                             :: lwork, liwork, m


        ddxx   = 1.0_dp / (dx*dx)

        !diagonal
        select case(first)
        case (dirichlet)
            d(1)         = -2.0_dp * ddxx
        case (neumann)
            d(1)         = -ddxx
        case default
            stop "Unknown boundary condition, please select 'dirichlet' or 'neumann'."
        end select

        select case(last)
        case (dirichlet)
            d(n)      = -2.0_dp * ddxx
        case (neumann)
            d(n)      = -ddxx
        case default
            stop "Unknown boundary condition, please select 'dirichlet' or 'neumann'."
        end select

        d(2:n-1)  = -2.0_dp * ddxx

        !subdiagonal
        e = ddxx

        !========================
        ! call lapack subroutine
        !========================
        !first calling is getting the optimal workspace
        lwork  = -1
        liwork = -1
        allocate(work(1), iwork(1))
        call dstevr('v','a',n,d,e,0.0_dp,0.0_dp,0,0,0.0_dp,m,ev,ee,n,isuppz,work,lwork,iwork,liwork,info)

        lwork  = int(work(1))
        liwork = iwork(1)
        deallocate(work, iwork)

        !second calling is calculating eigenvalues and eigenvectors
        allocate(work(lwork), iwork(liwork))
        call dstevr('v','a',n,d,e,0.0_dp,0.0_dp,0,0,0.0_dp,m,ev,ee,n,isuppz,work,lwork,iwork,liwork,info)
    end subroutine eigen_dp

end module poisson_direct_dp

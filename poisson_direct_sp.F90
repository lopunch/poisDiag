module poisson_direct_sp

    use param    
    use iso_fortran_env, only : sp => real32
                                   
    implicit none

    private

    public  :: poisson_solver_2d_sp, &
               poisson_solver_3d_sp

    !====================================

    type poisson_solver_2d_sp

        private

#ifdef OPENMP
        integer                                :: nproc, np, resi
#endif
        integer                                :: order
        integer                                ::   nx,  ny
        real(sp), dimension(:),   allocatable  ::  evx, evy
        real(sp), dimension(:,:), allocatable  ::   ex,  ey

    contains

        procedure, public, pass(self) :: init => poisson_2d_init_sp
        procedure, public, pass(self) :: solve => poisson_2d_solve_sp
        procedure, public, pass(self) :: destroy => poisson_2d_destroy_sp

    end type poisson_solver_2d_sp

    !====================================

    type poisson_solver_3d_sp

        private

#ifdef OPENMP
        integer                                :: nproc, np, resi
#endif
        integer                                :: order
        integer                                ::   nx,  ny,  nz
        real(sp), dimension(:),   allocatable  ::  evx, evy, evz
        real(sp), dimension(:,:), allocatable  ::   ex,  ey,  ez

    contains

        procedure, public, pass(self) :: init => poisson_3d_init_sp
        procedure, public, pass(self) :: solve => poisson_3d_solve_sp
        procedure, public, pass(self) :: destroy => poisson_3d_destroy_sp

    end type poisson_solver_3d_sp


    real(sp), parameter               :: eps = epsilon(1.0_sp)

contains

    !2d section double precision

    subroutine poisson_2d_init_sp(self, np, step, bdy)
#ifdef OPENMP
        use omp_lib
#endif
        implicit none
        class(poisson_solver_2d_sp), intent(inout)  :: self
        integer,  dimension(2),      intent(in)     :: np
        real(sp), dimension(2),      intent(in)     :: step
        integer,  dimension(4),      intent(in)     :: bdy  !boundary condition

        !local
        real(sp), dimension(2)                      :: tmp

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

        call eigen_sp(np(1), tmp(1), self%evx, self%ex, bdy(1), bdy(2))
        call eigen_sp(np(2), tmp(2), self%evy, self%ey, bdy(3), bdy(4))
    end subroutine poisson_2d_init_sp

    subroutine poisson_2d_solve_sp(self, f, p)
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
        class(poisson_solver_2d_sp), intent(in)  :: self
        real(sp), dimension(:, :),   intent(in)  :: f
        real(sp), dimension(:, :),   intent(out) :: p       


        !local variables
        integer                                :: i, j
        integer                                :: nx, ny
        real(sp), dimension(:, :), allocatable :: f1, f2
#ifdef OPENMP
        integer                                :: nproc, np, resi, myid, jds, jde
        real(sp), dimension(:, :), allocatable :: ftmp

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

        call sgemm('n','n',nx,np,ny,1.0_sp,      f,nx,self%ey(:,jds:jde),ny,0.0_sp,f1,nx)
        call sgemm('t','n',nx,np,nx,1.0_sp,self%ex,nx,                f1,nx,0.0_sp,f2,nx)

        do concurrent (i=1:nx,j=jds:jde, abs(self%evx(i) + self%evy(j)) > eps)
            f2(i,j) = f2(i,j) / (self%evx(i) + self%evy(j))
        end do

        call sgemm('n','n',nx,np,nx,1.0_sp,self%ex,nx,                f2,nx,0.0_sp,ftmp(:,jds:jde),nx)
        !$omp barrier
        call sgemm('n','t',nx,np,ny,1.0_sp,   ftmp,nx,self%ey(jds:jde,:),np,0.0_sp,   p(:,jds:jde),nx)

        deallocate(f1, f2)

        p(:,jds:jde) = p(:,jds:jde) * 10 ** (2 * self % order)

        !$omp end parallel

        deallocate(ftmp)
#else
        allocate(f1(nx, ny), f2(nx, ny))

        call sgemm('n','n',nx,ny,ny,1.0_sp,      f,nx,self%ey,ny,0.0_sp,f1,nx)
        call sgemm('t','n',nx,ny,nx,1.0_sp,self%ex,nx,     f1,nx,0.0_sp,f2,nx)

        do concurrent (j=1:ny,i=1:nx, abs(self%evx(i) + self%evy(j)) > eps)
            f2(i,j) = f2(i,j) / (self%evx(i) + self%evy(j))
        end do

        call sgemm('n','n',nx,ny,nx,1.0_sp,self%ex,nx,     f2,nx,0.0_sp,f1,nx)
        call sgemm('n','t',nx,ny,ny,1.0_sp,     f1,nx,self%ey,ny,0.0_sp, p,nx)

        deallocate(f1, f2)      

        p = p * 10 ** (2 * self % order)
#endif
    end subroutine poisson_2d_solve_sp

    subroutine poisson_2d_destroy_sp(self)
        implicit none
        class(poisson_solver_2d_sp),  intent(inout)  :: self

        deallocate( self % evx, &
                    self % evy, &
                    self %  ex, &
                    self %  ey )
    end subroutine poisson_2d_destroy_sp


    !3d section double precision

    subroutine poisson_3d_init_sp(self, np, step, bdy)
#ifdef OPENMP
        use omp_lib
#endif
        implicit none
        class(poisson_solver_3d_sp), intent(inout)  :: self
        integer,  dimension(3),      intent(in)     :: np
        real(sp), dimension(3),      intent(in)     :: step
        integer,  dimension(6),      intent(in)     :: bdy  !boundary condition

        !local
        real(sp), dimension(3)                      :: tmp

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

        call eigen_sp(np(1), tmp(1), self%evx, self%ex, bdy(1), bdy(2))
        call eigen_sp(np(2), tmp(2), self%evy, self%ey, bdy(3), bdy(4))
        call eigen_sp(np(3), tmp(3), self%evz, self%ez, bdy(5), bdy(6))
    end subroutine poisson_3d_init_sp

    subroutine poisson_3d_solve_sp(self, f, p)
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
        class(poisson_solver_3d_sp),  intent(in)  :: self
        real(sp), dimension(:, :, :), intent(in)  :: f
        real(sp), dimension(:, :, :), intent(out) :: p       


        !local variables
        integer                                   :: i, j, k
        integer                                   :: nx, ny, nz, nxny
        integer                                   :: nxnz, nynz
        real(sp), dimension(:, :, :), allocatable :: f1, f2, f3, f4
#ifdef OPENMP
        integer                                   :: nproc, np, resi, myid, kds, kde
        real(sp), dimension(:, :, :), allocatable :: ftmp
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

        call sgemm('n','n',nxny,np,nz,1.0_sp,      f,nxny,self%ez(:,kds:kde),nz,0.0_sp,f1,nxny)
        call sgemm('t','n',nx,nynz,nx,1.0_sp,self%ex,  nx,                f1,nx,0.0_sp,f2,  nx)

        do concurrent (k=kds:kde)
            f3(:,:,k) = transpose(f2(:,:,k))
        end do

        call sgemm('t','n',ny,nxnz,ny,1.0_sp,self%ey,  ny,     f3,ny,0.0_sp,f4,  ny)

        do concurrent (j=1:ny,i=1:nx,k=kds:kde, abs(self%evx(i) + self%evy(j) + self%evz(k)) > eps)
            f4(j,i,k) = f4(j,i,k) / (self%evx(i) + self%evy(j) + self%evz(k))
        end do

        call sgemm('n','n',ny,nxnz,ny,1.0_sp,self%ey,  ny,     f4,ny,0.0_sp,f3,  ny)

        do concurrent (k=kds:kde)
            f2(:,:,k) = transpose(f3(:,:,k))
        end do

        call sgemm('n','n',nx,nynz,nx,1.0_sp,self%ex,  nx,                f2,nx,0.0_sp, ftmp(:,:,kds:kde),  nx)
        !$omp barrier
        call sgemm('n','t',nxny,np,nz,1.0_sp,   ftmp,nxny,self%ez(kds:kde,:),np,0.0_sp,    p(:,:,kds:kde),nxny)

        deallocate(f1, f2, f3, f4)

        p(:,:,kds:kde) = p(:,:,kds:kde) * 10 ** (2 * self % order)

        !$omp end parallel

        deallocate(ftmp)
#else
        !dimensions are Ok, then allocate work arrays
        allocate( f1(nx, ny, nz), f2(nx, ny, nz), &
                  f3(ny, nx, nz), f4(ny, nx, nz) )

        call sgemm('n','n',nxny,nz,nz,1.0_sp,      f,nxny,self%ez,nz,0.0_sp,f1,nxny)
        call sgemm('t','n',nx,nynz,nx,1.0_sp,self%ex,  nx,     f1,nx,0.0_sp,f2,  nx)

        do concurrent (k=1:nz)
            f3(:,:,k) = transpose(f2(:,:,k))
        end do

        call sgemm('t','n',ny,nxnz,ny,1.0_sp,self%ey,  ny,     f3,ny,0.0_sp,f4,  ny)

        do concurrent (j=1:ny,i=1:nx,k=1:nz, abs(self%evx(i) + self%evy(j) + self%evz(k)) > eps)
            f4(j,i,k) = f4(j,i,k) / (self%evx(i) + self%evy(j) + self%evz(k))
        end do

        call sgemm('n','n',ny,nxnz,ny,1.0_sp,self%ey,  ny,     f4,ny,0.0_sp,f3,  ny)

        do concurrent (k=1:nz)
            f2(:,:,k) = transpose(f3(:,:,k))
        end do

        call sgemm('n','n',nx,nynz,nx,1.0_sp,self%ex,  nx,     f2,nx,0.0_sp,f1,  nx)
        call sgemm('n','t',nxny,nz,nz,1.0_sp,     f1,nxny,self%ez,nz,0.0_sp, p,nxny)

        deallocate(f1, f2, f3, f4)

        p = p * 10 ** (2 * self % order)
#endif
    end subroutine poisson_3d_solve_sp

    subroutine poisson_3d_destroy_sp(self)
        implicit none
        class(poisson_solver_3d_sp),  intent(inout)  :: self
        deallocate( self % evx, &
                    self % evy, &
                    self % evz, &
                    self %  ex, &
                    self %  ey, &
                    self %  ez )
    end subroutine poisson_3d_destroy_sp

    subroutine eigen_sp(n, dx, ev, ee, first, last)
        implicit none
        integer,                   intent(in)  :: n, first, last
        real(sp),                  intent(in)  :: dx
        real(sp), dimension(n),    intent(out) :: ev
        real(sp), dimension(n, n), intent(out) :: ee

        !local variables
        integer                             :: info
        real(sp)                            :: ddxx
        real(sp), dimension(n)              :: d
        real(sp), dimension(n-1)            :: e
        real(sp), dimension(:), allocatable :: work
        integer,  dimension(:), allocatable :: iwork
        integer,  dimension(2*n)            :: isuppz
        integer                             :: lwork, liwork, m


        ddxx   = 1.0_sp / (dx*dx)

        !diagonal
        select case(first)
        case (dirichlet)
            d(1)         = -2.0_sp * ddxx
        case (neumann)
            d(1)         = -ddxx
        case default
            stop "Unknown boundary condition, please select 'dirichlet' or 'neumann'."
        end select

        select case(last)
        case (dirichlet)
            d(n)      = -2.0_sp * ddxx
        case (neumann)
            d(n)      = -ddxx
        case default
            stop "Unknown boundary condition, please select 'dirichlet' or 'neumann'."
        end select

        d(2:n-1)  = -2.0_sp * ddxx

        !subdiagonal
        e = ddxx

        !========================
        ! call lapack subroutine
        !========================
        !first calling is getting the optimal workspace
        lwork  = -1
        liwork = -1
        allocate(work(1), iwork(1))
        call sstevr('v','a',n,d,e,0.0_sp,0.0_sp,0,0,0.0_sp,m,ev,ee,n,isuppz,work,lwork,iwork,liwork,info)

        lwork  = int(work(1))
        liwork = iwork(1)
        deallocate(work, iwork)

        !second calling is calculating eigenvalues and eigenvectors
        allocate(work(lwork), iwork(liwork))
        call sstevr('v','a',n,d,e,0.0_sp,0.0_sp,0,0,0.0_sp,m,ev,ee,n,isuppz,work,lwork,iwork,liwork,info)
    end subroutine eigen_sp

end module poisson_direct_sp

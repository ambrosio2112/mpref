!Modelo de preferencias
!para el ECC de Cerveceria
!
! Usa un algoritmo gibbs - metropolis
! para estimar los coeficientes del modelo
!

module ModPreferencias

  use randomlib
  implicit none
  integer N ! Numero total de individuos
  integer nbeta ! Numero de coeficientes
  integer nmarcas ! Numero de marcas en la escogencia (sin incluir otras)
  integer nclass ! Numero de conjuntos de beta segun la clasificacion del individuo

  real(8), allocatable :: p(:,:) ! probabilidad de seleccion de las marcas N x nmarcas+1
  real(8), allocatable :: Xe(:,:,:) ! evaluaciones de las marcas N x nmarcas x nbeta
  integer, allocatable :: y(:) ! preferecias dim=N
  integer, allocatable :: missing(:) ! valores perdidos dim=N
  integer, allocatable :: class(:) ! tipo de individuo dim=N (1...a nclass)
  real(8), allocatable :: beta( :, :, : ) !betas dimensionados por nclass x nmarcas x nbeta
  real(8), allocatable :: Sigma(:,:), chSigma(:,:), diagSigma(:) ! Matriz de covarianzas de los betas

  integer nmissing

  real(8), allocatable :: S(:,:) ! matriz de covarianza previas
  real(8) nu ! grados de libertad de Sigma

  real(8), allocatable :: mp(:,:)
  real(8), allocatable :: mSigma(:,:)
  real(8), allocatable :: mbeta( :, :, : )

contains

  !---------------------------------------------------------------------------------------------------
  subroutine genSigma()
    real(8) Sx(nbeta,nbeta), diagSx(nbeta), chS(nbeta,nbeta), dS(nbeta), Sg(nbeta,nbeta)
    real(8) nux

    logical warn
    integer i,j

    nux = nclass*nmarcas + nu
    call rIWishart( nux, chS, dS) !comes back with chol decomp of Sigma ~ IW(df,I)

    Sx = S
    do i=1,nclass
       do j=1,nmarcas
          Sx = Sx + matmul( reshape(beta(i,j,:), (/nbeta,1/)), reshape(beta(i,j,:),(/1,nbeta/)) )
       end do
    end do


    !rescale
    warn=.false.
    call chol(Sx,diagSx,warn)
    if(warn) then
       return
    end if
    ! Sx x Sigma tiene cholesky dada por la multiplicacion chSx x chSigma
    do i=1,nbeta
       diagSigma(i) = diagSx(i) * dS(i)
       do j=1,i-1
          chSigma(i,j) = dot_product(Sx(i, (j+1):(i-1) ), chS( (j+1):(i-1) ,j)) + &
               diagSx(i)*chS(i ,j) + Sx(i,j)*dS(j)
       end do
    end do

    Sigma=chSigma
    call recoverchol(Sigma, diagSigma) !solo recupera la parte de arriba
    do i=1,nbeta-1 ! symmetrize
       Sigma((i+1):nbeta,i) = Sigma(i, (i+1):nbeta )
    end do

  end subroutine genSigma

  !---------------------------------------------------------------------------------------------------
  subroutine genBeta( nc, m )
    integer nc ! clase de la beta a generar
    integer m ! marca a generar

    real(8) :: b(nbeta)
    real(8) :: alpha
    real(8) :: x1(nbeta), x2(nbeta)
    real(8) :: q(N,nmarcas+1), sq, sp

    integer i,j

    b = rmvnorm( beta(nc,m,:), 2*chSigma, 2*diagSigma)
    !b = rnorm(nbeta)/sqrt(diagSigma) + beta(nc,m,:)

    ! calcula probabilidades con las nuevas betas
    sq=0
    sp=0
    do i=1,N
       if (class(i) == nc) then
          do j=1,nmarcas
             if (j/=m) then
                q(i,j) = exp( sum( beta(class(i),j,:) * Xe(i,j,:)) )
             else
                q(i,j) = exp( sum( b * Xe(i,j,:)) )
             end if
          end do
          q(i,nmarcas+1) = 1.
          q(i,:) = q(i,:)/sum(q(i,:))
          sq = sq + log(q(i,y(i)))
          sp = sp + log(p(i,y(i)))
       end if
    end do


    ! aceptance
    call SolveChol( chSigma,diagSigma, b, x1)
    call SolveChol( chSigma,diagSigma, beta(nc,m,:),x2)

    alpha =  -0.5 * sum( b*x1 ) + 0.5 * sum( beta(nc,m,:)*x2 ) + sq - sp

    if ( log(runif()) < alpha ) then
       ! propuesta aceptada
       beta(nc,m,:) = b
       do i=1,N
          if (class(i)==nc) then
             p(i,:) = q(i,:)
          end if
       end do
    end if

  end subroutine genBeta

  !---------------------------------------------------------------------------------------------------
  subroutine genMissing()
    integer i
    do i=1,N
       if( missing(i)==1 ) y(i) = sample( p(i,:) )
    end do
  end subroutine genMissing

  !---------------------------------------------------------------------------------------------------
  subroutine saveState(is)
    integer i,j,is
    open(unit=10,file="betas.txt",access='APPEND')
    do i=1,nclass
       do j=1,nmarcas
          write(10,*) beta(i,j,:)
       end do
    end do
    close(10)
    open(unit=10,file="Sigmas.txt",access='APPEND')
    do i=1,nbeta
       write(10,*)Sigma(i,:)
    end do
    close(10)
    open(unit=10,file="prefs.txt",access='APPEND')
    do i=1,N
       write(10,*) p(i,:)
    end do
    close(10)
    open(unit=10,file="status.txt",access='APPEND')
    write(10,*) is
    close(10)

    mbeta = mbeta + beta
    mp = mp + p
    msigma = msigma + sigma
  end subroutine saveState

  !---------------------------------------------------------------------------------------------------

  subroutine saveStats()
    integer i,j
    open(unit=10,file="betas.stats.txt",access='APPEND')
    do i=1,nclass
       do j=1,nmarcas
          write(10,*) mbeta(i,j,:)
       end do
    end do
    close(10)
    open(unit=10,file="Sigmas.stats.txt",access='APPEND')
    do i=1,nbeta
       write(10,*) mSigma(i,:)
    end do
    close(10)
    open(unit=10,file="prefs.stats.txt",access='APPEND')
    do i=1,N
       write(10,*) mp(i,:)
    end do
    close(10)
  end subroutine saveStats

  !---------------------------------------------------------------------------------------------------
  subroutine iterar( niter, nburn, nstep )
    integer niter,nburn,nstep ! numero de iteraciones a realizar
    integer it, i,j, is

    is = 0
    do it=1,niter
       if (nmissing>0) then
          call genMissing()
       end if

       do i=1,nclass
          do j=1,nmarcas
             call genBeta(i,j)
          end do
       end do

       call genSigma()
       if( it >= nburn .and. mod(it-nburn,nstep)==0) then
          is = is + 1
          call saveState(is)
       end if
    end do

    mbeta = mbeta/is
    mp = mp/is
    msigma = msigma/is
    call saveStats()
  end subroutine iterar

end module ModPreferencias


subroutine preferencias(niter,nburn,nstep, N_,nbeta_,nmarcas_,nclass_,Xe_,y_, missing_, class_, S_,nu_ )
!dec$attributes dllexport :: preferencias
use ModPreferencias
implicit none
integer niter,nburn,nstep
integer N_ ! Numero total de individuos
integer nbeta_ ! Numero de coeficientes
integer nmarcas_ ! Numero de marcas en la escogencia (sin incluir otras)
integer nclass_ ! Numero de conjuntos de beta segun la clasificacion del individuo

real(8) Xe_(N_,nmarcas_,nbeta_) ! evaluaciones de las marcas N x nmarcas x nbeta
integer y_(N_) ! preferecias dim=N
integer missing_(N_) ! valores perdidos (0 si no, 1 si esta perdido)
integer class_(N_) ! tipo de individuo dim=N (1...a nclass)
real(8) S_(nbeta_,nbeta_) ! matriz de covarianza previas
real(8) nu_ ! grados de libertad de Sigma
logical warn

integer i,j
real(8) :: zero(nbeta_)

	N=N_
	nbeta=nbeta_
	nmarcas=nmarcas_
	nclass=nclass_
	nu=nu_
	allocate (p(N,nmarcas+1),Xe(N,nmarcas,nbeta), y(N), missing(N), class(N), beta(nclass,nmarcas,nbeta), &
		      Sigma(nbeta,nbeta), chSigma(nbeta,nbeta), diagSigma(nbeta), S(nbeta,nbeta) )
	Xe=Xe_
	y=y_
	class=class_
	missing=missing_

	nmissing=sum(missing)


	!inicializa Sigma en S
	S=nu*S_
	Sigma=S
	chSigma=Sigma
	warn=.false.
	call chol(chSigma,diagSigma,warn)

	! incializa los beta en 0
	zero=0
	do i=1,nclass
		do j=1,nmarcas
			beta(i,j,:)= rmvnorm( zero , chSigma,diagSigma)
		end do
	end do

	!calcula las probabilidades
	do i=1,N
		do j=1,nmarcas
			p(i,j) = exp( sum( beta(class(i),j, :)*Xe(i,j,:) ) )
		end do
		p(i,nmarcas+1)=1.
		p(i,:) = p(i,:)/sum(p(i,:))
	end do

	!inicializa las medias
	allocate (mp(N,nmarcas+1),mbeta(nclass,nmarcas,nbeta),mSigma(nbeta,nbeta))
	mp=0
	mbeta=0
	mSigma=0
	call iterar(niter,nburn,nstep)

	deallocate (p,Xe, y, missing, class, beta, &
		      Sigma, chSigma, diagSigma, S )
	deallocate (mp,mbeta,mSigma)

end subroutine preferencias

subroutine preferencias0(niter,nburn,nstep, N_,nbeta_,nmarcas_,nclass_,Xe_,y_, missing_, class_, S_,nu_, &
                         beta_,Sigma_)
!dec$attributes dllexport :: preferencias0
use ModPreferencias
implicit none
integer niter,nburn,nstep
integer N_ ! Numero total de individuos
integer nbeta_ ! Numero de coeficientes
integer nmarcas_ ! Numero de marcas en la escogencia (sin incluir otras)
integer nclass_ ! Numero de conjuntos de beta segun la clasificacion del individuo

real(8) Xe_(N_,nmarcas_,nbeta_) ! evaluaciones de las marcas N x nmarcas x nbeta
integer y_(N_) ! preferecias dim=N
integer missing_(N_) ! valores perdidos (0 si no, 1 si esta perdido)
integer class_(N_) ! tipo de individuo dim=N (1...a nclass)
real(8) S_(nbeta_,nbeta_) ! matriz de covarianza previas
real(8) nu_ ! grados de libertad de Sigma

real(8) beta_(nclass_,nmarcas_,nbeta_), Sigma_(nbeta_,nbeta_)
logical warn

integer i,j

	N=N_
	nbeta=nbeta_
	nmarcas=nmarcas_
	nclass=nclass_
	nu=nu_
	allocate (p(N,nmarcas+1),Xe(N,nmarcas,nbeta), y(N), missing(N), class(N), beta(nclass,nmarcas,nbeta), &
		      Sigma(nbeta,nbeta), chSigma(nbeta,nbeta), diagSigma(nbeta), S(nbeta,nbeta) )
	Xe=Xe_
	y=y_
	class=class_
	missing=missing_

	nmissing=sum(missing)

	! incializa los beta en 0
	beta=beta_

	!inicializa Sigma en S
	S=nu*S_
	Sigma=Sigma_
	chSigma=Sigma
	warn=.false.
	call chol(chSigma,diagSigma,warn)

	!calcula las probabilidades
	do i=1,N
		do j=1,nmarcas
			p(i,j) = exp( sum( beta(class(i),j, :)*Xe(i,j,:) ) )
		end do
		p(i,nmarcas+1)=1.
		p(i,:) = p(i,:)/sum(p(i,:))
	end do

	!inicializa las medias
	allocate (mp(N,nmarcas+1),mbeta(nclass,nmarcas,nbeta),mSigma(nbeta,nbeta))
	mp=0
	mbeta=0
	mSigma=0
	call iterar(niter,nburn,nstep)

	deallocate (p,Xe, y, missing, class, beta, &
		      Sigma, chSigma, diagSigma, S )
	deallocate (mp,mbeta,mSigma)

end subroutine preferencias0


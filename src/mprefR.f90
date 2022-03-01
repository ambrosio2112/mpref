! Modelo de Preferencia:
! Modela la preferencia hacia una marca en particular
! versus otras marcas en funcion de las evaluaciones dadas por cada individuo
! para las marcas (no necesariamente evaluando las mismas marcas)
!
! logit(p_{i,M_{ij}}) = beta_i X_{i,M_{ij}}
! donde p es la probabilidad de preferir la marca particular contra la marca
! M_{ij} evaluada por el individuo i. beta_i es el vector de coeficientes,
! modelados a priori independientemente como una gamma(nu, lambda_k) para
! cada atributo k (un coeficiente por atributo). X contiene la diferencia
! entre las evaluaciones de la marca particular y la marca M_{ij} para
! cada una de los K atributos. lambda a su vez es modelado como una
! iid igamma (alpha, lambda0)
!
! El codigo suministra el log de la posterior (excepto por constantes)
! y el gradiente del mismo para los coeficientes beta
! y evalua el maximo para los hiperparametros lambda de los coeficientes


! mpreflpost devuelve el log de la posterior (menos constantes)
! N=numero de individuos
! K=numero de atributos evaluados
! NM= numero de marcas evaluadas por cada individuo
! maxNM: total de marcas evaluadas por los individuos
! M: Otras marcas evaluadas por los individuos
! y: 1 si la marca particular es preferida contra la otra marca, 0 si no
! X: diferencia entre los atributos entre la marca particular y cada una de las otras marcas

subroutine mpreflpost( res, N,K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
!dec$attributes dllexport:: mpreflpost
implicit none
real(8) res
integer N,K,NM(N), maxNM, M(N,maxNM), y(N,maxNM+1)
real(8) X(N,maxNM+1,K), beta(N,K), lambda(K), nu, alpha, lambda0

integer i,j
real(8) xb, sl, sm

	res = 0
	sl = sum(log(lambda(1:K)))
	do i=1,N
		do j=1, NM(i)
			sm = min(500.,max(-500.,sum( beta(i,1:K)*X(i,M(i,j),1:K) )))
			xb = exp(sm)
			if( y(i,M(i,j)) == 1 ) then
				res = res + log( xb/(1.+xb) )
			else
				res = res - log( 1.+ xb )
			end if
		end do
		res = res + sum( (nu-1.)*log(beta(i,1:K)) - beta(i,1:K)/lambda(1:K)) - nu*sl
	end do
	res = res - (alpha+1.)*sl - 1./sum(lambda(1:K))/lambda0
end subroutine mpreflpost


! esta funcion es auxiliar calcula el log post excluyendo i0
subroutine mpreflpostexi( res, i0, N,K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
!dec$attributes dllexport:: mpreflpostexi
implicit none
real(8) res
integer N,K,NM(N), maxNM, M(N,maxNM), y(N,maxNM+1)
real(8) X(N,maxNM+1,K), beta(N,K), lambda(K), nu, alpha, lambda0
integer i0

integer i,j
real(8) xb, sl, sm

	res = 0
	sl = sum(log(lambda(1:K)))
	do i=1,N
		if (i==i0) continue
		do j=1, NM(i)
			sm = min(500.,max(-500.,sum( beta(i,1:K)*X(i,M(i,j),1:K) )))
			xb = exp(sm)
			if( y(i,M(i,j)) == 1 ) then
				res = res + log( xb/(1.+xb) )
			else
				res = res - log( 1.+ xb )
			end if
		end do
		res = res + sum( (nu-1.)*log(beta(i,1:K)) - beta(i,1:K)/lambda(1:K)) - nu*sl
	end do
	res = res - (alpha+1.)*sl - 1./sum(lambda(1:K))/lambda0
end subroutine mpreflpostexi


!esta es auxiliar, calcula log post a partir de res0 mas la parte de i
subroutine mpreflposti( res, res0, i, N,K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
!dec$attributes dllexport:: mpreflposti
implicit none
real(8) res,res0
integer N,K,NM(N), maxNM, M(N,maxNM), y(N,maxNM+1)
real(8) X(N,maxNM+1,K), beta(N,K), lambda(K), nu, alpha, lambda0
integer i

integer j
real(8) xb, sm

	res = res0

	do j=1, NM(i)
		sm = min(500.,max(-500.,sum( beta(i,1:K)*X(i,M(i,j),1:K) )))
		xb = exp(sm)
		if( y(i,M(i,j)) == 1 ) then
			res = res + log( xb/(1.+xb) )
		else
			res = res - log( 1.+ xb )
		end if
	end do
	res = res + sum( (nu-1.)*log(beta(i,1:K)) - beta(i,1:K)/lambda(1:K))

end subroutine mpreflposti



! mpreflambda devuelve los MAP de lambda dados
! los beta.
subroutine mpreflambda( lambda, N,K,beta,nu,alpha,lambda0 )
!dec$attributes dllexport:: mpreflambda
implicit none
integer N,K
real(8) lambda(K), beta(N,K), nu, alpha, lambda0

	lambda = ( reshape(sum(beta,1),(/K/)) + 1./lambda0 )/(N*nu + alpha+1.)

end subroutine mpreflambda


!mprefgrlpostbetai devuelve el vector gradiente para los coeficientes
! beta_i para un i fijo dados los lambda.

subroutine mprefgrlpostbetai( res, i, N,K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
!dec$attributes dllexport:: mprefgrlpostbetai
implicit none
integer i, N,K,NM(N), maxNM, M(N,maxNM), y(N,maxNM+1)
real(8) X(N,maxNM+1,K), beta(N,K), lambda(K), nu, alpha, lambda0
real(8) res(K)

integer j
real(8) xb,sm

	res= (nu-1.)/beta(i,1:K) - 1./lambda(1:K)
	do j=1,NM(i)
		sm=min(500.,max(-500.,sum( beta(i,1:K)*X(i,M(i,j),1:K) ) ))
		xb = exp(sm)
		if( y(i,M(i,j)) == 1 ) then
			res = res + X(i,M(i,j),1:K)/(1.+ xb)
		else
			res = res - X(i,M(i,j),1:K)*xb/(1.+xb)
		end if
	end do

end subroutine mprefgrlpostbetai

!======================================================================================
!!$$ module imsl
!!$$   interface
!!$$      subroutine dbcong (fcn, grad, n, xguess, ibtype, xlb,&
!!$$           xub, xscale, fscale, iparam, rparam, x, fvalue)
!!$$ !dec$attributes stdcall, alias:'_DBCONG@52' :: dbcong
!!$$ !dec$attributes reference :: fcn, grad, n, xguess, ibtype, xlb, xub, xscale, fscale, iparam, rparam, x, fvalue
!!$$        integer    n, ibtype, iparam(*)
!!$$        double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
!!$$             xscale(*), rparam(*), x(*)
!!$$        external   fcn, grad
!!$$      end subroutine dbcong

!!$$      subroutine umach(id,un)
!!$$ !dec$attributes stdcall, alias:'_UMACH@8' ::umach
!!$$ !dec$attributes reference :: id,un
!!$$        integer id,un
!!$$      end subroutine umach

!!$$      subroutine du4inf(ip,rp)
!!$$ !dec$attributes stdcall, alias:'_DU4INF@8' :: du4inf
!!$$ !dec$attributes reference :: ip,rp
!!$$        integer ip(*)
!!$$        real(8) rp(*)
!!$$      end subroutine du4inf

!!$$      subroutine dlslsf (n, a, lda, b, x)
!!$$ !dec$attributes stdcall, alias:'_DLSLSF@20' :: dlslsf
!!$$ !dec$attributes reference :: n, a, lda, b, x
!!$$        integer    n, lda
!!$$        double precision a(lda,*), b(*), x(*)
!!$$      end subroutine dlslsf

!!$$ 	 subroutine erset(iersvr,ipact,isact)
!!$$ !dec$attributes stdcall, alias:'_ERSET@12' :: erset
!!$$ !dec$attributes reference :: iersvr,ipact,isact
!!$$ 	   integer iersvr,ipact,isact
!!$$ 	 end subroutine erset
!!$$   end interface
!!$$ end module imsl

!======================================================================================

module modmprefoptim

integer, allocatable :: NM(:),  M(:,:), y(:,:)
real(8) res0,  nu, alpha, lambda0
real(8), allocatable ::  X(:,:,:), lambda(:)
integer i0

end module

!======================================================================================

subroutine mprefoptimFCN (K, nbeta, F)
!dec$attributes stdcall :: mprefoptimFCN
!dec$attributes reference :: K, nbeta, F
use modmprefoptim
implicit none

integer K
real(8) nbeta(K)
real(8) F

integer j
real(8) xb, sm

	F = res0

	do j=1, NM(i0)
		sm = min(500.,max(-500.,sum( nbeta(1:K)*X(i0,M(i0,j),1:K) )))
		xb = exp(sm)
		if( y(i0,M(i0,j)) == 1 ) then
			F = F + log( xb/(1.+xb) )
		else
			F = F - log( 1.+ xb )
		end if
	end do
	F = F + sum( (nu-1.)*log(nbeta(1:K)) - nbeta(1:K)/lambda(1:K))

	F=-F
end subroutine mprefoptimFCN


subroutine mprefoptimGR (K, nbeta, G)
!dec$attributes stdcall :: mprefoptimGR
!dec$attributes reference :: K, nbeta, G
use modmprefoptim
implicit none

integer K
real(8) nbeta(K)
real(8) G(K)

integer j
real(8) xb,sm

	G = (nu-1.)/nbeta(1:K) - 1./lambda(1:K)
	do j=1,NM(i0)
		sm=min(500.,max(-500.,sum( nbeta(1:K)*X(i0,M(i0,j),1:K) ) ))
		xb = exp(sm)
		if( y(i0,M(i0,j)) == 1 ) then
			G = G + X(i0,M(i0,j),1:K)/(1.+ xb)
		else
			G = G - X(i0,M(i0,j),1:K)*xb/(1.+xb)
		end if
	end do
	G=-G
end subroutine mprefoptimGR


subroutine mprefmapbeta(N, K, oNM, maxNM, oM, oy, oX, beta, olambda, onu, oalpha, olambda0)
!dec$attributes dllexport :: mprefmapbeta
  use modmprefoptim
!!  use imsl
  implicit none

  integer N,K,oNM(N), maxNM, oM(N,maxNM), oy(N,maxNM+1)
  real(8) oX(N,maxNM+1,K), beta(N,K), olambda(K), onu, oalpha, olambda0

  integer i

  real(8)  xlb(K), xub(K),  factr, pgtol
  integer ifail, nbd(K), fncount, grcount, maxit
  real(8) rbeta(K), fvalue

  interface
     subroutine mprefoptimFCN (K, nbeta, F)
!dec$attributes stdcall :: mprefoptimFCN
!dec$attributes reference :: K, nbeta, F
       integer K
       real(8) nbeta(K)
       real(8) F
     end subroutine mprefoptimFCN

     subroutine mprefoptimGR (K, nbeta, G)
!dec$attributes stdcall :: mprefoptimGR
!dec$attributes reference :: K, nbeta, G
       integer K
       real(8) nbeta(K)
       real(8) G(K)
     end subroutine mprefoptimGR
  end interface


  allocate(NM(N), M(N,maxNM), y(N,maxNM+1), X(N,maxNM+1,K), lambda(K) )
  NM=oNM
  M=oM
  y=oy
  X=oX
  lambda=olambda
  nu=onu
  alpha=oalpha
  lambda0=olambda0

  xlb=0.0001
  xub=10000.0
  nbd=2
  factr=1.e+7
  pgtol=0.
  maxit=100
  fncount=1000
  grcount=1000
  fvalue=0.

!!$$  CALL UMACH (-3, 9)
!!$$  CALL UMACH (-2, 9)
!!$$  OPEN (UNIT=9,FILE='CHECKERR')
!!$$  call ERSET(0,1,0)
!!$$  call DU4INF(iparam,rparam)
!!$$  iparam(3)=100 ! Maximum number of iterations.
!!$$  iparam(4)=1000 ! Maximum number of function evaluations.
!!$$  iparam(5)=1000 ! Maximum number of gradient evaluations.

  do i=1,N
     i0=i
     call mpreflpostexi( res0, i, N,K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)

     rbeta=beta(i,1:K)

     call lbfgsb_f ( &
          K,& ! n=K n is the number of parameters
          5, & ! lmm = 5 is an integer giving the number of BFGS updates retained in
               ! the "L-BFGS-B" method, It defaults to 5
          rbeta, & ! x =  is the starting parameters on entry and x the final parameters on exit,
          XLB, XUB, & ! bounds
          nbd, & ! type of bound per variable (0 unbounded, 1, low, 2 both,3 upper
          FVALUE, & ! Fmin = function value at minimum
          mprefoptimFCN, & !fn
          mprefoptimGR, & ! grad
          ifail, & ! fail status
          0, & ! ex pointer to additional data passed to functions
          factr, & ! controls the convergence of the "L-BFGS-B" method.
          pgtol, & ! helps control the convergence of the "L-BFGS-B" method.
                   !It is a tolerance on the projected gradient in the current
                   !search direction. This defaults to zero
          fncount, grcount, & !  fn and grad max evaluation count
          maxit & ! max iterations
          )
     beta(i,1:K)=rbeta
  end do

  !! close(9)
  deallocate(NM, M, y, X, lambda )
end subroutine mprefmapbeta

!======================================================================

subroutine mprefiterbetai(i0, N, K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
!!$$  use imsl
  implicit none

  interface !! from linpack
     SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
       CHARACTER      UPLO
       INTEGER        INFO, LDA, LWORK, N
       INTEGER        IPIV( * )
       real(8)        A( LDA, * ), WORK( * )
     end subroutine DSYTRF

     SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       CHARACTER      UPLO
       INTEGER        INFO, LDA, LDB, N, NRHS
       INTEGER        IPIV( * )
       real(8)        A( LDA, * ), B( LDB, * )
     end SUBROUTINE DSYTRS
  end interface


  integer N,K,NM(N), maxNM, M(N,maxNM), y(N,maxNM+1)
  real(8) :: X(N,maxNM+1,K), beta(N,K), lambda(K), nu, alpha, lambda0
  integer i0
  real(8) :: H(K,K), G(K), work(K*K)
  integer :: ipiv(K),info

  real(8) :: p(K),sm
  real(8) :: x0(K), beta0(K)
  integer j,ik
  !construye la hessiana

  integer,parameter :: nmaxiter=100
  real(8), parameter ::  tol=1.E-4
  logical done
  integer iter

  beta0 = beta(i0,1:K)
  done=.false.
  iter=0
  do while( .not. done )
     H=0.
     G= (nu-1.)/beta0 - 1./lambda(1:K)
     do j=1,NM(i0)
        sm=min(500.,max(-500.,sum( beta0(1:K)*X(i0,M(i0,j),1:K) ) ))
        p = 1./(1. + exp(-sm))
        H= H - matmul( reshape(p*(1.-p)*X(i0,M(i0,j),1:K), (/K ,1 /)), reshape(X(i0,M(i0,j),1:K), (/1,K/) ) )
        if (y(i0,M(i0,j)) == 1 ) then
           G = G - p*X(i0,M(i0,j),1:K) + X(i0,M(i0,j),1:K)
           do ik=1,K
              H(ik,ik) = H(ik,ik) -  (nu-1.)/beta0(ik)**2
           end do
        else
           G = G - p*X(i0,M(i0,j),1:K)
        end if
     end do

!!$$ call dlslsf( K, H, K, G, x0 ) ! N=K number of eq, lda=dim(A,1), Ax=b A symmetric

     call dsytrf('U',K,H,K,ipiv,work,K*K,info)
     x0=g
     call dsytrs('U',K,1,H,K,ipiv,x0,K,info)

     beta0 = beta0 - x0

     do ik=1,K
        beta0(ik)=max(beta0(ik),.0001)
     end do

     iter=iter+1
     done = (iter>nmaxiter) .or. ( sum( (beta0 - beta(i0,1:K) )**2 )/ sum( beta(i0,1:K)**2 ) < tol )
     beta(i0,1:K)=beta0
  end do

end subroutine mprefiterbetai


subroutine mprefitermapbeta(N, K, NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
!dec$attributes dllexport :: mprefitermapbeta
  !!use imsl
  !!use dfport
  implicit none

  integer N,K, NM(N), maxNM, M(N,maxNM), y(N,maxNM+1)
  real(8) X(N,maxNM+1,K), beta(N,K), lambda(K), nu, alpha, lambda0

  integer i

!!$$  CALL UMACH (-3, 9)
!!$$  CALL UMACH (-2, 9)
!!$$  OPEN (UNIT=9,FILE='CHECKERR')
!!$$  call ERSET(0,1,0)

  do i=1,N
     call mprefiterbetai(i,N,K,NM, maxNM, M, y, X, beta, lambda, nu, alpha, lambda0)
	 call mpreflambda( lambda, N,K,beta,nu,alpha,lambda0 )
	 !!write(9,*) i, beta(i,1:K)
	 !!call flush(9)
  end do
  !!close(9)
end subroutine mprefitermapbeta

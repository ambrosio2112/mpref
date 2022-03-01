#' @title Modelo de Preferencias StatConsult
#'
#' @description
#' Modela la preferencia hacia una marca en particular
#' versus otras marcas en funcion de las evaluaciones dadas por cada individuo
#' para las marcas (no necesariamente evaluando las mismas marcas)
#'
#' logit(p_{i,M_{ij}}) = beta_i X_{i,M_{ij}}
#' donde p es la probabilidad de preferir la marca particular contra la marca
#' M_{ij} evaluada por el individuo i. beta_i es el vector de coeficientes,
#' modelados a priori independientemente como una gamma(nu, lambda_k) para
#' cada atributo k (un coeficiente por atributo). X contiene la diferencia
#' entre las evaluaciones de la marca particular y la marca M_{ij} para
#' cada una de los K atributos. lambda a su vez es modelado como una
#' iid igamma (alpha, lambda0)
#'
#' Modela la preferencia hacia una marca en particular
#' versus otras marcas en funcion de las evaluaciones dadas por cada individuo
#' para las marcas (no necesariamente evaluando las mismas marcas)
#'
#' @param N Individuos
#' @param nbeta Cantidad de coeficientes
#' @param nmarcas Cantidad de objetos a preferir
#' @param nclase Cantidad de posibles clases de coeficientes (nclases=N indica uno por individuo)
#' @param clases Mapa entre individuos y clases
#' @param X Covariables que explican la preferencia
#' @param y Vector con las preferidas (1...nmarcas, 0 para ninguna, NA para valor perdido)
#' @param missing Indicador con cuales observaciones faltan
#' @param niter Iteraciones del MCMC
#' @param nburn Iteraciones a quemar por el MCMC
#' @param nstep Pasos de iteraciones para tomar una muestra
#' @param nu Parametros de la previa
#' @param S Parametros de la previa
#'
#' @return La funcion devuelve sus resultados en unos archivos de texto
#' @export
#' @useDynLib mpref
#' @examples
#' ## buscar los ejemplos de los proyectos de modelos de marca
#'
preferencias <- function(
                         N,nbeta,nmarcas,
                         nclase,clases,
                         X,y,missing,
                         niter=15000, nburn=100,nstep=10,
                         nu=NULL,S=NULL)
  {
    if(!is.loaded("preferencias_"))
	{
      		dyn.load("/usr/local/lib64/libmpref.so")
	}

    if(is.null(nu)) nu <- as.double(nbeta + 0.5)
    if(is.null(S)) S <- diag(nbeta)*1

    z <- .Fortran("preferencias",
                  niter=as.integer(niter),
                  nburn=as.integer(nburn),
                  nstep=as.integer(nstep),
                  N=as.integer(N),
                  nbeta=as.integer(nbeta),
                  nmarcas=as.integer(nmarcas), #numero de marcas sin otras
                  nclass=as.integer(nclase),
                  Xe=as.double(X),
                  y=as.integer(y),
                  missing=as.integer(missing),
                  class=as.integer(clases),
                  S=as.double(S),
                  nu=as.double(nu),
                  NAOK=TRUE
                  )
    NULL
  }


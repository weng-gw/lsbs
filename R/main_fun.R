#' Calculate the optiaml bandwidth matrix for highest density region
#' estimation
#'
#' This function allow you to compute the optiaml bandwidth matrix
#' for highest density region estimation by using a plug-in strategy.
#' @param X a matrix with two columns containing the data from the
#' density function.
#' @param tau a probability value between 0 and 1
#' @param xrange a vector of of length 2, e.g., \code{c(xmin, xmax)},
#' indicating the range the grid points to be generated on x-axis
#' @param yrange a vector of of length 2, e.g., \code{c(ymin, ymax)},
#' indicating the range the grid points to be generated on y-axis
#' @param gridwidth width between grid points.
#' @param init starting value of the bandwidth matrix for
#' optimization. If not specified, use direct-plug estimator from
#' \code{ks} package as starting value
#' @param maxit maximum number of iterations for optimization
#' @param tol  tolerance value for stopping the optimization algorithm
#' @param print_obj a flag (boolean type) indicates printing the loss
#' function values during optimizatin or not.
#' @return the optimal bandwidth matrix.
#' @references  Doss, C.R. and Weng, G., 2018. \emph{Bandwidth selection for
#' kernel density estimators of multivariate level sets and highest
#' density regions}. arXiv preprint arXiv:1806.00731.
#' @examples 
#' X <- matrix(rnorm(100),ncol=2)
#' xrange <- c(-2.5,2.5)
#' yrange <- c(-2.5,2.5)
#' hdrbs(X,0.1,xrange,yrange,0.1)
#' 
#' @export


hdrbs <- function(X, tau, xrange, yrange, gridwidth, init=NULL,
                  maxit=200,tol=1e-6,print_obj=FALSE){
    if(is.null(init)){
        init <- ks::Hpi(X,binned=TRUE)
    }
    xxgrid <- seq(xrange[1],xrange[2],gridwidth)
    yygrid <- seq(yrange[1],yrange[2],gridwidth)
    ep <- cbind(rep(xxgrid,each=length(yygrid)),rep(yygrid,length(xxgrid)))
    mykde_fval <- ks::kde(X,eval.points=ep,binned=TRUE)
    zzest_dpi <- matrix(mykde_fval$estimate,nrow=length(xxgrid),byrow=TRUE)
    mykde_hess <- ks::kdde(X, deriv.order=2, eval.points=ep)
    hess_grid_est <-
        list(dxx=matrix(mykde_hess$estimate[,1],nrow=length(xxgrid),byrow=TRUE),
             dxy=matrix(mykde_hess$estimate[,2],nrow=length(xxgrid),byrow=TRUE),
             dyy= matrix(mykde_hess$estimate[,4],nrow=length(xxgrid),byrow=TRUE))
    ftau_dpi <- biftau(zzest_dpi,tau,gridwidth,tol)

    cl_est <- cp_comp(xxgrid,yygrid,zzest_dpi,ftau_dpi)
    contpoint_est <- cl_est$contpoint
    contpiece_est <- cl_est$contpiece

    mykde_grad <- ks::kdde(X,deriv.order=1,eval.points=contpoint_est)
    mykde_betahess <- ks::kdde(X,deriv.order=2,eval.points=contpoint_est)
    fx_est <- mykde_grad$estimate[,1]
    fy_est <- mykde_grad$estimate[,2]
    grad_norm_est <- sqrt(fx_est^2+fy_est^2)
    hess_est <- list(dxx=mykde_betahess$estimate[,1],
                 dxy=mykde_betahess$estimate[,2],
                 dyy=mykde_betahess$estimate[,4])
    h_init <- c(init[1,1],init[1,2],init[2,2])
    optH <- hdr_opt(n=nrow(X), ftau=ftau_dpi,grad_norm=grad_norm_est,
                    hess=hess_est,hess_grid=hess_grid_est,
                    contpiece=contpiece_est, zztrue=zzest_dpi,
                    initial=h_init,gridwidth=gridwidth,
                    maxit=maxit,tol=tol, print_obj = print_obj)
    return(optH)
}

#' Calculate the optiaml bandwidth matrix for level set
#' estimation
#'
#' This function allow you to compute the optiaml bandwidth matrix
#' for level set  estimation by using a plug-in strategy.
#' @param X a matrix with two columns containing the data from the
#' density function.
#' @param levelc a positive value indicating the height of the level
#' set
#' @param xrange a vector of of length 2, e.g., \code{c(xmin, xmax)},
#' indicating the range the grid points to be generated on x-axis
#' @param yrange a vector of of length 2, e.g., \code{c(ymin, ymax)},
#' indicating the range the grid points to be generated on y-axis
#' @param gridwidth width between grid points.
#' @param init starting value of the bandwidth matrix for
#' optimization. If not specified, use direct-plug estimator from
#' \code{ks} package as  starting value
#' @param maxit maximum number of iterations for optimization
#' @param tol  tolerance value for stopping the optimization algorithm
#' @param print_obj a flag (boolean type) indicates printing the loss
#' function values during optimizatin or not.
#' @return the optimal bandwidth matrix.
#' @references  Doss, C.R. and Weng, G., 2018. \emph{Bandwidth selection for
#' kernel density estimators of multivariate level sets and highest
#' density regions}. arXiv preprint arXiv:1806.00731.
#' @examples
#' X <- matrix(rnorm(100),ncol=2)
#' xrange <- c(-3,3)
#' yrange <- c(-3,3)
#' lsbs(X,0.1,xrange,yrange,0.05)
#' @export


lsbs <- function(X, levelc, xrange, yrange, gridwidth,
                 init=NULL, maxit=200,tol=1e-6,print_obj=FALSE){
    if(is.null(init)){
        init <- ks::Hpi(X,binned=TRUE)
    }
    xxgrid <- seq(xrange[1],xrange[2],gridwidth)
    yygrid <- seq(yrange[1],yrange[2],gridwidth)
    ep <- cbind(rep(xxgrid,each=length(yygrid)),rep(yygrid,length(xxgrid)))
    mykde_fval <- ks::kde(X,eval.points=ep,binned=TRUE)
    zzest_dpi <- matrix(mykde_fval$estimate,nrow=length(xxgrid),byrow=TRUE)

    cl_est <- cp_comp(xxgrid,yygrid,zzest_dpi,levelc)
    contpoint_est <- cl_est$contpoint
    contpiece_est <- cl_est$contpiece

    mykde_grad <- ks::kdde(X,deriv.order=1,eval.points=contpoint_est)
    mykde_betahess <- ks::kdde(X,deriv.order=2,eval.points=contpoint_est)
    fx_est <- mykde_grad$estimate[,1]
    fy_est <- mykde_grad$estimate[,2]
    grad_norm_est <- sqrt(fx_est^2+fy_est^2)
    hess_est <- list(dxx=mykde_betahess$estimate[,1],
                 dxy=mykde_betahess$estimate[,2],
                 dyy=mykde_betahess$estimate[,4])
    h_init <- c(init[1,1],init[1,2],init[2,2])
    optH <- ls_opt(nrow(X), levelc,grad_norm_est,hess_est,
                   contpiece_est,initial=h_init,tol=tol,
                   gridwidth=gridwidth, maxit=maxit,print_obj = print_obj)
    return(optH)
}

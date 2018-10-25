cp_comp <- function(xxgrid,yygrid,zztrue,level){
    cl <- contourLines(xxgrid, yygrid, zztrue, levels=level)
    ## Careful, if the level set is several connects sets, then cl
    ## could be list of several lists, each standards for one
    ## connected set
    ncp <- length(cl)
    if(ncp == 1){
        contpoint <- cbind(cl[[1]]$x, cl[[1]]$y)
        contdiff <- contpoint[-nrow(contpoint),]-contpoint[-1,]
        contpiece <- sqrt(contdiff[,1]^2+contdiff[,2]^2)
        re <- list(contpoint = contpoint[-1,],contpiece=contpiece)
    }else{
        set_len <- unlist(lapply(cl, function(l) length(l$x)-1))
        contpoint <- matrix(0.,nrow = sum(set_len), ncol=2)
        contpiece <- rep(0,nrow = sum(set_len))
        pos <- cumsum(c(0,set_len))
        for(i in 1:ncp){
            start <- pos[i]+1
            end <- pos[i+1]
            contpoint_tmp <- cbind(cl[[i]]$x, cl[[i]]$y)
            contdiff_tmp <-
                contpoint_tmp[-nrow(contpoint_tmp),]-contpoint_tmp[-1,]
            contpiece_tmp <-
                sqrt(contdiff_tmp[,1]^2+contdiff_tmp[,2]^2)
            contpoint[start:end,] <- contpoint_tmp[-1,]
            contpiece[start:end] <- contpiece_tmp
        }
        re <- list(contpoint = contpoint,contpiece=contpiece)
    }
    return(re)
}


biftau <- function(zztrue,tau,gridwidth=0.01,tol=1e-6){
    high <- max(zztrue)
    low <- min(zztrue)
    while(high-low>=tol){
        mid <- (high+low)/2
        int <- sum(zztrue*(zztrue>=mid))*gridwidth^2
        if(int > 1-tau){
            low <- mid
        }else{
            high <- mid
        }
    }
    return(mid)
}


hess_haus_int <- function(hess, grad_norm, contpiece){
    intxx <- sum(hess$dxx/2/grad_norm*contpiece)
    intxy <- sum(hess$dxy/2/grad_norm*contpiece)
    intyy <- sum(hess$dyy/2/grad_norm*contpiece)
    return(list(intxx=intxx, intxy=intxy,intyy=intyy))
}


hess_leb_int <- function(ftau, zztrue, hess_grid,  gridwidth){
    intxx <- sum(hess_grid$dxx*(zztrue>=ftau))*gridwidth^2/2
    intxy <- sum(hess_grid$dxy*(zztrue>=ftau))*gridwidth^2/2
    intyy <- sum(hess_grid$dyy*(zztrue>=ftau))*gridwidth^2/2
    return(list(intxx=intxx, intxy=intxy,intyy=intyy))
}


hdr_opt <- function(n, ftau,grad_norm, hess, hess_grid, contpiece, zztrue,
                   initial,gridwidth, maxit, tol,print_obj){
    #preparing for risk function
    v1_integrals <- hess_haus_int(hess,grad_norm,contpiece)
    v2_integrals <- hess_leb_int(ftau,zztrue,hess_grid,gridwidth)
    rk <- 1/(4*pi)
    aterm <- -grad_norm/sqrt(rk*ftau)
    w0 <- 1/sum(1/grad_norm*contpiece)
    # define objective function
    wrapper <- function(h){
        nhvar_term <- sqrt(n*sqrt(h[1]*h[3]-h[2]^2)/rk/ftau)
        bterm <- -nhvar_term*(h[1]*hess$dxx+h[3]*hess$dyy+2*h[2]*hess$dxy)/2
        v1 <- h[1]*v1_integrals$intxx+2*h[2]*v1_integrals$intxy+
            h[3]*v1_integrals$intyy
        v2 <- (h[1]*v2_integrals$intxx+2*h[2]*v2_integrals$intxy+
               h[3]*v2_integrals$intyy)/ftau
        cterm <- bterm+nhvar_term*w0*(v1+v2)
        integrand <- (2*dnorm(cterm)+2*pnorm(cterm)*cterm-cterm)/(-aterm)
        hauss_int <- sum(integrand*contpiece)
        fvalue <- hauss_int*ftau/sqrt(n*sqrt(h[1]*h[3]-h[2]^2))
        return(fvalue)
    }

    #begin optimization
    hold <- c(0,0,0)
    hnew <- initial
    ite <- 0
    while (sum(abs(hold-hnew))>tol && ite<maxit){
        ite <- ite+1
        hold <- hnew
        g <- numDeriv::grad(wrapper,hold)
        H <- numDeriv::hessian(wrapper,hold)
        hnew <- hold-as.numeric(qr.solve(H,g))
        hmat <- matrix(c(hnew[1],hnew[2],hnew[2],hnew[3]),nrow=2)
        h.ei <- eigen(hmat)
        if(any(h.ei$values<=0)){
            values.new <- h.ei$values
            values.new[values.new<=0] <- 0.001
            hmatnew <-
                h.ei$vectors%*%diag(values.new)%*%t(h.ei$vectors)
            hnew <- c(hmatnew[1,1],hmatnew[1,2],hmatnew[2,2])
        }        
        if(print_obj) cat("The objective function at iteration ",
                          ite," is equal to", wrapper(hnew),"\n")
    }
    if (ite >= maxit) cat("Warning: reach maximum iterations\n")
    return(matrix(c(hnew[1],hnew[2],hnew[2],hnew[3]),nrow=2,ncol=2))

}


ls_opt <- function(n, levelc, grad_norm, hess, contpiece, 
                   initial, gridwidth, maxit, tol, print_obj){
    #preparing for risk function
    rk <- 1/(4*pi)
    aterm <- -grad_norm/sqrt(rk*levelc)
    # define objective function
    wrapper <- function(h){
        nhvar_term <- sqrt(n*sqrt(h[1]*h[3]-h[2]^2)/rk/levelc)
        bterm <- -nhvar_term*(h[1]*hess$dxx+h[3]*hess$dyy+2*h[2]*hess$dxy)/2
        integrand <- (2*dnorm(bterm)+2*pnorm(bterm)*bterm-bterm)/(-aterm)
        hauss_int <- sum(integrand*contpiece)
        fvalue <- hauss_int*levelc/sqrt(n*sqrt(h[1]*h[3]-h[2]^2))
        return(fvalue)
    }

    #begin optimization
    hold <- c(0,0,0)
    hnew <- initial
    ite <- 0
    while(sum(abs(hold-hnew))>tol && ite<maxit){
        ite <- ite+1
        hold <- hnew
        g <- numDeriv::grad(wrapper,hold)
        H <- numDeriv::hessian(wrapper,hold)
        hnew <- hold-as.numeric(qr.solve(H,g))
        hmat <- matrix(c(hnew[1],hnew[2],hnew[2],hnew[3]),nrow=2)
        h.ei <- eigen(hmat)
        if(any(h.ei$values<=0)){
            values.new <- h.ei$values
            values.new[values.new<=0] <- 0.001
            hmatnew <-
                h.ei$vectors%*%diag(values.new)%*%t(h.ei$vectors)
            hnew <- c(hmatnew[1,1],hmatnew[1,2],hmatnew[2,2])
        }        
        if(print_obj) cat("The objective function at iteration ",
                          ite," is equal to", wrapper(hnew),"\n")
    }
    if (ite >=maxit) cat("Warning: reach maximum iterations\n")
    return(matrix(c(hnew[1],hnew[2],hnew[2],hnew[3]),nrow=2,ncol=2))
}

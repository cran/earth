# earth.leaps.R:
#
# Copied from Thomas Lumley's leaps 2.9 package for earth 3.2-6 to avoid
# use of leaps::: in the earth code, to prevent complaints from CRAN check.


# leaps.setup is modified to handle linear dependencies in x properly.  I think
# this fix is needed only if leaps.setup is called with intercept=FALSE.
#
# The fix is needed because if there are linear dependencies in the matrix
# x passed to the original leaps.setup, it gives an incorrect error
# message: "missing value where TRUE/FALSE needed".
#
# In the earth context, this happened if you call earth.update with new
# data and the bx generated from that data has linear dependencies (which
# is actually ok).
#
# I also touched up the warning messages to be more informative,
# and changed most warning messages to stop messages.

leaps.setup<-function(x,y,wt=rep(1,length(y)),force.in=NULL,
                      force.out=NULL,intercept=TRUE,
                      nvmax=8,nbest=1,warn.dep=TRUE){
  make.names<-function(np){
    if (np<27) letters[1:np] else as.character(1:np)
  }
  np<-NCOL(x)
  nn<-NROW(x)
  if (length(y)!=nn) stop("y and x different lengths")
  if (length(wt)!=nn) stop("wt and x different lengths")
  if (is.null(colnames(x))) colnames(x)<-make.names(np)
  index<-rep(0,np)
  names(index)<-colnames(x)
  index[force.in]<--1
  if (any(index[force.out]==-1)) stop("Can't force the same variable in and out")
  index[force.out]<-1
  force.in<-(index==-1) ## make force.in, force.out logical vectors
  force.out<-(index==1)
  ii<-order(index)
  xx<-x[,ii]
  force.in<-force.in[ii]
  force.out<-force.out[ii]
  ones<-rep(1,np)
  names(ones)<-colnames(x)
  first<-1+sum(ones[force.in])
  last<-np-sum(ones[force.out])
  nvmax<-min(nvmax,np)
  if (intercept){
    np<-np+1
    xnames<-c("(Intercept)",colnames(xx))
    xx<-cbind(1,xx)
    colnames(xx)<-xnames
    first<-first+1
    last<-last+1
    nvmax<-nvmax+1
    index<-c(-1,index)
  }
  vorder<-1:np
  il<-nvmax*(nvmax+1)/2
  nrbar<-np*(np-1)/2
  qrleaps<-.Fortran("makeqr",np=as.integer(np),nn=as.integer(nn),
                    wt=as.double(wt),tx=t(xx),y=as.double(y),d=numeric(np),
                    rbar=numeric(nrbar),
                    thetab=numeric(np),sserr=numeric(1),ier=as.integer(0),
                    PACKAGE="earth")
  if (qrleaps$ier!=0)
      stopf("Fortran routine MAKEQR returned error code 0x%4.4X", qrleaps$ier)
  qrleaps$tx<-NULL
  qrleaps$wt<-NULL
  tolset<-.Fortran("tolset",as.integer(np),
                   as.integer(nrbar),qrleaps$d,qrleaps$rbar,
                   tol=numeric(np),numeric(np),ier=as.integer(0), PACKAGE="earth")
  if (tolset$ier!=0)
      stopf("Fortran routine TOLSET returned error code 0x%4.4X", tolset$ier)
  ss<-.Fortran("ssleaps",as.integer(np),qrleaps$d,
               qrleaps$thetab,qrleaps$sserr,rss=numeric(np),
               ier=as.integer(0), PACKAGE="earth")
  if (ss$ier!=0)
      stopf("Fortran routine SSLEAPS returned error code 0x%4.4X", ss$ier)
  sing<-.Fortran("sing",np=as.integer(qrleaps$np),nrbar=as.integer(nrbar),
                d=qrleaps$d,rbar=qrleaps$rbar,thetab=qrleaps$thetab,
                sserr=qrleaps$sserr,tol=tolset$tol,lindep=integer(qrleaps$np),
                work=numeric(qrleaps$np),ier=as.integer(0), PACKAGE="earth")
  if (sing$ier>0)
       stopf("Fortran routine SING returned error code 0x%4.4X", sing$ier)
  sing$lindep <- as.logical(sing$lindep) # from integer (0 or 1) to logical
  sing$work<-NULL
  if(any(sing$lindep)) { # linear dependencies in x? (should never happen with earth bx)
      if (intercept) {
          new.force.out <- sing$lindep | c(FALSE,force.out)
          reordered.col.nbrs <- order(new.force.out[-1]) # put lin dep cols at end
          try.again <- any((c(new.force.out,1) - c(0,new.force.out)) < 0) # huh?
          lindep.in.force.in <- any(sing$lindep[-1] & force.in)
          colnames.with.intercept <- c("(Intercept)", colnames(x))
      } else {
          new.force.out <- sing$lindep | force.out
          reordered.col.nbrs <- order(new.force.out)
          try.again <- any((c(new.force.out,1) - c(0,new.force.out)) < 0)
          lindep.in.force.in <- any(sing$lindep & force.in)
          colnames.with.intercept <- colnames(x)
      }
      if (warn.dep) {
          nsingular <- sum(sing$lindep)
          warning0("In leaps.setup, ",
                   if(try.again) "discarding " else "",
                   if(nsingular > 1) paste0(nsingular, " linearly dependent variables: ")
                   else "linearly dependent variable: ",
                   paste(colnames.with.intercept[sing$lindep], collapse=", "))
      }
      if (lindep.in.force.in)
          stop("Linear dependency in force.in variable(s)")
      if (try.again) { # recursive call
          rval<-leaps.setup(x[,ii[reordered.col.nbrs],drop=FALSE], y, wt,
                            force.in[reordered.col.nbrs], force.out[reordered.col.nbrs],
                            intercept, nvmax, nbest, warn.dep=FALSE)
          rval$reorder<-ii[reordered.col.nbrs]
          return(rval)
      }
      lastsafe<-max((1:np)[!new.force.out])
      if (lastsafe<min(nvmax,last)) {
          if (warn.dep)
              warning0("nvmax reduced from ", nvmax, " to ",
                       lastsafe-intercept, " because of linear dependencies")
          nvmax<-lastsafe
      }
  }
  if (any(sing$lindep)){
      ss<-.Fortran("ssleaps",as.integer(np),sing$d,sing$thetab,
                   sing$sserr,rss=numeric(np),ier=as.integer(0),PACKAGE="earth")
    if (ss$ier!=0)
            stopf("Fortran routine SSLEAPS returned error code 0x%4.4X", ss$ier)
  }
  initr<-.Fortran("initr",as.integer(np),as.integer(nvmax),as.integer(nbest),
                  bound=numeric(np),ress=numeric(nbest*nvmax),as.integer(nvmax),
                  lopt=integer(nbest*il),as.integer(il),vorder=as.integer(vorder),
                  ss$rss,ier=as.integer(0), PACKAGE="earth")
  if (initr$ier!=0)
      stopf("Fortran routine INITR returned error code 0x%4.4X", initr$ier)
  nullrss<-if (intercept) ss$rss[1] else sum(y^2)
  rval<-c(sing,list(nn=qrleaps$nn,rss=ss$rss,bound=initr$bound,
                    ress=matrix(initr$ress,ncol=nbest),
                    lopt=matrix(initr$lopt,ncol=nbest),
                    nvmax=nvmax,nbest=nbest,nrbar=nrbar,il=il,
                    ir=nvmax,vorder=initr$vorder,
                    first=first,last=last,xnames=colnames(xx),
                    force.in=(index==-1),force.out=(index==1),
                    intercept=intercept,nullrss=nullrss))
  class(rval)<-"regsubsets"
  invisible(rval)
}

leaps.exhaustive<-function(leaps.obj){
    if (!inherits(leaps.obj,"regsubsets"))
        stop("Not a regsubsets object -- must run leaps.setup")
    nbest<-leaps.obj$nbest
    dimwk<-3*leaps.obj$last
    dimiwk<-leaps.obj$nvmax
    rval<-.Fortran("xhaust",
                   np=as.integer(leaps.obj$np),
                   nrbar=as.integer(leaps.obj$nrbar),
                   d=leaps.obj$d,
                   rbar=leaps.obj$rbar,
                   thetab=leaps.obj$thetab,
                   first=as.integer(leaps.obj$first),
                   last=as.integer(leaps.obj$last),
                   vorder=as.integer(leaps.obj$vorder),
                   tol=leaps.obj$tol,
                   rss=leaps.obj$rss,
                   bound=leaps.obj$bound,
                   nvmax=as.integer(leaps.obj$nvmax),
                   ress=leaps.obj$ress,
                   ir=as.integer(leaps.obj$ir),
                   nbest=as.integer(leaps.obj$nbest),
                   lopt=matrix(as.integer(leaps.obj$lopt),ncol=nbest),
                   il=as.integer(leaps.obj$il),
                   wk=numeric(dimwk),
                   dimwk=as.integer(dimwk),
                   iwk=integer(dimiwk),
                   dimiwk=as.integer(dimiwk),
                   ier=as.integer(0), PACKAGE="earth")
    if(rval$ier!=0) stopf("Fortran routine XHAUST returned error code 0x%4.4X", rval$ier)
    rval$dimwk<-rval$dimiwk<-rval$iwk<-rval$wk<-NULL
    rval$xnames<-leaps.obj$xnames
    rval$method<-c("exhaustive",leaps.obj$method)
    rval$force.in<-leaps.obj$force.in
    rval$force.out<-leaps.obj$force.out
    rval$sserr<-leaps.obj$sserr
    rval$intercept<-leaps.obj$intercept
    rval$lindep<-leaps.obj$lindep
    rval$reorder<-leaps.obj$reorder
    rval$nullrss<-leaps.obj$nullrss
    rval$nn<-leaps.obj$nn
    class(rval)<-"regsubsets"
    rval
}


leaps.backward<-function(leaps.obj){
  if (!inherits(leaps.obj,"regsubsets"))
      stop("Not a regsubsets object -- must run leaps.setup")
  nbest<-leaps.obj$nbest
  dimwk<-2*leaps.obj$last
  rval<-.Fortran("bakwrd",np=as.integer(leaps.obj$np),
                 nrbar=as.integer(leaps.obj$nrbar),d=leaps.obj$d,
                 rbar=leaps.obj$rbar,thetab=leaps.obj$thetab,
                 first=as.integer(leaps.obj$first),
                 last=as.integer(leaps.obj$last),
                 vorder=as.integer(leaps.obj$vorder),tol=leaps.obj$tol,
                 rss=leaps.obj$rss,bound=leaps.obj$bound,
                 nvmax=as.integer(leaps.obj$nvmax),
                 ress=leaps.obj$ress,ir=as.integer(leaps.obj$ir),
                 nbest=as.integer(leaps.obj$nbest),
                 lopt=matrix(as.integer(leaps.obj$lopt),ncol=nbest),
                 il=as.integer(leaps.obj$il),wk=numeric(dimwk),
                 iwk=as.integer(dimwk),ier=as.integer(0), PACKAGE="earth")
  if(rval$ier!=0)
      stopf("Fortran routine BAKWRD returned error code 0x%4.4X", rval$ier)
  rval$dimwk<-rval$wk<-NULL
  rval$xnames<-leaps.obj$xnames
  rval$method<-c("backward",leaps.obj$method)
  rval$force.in<-leaps.obj$force.in
  rval$force.out<-leaps.obj$force.out
  rval$sserr<-leaps.obj$sserr
  rval$intercept<-leaps.obj$intercept
  rval$lindep<-leaps.obj$lindep
  rval$reorder<-leaps.obj$reorder
  rval$nullrss<-leaps.obj$nullrss
  rval$nn<-leaps.obj$nn
  class(rval)<-"regsubsets"
  rval
}


leaps.forward<-function(leaps.obj){
  if (!inherits(leaps.obj,"regsubsets"))
      stop("Not a regsubsets object -- must run leaps.setup")
  nbest<-leaps.obj$nbest
  dimwk<-3*leaps.obj$last
  rval<-.Fortran("forwrd",np=as.integer(leaps.obj$np),
                 nrbar=as.integer(leaps.obj$nrbar),
                 d=leaps.obj$d,rbar=leaps.obj$rbar,
                 thetab=leaps.obj$thetab,first=as.integer(leaps.obj$first),
                 last=as.integer(leaps.obj$last),
                 vorder=as.integer(leaps.obj$vorder),
                 tol=leaps.obj$tol,rss=leaps.obj$rss,
                 bound=leaps.obj$bound,nvmax=as.integer(leaps.obj$nvmax),
                 ress=leaps.obj$ress,ir=as.integer(leaps.obj$ir),
                 nbest=as.integer(leaps.obj$nbest),
                 lopt=matrix(as.integer(leaps.obj$lopt),ncol=nbest),
                 il=as.integer(leaps.obj$il),wk=numeric(dimwk),
                 iwk=as.integer(dimwk),ier=as.integer(0), PACKAGE="earth")
  if(rval$ier!=0)
      stopf("Fortran routine FORWARD returned error code 0x%4.4X", rval$ier)
  rval$dimwk<-rval$wk<-NULL
  rval$xnames<-leaps.obj$xnames
  rval$method<-c("forward",leaps.obj$method)
  rval$force.in<-leaps.obj$force.in
  rval$force.out<-leaps.obj$force.out
  rval$sserr<-leaps.obj$sserr
  rval$intercept<-leaps.obj$intercept
  rval$lindep<-leaps.obj$lindep
  rval$reorder<-leaps.obj$reorder
  rval$nullrss<-leaps.obj$nullrss
  rval$nn<-leaps.obj$nn
  class(rval)<-"regsubsets"
  rval
}


leaps.seqrep<-function(leaps.obj){
    if (!inherits(leaps.obj,"regsubsets"))
        stop("Not a regsubsets object -- must run leaps.setup")
    nbest<-leaps.obj$nbest
    dimwk<-3*leaps.obj$last
    rval<-.Fortran("seqrep",np=as.integer(leaps.obj$np),
                   nrbar=as.integer(leaps.obj$nrbar),
                   d=leaps.obj$d,rbar=leaps.obj$rbar,
                   thetab=leaps.obj$thetab,
                   first=as.integer(leaps.obj$first),
                   last=as.integer(leaps.obj$last),
                   vorder=as.integer(leaps.obj$vorder),
                   tol=leaps.obj$tol,rss=leaps.obj$rss,
                   bound=leaps.obj$bound,nvmax=as.integer(leaps.obj$nvmax),
                   ress=leaps.obj$ress,ir=as.integer(leaps.obj$ir),
                   nbest=as.integer(leaps.obj$nbest),
                   lopt=matrix(as.integer(leaps.obj$lopt),
                   ncol=nbest),il=as.integer(leaps.obj$il),
                   wk=numeric(dimwk),iwk=as.integer(dimwk),
                   ier=as.integer(0), PACKAGE="earth")
  if(rval$ier!=0)
      stopf("Fortran routine SEQREP returned error code 0x%4.4X", rval$ier)
  rval$dimwk<-rval$wk<-NULL
  rval$xnames<-leaps.obj$xnames
  rval$method<-c("seqrep",leaps.obj$method)
  rval$force.in<-leaps.obj$force.in
  rval$force.out<-leaps.obj$force.out
  rval$sserr<-leaps.obj$sserr
  rval$intercept<-leaps.obj$intercept
  rval$lindep<-leaps.obj$lindep
  rval$reorder<-leaps.obj$reorder
  rval$nullrss<-leaps.obj$nullrss
  rval$nn<-leaps.obj$nn
  class(rval)<-"regsubsets"
  rval
}

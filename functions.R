################



######
# Estimator of inverse density weighted expectation
# over cube [a,b]

expec.bias.cor = function(Y,Z,K=6,L=1,a,b,m=10000, tol=0.2)
{
  if (!require(dbscan)) {
    
    stop("dbscan not installed")
    
  } 
  else {
    K=as.integer(K)
    L=as.integer(L)
    
    if((is.vector(a)==FALSE) || (is.vector(b)==FALSE) || (is.matrix(Z)==FALSE) || 
       (is.vector(Y)==FALSE) || length(Y) != nrow(Z) || ncol(Z) != length(a) || 
       length(b) != length(a) || prod(b>a)==0 || L<0 || K < 1 || L>1 || (K <= ncol(Z) +1 & L==1))
      
      stop("Wrong Parameters")
    
    else
    {
   
      Z_max_tol=apply(Z,2,max)+tol
      Z_min_tol=apply(Z,2,min)-tol
      
      if(prod(a>Z_min_tol)*prod(b<Z_max_tol)==0) stop("Region [a,b] not suitable for this tolerance")
      
      grid=runif(m,min=a[1],max=b[1])
      if(ncol(Z)>=2)
      {
        for(i in 2:ncol(Z)) grid=cbind(grid,runif(m,min=a[i],max=b[i]))
      }
      
      
      nn=kNN(Z,K,grid)
      
      if(L==0) 
      {
        return(mean(Y[nn$id])*prod((b-a)))
      }  
      else
      {
        estval=0
        for(i in 1:m)
        {
          I=nn$id[i,]
          M=cbind(rep(1,K),t(apply(Z[I,],1,function(x){x-grid[i,]})))
          estval=estval+(solve(t(M)%*%M)%*%t(M)%*%Y[I])[1]
        }
        return(prod((b-a))*estval/m)
      }
    }    
  }
}


######
# Matching estimator, additional sample X required

match.bias.cor = function(Y,Z,X,K=6,L=1,del=0)
{
  if (!require(dbscan)) {
    
    stop("dbscan not installed")
    
  } 
  else {
    K=as.integer(K)
    L=as.integer(L)
    
    if((is.matrix(X)==FALSE) || (is.matrix(Z)==FALSE) || (is.vector(Y)==FALSE) || 
       length(Y) != nrow(Z) || ncol(Z) != ncol(X) || del < 0 
       || L<0 || K < 1 || L>1 || (K <= ncol(Z) +1 & L==1))
      
      stop("Wrong Parameters")
    
    else
    {
      nn=kNN(Z,K,X)
      Z_max_del=apply(Z,2,max)-del
      Z_min_del=apply(Z,2,min)+del
      ind_max=X <=Z_max_del
      ind_max=apply(ind_max,1,prod)
      ind_min=X >=Z_min_del
      ind_min=apply(ind_min,1,prod)
      
      if(L==0) 
      {
        y_res=Y[nn$id[,1]]
        if(K >=2)
        {
          for(i in 2:K) y_res = cbind(y_res,Y[nn$id[,i]])
          y_res=apply(y_res,1,mean)
        }
        return(mean(y_res*ind_min*ind_max))
      }  
      else
      {
        estval=0
        for(i in 1:nrow(X))
        {
          I=nn$id[i,]
          M=cbind(rep(1,K),t(apply(Z[I,],1,function(x){x-X[i,]})))
          estval=estval+(solve(t(M)%*%M)%*%t(M)%*%Y[I])[1]*ind_min[i]*ind_max[i]
        }
        return(estval/nrow(X))
      }
    }    
  }
}


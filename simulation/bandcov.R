P <-
function(m,n){
  if(m>n){
    stop("m cannot be larger than n")
  }
  perm<-.C("P",as.integer(m),as.integer(n),as.double(1))
  perm[[3]]
}
sscov <-
function(x){
  if(nrow(x)<4){
    stop("X must have more than 4 rows")
  }
  if(typeof(x)!="double"){
    stop("X must be numeric")
  }
  M<-c(nrow(x),ncol(x))
  a<-matrix(1,nrow=1,ncol=M[1])
  A<-t(a%*%x)
  B<-sum(diag(x%*%t(x)))
  mf22<-t(A)%*%A-B
  mf1<-1/P(2,M[1])+2/P(3,M[1])+2/P(4,M[1])
  mf2<-1/P(4,M[1])
  mf3<-2/P(3,M[1])+4/P(4,M[1])
  res<-.C("adT",as.double(x),as.double(M[1]),as.double(M[2]),as.integer(a),as.double(A),
          as.double(B),as.double(mf22),as.double(mf1),as.double(mf2),as.double(mf3),as.double(0))
  as.numeric(res[11])
}
ssbandcov <-
function(x,k){
  if(typeof(x)!="double"){
    stop("X must be numeric")
  }
  if(k>=ncol(x)){
    stop("k must be less than p")
  }
  res<-.C("adb",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),as.integer(k),as.double(0))[5]
  as.numeric(res)
}
bandtest.stat <-
function(x,k){
  if(typeof(x)!="double"){
    stop("X must be numeric")
  }
  M<-c(nrow(x),ncol(x))
  ts<-M[1]*(sscov(x)/ssbandcov(x,k)-1)
  ts
}
bandpen2 <-
function(x,q){
  if(typeof(x)!="double"){
    stop("X must be numeric")
  }
  if(q>=ncol(x)){
    stop("q must be less than p")
  }
  e1<-as.numeric(.C("sbandpen2call",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),as.integer(q),as.double(0))[5])
  e2<-as.numeric(.C("pbandpen2call",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),as.integer(q),as.double(0))[5])
  c(e1,e2)
  
}
bandpen <-
function(x,K){
  if(typeof(x)!="double"){
    stop("X must be numeric")
  }
  if(K>=ncol(x)){
    stop("K must be less than p")
  }
  store = matrix(0, K + 1, 2)
  res<-.C("bandpen",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),as.integer(K),as.double(store))[[5]]
  resM<-matrix(res,nrow=K+1,ncol=2,byrow=T)
  resM[1,]<-bandpen2(x,0)
  resM
}
bandloss <-
function(x,K){
  adt<-sscov(x)
  tempres<-bandpen(x,K)
  res.b<-rep(0,times=K+1)
  bl<-.C("bandloss",as.integer(nrow(x)),as.integer(ncol(x)),as.integer(K),
         as.double(adt),as.double(tempres),as.double(res.b))[[6]]
  bl
}
tapeloss <-
function(x,K){
  adt<-sscov(x)
  tempres<-bandpen(x,K)
  res.t<-rep(0,times=floor(K/2)+1)
  tl<-.C("tapeloss",as.integer(nrow(x)),as.integer(ncol(x)),as.integer(K),
         as.double(adt),as.double(tempres),as.double(res.t))[[6]]
  tl
}
bandmin <-
function(x, K){
  bandprop = bandloss(x, K)
  IXband = order(bandprop)[1]-1
  IXband
}
tapemin <-
function(x, K){
  tapeprop = tapeloss(x, K)
  IXtape = order(tapeprop)[1]-1
  IXtape
}
btplot <-
function(x,K,plot="S"){
  K.vec<-0:K
  adt<-sscov(x)
  tempres<-bandpen(x,K)
  res.b<-rep(0,times=K+1)
  res.t<-rep(0,times=floor(K/2)+1)
  if(plot=="B"||plot=="Band"||plot=="Bandloss"){
    band.plot<-.C("bandloss",as.integer(nrow(x)),as.integer(ncol(x)),as.integer(K),
                  as.double(adt),as.double(tempres),as.double(res.b))[[6]] 
    plot(x=K.vec,y=band.plot,xlab="Bandwidth",ylab="Loss",main="Banding Loss")
    lines(K.vec[order(K.vec)],band.plot[order(K.vec)],xlim=range(K.vec),ylim=range(band.plot))
  }
  
  if(plot=="T"||plot=="Tape"||plot=="Tapeloss"){
    f<-floor(K/2)
    K.vech<-0:f
    tape.plot<-.C("tapeloss",as.integer(nrow(x)),as.integer(ncol(x)),as.integer(K),
                  as.double(adt),as.double(tempres),as.double(res.t))[[6]]
    plot(x=K.vech,y=tape.plot,xlab="Bandwidth",ylab="Loss",main="Tapering Loss")
    lines(K.vech[order(K.vech)],tape.plot[order(K.vech)],xlim=range(K.vech),ylim=range(tape.plot))
    
  }
  if(plot=="S"||plot=="SS"||plot=="Squares"){
    ss.plot<-c()
    
    for(i in K.vec){
      ss.plot[i+1]<-as.numeric(.C("sbandpen2call",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),as.integer(i),as.double(0))[5])
    }
    plot(x=K.vec,y=ss.plot,xlab="Bandwidth",ylab="Sum of Squares",main="Sum of Squares for Off Diagonal")
    lines(K.vec[order(K.vec)],ss.plot[order(K.vec)],xlim=range(K.vec),ylim=range(ss.plot))
  }
}

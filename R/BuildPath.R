BuildPath=function(p, K, propOverl=0){
  xx=round(p/(K-(K-1)*propOverl/100));
  path=matrix(0,K,p)
  for (k in seq(0,K-1)){
    x1=round((xx-propOverl*xx/100)*k);
    for (j in (x1):(x1+xx-1)){
      if (j<p){ path[k+1,j+1]=1;}
    }
  }
  for (j in seq(0,p-1)){
    x=1;
    for (k in seq(0,K-1)){
      x=x*(1-path[k+1,j+1]);
    }
    if (x==1) {path[K,j+1]=1;}
  }
  return (path);
}

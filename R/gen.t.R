gen.t = function(x, df = 4, mu, sig, log = TRUE){
  tmp = (gamma((df+1)/2)/gamma(df/2)/sqrt(pi*df)/sig)*(1+(x-mu)^2/sig^2/df)^(-(df+1)/2);
  if(log){tmp = log(tmp);}
  return(tmp);
}

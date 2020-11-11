#############################################################################
#                                                                           #     
#       Transfer functions between three (or four) functional traits        #
#               and the demographic rates                                   #
#       Loïc Chalmandrier, Florian Hartig, Loïc Pellissier                  #
#                                                                           #
#############################################################################

Trait2Demo <- function(t1, t2, t3, phi1, phi2,a,b, exp = F){
  x1 = cos(phi1)
  x2 = sin(phi1)*cos(phi2)
  x3  = sin(phi1)*sin(phi2)
  
  if (exp){
    theo <- exp(a*(x1*t1 + x2*t2 + x3*t3)+b)
  }else{
    theo <- a*(x1*t1 + x2*t2 + x3*t3)+b
  }
  return(theo)
}

Trait2Demo_n4 <- function(t1, t2, t3, t4, phi1, phi2,phi3, a,b, exp = F){
  x1 = cos(phi1)
  x2 = sin(phi1)*cos(phi2)
  x3  = sin(phi1)*sin(phi2)*cos(phi3)
  x4  = sin(phi1)*sin(phi2)*sin(phi3) 
  
  if (exp){
    theo <- exp(a*(x1*t1 + x2*t2 + x3*t3+ x4*t4)+b)
  }else{
    theo <- a*(x1*t1 + x2*t2 + x3*t3 + x4*t4)+b
  }
  return(theo)
}
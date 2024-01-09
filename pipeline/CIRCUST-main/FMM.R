# Funcon to simply represent an FMM
# Receives as arguments:
# - M --> Value for M
# - A --> value for A
# - a --> value for alpha
# - b --> value for beta
# - w --> value for omega
# - from --> starting time. By default it is 0.
# - to --> end time instant By default it is 2*pi.
# - timePoints --> It is also possible to define preset moments of time
# - plot --> if yes, it will be available through the standard graphics output of a
# graphical representation of the FMM
# - outvalues ​​--> if yes, a list with the time instants will be returned
# and the values ​​produced.
# - length.out --> simulation length
# - sigmaNoise --> if yes, normal noise is added with s.d. specified.
# Instead of setting individual values ​​for M, A, alpha, beta, and omega,
# an array can be established. In this case a linear combination of FMMs is generated.
# The missing ones are added due to replication

FMM <- function(M,A,a,b,w,from=0,to=2*pi,length.out=100,timePoints=seq(from,to,length=length.out),
                plot=T,outvalues=F,sigmaNoise=0){
  
  narg <- max(length(M),length(A),length(a),length(b),length(w))
  M <- rep(M,length.out=narg)
  A <- rep(A,length.out=narg)
  a <- rep(a,length.out=narg)
  b <- rep(b,length.out=narg)
  w <- rep(w,length.out=narg)
  
  t <- timePoints
  
  phi <- list()
  for(i in 1:narg){
    phi[[i]] <- b[i]+2*atan(w[i]*tan((t-a[i])/2))
  }
  
  ym <- list()
  for(i in 1:narg){
    ym[[i]] <- M[i]+A[i]*cos(phi[[i]])
  }
  
  y <- rep(0,length(t))
  for(i in 1:narg){
    y <- y + ym[[i]]
  }
  
  if (sigmaNoise > 0) y <- y + rnorm(length.out,0,sigmaNoise)
  
  if(plot) plot(t,y,type="l",lwd=2,col=2,xlab="tiempo",ylab="Respuesta",
       main=paste("FMM:   M=",M,"   A=",A,"   a=",a,"   b=",b,"   w=",w,sep=""))
  
  if(outvalues) return(list(t=t,y=y))
}


predict_FMM <- function(M,A,alpha,beta,omega,timePoints){
  return(FMM(M,a,alpha,beta,omega,timePoints = timePoints,plot=F,outvalues = T)$y)
}

# Function to fit an FMM model on data.
# Receives as input the data and the moments of time where they are observed:
# vData --> values ​​for tuning
# timePoints --> instants of time in which the adjustment occurs.
# By default they are equally spaced
# lengthAlphaGrid --> Length of the GRID used for alpha. It must be indicated with an integer.
# It is not necessary to indicate this if a custom GRID is to be used through
# the alphaGrid argument. By default it is 24.
# lengthOmegaGrid --> Length of the GRID used for omega Must be indicated with an integer.
# It is not necessary to indicate this if a custom GRID is to be used through
# the omegaGrid argument. By default it is 24.
# alphaGrid --> Alpha values ​​used for initial estimation
# If the approximate value of alpha is known, it is advisable to indicate it here
# By default, equally spaced between 0 and 2pi, using 24 elements
# omegaMax --> maximum value that the omega parameter can take. In no case is it returned
# an estimated omega greater than specified.
# If the approximate value of omega is known, it is advisable to indicate a dimension here
# omegaGrid --> Omega values ​​used for initial estimation
# By default it takes values ​​between 0.0001 and 1 on a logarithmic scale
# numReps --> number of times step 1 + 2 is repeated. The output of step 2 is taken as a new input
# around which to build a new GRID. The estimation is repeated for a total of numReps.
# By default it is 3, that is, it is done three times.
# REWRITTEN WITH RESPECT TO THE ORIGINAL
# PROBABLY NEEDS TO CHANGE THE NAME
fitFMM_Par<-function(vData, timePoints = seq(0,2*pi,length.out = length(vData)),
                     lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                     alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid), omegaMax = 1,
                     omegaGrid = exp(seq(log(0.0001),log(omegaMax),length.out=lengthOmegaGrid)),
                     numReps = 3){
  
  n <- length(vData)

  ## Step 1: Calculate the initial estimators of M, A, alpha, beta and omega
  #To do this, alpha and omega are set and the others are calculated using a Cosinor model
  grid <- expand.grid(alphaGrid,omegaGrid)
  
  #To avoid an inefficient double loop in R, the apply function is used with a helper function
  #This auxiliary function is responsible for estimating beta, M, A, and also returns the RSS as an aggregate
  #to facilitate the search for the best one after the fact
    
  step1 <- t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints))
  colnames(step1) <- c("M","A","alpha","beta","omega","RSS")
  
  #Of all the adjusted ones, we look for the one that has a lower RSS but that meets certain
  #stability conditions. For this the bestStep1 function is used.
  bestPar <- bestStep1(vData,step1)
  #The initial estimators for M, A, alpha, beta and omega have already been found
  
  #Step-2: Nelder-Mead Optimization for RSS
  nelderMead <- optim(par = bestPar[1:5], fn = step2FMM, vData = vData, timePoints = timePoints, omegaMax = omegaMax)
  parFinal <- nelderMead$par
  SSE <- nelderMead$value*n
  
  #Take alpha and beta to the interval [0.2pi]
  parFinal[3] <- parFinal[3]%%(2*pi)
  parFinal[4] <- parFinal[4]%%(2*pi)
  
  #The GRID is repeated the specified number of times (by default this option is NOT applied)
  numReps <- numReps - 1
  while(numReps > 0){
    
    #New GRID for alpha, the same length as the previous one. It is ensured that it is between 0 and 2pi.
    nAlphaGrid <- length(alphaGrid)
    amplitudeAlphaGrid <- 1.5*mean(diff(alphaGrid))
    alphaGrid <- seq(parFinal[3]-amplitudeAlphaGrid,parFinal[3]+amplitudeAlphaGrid,length.out = nAlphaGrid)
    alphaGrid <- alphaGrid%%(2*pi)
    
    #New GRID for omega, the same length as the previous one. It is ensured that it is between 0 and omegaMax
    nOmegaGrid <- length(omegaGrid)
    amplitudeOmegaGrid <- 1.5*mean(diff(omegaGrid))
    omegaGrid <- seq(max(parFinal[5]-amplitudeOmegaGrid,0),
                     min(omegaMax,parFinal[5]+amplitudeOmegaGrid),
                     length.out = nOmegaGrid)
    
    
    ## Step 1
    grid <- as.matrix(expand.grid(alphaGrid,omegaGrid))
    step1 <- t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints))
    colnames(step1) <- c("M","A","alpha","beta","omega","RSS")
    antBestPar <- bestPar
    bestPar <- bestStep1(vData,step1)
    
    #There may not be any that satisfy the conditions
    if(is.null(bestPar)){
      bestPar <- antBestPar
      numReps <- 0
      warning("Es posible que el FMM no sea adecuado en este caso")
    }
    
    
    ## Step 2
    nelderMead <- optim(par = bestPar[1:5], fn = step2FMM, vData = vData, timePoints = timePoints, omegaMax = omegaMax)
    parFinal <- nelderMead$par
    
    #Take alpha and beta to the interval [0.2pi]
    parFinal[3] <- parFinal[3]%%(2*pi)
    parFinal[4] <- parFinal[4]%%(2*pi)
    
    
    numReps <- numReps - 1
  }
  
  names(parFinal) <- c("M","A","alpha","beta","omega")
  
  #Returns, in addition to the parameters, an adjustment at the specified time intervals
  #and the SSE (sum squared error).
    
  adjMob <- parFinal["M"] + parFinal["A"]*cos(parFinal["beta"] + 2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((adjMob-vData)^2)
  
  outMobius<-list(FMM=adjMob,alpha=parFinal[3],beta=parFinal[4],omega=parFinal[5],M=parFinal["M"],
                  A=parFinal["A"], SSE = SSE)
  return(outMobius)
}



# Function to fit a multiple FMM model (backfitting) on ​​data.
# Receives as input the data and the moments of time where they are observed:
# vData --> values ​​for tuning
# timePoints --> instants of time in which the adjustment occurs.
# By default they are equally spaced
# nback --> number of components to be adjusted by backfitting
# maxiter --> maximum number of iterations to perform. By default it is the same as nback.
# stopFunction --> optional parameter. You must specify a function that takes three arguments:
# the actual values ​​(vData), the value adjusted by the FMM in the current iteration,
# and the value adjusted by the FMM in the previous iteration. Should return TRUE if desired
# stop backfitting, or FALSE if it should continue because convergence has not yet been reached.
# By default it always returns FALSE: stopFunction = function(vData,pred,prevPred){return(FALSE)}
# objectFMM --> an FMM object, with at most nback components, to continue iterating over.
# This option allows you to continue a setting from a previous FMM object
# saving time and resources by not starting from the beginning.
# If the number of components is not the same, the other parameters are padded with zero.
# lengthAlphaGrid --> Length of the GRID used for alpha. It must be indicated with an integer.
# It is not necessary to indicate this if a custom GRID is to be used through
# the alphaGrid argument. By default it is 24.
# lengthOmegaGrid --> Length of the GRID used for omega Must be indicated with an integer.
# It is not necessary to indicate this if a custom GRID is to be used through
# the omegaGrid argument. By default it is 24.
# alphaGrid --> Alpha values ​​used for initial estimation
# If the approximate value of alpha is known, it is advisable to indicate it here
# By default, equally spaced between 0 and 2pi, using 24 elements
# omegaMax --> maximum value that the omega parameter can take. In no case is it returned
# an estimated omega greater than specified.
# If the approximate value of omega is known, it is advisable to indicate a dimension here
# omegaGrid --> Omega values ​​used for initial estimation
# By default it takes values ​​between 0.0001 and 1 on a logarithmic scale
# numReps --> number of times step 1 + 2 is repeated. The output of step 2 is taken as a new input
# around which to build a new GRID. The estimation is repeated for a total of numReps.
# By default it is 3, that is, it is done three times.
# REWRITTEN WITH RESPECT TO THE ORIGINAL
# PROBABLY NEEDS TO CHANGE THE NAME

fitFMM_back<-function(vData, timePoints = seq(0,2*pi,length.out = length(vData)), nback, maxiter=nback,
                      stopFunction = function(vData,pred,prevPred){return(FALSE)}, objectFMM = NULL,
                     lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                     alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid), omegaMax = 1,
                     omegaGrid = exp(seq(log(0.0001),log(omegaMax),length.out=lengthOmegaGrid)),
                     numReps = 3, showProgress = T){
  
  n <- length(vData)
  
  if(showProgress){
    marcasTotales <- 50
    granoInforme <- 2
    cat("|")
    for(m in 1:marcasTotales) cat("-")
    cat("|\n")
    cat("|")
    porcentajeCompletado <- 0.00001
  }
  
  #Initialization of structures to store components.
  predichosComponente <- list()
  ajusteComponente <- list()
  
  #En caso de comenzar desde el principio (no se aporta un objectFMM)
  if(is.null(objectFMM)){
    for(i in 1:nback){
      predichosComponente[[i]] <- rep(0,n)
    }
    prevAdjMob <- NULL
    
  #In case of starting from the beginning (an objectFMM is not provided)
  } else {
    prevAdjMob <- objectFMM$FMM
    nbackAnterior <- length(objectFMM$alpha)
    if(nbackAnterior > nback){
      stop("Impossible to reduce dimensions from input objectFMM")
    }
    for(i in 1:nback){
      if(i <= nbackAnterior){
        predichosComponente[[i]] <- objectFMM$M/nbackAnterior + objectFMM$A[i]*cos(objectFMM$beta[i] + 
                                      2*atan(objectFMM$omega[i]*tan((timePoints-objectFMM$alpha[i])/2)))
      } else {
        predichosComponente[[i]] <- rep(0,n)
      }
    }
  }
  
  
  #Backfitting up to a maximum of maxiter
  for(i in 1:maxiter){
    
    #Backfitting por componentes
    for(j in 1:nback){
      
      #Obtenemos la componente j a ajustar como la diferencia de los datos originales con respecto
      #al ajuste de todas las demás componentes
      vDataAjuste <- vData
      for(k in 1:nback){
        if(j != k){
          vDataAjuste <- vDataAjuste - predichosComponente[[k]]
        }
      }
      
      #Realizamos el ajuste de la componente j con las opciones especificadas
      ajusteComponente[[j]] <- fitFMM_Par(vDataAjuste,timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                          lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid, omegaMax = omegaMax,
                                          omegaGrid = omegaGrid, numReps = numReps)
      predichosComponente[[j]] <- ajusteComponente[[j]]$FMM
      
      #Esto se utiliza para mostrar progreso visible
      if(showProgress){
        porcentajeAntes <- porcentajeCompletado
        porcentajeCompletado <- porcentajeCompletado + 100/(nback*maxiter)
        numMarcas <- sum((seq(ceiling(porcentajeAntes),floor(porcentajeCompletado),by=1)%%granoInforme == 0))
        if (numMarcas > 0) for(m in 1:numMarcas) cat("=")
      }
      
    }
    
    #Comprobación de la condición de parada
    #Calculamos los valores predichos por el modelo como la suma de todas las componentes
    adjMob <- rep(0,n)
    for(j in 1:nback){
      adjMob <- adjMob + predichosComponente[[j]]
    }
    if(!is.null(prevAdjMob)){
      if(stopFunction(vData,adjMob,prevAdjMob)){
        break
      }
    }
    prevAdjMob <- adjMob
    
    
  }
  
  #Esto se utiliza para mostrar progreso visible
  if(showProgress){
    if(porcentajeCompletado < 100){
      porcentajeAntes <- porcentajeCompletado
      porcentajeCompletado <- 100
      numMarcas <- sum((seq(ceiling(porcentajeAntes),floor(porcentajeCompletado),by=1)%%granoInforme == 0))
      for(m in 1:numMarcas) cat("=")
    }
    cat("|\n")
    if(i == maxiter){
      cat("Stopped by reaching maximum iterations\n")
    } else {
      cat("Stopped by the stopFunction\n")
    }
  }

  #El parámetro M es común a todas las componentes (indica el desplazamiento)
  #Debe ser así para evitar singularidades
  #Se toma el parámetro M como la suma de todos los M de cada componente
  #Los parámetros A, alpha, beta y omega son únicos de cada componente
  alpha <- rep(0,nback)
  beta <- rep(0,nback)
  omega <- rep(0,nback)
  for(j in 1:nback){
    alpha[j] <- ajusteComponente[[j]]$alpha
    beta[j] <- ajusteComponente[[j]]$beta
    omega[j] <- ajusteComponente[[j]]$omega
  }
  
  #Se recalculan A y M con una regresion lineal
  phi <- list()
  for(j in 1:nback){
    phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
  }
  M <- matrix(unlist(phi),ncol=nback)
  regresion <- lm(vData ~ M)
  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]
  
  #Calculamos los valores predichos por el modelo
  adjMob <- predict(regresion)

  #Calculo del SSE
  SSE <- sum((adjMob-vData)^2)

  
  #Se devuelve una salida similar al caso sin backfitting
  outMobius<-list(FMM=adjMob,alpha=alpha,beta=beta,omega=omega,M=M,A=A,SSE = SSE)
  return(outMobius)

}


#dropComponents_FMM <- function(objectFMM,numOfComponents=0,components=NULL,revise=T){
  
  
#  FMM(20,2,c(1.5,3.4,5),c(0.2,2.3,4.7),c(0.1,0.2,0.05))
#  res <- FMM(20,2,c(1.5,3.4,5),c(0.2,2.3,4.7),c(0.1,0.2,0.05),outvalues = T,plot=F,sigmaNoise = 0.1)
#  fit <- fitFMM_back(res$y,res$t,nback = 4)
#  plot(res$t,res$y,xlab="Tiempo",ylab="Respuesta")
#  points(res$t,fit$FMM,col=2,type="l",lwd=2)
#  plot_FMM(fit,res$y,res$t,components = T)
  
#  nback <- length(objectFMM$alpha)
#  if(is.null(components)){
#    components <- seq(len)
#  }
  
#}


# Función auxiliar que calcula el porcentaje de variabilidad recogido
# Recibe como entradas:
#   vData --> datos originales
#   pred --> datos predichos por el modelo
# FUNCTION OCULTA
PV <- function(vData,pred){
  meanVData <- mean(vData)
  return(1 - sum((vData-pred)^2)/sum((vData-meanVData)^2))
}

# Función posible para utilizar como criterio de parada en el backfitting.
# Devuelve TRUE (parada) si se cumple al menos una de estas dos condiciones:
# cond1 -> ((PV(i) - PV(i-1)) > 0.0001) & ((PV(i)-PV(i-1)) < 0.025*(1-PV(i-1)))
# cond2 -> (PV(i) - PV(i-1)) <= 0.0001
# Donde PV(i) es el porcentaje de variabilidad recogido en la iteración i
# Recibe como entradas:
#   vData --> datos originales
#   pred --> valores predichos en la iteración actual
#   prevPred --> valores predichos en la iteración anterior
stopFunction1 <- function(vData,pred,prevPred){
  PV_pred <- PV(vData,pred)
  PV_prevPred <- PV(vData,prevPred)
  cond1 <- ((PV_pred-PV_prevPred) > 0.0001) & ((PV_pred-PV_prevPred) < 0.025*(1-PV_prevPred))
  cond2 <- (PV_pred-PV_prevPred) <= 0.0001
  if(cond1 | cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# Función auxiliar para el primer paso de la estimación FMM
# Recibe como entradas:
#   alphaOmega --> un vector que contiene dos parámetros: alpha y omega
#   vData --> valores para el ajuste
#   timePoints --> instantes de tiempo en el que se produce el ajuste.
#La salida producida es un vector de seis componentes: M, A, alpha, beta, omega y RSS
# FUNCION OCULTA
step1FMM <- function(alphaOmega, vData, timePoints) {
  
  
  alpha <- alphaOmega[1]
  omega <- alphaOmega[2]
  
  parteMobius <- 2*atan(omega*tan((timePoints-alpha)/2))
  t_star <- alpha + parteMobius
  
  #Modelo cosinor suponiendo fijos alpha y omega
  xx<-cos(t_star)
  zz<-sin(t_star)
  fit<-lm((vData)~xx+zz)
  M <-fit$coefficients[1] #intercept
  bb<-fit$coefficients[2] #Coeficiente del coseno
  gg<-fit$coefficients[3] #Coeficiente del seno
  phiEst<-atan2(-gg,bb) #acrophase (phi)
  A<-sqrt(bb^2+gg^2) #Amplitud de la onda
  beta <- (phiEst+alpha)%%(2*pi)
  
  dataReg<-cos(beta+parteMobius) #Resultado de Mobius sin amplitud e intercept
  
  adj0<-M+A*dataReg #Mobius Reg
  RSS<-sum((vData-adj0)^2)/length(timePoints) #suma de cuadrados de residuos (medio)
  
  devolver <- c(M,A,alpha,beta,omega,RSS)
  return(devolver)
  
}




#Función para uso interno que permite encontrar el parámetro con menor RSE. Toma como entradas:
#  vData --> el vector que contiene los datos originales
#  step1 --> la salida de t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints)),
#            que es un data.frame que tiene por columnas estimaciones de M, A, alpha, beta, omega, RSS
#  Devuelve la fila de step1 que se haya elegido.
# FUNCION OCULTA
bestStep1 <- function(vData,step1){
  
  #De todos los ajustados, buscamos el que tiene un menor RSS pero que cumpla ciertas
  #condiciones de estabilidad. 
  ordenMinRSS <- order(step1[,"RSS"])
  maxVData <- max(vData)
  minVData <- min(vData)
  n <- length(vData)
  
  #Empezamos desde el que tiene un RSS menor y vamos buscando el primero que cumpla las condiciones
  condicionContinuar <- TRUE
  i <- 1
  while(condicionContinuar){
    
    #Recuperar los parámetros y el RSS
    M <- step1[ordenMinRSS[i],"M"]
    A <- step1[ordenMinRSS[i],"A"]
    alpha <- step1[ordenMinRSS[i],"alpha"]
    beta <- step1[ordenMinRSS[i],"beta"]
    omega <- step1[ordenMinRSS[i],"omega"]
    sigma <- sqrt(step1[ordenMinRSS[i],"RSS"]*n/(n-5))
    
    #Restricciones de estabilidad
    maxi <- M + A
    mini <- M - A
    rest1 <- maxi <= maxVData+0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
    rest2 <- mini >= minVData-0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
    
    #Si se cumple, parar aquí. Hay que comprobar que todos son distintos de NA, porque puede ser un ajuste extremo
    if(is.na(rest1)) rest1 <- FALSE
    if(is.na(rest2)) rest2 <- FALSE
    
    if(rest1 & rest2){
      condicionContinuar <- FALSE
    } else {
      i <- i+1
    }
    
    if(i > nrow(step1)){
      return(NULL)
    }
  }
  
  return(step1[ordenMinRSS[i],])
  
}

# Función que realiza el segundo paso del FMM. En este caso, es una implementación
# de la función objetivo a minimizar con Nelder-Mead. Calcula el RSS dados los parámetros
# y los datos. Los argumentos que recibe son:
#   param --> un vector que contiene, en orden, M, A, alpha, beta, omega
#   vData --> los datos para ser ajustados
#   timePoints --> instantes de tiempo donde se ajustan los datos
# FUNCION OCULTA
step2FMM <- function(param, vData, timePoints, omegaMax){
  
  n <- length(timePoints)
  
  #Ajuste FMM
  ffMob<-param[1] + param[2] * cos(param[4]+2*atan2(param[5]*sin((timePoints-param[3])/2),cos((timePoints-param[3])/2)))
  
  #La media de la suma de residuos al cuadrado
  RSS<-sum((ffMob - vData)^2)/n
  sigma <- sqrt(RSS*n/(n-5))
  
  #En caso de que se cumpla la condición de la amplitud, se devuelve el RSS.
  #En caso contrario, se devuelve infinito.
  maxi<-param[1]+param[2]
  mini<-param[1]-param[2]
  rest1 <- maxi <= max(vData)+0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
  rest2 <- mini >= min(vData)-0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
  
  #Otras condiciones de integridad que deben cumplirse
  rest3 <- param[2] > 0 #A > 0
  rest4 <- param[5] > 0 #omega > 0
  rest5 <- param[5] <= omegaMax #omega <= omegaMax
  if(rest1 & rest2 & rest3 & rest4 & rest5){
    return(RSS)
  }else{
    return(Inf)
  }
  
}


#Dibuja el FMM. Si components = F, dibuja las componentes centradas sin los datos.
plot_FMM <- function(object_FMM,vData,timePoints = seq(0,2*pi,length.out = length(vData)),components=F){
  
  if(components){
    
    nc <- length(object_FMM$alpha)
    predicted <- list()
    
    for(i in 1:nc){
      predicted[[i]] <- object_FMM$A[i]*cos(object_FMM$beta[i] + 2*atan(object_FMM$omega[i]*tan((timePoints-object_FMM$alpha[i])/2)))
    }
    
    mini <- min(vData)
    maxi <- max(vData)
    for(i in 1:nc){
      #predicted[[i]] <- predicted[[i]] - median(predicted[[i]])
      predicted[[i]] <- predicted[[i]] - predicted[[i]][1] + vData[1]
      mini <- min(mini,min(predicted[[i]]))
      maxi <- max(maxi,max(predicted[[i]]))
    }
    
    plot(timePoints,vData,ylim=c(mini,maxi),xlab="Time",ylab="Response",main="Components FMM",type="n")
    for(i in 1:nc){
      points(timePoints,predicted[[i]],type="l",lwd=2,col=i+1)
    }
    
  } else {
    
    plot(timePoints,vData,xlab="Time",ylab="Response",main="Adjusted FMM")
    points(timePoints,object_FMM$FMM,type="l",col=2,lwd=2)
    
  }
  
}


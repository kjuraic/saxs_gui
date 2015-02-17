
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# Komentar

library(shiny)



# START Functions ---------------------------------------------------------------
calcGuinier <- function(dat, bnd){
  interval <- bnd[1]:bnd[2]
  q <- dat[interval,1]
  intens <- dat[interval,2]
  y <- log(intens)
  x <- q^2
  fit <- lm(y ~ x)
  fitCoef <- summary(fit)$coefficients
  Rg  <- sqrt(-3 * as.numeric(fitCoef[2,1]))
  RgErr <- .5 * sqrt(-3./as.numeric(fitCoef[2,1])) * as.numeric(fitCoef[2,2])
  I0  <- exp(as.numeric(fitCoef[1,1]))
  I0Err  <- exp(as.numeric(fitCoef[1,1])) * as.numeric(fitCoef[1,2])
  qRg <- q[bnd[2]] * Rg
  out <- list(Rg = Rg, RgErr = RgErr, I0 = I0, I0Err = I0Err, qRg = qRg, fit = fit)
  return(out)
}

calcPorod <- function(dat, bnd, n = 4){
  interval <- bnd[1]:bnd[2]
  q <- dat[interval,1]
  intens <- dat[interval,2]
  y <- intens * q^n
  x <- q^n
  fit <- lm(y ~ x)
  k1 <- as.numeric(coef(fit)[2])
  k2 <- as.numeric(coef(fit)[1])
  k1Err <- as.numeric(coef(summary(fit))[2,2])
  k2Err <- as.numeric(coef(summary(fit))[1,2])
  out <- list(k1 = k1, k1Err = k1Err, k2 = k2, k2Err = k2Err, n =n, fit = fit)
  return(out)
}


# END Function ------------------------------------------------------------



readNika1D <- function(fname)
{
  dat <- read.table(fname,skip=23,col.names=c("q","intens","err","unknown"))
  return(dat)
}


shinyServer(function(input, output, clientData, session) {
  
  saxsData <- reactive({
    saxsFile <- input$saxsFile
    if (is.null(saxsFile)) {
      out <- NULL
    } else {
      #saxsFile <- "/home/kruno/TU/experiments/SAXS/141103_Aden/Pil1M/CSI_rot_stat/1D/cap_r_1M_C.dat"
      out <- readNika1D(saxsFile$datapath)
    }
    out
    
  })
  
  fitGun <- reactive({
    if (is.null(saxsData()))
      return(NULL)
    else {
      dat <- saxsData()
      fit <- calcGuinier(dat,input$rangeGun)
      fit
    }
  })
  
  fitPor <- reactive({
    if (is.null(saxsData()))
      return(NULL)
    else {
      dat <- saxsData()
      fit <- calcPorod(dat,input$rangePor,n = 4)
      fit
    }
  })
  
  output$saxsPlot <- renderPlot({
    if (is.null(saxsData()))
      return(NULL)
    else {
      dat <- saxsData()
      plot(dat[,1],dat[,2],log='y',col=2,pch=16,cex=.5,xlab="Q [nm-1]", ylab = "Intensity [a.u.]")
    }
  })
  
  output$plotGun <- renderPlot({
    if (is.null(saxsData()) | is.null(fitGun()))
      return(NULL)
    else {
      dat <- saxsData()
      plot(dat[,1]^2,log(dat[,2]),col='black', 
           xlim=c(dat[input$rangePlotGun[1],1]^2,dat[input$rangePlotGun[2],1]^2),
           xlab = 'Q^2', ylab = "log(intensity)")
      abline(v=dat[input$rangeGun[1],1]^2,col=4)
      abline(v=dat[input$rangeGun[2],1]^2,col=4)
      fit <- fitGun()
      abline(fit$fit,col='red')
    }
  })
  
  output$summaryGun <- renderTable({
    if (is.null(fitGun()))
      return(NULL)
    fit <- fitGun()
    data.frame(paramter = c("Rg","Rg Err","I0","I0 Err","Q*Rg"),
               value = c(fit$Rg,fit$RgErr,fit$I0,fit$I0Err, fit$qRg))
    #dat <- saxsData()
    #dat
  })
  
  output$plotPor <- renderPlot({
    if (is.null(saxsData()) )
      return(NULL)
    else {
      dat <- saxsData()
      plot(dat[,1],dat[,2],col='black',
           log='xy',
           xlim=c(dat[input$rangePlotPor[1],1],dat[input$rangePlotPor[2],1]),
           xlab = 'Q', ylab = "Intensity")
      abline(v=dat[input$rangePor[1],1],col=4)
      abline(v=dat[input$rangePor[2],1],col=4)
      fit <- fitPor()
      points(dat[,1],fit$k1+fit$k2/dat[,1]^fit$n, type='l', col = 'red')
    }
  })
  
  output$summaryPor <- renderTable({
    if (is.null(fitPor()))
      return(NULL)
    fit <- fitPor()
    data.frame(paramter = c("K1","K1 Err","K2","K2 Err","n"),
               value = c(fit$k1,fit$k1Err,fit$k2,fit$k2Err, fit$n))
    #dat <- saxsData()
    #dat
  })
  
  output$plotInv <- renderPlot({
    if (is.null(saxsData()) )
      return(NULL)
    else {
      dat <- saxsData()
      plot(dat[,1],dat[,1]^2*dat[,2],col='black',
           log='y',
           xlim=c(0,dat[input$rangePlotInv[2],1]),
           xlab = 'Q', ylab = "Intensity * Q^2")
      abline(v=dat[input$rangeInv[1],1],col=4)
      abline(v=dat[input$rangeInv[2],1],col=4)
      seqGun <- seq(0,dat[input$rangeInv[1],1],length.out = 1000)
      seqPor <- seq(dat[input$rangeInv[2],1],max(dat[,1]),length.out = 1000)
      fitGun <- fitGun()
      fitPor <- fitPor()
      points(seqGun,seqGun^2*fitGun$I0*exp(-1*fitGun$Rg^2*seqGun^2/3.), type='l', col = 'red')
      points(seqPor,seqPor^2*(fitPor$k2/seqPor^fitPor$n), type='l', col = 'red')
    }
  })
  
  observe({
    numRows <- nrow(saxsData())
    updateSliderInput(session, "rangeGun", max = numRows)
    updateSliderInput(session, "rangePor", max = numRows)
    updateSliderInput(session, "rangeInv", max = numRows)
    updateSliderInput(session, "rangePlotGun", max = numRows)
    updateSliderInput(session, "rangePlotPor", max = numRows)
    updateSliderInput(session, "rangePlotInv", max = numRows)
  })
})

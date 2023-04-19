# TrainingLoad
#
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#'This function estimates the Banister Training Impulse
#'
#' \code{bTRIM} Estimates the Banister Training Impulse
#'
#'
#' @param tcxTable Heart Rate data from tcx file
#' @param restHR Resting heart rate
#' @param maxHR Maximum Heart rate
#' @param gender Gender of the subject 0: Male 1: Female
#' @return tl Training Load
#' @export
#' @seealso \code{\link{eTRIM}}
#'
bTRIM <- function(tcxTable, restHR, maxHR, gender) {

  tl = 0;
  D = (as.numeric(tail(tcxTable$time,1))-as.numeric(tcxTable$time[1]))
  mhr = mean(tcxTable$heart_rate, na.rm = TRUE)

  deltaHR = (mhr-restHR)/(maxHR-restHR)

  b = 1.92
  if (gender==1) {
    b = 1.67
  }

  Y = exp(b*deltaHR)

  tl =  D*deltaHR*Y

  return(tl)
}


#'This function estimates the intensity distribution Training Impulse
#'
#' \code{bTRIM} Estimates the intensity distribution Training Impulse
#'
#'
#' @param tcxTable Heart Rate data from tcx file
#' @param restHR Resting heart rate
#' @param maxHR Maximum Heart rate
#' @param gender Gender of the subject 0: Male 1: Female
#' @return tl Training Load
#' @importFrom pracma interp1
#' @export
#' @seealso \code{\link{bTRIM}}
#'
eTRIM <- function(tcxTable, restHR, maxHR, gender) {

  #Time abs from start
  time = as.numeric(tcxTable$time)-as.numeric(tcxTable$time[1])
  #New time
  time2 = seq(time(1),tail(time,1))

  defaultW <- getOption("warn")
  options(warn = -1)
  #Lets make sure the data is evenly distributed in seconds
  yi = interp1(time, tcxTable$heart_rate, xi = time2, method = "nearest")

  options(warn = defaultW)

  #
  yirel = (yi-restHR)/(maxHR-restHR)*100

  #Categorizing in intensity zones
  res = hist(yirel, c(0,50,60,70,80,90,1000), plot = FALSE)

  #Estimate the training load
  tl = sum(res$counts[2:6]*c(1,2,3,4,5))

}

#'This function estimates the external training load from accelerometry
#'
#' \code{playerLoad} estimates the external training load from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @export
#' @seealso \code{\link{accelRate}}
#'
playerLoad <- function(ax3) {
  dax = diff(ax3$data$x)^2
  day = diff(ax3$data$y)^2
  daz = diff(ax3$data$z)^2

  tl = sum(sqrt((dax+day+daz)/100))

  return(tl)
}

#'This function estimates the external training load from accelerometry
#'
#' \code{playerLoad2} estimates the external training load from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @export
#' @seealso \code{\link{accelRate}}
#'
playerLoad2 <- function(ax3) {
  dax = diff(ax3$data$x)^2
  day = diff(ax3$data$y)^2
  daz = diff(ax3$data$z)^2

  tl = sum(sqrt(dax+day+daz)/100)

  return(tl)
}

#'This function estimates the external training load from accelerometry
#'
#' \code{playerLoadMAG} estimates the external training load from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @export
#' @seealso \code{\link{accelRate}}
#'
playerLoadMAG <- function(ax3) {
  dax = ax3$data$x^2
  day = ax3$data$y^2
  daz = ax3$data$z^2

  tl = sum(sqrt(dax+day+daz))

  return(tl)
}

#'This function estimates the external training load from accelerometry
#'
#' \code{playerLoadSUM} estimates the external training load from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @export
#' @seealso \code{\link{accelRate}}
#'
playerLoadSUM <- function(ax3) {
  dax = sqrt(ax3$data$x^2)
  day = sqrt(ax3$data$y^2)
  daz = sqrt(ax3$data$z^2)

  tl = sum(dax+day+daz)

  return(tl)
}

#'This function estimates the external training load from accelerometry
#'
#' \code{impulseLoad} estimates the external training load from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @export
#' @seealso \code{\link{accelRate}}
#'
impulseLoad <- function(ax3) {
  dax = ax3$data$x^2
  day = ax3$data$y^2
  daz = ax3$data$z^2

  vm = sqrt(dax+day+daz)

  #removing all acceleration below walking intensity
  vm = vm[which(vm>0.1)]

  tl = sum(vm/9.8067)

  return(tl)
}

#'This function estimates the external training load from accelerometry
#'
#' \code{accelRate} estimates the external training load from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @export
#' @seealso \code{\link{playerLoad}}
#'
accelRate <- function(ax3) {
  vm = sqrt(rowSums(ax3$data[,2:4]^2))

  tl = sum(abs(diff(vm)))

  return(tl)
}

#'This function estimates the external training load from accelerometry
#'
#' \code{velocityLoad} estimates the external training load from accelerometry using ENMO and filtering the data
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @return tl Training Load
#' @importFrom signal butter filtfilt
#' @importFrom caTools runmean
#' @export
#' @seealso \code{\link{playerLoad}}
#'
velocityLoad <- function(ax3) {
  vm = sqrt(rowSums(ax3$data[,2:4]^2))

  wnds = ax3$header$frequency * 3
  #Lets remove gravity and high frequencies
  #Might want to reduce the high pass filter to 2Hz
  bf = signal::butter(4,c(0.01/(ax3$header$frequency/2),10/(ax3$header$frequency/2)), type = "pass")

  vmf = caTools::runmean(abs(signal::filtfilt(bf,vm-1)),wnds,align = 'left')

  #Lets estimate MAD in 3 seconds epoch

  vmf = vmf[seq(1,length(vmf),wnds)]
  tl = sum(vmf)

  return(tl)
}

#'This function estimates Mean Absolut deviation from accelerometry
#'
#' \code{mad} estimates MAD from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @param epoch Epoch length in seconds - Default is 5 seconds
#' @importFrom caTools runmad
#' @return mad data
#' @export
#' @seealso \code{\link{playerLoad}}
#'
mad <- function(ax3,epoch = 5) {
  vm = sqrt(rowSums(ax3$data[,2:4]^2))

  wnds = ax3$header$frequency * epoch + 1

  vmf = caTools::runmad(vm,wnds,align = 'left')

  #Lets estimate MAD in 3 seconds epoch

  vmf = vmf[seq(1,length(vmf),wnds)]

  return(vmf)
}

#'This function estimates Euclidean Norm minus one from accelerometry
#'
#' \code{enmo} estimates ENMO from accelerometry
#'
#'
#' @param ax3 AX3 Accelerometry data
#' @param epoch Epoch length in seconds - Default is 5 seconds
#' @importFrom caTools runmean
#' @return enmo data
#' @export
#' @seealso \code{\link{playerLoad}}
#'
enmo <- function(ax3,epoch = 5) {

  vm = sqrt(rowSums(ax3$data[,2:4]^2))-1

  vm[which(vm<0)] = 0

  wnds = ax3$header$frequency * epoch + 1

  vmf = caTools::runmean(vm,wnds,align = 'left')

  #Lets estimate ENMO in 5 seconds epoch

  vmf = vmf[seq(1,length(vmf),wnds)]

  return(vmf)
}

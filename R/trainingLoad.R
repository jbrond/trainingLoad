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
#' @param gender Gender of the subject
#' @return tl Training Load
#' @export
#' @seealso \code{\link{eTRIM}}
#'
bTRIM <- function(tcxTable, restHR, maxHR, gender) {

  tl = 0;
  D = as.numeric(tcxTable$time[362]-tcxTable$time[1])
  mhr = mean(run2$heart_rate)

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
#' @param gender Gender of the subject
#' @return tl Training Load
#' @export
#' @seealso \code{\link{bTRIM}}
#'
eTRIM <- function(tcxTable, restHR, maxHR, gender) {

  #Time abs from start
  time = as.numeric(run2$time-run2$time[1])
  #New time
  time2 = seq(time(1),tail(time,1))

  #Lets make sure the data is evenly distributed in seconds
  yi = interp1(time, run2$heart_rate, xi = time2)

  #
  yirel = (yi-restHR)/(maxHR-restHR)*100

  #Categorizing in intensity zones
  res = hist(yirel, c(0,50,60,70,80,90,1000))

  #Estimate the training load
  tl = sum(res$counts[2:6]/60*c(1,2,3,4,5))

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

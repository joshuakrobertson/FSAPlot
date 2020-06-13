#' FSAPlot: read_fsa
#'
#' @description
#' A function to read in .fsa files from an ABI analyser and produce a data.frame with relative fluoresences according to channel and basepair length. 
#' @details The default ladder for fragment length calibration is GS500.
#' @param File A path to a .fsa file to be opened.
#' @param Lad_Chan The loading channel in which the ladder may be found. Default is 4.
#' @param Sig_Chan Loading channels wherein sample reads lie. Default is 1 to 3.
#' @param Range Basepair range wherein peaks will be extracted. Default is 100 to 500.
#' @keywords fsa, microsatellites, ABI
#' @import seqinr caTools dplyr
#' @examples
#' source = "/home/fish332_microsat.fsa"
#' Peak_Dat = read_fsa(File = source, Lad_Chan = 4, Range = c(50:80))
#' @export

read_fsa <- function(File, Lad_Chan = 4, Sig_Chan = 1:3, Range = c(100:500)) {

  require("seqinr")
  require('caTools')
  require('dplyr')
  
  set.ladder <- function(lad.dat, ladder, SNR, ladder.check = NULL) {
    n <- 10 * length(ladder)
    bp <- rep(NA, length(lad.dat))
    Peak <- cumsum(abs(c(0, diff(lad.dat > quantile(lad.dat, 1 - n / length(lad.dat))))))
    PeakInd <- Peak %% 2 == 1
    Index <-
      seq_len(length(lad.dat))[PeakInd][lad.dat[PeakInd] == 
                                          ave(lad.dat[PeakInd],
                                              Peak[PeakInd], FUN = max)]    
    if(max(lad.dat) > SNR){
      primer.peak <-
        tail(which(lad.dat[Index] == max(lad.dat)), 1) 
      Index <- Index[-1 * (1:primer.peak)] # exclude the primer peak and
      # everything smaller 
    }
    
    while(length(Index) < length(ladder)){
      n <- n + 10
      Peak <- cumsum(abs(c(0, diff(lad.dat > quantile(lad.dat, 1 - n /
                                                        length(lad.dat))))))
      PeakInd <- Peak %% 2 == 1
      Index <-
        seq_len(length(lad.dat))[PeakInd][lad.dat[PeakInd]
                                          == ave(lad.dat[PeakInd],
                                                 Peak[PeakInd], FUN = max)]  
      if(max(lad.dat) > SNR){
        primer.peak <- which(lad.dat[Index] == max(lad.dat))
        Index <- Index[-1 * (1:primer.peak)] # exclude the primer peak and everything smaller
      }
    }
    if(length(Index) > length(ladder)){
      
      if(length(Index) - length(ladder) == 1){
        Index <- Index[-which.max(sapply(seq_along(Index), function(i){
          summary(lm(ladder ~ stats::poly(Index[-i], 2)))$r.squared
        }))]
      } else {
        if (length(Index) - length(ladder) >= 10){
            return(paste("To many ladder peaks"))
        } else { 
          toTry <- combn(length(Index), length(Index) - length(ladder))
        }
                
        if(ncol(toTry) > 4000){
          toTry <- toTry[, sample(ncol(toTry), 4000)]
        }
        
        sap <- sapply(seq_len(ncol(toTry)),
                      function(i){
                        summary(
                          lm(ladder ~ stats::poly(Index[-toTry[, i]], 2))
                        )$r.squared
                      })
        Index <- Index[-toTry[, which.max(sap)]]
      }
    }
    
    bp[Index] <- ladder
    
    if(! is.null(ladder.check)) {
      lad.ind <- which(ladder == ladder.check)
      ladder <- ladder[-lad.ind]
      Index <- Index[-lad.ind]
    }
    
    fit.value <- summary(lm(ladder ~ stats::poly(Index, 2)))$r.squared
    message("  ladder fit r2: ", round(fit.value, 5))
    return(list(bp = bp, val = fit.value))
  }
  
  File <- File
  EPG <- read.abif(File)
  tag <- gsub("\\.fsa", "", basename(File))
  baseline_width = 51
  smoothing = 3
  ladder.check = 250
  ladder = c(35, 50, 75, 100, 139, 150, 160, 200, 250,
             300, 340, 350, 400, 450, 490, 500)
  SNR = 5000
  
  # Pulling ladder data

  lad.dat <- EPG$Data[[paste("DATA.", Lad_Chan, sep = "")]]
  lad.dat <- lad.dat - caTools::runmin(lad.dat, baseline_width)
  lad.dat <- caTools::runmean(lad.dat, k = smoothing)
  
  scans <- data.frame(standard = lad.dat)
  tmp <- set.ladder(lad.dat, ladder, SNR, ladder.check)
  if (is.character(tmp)) {
    return("Too many peaks in ladder")
  } else {
  scans$bp <- tmp$bp
  lad_mat <- as.data.frame(cbind(bp = scans$bp[!is.na(scans$bp)], time = which(!is.na(scans$bp))))
  pred_mod = lm(bp ~ time, data = lad_mat)
  scans = scans %>% mutate("time" = row_number())
  scans$bp_adj = round(predict(pred_mod, newdata = scans, na.action = na.pass, type = "response"))
  scans = scans %>% dplyr::select(-bp, "bp" = bp_adj)
  
  EPG_Dat <- data.frame("bp" = scans$bp,
                        "Peak" = scans$standard,
                        "Channel" = "Standard")
  EPG_Dat$Channel = as.character(EPG_Dat$Channel)
  
  for (i in Sig_Chan) {
    chan.dat <- EPG$Data[[paste("DATA.", Sig_Chan[i], sep = "")]]
    chan.dat <- chan.dat - caTools::runmin(chan.dat, baseline_width)
    chan.dat <- caTools::runmean(chan.dat, k = smoothing)
    chan.dat = as.data.frame(chan.dat) %>% mutate("time" = row_number())
    chan.dat$bp = round(predict(pred_mod, newdata = chan.dat, na.action = na.pass, type = "response"))
    chan.dat = chan.dat %>% dplyr::rename("standard" = chan.dat)
    
    EPG_Dat <- rbind(EPG_Dat, data.frame("bp" = chan.dat$bp,
                          "Peak" = chan.dat$standard,
                          "Channel" = as.character(i)))
  }
  EPG_Dat <- subset(EPG_Dat, bp %in% Range)

  # Appending file name

  comment(EPG_Dat) = paste(basename(File))
  return(EPG_Dat)
  }
}
#' FSAPlot: plot_fsa
#'
#' @description
#' A function to plot fsa files produce from an ABI analyser and rendered into a data.frame using the read_fsa function. Plots are produced using ggplot2 by Hadley Wickham.
#' @details The default ladder for fragment length calibration is GS500. Note that only data.frames produced by read_fsa may be plotted with this function.
#' @param fsa A data.frame containing microsatellite fragment lengths and relative fluoresence produced using read_fsa.
#' @param dye Fluorescent dye bound to loci that one wishes to view. Only 6FAM, HEX, NED, or ROX are supported. This call can either be a single dye or a vector of multiple dyes. No default is provided.
#' @param title An argument defining whether the title of the fsa file is to be given at the head of a given plot. This argument must by TRUE or FALSE.
#' @param x_limits A vector of length 2 defining the basepair range to be viewed. Default is 100 to 500
#' @keywords fsa, microsatellites, ABI
#' @import ggplot2 dplyr
#' @examples
#' source = "/home/fish332_microsat.fsa"
#' Peak_Dat = read_fsa(File = source, Lad_Chan = 4, Range = c(50:80))
#' plot_fsa(Peak_Dat, dye = c("6FAM", "NED"), title = TRUE, x_limits = c(50,80))
#' @export

plot_fsa <- function(fsa, dye, title = TRUE, x_limits){
  require('ggplot2')
  require('dplyr')

  Legend = data.frame("Dye" = c("6FAM", "HEX", "NED", "ROX"), 
      "Channel" = c("1","2","3","Standard"), "Colour" = 
      c("royalblue3", "mediumseagreen", "black", "firebrick"))

  if (length(dye) == 1){
    Plot = fsa %>% filter(Channel == Legend[which(Legend$Dye == dye),2]) %>%
    ggplot(aes(x = bp, y = Peak)) + geom_line(size = 0.75, 
    colour = Legend[which(Legend$Dye == dye),3]) + theme_bw() + ylim(c(-20, max(subset(fsa, Channel == Legend[which(Legend$Dye == dye),2])$Peak))) + scale_x_continuous(limits = c(x_limits[1], x_limits[2]), breaks = seq(x_limits[1],x_limits[2],1)) + ylab("Signal\n(Relative Fluorescence)") + theme(axis.text.x = element_text(angle = 90))
  } else if (length(dye) > 1){
    Drop = Legend$Channel[c(which(!(Legend$Dye %in% dye)))]
    Plot = fsa %>% filter(!(Channel %in% Drop)) %>% 
    ggplot(aes(x = bp, y = Peak, colour = Channel)) + geom_line(size = 0.75) + theme_bw() + ylim(c(-20, max(fsa %>% filter(!(Channel %in% Drop)) %>% dplyr::select(Peak)))) + scale_x_continuous(limits = c(x_limits[1], x_limits[2]), breaks = seq(x_limits[1],x_limits[2],1)) + ylab("Signal\n(Relative Fluorescence)") + theme(axis.text.x = element_text(angle = 90)) + 
    scale_colour_manual(values = c(Legend$Colour[c(which(Legend$Channel %in% unique(as.data.frame(fsa %>% filter(!(Channel %in% Drop)))$Channel)))]))  
    }

  if (title == TRUE){
   Plot = Plot + 
      ggtitle(paste(attributes(fsa)$comment)) + 
      theme(plot.title = element_text(size = 8))
  } else if (title == FALSE){
    Plot = Plot
  } else {
    return(paste("Title argument must be TRUE or FALSE"))
  }

  return(Plot)
}
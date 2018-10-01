#' Plot molecules 
#' 
#' Plot molecules in R plot window instead of separate Java window
#' 
#' @param molecule an RCDK molecule 
#' @param width an integer width of the molecule
#' @param height an integer height of the molecule
#' @param marg margin for all side of the plot (default: 0)
#' @param main a string the main title of the figure
#' @return None
#' 
#' @details Code taken from: 
#' http://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/
#' 
#' @examples 
#' tmp <- rcdk::parse.smiles("C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)(F)F")
#' rcdkplot(tmp[[1]], width=300, height=300, main="Gemcitabine")
#'
#' @concept rcellminer
#' @export 
#' 
#' @importFrom rcdk view.image.2d get.depictor 
rcdkplot <- function(molecule,width=300,height=300,marg=0,main='') {
	par(mar=c(marg,marg,2.2,marg)) # set margins to zero since this isn't a real plot; top is 1.5 to accomodate title
  depictor <- get.depictor(width=width, height=height) # Necessary as of Summer 2017
	temp1 <- view.image.2d(molecule,depictor) # get Java representation into an image matrix. set number of pixels you want horiz and vertical
	plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='',main=main) # create an empty plot
	rasterImage(temp1,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
}
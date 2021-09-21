chromatogram_data <- data.frame(cbind(test_spectra[1,]@rtime, test_spectra[1,]@intensity))

test_chromatogram <- DrawChromatogram(time = chromatogram_data$X1, 
                                      intensity = chromatogram_data$X2, 
                                      range = list(start = 1.8, stop = 2.2),
                                      color = "blue",
                                      xlab = "retention time", 
                                      ylab = "intensity")


smoothingSpline = smooth.spline(chromatogram_data$X1, chromatogram_data$X2, spar=0.5)
plot(chromatogram_data$X1,chromatogram_data$X2)
lines(smoothingSpline)


peaksWithCentWave(peakwidth = c(2, 100),
                  snthresh = 1,
                  prefilter = c(3, 5),
                  integrate = 1,
                  fitgauss = FALSE,
                  noise = 0,
  verboseColumns = FALSE,
  firstBaselineCheck = TRUE,
  extendLengthMSW = FALSE)

pks <- peaksWithMatchedFilter(rt= chromatogram_data$X1,
                       int=chromatogram_data$X2,
  fwhm = 50,
  sigma = 50/2.3548,
  max = 20,
  snthresh = 1
);pks

plot(chromatogram_data$X1, chromatogram_data$X2, type = "h")
rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"], ybottom = c(0, 0),
     ytop = pks[, "maxo"], border = "red")



idx_apex <- which(chromatogram_data$X2 == max(chromatogram_data$X2))



chromatogram_list <- list()

print(Sys.time())
for(idx_lipid in 1:length(test_spectra)){
  
chromatogram_data <- data.frame(cbind(test_spectra[idx_lipid,]@rtime, test_spectra[idx_lipid,]@intensity))

# plot(x = chromatogram_data$X1, 
#      y = chromatogram_data$X2, 
#      type = "l", 
#      col = "gray")
# #pk.pos <- findpeaks(chromatogram_data$X2, span = 50)

if(max(chromatogram_data$X2) > 0){
span_number <- nrow(chromatogram_data)-1
pk.pos = NULL
while(length(pk.pos) != 1 & span_number > 1){
  pk.pos <- findpeaks(chromatogram_data$X2, span = span_number)
  span_number <- span_number - 1
}


# abline(v = chromatogram_data$X1[pk.pos], col = 4)

pks <- fitpeaks(y = chromatogram_data$X2, 
                pos = pk.pos)
# apply(pks, 1,
#       function(pkmodel) {
#       lines(chromatogram_data$X1,
#               dnorm(1:length(chromatogram_data$X1), 
#                     pkmodel["rt"], 
#                     pkmodel["sd"]) * pkmodel["area"],
#               col = 2)
#         invisible()
#       }) 
# 
# grid.echo()
# fig_1 <- grid.grab()
 }

#idx_lipid <- idx_lipid+1


chromatogram_list[idx_lipid] <- list(pks)
                                     #grid.arrange(fig_1))
}

print(Sys.time())









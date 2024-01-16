

library(enmSdmX)
library(scales)
library(terra)

ag <- vect('C:/Ecology/Drive/Research/Andropogon/Andropogon/data/occurrence_data/andropogon_gerardi_occurrences_with_environment.gpkg')

xy <- read.csv('C:/!Scratch/Spatial Sampling Site Information.csv')
traits <- read.csv('C:/Ecology/Drive/Research/Andropogon/Andropogon/data/data_from_loretta/jack_sytsma_2023_11_09a/Spatial Sampling Phenotype Data 2023.csv')

traits <- aggregate(traits, by = list(traits$Site), FUN = mean, na.rm = TRUE)
names(traits)[1] <- 'site'


traits$lat <- traits$long <- NA_real_
traits$lat <- xy$lat[match(xy$site, traits$site)] 
traits$long <- xy$long[match(xy$site, traits$site)] 

traits <- vect(traits, c('long', 'lat'), crs = getCRS('wgs84'))

trs <- c('Nitrogen.Content....', 'C.N.Ratio', 'Height..cm.', 'Blade.Width..cm.', 'Leaf.Thickness..mm.', 'SPAD', 'Canopy.Diameter..cm.', 'Midday.Water.Potential..bars.', 'Photosynthetic.Rate', 'Stomatal.Conductance', 'Internal.Carbon.Dioxide', 'Vegetative.Biomass..g.', 'Reproductive.Biomass..g.', 'Total.Biomass..g..2', 'Veg.Rep.Ratio', 'Transpiration.Rate')

ag <- ag[!is.na(ag$any_ag_quality1to3)]

col <- ag$any_ag_quality1to3
col <- log10(col + 1)
col <- col - min(col, na.rm = TRUE)
col <- col / max(col, na.rm = TRUE)

for (tr in trs) {

	say(tr)

	png(paste0('C:/!Scratch/AG/', tr, '.png'), width = 2000, height = 1200)

	plot(ag, col = alpha('forestgreen', col), border = 'gray60', main = tr, cex.main = 3)
	
	cex <- traits[ , tr]
	cex <- as.data.frame(cex)
	cex <- cex[ , 1]
	cex <- cex - min(cex, na.rm = TRUE)
	cex <- cex / max(cex, na.rm = TRUE)
	cex <- cex + 0.5
	cex <- 4 * cex
	
	plot(traits, cex = cex, add = TRUE)
	
	dev.off()

}

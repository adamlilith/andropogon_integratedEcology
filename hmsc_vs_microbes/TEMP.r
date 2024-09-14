# source('C:/Ecology/R/andropogon_integratedEcology/hmsc_vs_microbes/TEMP.r')

	# plotting extent
	maps_mean <- rast(paste0(output_subfolder, '/prediction_rasters_mean.tif'))

	extent <- c(-754404, 2156962, -2232910, 836303)
	extent <- ext(extent)
	extent <- as.polygons(extent, crs = crs(maps_mean))
	extent <- buffer(extent, 100000)
	maps_mean <- crop(maps_mean, extent)
	
	# ...but actually plot a smaller area (looks nicer)
	extent <- buffer(extent, -100000)
	extent <- ext(extent)
	extent_for_plot <- as.vector(extent)

		taxon <- as.character(taxa$taxon[i])
		say(taxon)

		# get map and project		
		this_map <- maps_mean[[taxon]]

		# to define colors, get min/max abundance across sites and the map
		this_abundance <- as.data.frame(sites_abundances[ , taxon])[ , 1, drop = TRUE]
		this_abundance <- log10(this_abundance + 1)
		map_min <- minmax(this_map)[1, ]
		map_max <- minmax(this_map)[2, ]
		min_abund <- min(map_min, this_abundance)
		
		# remove extreme abundances
		max_abund <- globalx(this_map, quantile, probs = 0.99, na.rm = TRUE)
		this_map <- app(this_map, fun = function(x, max_abund) ifelse(x > max_abund, 0.999 * max_abund, x), max_abund = max_abund)

		if (is.infinite(min_abund)) {
		
			neg_inf_abund_masked <- app(this_map, fun = function(x) ifelse(is.infinite(x), max_abund, x))
			min_abund <- globalx(neg_inf_abund_masked, min)
			this_map <- app(this_map, fun = function(x, min_abund) ifelse(x < min_abund, min_abund, x), min_abund = min_abund)
			
		}

		abund_seq <- seq(min_abund, max_abund, length.out = 100)
		this_abundance_scaled <- (this_abundance - min_abund) / (max_abund - min_abund)
		this_abundance_scaled <- round(100 * this_abundance_scaled)

		pallette <- viridis(n = 100)
		site_colors <- pallette[this_abundance_scaled]

		plot(nam, col = 'gray30', ext = extent_for_plot, lwd = 0.1, main = taxon)
		plot(this_map, range = c(min_abund, max_abund), add = TRUE)
		plot(nam, border = 'gray20', lwd = 0.1, add = TRUE)

		plot(sites_abundances, pch = 21, bg = site_colors, cex = 1.6, add = TRUE)
	
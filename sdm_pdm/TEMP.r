# source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/TEMP.r')

		focal_pzero <- unlist(pred_vect_conus[[paste0('pzero_occs_county_', period)]])
		delta_vect$delta_pzero <- delta_vect$pzero_occs_county_sq - focal_pzero
		# delta_vect$delta_pzero <- focal_pzero / delta_vect$pzero_occs_county_sq

		delta_vect_display <- crop(delta_vect, extent_display)

		if (period == 'thirties') {
			title <- 'B'
		} else if (period == 'fifties' ) {
			title <- 'D'
		}
		title <- paste0(title, ') Change in probability of absence: Present vs. ', decade)

		map_change_pzero <- ggplot() +
			layer_spatial(nam, fill = 'gainsboro', color = NA) +
			layer_spatial(
				delta_vect_display,
				aes(fill = delta_pzero),
				color = NA
			) +
			scale_fill_gradient2(
				name = paste0('present - ', decade),
				low = '#b2182b',
				mid = 'beige',
				high = '#2166ac',
				midpoint = 0,
				# trans = 'log10',
				limits = rep_limits_pzero
			) +
			layer_spatial(nam, color = 'gray40', fill = NA) +
			layer_spatial(site_vect, pch = 3) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

print(map_change_pzero)

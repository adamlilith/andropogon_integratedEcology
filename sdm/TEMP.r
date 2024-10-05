# source('E:/Adam/R/andropogon_integratedEcology/sdm/TEMP.r')
	fill_scale <- c('[0 - 0.25)' = 'gray90', '[0.25 - 0.5)' = alpha('forestgreen', 0.125), '[0.5 - 0.75)' = alpha('forestgreen', 0.25), '[0.75 - 0.95)' = alpha('forestgreen', 0.5), '[0.95-1]' = 'forestgreen')

	map <- ggplot() +
			layer_spatial(this_ag_vect, aes(fill = quant_col), color = NA) +
			# layer_spatial(nam, color = 'gray50', fill = NA, linewidth = 0.1) +
			scale_fill_manual(
				name = 'Quantile\n of λ',
				values = fill_scale
			) +
			# layer_spatial(cents_with_ag, pch = 16, col = 'black', size = 0.2) +
			guides(fill = guide_legend(reverse = TRUE)) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(expression('Future equilibrial distribution of ' * italic('Andropogon gerardi')), subtitle = pretty_title)

	ggsave(plot = map, filename = paste0(out_dir, '/nonintegrated_sdm_lambda_', fut, '.png'), width = 12, height = 10, dpi = 600)

	print(map)

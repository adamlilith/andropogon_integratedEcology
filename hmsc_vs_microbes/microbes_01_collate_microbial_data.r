### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles abundance data for microbes associated with the rhizobiome of Andropogon gerardi and bulk soil samples taken in the plants' vicinities.
###
### source('C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_01_collate_microbial_data.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_01_collate_microbial_data.r')
###
### CONTENTS ###
### setup ###
### collate rhizobiome & bulk soil microbial abundances

#############
### setup ###
#############

rm(list = ls())

# drive <- 'C:/Ecology/'
drive <- 'E:/Adam/'

setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

library(data.table)
library(read_xl)
library(omnibus)

say('###########################################################')
say('### collate rhizobiome & bulk soil microbial abundances ###')
say('###########################################################')

	### collate rhizobiome & bulk soil microbial abundances
	#######################################################

	# We want a matrix with cell values representing abundances
	# each column representing a taxon
	# each row a site
	# column 1: site
	# column 2: sample ID

	rhizo <- fread('./data/data_from_sonny/ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Rhizomicrobiome_Level_7.csv')
	bulk <- fread('./data/data_from_sonny/DATA_Microbiome_ASVs_Bulk_Soil_Level_7_R-Friendly.csv')

	# remove control columns/rows
	bulk <- bulk[ , `NC-11292023-16s` := NULL]
	rhizo <- rhizo[index != 'NC_S8']

	bulk <- bulk[Kingdom != ''] # reads empty rows for some reason

	microbes <- cbind(
		data.table(
			site = rhizo$location,
			sample = rhizo$plant,
			type = 'rhizobiome'
		),
		rhizo
	)
	colnames(microbes)[colnames(microbes) == 'index'] <- 'site_sample'
	microbes[ , c('location', 'plant') := NULL]

	# add bulk samples to rhizobiome
	bulk_meta_cols <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Samples')
	sample_cols <- colnames(bulk)[colnames(bulk) %notin% bulk_meta_cols]

	microbes_meta_cols <- c('site', 'sample', 'type', 'site_sample')
	microbes_abund_cols <- which(colnames(microbes) %notin% microbes_meta_cols)
	for (i in seq_along(sample_cols))
	
		site_sample <- sample_cols[i]
		site <- substr(site_sample, 1, 3)
		sample <- substr(site_sample, 5, nchar(site_sample))
	
		# grab existing row and use it as a template
		new_row <- microbes[1]
		new_row$site <- site
		new_row$sample <- sample
		new_row$type <- 'bulk'
		new_row$site_sample <- site_sample
		new_row[ , (microbes_abund_cols) := 0]

		# populate cells with abundances	
		bulk_taxa <- apply(bulk[ , c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')], 1, paste, collapse = ';')
		microbe_taxa <- colnames(microbes[ , ..microbes_abund_cols])
		for (j in seq_along(bulk_taxa)) {

			bulk_taxon <- bulk_taxa[j]
			
			# this bulk taxon is already in microbe data
			if (bulk_taxon %in% microbe_taxa) {
				stop()
				this_col <- which(microbe_taxa == bulk_taxon)
				new_row[1, this_col] <- bulk[XYZ]

			# this bulk taxon is not yet in microbe data
			} else {

				# add taxon to microbe data
				microbes[ , DUMMY := 0]
				colnames(microbes)[ncol(microbes)] <- bulk_taxon

				# add taxon to new data
				new_row[1, ..bulk_taxon := bulk[j, ..site_sample]]
				colnames(new_row)[ncol(new_row)] <- bulk_taxon

				microbe_taxa <- c(microbe_taxa, bulk_taxon)

			}

			microbes <- rbind(microbes, new_row)

		}

	fwrite(microbes, './microbes_rhizobiome_bulk.csv', row.names = FALSE)

	# ### collate environmental data
	# ##############################

	# extract_env <- fread('./data/data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024_extracted_environment.csv')
	# env <- extract_env[ , c('my.id', 'site', 'x-coordinate', 'y-coordinate', 'summer_tmean_C', 'annual_ppt_mm', 'sampling_ppt_mm', 'ppt_cv')]

	# colnames(env)[colnames(env) == 'my.id'] <- 'site_sample'
	# colnames(env)[colnames(env) == 'x-coordinate'] <- 'longitude'
	# colnames(env)[colnames(env) == 'y-coordinate'] <- 'latitude'

	# ### soil texture measured in situ (from soil taken at sites)
	# soil_in_situ <- read_xl('./data/data_from_sonny_and_loretta/site_soil_data_Order_52761.xlsx', sheet = 'Test Results For Order', skip = 9)
	# soil_in_situ_site_names <- soil_in_situ$`Sample Name`
	# soil_in_situ_site_names <- strsplit(soil_in_situ_site_names, split = ' ')
	# soil_in_situ_site_names <- do.call(rbind, soil_in_situ_site_names)
	# soil_in_situ_site_names <- soil_in_situ_site_names[ , 2, drop = TRUE]
	# soil_in_situ_site_names <- trimws(soil_in_situ_site_names)
	# soil_in_situ_site_names <- gsub(soil_in_situ_site_names, pattern = '-', replacement = '')

	# env[ c('sand_in_situ', 'silt_in_situ', 'clay_in_situ') := NA_real_]
	# for (i in 1:nrow(env)) {

	# 	site <- env$site[i]
	# 	soil_row_index <- which(soil_in_situ_site_names == site)

	# 	env$sand_in_situ[i] <- soil_in_situ$`Sand %`[soil_row_index] / 100
	# 	env$silt_in_situ[i] <- soil_in_situ$`Silt %`[soil_row_index] / 100
	# 	env$clay_in_situ[i] <- soil_in_situ$`Clay %`[soil_row_index] / 100

	# }

	# ### soil chemistry measured in situ
	# soil_chem <- read_xl('./data/data_from_sonny/DATA_Biogeographical_Sampling.xlsx', sheet = 'Soil Chemistry', skip = 2)
	# soil_chem_sites <- soil_chem$`Sample Name`
	# soil_chem_sites <- gsub(soil_chem_sites, pattern = '\\*', replacement = '')
	# soil_chem_sites <- substr(soil_chem_sites, 1, 3)

	# env[ , c('ph_in_situ', 'k_in_situ', 'n_in_situ', 'c_in_situ') := NA_real_]
	# for (i in 1:nrow(env)) {
	
	# 	site <- env$site[i]
	# 	row_index <- which(soil_chem_sites == site)

	# 	env$ph_in_situ[i] <- soil_chem$`pH`[row_index]
	# 	env$k_in_situ[i] <- soil_chem$`K ppm`[row_index]
	# 	env$n_in_situ[i] <- soil_chem$`Total N %`[row_index]
	# 	env$c_in_situ[i] <- soil_chem$`Total C %`[row_index]
	
	# }





say('DONE!', level = 1, deco = '!')

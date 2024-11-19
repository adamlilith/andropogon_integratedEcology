### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles abundance data for microbes associated with the rhizobiome of Andropogon gerardi and bulk soil samples taken in the plants' vicinities.
###
### source('C:/Ecology/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_01_collate_microbial_and_environmental_data.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_01_collate_microbial_and_environmental_data.r')
###
### CONTENTS ###
### setup ###
### compile microbial abundance, environmental data, and study design into forms needed by HMSC ###
### construct raster template and data table of predictors for making spatial predictions ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'

	library(data.table)
	library(enmSdmX)
	library(omnibus)
	library(readxl)
	library(terra)

	devtools::load_all(paste0(drive, '/R/airUpThere'))

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	predictor_names <- c('aridity', 'bio7', 'bio12', 'bio15', 'soc', 'ph', 'sand', 'silt')
	MASTER_climate_predictor_names <- c('aridity', 'bio7', 'bio12', 'bio15')

say('###################################################################################################')
say('### compile microbial abundance, environmental data, and study design into forms needed by HMSC ###')
say('###################################################################################################')

	### user-defined
	################
	
		# collate at this taxonomic level (only "class" and "phylum" are allowed)
		# collapse_taxon_level <- 'class'
		collapse_taxon_level <- 'phylum'
	
		# folder with PRISM monthly values
		# pr_dir <- 'D:/ecology/PRISM/working/an81'
		pr_dir <- 'E:/ecology/PRISM/working/an81'
		# pr_dir <- 'F:/ecology/PRISM/working/an81'

		# names of control sites in rhizobiome and bulk soil
		rhizo_control_site <- 'NC_S8'
		bulk_control_site <- 'NC-11292023-16s_S92'
		
	# Want:
	# * An "abundances" matrix with each row = sample, each column = taxon abundance.
	# 	* Need to collapse the abundances by Domain/Phylum/Class because some are represented >1 time.
	#	* Need to remove "control" samples (positive and negative controls.)
	# * Three "taxonomy" matrices (rhizobiome, bulk soil, and combined) with each row = one taxon, where there are at least three columns, named "Domain", "Phylum", and "Class". All taxa in "abundances" must be represented. No duplicates allowed in these taxonomy matrices.
	# * An "environment" frame with one row per site/sample and with environmental data, including a flag indicating if a sample was from the bulk or rhizomicrobiome soil.
	#' * A "study design" matrix with factors for each site and each sample within a site.

	### save to
	###########
	
		out_dir <- paste0('./data_from_sonny/collated_for_hmsc_collapsed_to_', collapse_taxon_level)
		dirCreate(out_dir)

	### inputs
	##########
		
		raw_abund_sites_rhizo <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/brooke_2024_10_04/level-7_Final_Root.csv')
		
		raw_abund_sites_bulk <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/brooke_2024_10_04/level-7_Final_Soil.csv')
		
		soil_chemistry <- read_xlsx('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Biogeographical_Sampling.xlsx', sheet = 'Soil Chemistry', skip = 1)
		soil_chemistry <- as.data.table(soil_chemistry)
		
		sampling_dates <- read_xlsx('./data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024.xlsx', sheet = 'Sheet1')
		sampling_dates <- as.data.table(sampling_dates)

	### RHIZOBIOME & BULK: remove controls
	######################################
	
		say('RHIZOBIOME & BULK: remove controls')

		raw_abund_sites_rhizo <- raw_abund_sites_rhizo[index != rhizo_control_site]
		raw_abund_sites_bulk <- raw_abund_sites_bulk[index != bulk_control_site]

	### RHIZOBIOME: clean taxon names
	#################################
	
		say('RHIZOBIOME: clean taxon names')
		
		cols <- colnames(raw_abund_sites_rhizo)
		taxa <- cols[cols %notin% c('index', 'location', 'plant')]
		
		taxa_rhizo <- data.table()
		for (i in seq_along(taxa)) {
		
			taxon <- taxa[[i]]
			taxon_orig <- taxon
			taxon <- trimws(taxon)
			taxon <- strsplit(taxon, split = ';')[[1]]

			taxon <- gsub(taxon, pattern = 'k__', replacement = '')
			taxon <- gsub(taxon, pattern = 'd__', replacement = '')
			taxon <- gsub(taxon, pattern = 'p__', replacement = '')
			taxon <- gsub(taxon, pattern = 'c__', replacement = '')
			taxon <- gsub(taxon, pattern = 'o__', replacement = '')
			taxon <- gsub(taxon, pattern = 'f__', replacement = '')
			taxon <- gsub(taxon, pattern = 'g__', replacement = '')
			taxon <- gsub(taxon, pattern = 's__', replacement = '')
			taxon <- gsub(taxon, pattern = ';', replacement = '')
			taxon <- gsub(taxon, pattern = '__', replacement = 'unknown')
			taxon <- gsub(taxon, pattern = '-', replacement = '_')
			taxon <- gsub(taxon, pattern = '\\[', replacement = '')
			taxon <- gsub(taxon, pattern = '\\]', replacement = '')
			
			domain <- taxon[1]
			phylum <- taxon[2]
			class <- taxon[3]
			
			if (domain == 'unknown') stop('Unknown domain!')
			if (is.na(domain) || domain == '') stop('Empty domain!')
			if (is.na(phylum) || phylum == '') phylum <- 'unknown'
			if (is.na(class) || class == '') class <- 'unknown'
			
			taxon_combined <- if (collapse_taxon_level == 'class') {
				paste(c(domain, phylum, class), collapse = '_')
			} else {
				paste(c(domain, phylum), collapse = '_')
			}

			taxa_rhizo <- rbind(
				taxa_rhizo,
				data.table(
					taxon = taxon_combined,
					taxon_orig = taxon_orig,
					domain = domain,
					phylum = phylum,
					class = class			
				)
			)
		
		}

		names(raw_abund_sites_rhizo)[names(raw_abund_sites_rhizo) %notin% c('index', 'location', 'plant')] <- taxa_rhizo$taxon
		
		taxa_rhizo <- taxa_rhizo[!duplicated(taxa_rhizo$taxon)]

	### BULK SOIL: clean taxon names
	################################
	
		say('BULK SOIL: clean taxon names')
		
		cols <- colnames(raw_abund_sites_bulk)
		taxa <- cols[cols %notin% c('index', 'location', 'plant')]
		
		taxa_bulk <- data.table()
		for (i in seq_along(taxa)) {
		
			taxon <- taxa[[i]]
			taxon_orig <- taxon
			taxon <- trimws(taxon)
			taxon <- strsplit(taxon, split = ';')[[1]]

			taxon <- gsub(taxon, pattern = 'k__', replacement = '')
			taxon <- gsub(taxon, pattern = 'd__', replacement = '')
			taxon <- gsub(taxon, pattern = 'p__', replacement = '')
			taxon <- gsub(taxon, pattern = 'c__', replacement = '')
			taxon <- gsub(taxon, pattern = 'o__', replacement = '')
			taxon <- gsub(taxon, pattern = 'f__', replacement = '')
			taxon <- gsub(taxon, pattern = 'g__', replacement = '')
			taxon <- gsub(taxon, pattern = 's__', replacement = '')
			taxon <- gsub(taxon, pattern = ';', replacement = '')
			taxon <- gsub(taxon, pattern = '__', replacement = 'unknown')
			taxon <- gsub(taxon, pattern = '-', replacement = '_')
			taxon <- gsub(taxon, pattern = '\\[', replacement = '')
			taxon <- gsub(taxon, pattern = '\\]', replacement = '')
			
			domain <- taxon[1]
			phylum <- taxon[2]
			class <- taxon[3]
			
			if (domain == 'unknown') stop('Unknown domain!')
			if (is.na(domain) || domain == '') stop('Empty domain!')
			if (is.na(phylum) || phylum == '') phylum <- 'unknown'
			if (is.na(class) || class == '') class <- 'unknown'
			
			taxon_combined <- if (collapse_taxon_level == 'class') {
				paste(c(domain, phylum, class), collapse = '_')
			} else {
				paste(c(domain, phylum), collapse = '_')
			}

			taxa_bulk <- rbind(
				taxa_bulk,
				data.table(
					taxon = taxon_combined,
					taxon_orig = taxon_orig,
					domain = domain,
					phylum = phylum,
					class = class			
				)
			)
		
		}
		
		names(raw_abund_sites_bulk)[names(raw_abund_sites_bulk) %notin% c('index', 'location', 'plant')] <- taxa_bulk$taxon
		
		taxa_bulk <- taxa_bulk[!duplicated(taxa_bulk$taxon)]

	### RHIZOBIOME: sum abundances across columns with the same domain/phylum/class
	###############################################################################
	
		say('RHIZOBIOME: sum abundances across columns with the same domain/phylum/class')
		collapsed_abund_sites_rhizo <- raw_abund_sites_rhizo[ , c('index', 'location', 'plant')]
		for (taxon in taxa_rhizo$taxon) {
			
			these <- which(names(raw_abund_sites_rhizo) == taxon)
			these <- raw_abund_sites_rhizo[ , ..these]
			n <- rowSums(these)
			collapsed_abund_sites_rhizo[ , DUMMY := n]
			names(collapsed_abund_sites_rhizo)[ncol(collapsed_abund_sites_rhizo)] <- taxon

		}

		### check the collapsing process for some taxa
			
			raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/brooke_2024_10_04/level-7_Final_Root.csv')
			raw <- raw[index != rhizo_control_site]

			# taxon 1
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes;c__Clostridia' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Firmicutes' # >1 column has this kingdom/phylum
			}

			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(collapsed_abund_sites_rhizo[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

			# taxon 2
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Proteobacteria' # >1 column has this kingdom/phylum
			}

			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(collapsed_abund_sites_rhizo[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

	### BULK: sum abundances across columns with the same domain/phylum/class
	#########################################################################
	
		say('BULK: sum abundances across columns with the same domain/phylum/class')
		
		collapsed_abund_sites_bulk <- raw_abund_sites_bulk[ , c('index', 'location', 'plant')]
		for (taxon in taxa_bulk$taxon) {
			
			these <- which(names(raw_abund_sites_bulk) == taxon)
			these <- raw_abund_sites_bulk[ , ..these]
			n <- rowSums(these)
			collapsed_abund_sites_bulk[ , DUMMY := n]
			names(collapsed_abund_sites_bulk)[ncol(collapsed_abund_sites_bulk)] <- taxon

		}

		### check the collapsing process for some taxa
			
			raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/brooke_2024_10_04/level-7_Final_Soil.csv')
			raw <- raw[index != bulk_control_site]

			# taxon 1
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes;c__Clostridia' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Firmicutes' # >1 column has this kingdom/phylum
			}
			
			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(collapsed_abund_sites_bulk[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

			# taxon 2
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Proteobacteria' # >1 column has this kingdom/phylum/class
			}

			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(collapsed_abund_sites_bulk[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

	### COMBINED: taxon list
	########################
	
		say('COMBINED: taxon list')
		
		# flag for in raw_abund_sites_bulk or not
		tdpc <- c('taxon', 'domain', 'phylum', 'class')
		taxa_combined <- rbind(taxa_bulk[ , ..tdpc], taxa_rhizo[ , ..tdpc]) # order matters here!!!
		taxa_combined <- taxa_combined[!duplicated(taxa_combined)]
		taxa_combined[ , in_bulk := taxa_combined$taxon %in% taxa_bulk$taxon]
		taxa_combined[ , in_rhizobiome := taxa_combined$taxon %in% taxa_rhizo$taxon]

		if (collapse_taxon_level == 'phylum') {
			taxa_combined[ , class := NULL]
		}

		write.csv(taxa_combined, paste0(out_dir, '/taxa_combined.csv'), row.names = FALSE)

	### RHIZOBIOME: collate abundances
	##################################
	
		say('RHIZOBIOME: collate abundances')

		# site names
		sites_rhizo <- collapsed_abund_sites_rhizo[ , c('index', 'location', 'plant')]
		sites_rhizo$sample <- NA_character_

		### compile site-by-taxon rhizobiome abundance matrix
		# make one column for each taxon found in either rhizobiome or bulk soil (we need to have columns for taxa only in bulk so we can later stack it onto the bulk abundance matrix)
		abund_rhizo <- sites_rhizo
		abund_rhizo$rhizobiome_or_bulk <- 'rhizobiome'
		taxa <- taxa_combined$taxon
		for (taxon in taxa) {
			abund_rhizo[ , DUMMY := 0]
			names(abund_rhizo)[ncol(abund_rhizo)] <- taxon
		}

		for (i in 1:nrow(abund_rhizo)) {
		
			site <- abund_rhizo$index[i]

			for (taxon in taxa_rhizo$taxon) {
			
				n <- collapsed_abund_sites_rhizo[collapsed_abund_sites_rhizo$index == site, ..taxon]
				n <- unlist(n)

				abund_rhizo[i, taxon] <- n
			
			}
		
		}

		### check assignment process for some taxa
		
			raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/brooke_2024_10_04/level-7_Final_Root.csv')
			raw <- raw[index != rhizo_control_site]

			# taxon 1
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes;c__Clostridia' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Firmicutes' # >1 column has this kingdom/phylum
			}

			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(abund_rhizo[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

			# taxon 2
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Proteobacteria' # >1 column has this kingdom/phylum
			}

			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(abund_rhizo[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

	### BULK: collate abundances
	############################
	
		say('BULK: collate abundances')

		# site names
		sites_bulk <- collapsed_abund_sites_bulk[ , c('index', 'location', 'plant')]
		sites_bulk$sample <- NA_character_

		### compile site-by-taxon rhizobiome abundance matrix
		# make one column for each taxon found in either rhizobiome or bulk soil (we need to have columns for taxa only in bulk so we can later stack it onto the bulk abundance matrix)
		abund_bulk <- sites_bulk
		abund_bulk$rhizobiome_or_bulk <- 'bulk'
		taxa <- taxa_combined$taxon
		for (taxon in taxa) {
			abund_bulk[ , DUMMY := 0]
			names(abund_bulk)[ncol(abund_bulk)] <- taxon
		}

		for (i in 1:nrow(abund_bulk)) {
		
			site <- abund_bulk$index[i]

			for (taxon in taxa_bulk$taxon) {
			
				n <- collapsed_abund_sites_bulk[collapsed_abund_sites_bulk$index == site, ..taxon]
				n <- unlist(n)

				abund_bulk[i, taxon] <- n
			
			}
		
		}

		### check assignment process for some taxa
		
			raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/brooke_2024_10_04/level-7_Final_Soil.csv')
			raw <- raw[index != bulk_control_site]

			# taxon 1
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes;c__Clostridia' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Firmicutes' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Firmicutes' # >1 column has this kingdom/phylum
			}
			

			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(abund_bulk[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

			# taxon 2
			if (collapse_taxon_level == 'class') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria' # >1 column has this kingdom/phylum/class
				clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class
			} else if (collapse_taxon_level == 'phylum') {
				raw_test_taxon <- 'd__Bacteria;p__Proteobacteria' # >1 column has this kingdom/phylum
				clean_test_taxon <- 'Bacteria_Proteobacteria' # >1 column has this kingdom/phylum
			}
			
			these <- grepl(names(raw), pattern = raw_test_taxon)
			these <- raw[ , ..these]
			raw_n <- rowSums(these)

			clean_n <- unlist(abund_bulk[ , ..clean_test_taxon])

			stopifnot(all(raw_n == clean_n))

	### COMBINED: combined abundances
	#################################

		abund_combined <- rbind(abund_rhizo, abund_bulk)
		write.csv(abund_combined, paste0(out_dir, '/abundances_site_by_taxon_combined.csv'), row.names = FALSE)

		# check that the order of taxa in the taxonomic table is the same as the order in the abundance table

		stopifnot(all(colnames(abund_combined)[colnames(abund_combined) %notin% c('index', 'location', 'plant', 'sample', 'rhizobiome_or_bulk')] == taxa_combined$taxon))

	### RHIZOBIOME & BULK: set up environmental data frame
	######################################################
	say('RHIZOBIOME & BULK: set up environmental data frame')

		env_rhizo <- abund_rhizo[ , c('index', 'location', 'plant', 'sample', 'rhizobiome_or_bulk')]
		env_bulk <- abund_bulk[ , c('index', 'location', 'plant', 'sample', 'rhizobiome_or_bulk')]

		### add coordinates
		xy <- fread('./data_from_sonny/erica_harmonized_taxa_between_files_2024_07_26/environment_rhz_26JUL2024.csv')
		xy <- xy[ , c('site', 'latitude', 'longitude')]
		xy <- aggregate(xy, by = list(xy$site), FUN = mean)
		xy$site <- NULL
		names(xy)[1] <- 'site'
		
		# RHIZOBIOME
		env_rhizo[ , longitude := NA_real_]
		env_rhizo[ , latitude := NA_real_]

		for (i in 1:nrow(env_rhizo)) {
		
			index <- xy$site == env_rhizo$location[i]
			env_rhizo$latitude[i] <- xy$latitude[index]
			env_rhizo$longitude[i] <- xy$longitude[index]
		
		}

		# BULK
		env_bulk[ , longitude := NA_real_]
		env_bulk[ , latitude := NA_real_]

		for (i in 1:nrow(env_bulk)) {
		
			index <- xy$site == env_bulk$location[i]
			env_bulk$latitude[i] <- xy$latitude[index]
			env_bulk$longitude[i] <- xy$longitude[index]
		
		}

		env_combined <- rbind(env_rhizo, env_bulk)

	### COMBINED: add sampling date
	###############################

		sampling_dates$SITE_ID <- gsub(sampling_dates$SITE_ID, pattern = '_', replacement = '')

		env_combined[ , sampling_date := sampling_dates$SAMPLING_DATE[match(env_combined$location, sampling_dates$SITE_ID)]]

	### COMBINED: extract CHELSA climate to sites
	#############################################
	say('COMBINED: extract CHELSA climate to sites')

		env_combined_nad83 <- vect(env_combined, geom = c('longitude', 'latitude'), crs = getCRS('NAD83'))
		
		# elevation
		elevation <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/elevation.tif'))
		names(elevation) <- 'elevation_m'

		# BIOCLIMs
		bc <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))

		aridity <- bc$bio1 / (bc$bio12 + 1)
		names(aridity) <- 'aridity'
		
		env <- c(bc, aridity, elevation)
		env <- env[[c(MASTER_climate_predictor_names, 'elevation_m')]]

		env_combined_lambert <- project(env_combined_nad83, env)
		climate_extract <- extract(env, env_combined_lambert, ID = FALSE)

		env_combined <- cbind(env_combined, climate_extract)

	### COMBINED: extract PRISM weather to sites
	############################################

		say('COMBINED: extract PRISM weather to sites')

		# precipitation before sampling
		x <- prExtractRelativeDaily(
			pr_dir,
			x = env_combined_nad83,
			vars = 'ppt',
			date = 'sampling_date',
			res = 800,
			rastSuffix = 'tif',
			windowYears = 0,
			windowDays = 6,
			verbose = FALSE
		)
		
		x <- rowSums(x)
		
		env_combined$sampling_ppt_mm <- x

		# temperature before sampling
		x <- prExtractRelativeDaily(
			pr_dir,
			x = env_combined_nad83,
			vars = 'tmean',
			date = 'sampling_date',
			res = 800,
			rastSuffix = 'tif',
			windowYears = 0,
			windowDays = 6,
			verbose = FALSE
		)
		
		x <- rowMeans(x)
		
		env_combined$sampling_tmean_c <- x

	### COMBINED: extract FIELD SOIL TEXTURE data to sites
	######################################################

		say('COMBINED: extract FIELD SOIL TEXTURE data to sites')

		soil_texture <- read_xlsx('./data_from_sonny_and_loretta/site_soil_data_Order_52761.xlsx', sheet = 'Test Results For Order', skip = 9)
		soil_texture <- as.data.table(soil_texture)
		sample_name <- soil_texture$`Sample Name`
		sample_name <- strsplit(sample_name, split = ' ')
		sample_name <- do.call(rbind, sample_name)
		sample_name <- sample_name[ , 2]
		sample_name <- gsub(sample_name, pattern = '-', replacement = '')
		soil_texture$sample_name <- sample_name

		x <- soil_texture$`Sand\r\n %`[match(env_combined$location, soil_texture$sample_name)]
		x <- as.numeric(x)
		x <- x / 100
		env_combined$sand_field <- x

		x <- soil_texture$`Silt\r\n %`[match(env_combined$location, soil_texture$sample_name)]
		x <- as.numeric(x)
		x <- x / 100
		env_combined$silt_field <- x

		x <- soil_texture$`Clay\r\n %`[match(env_combined$location, soil_texture$sample_name)]
		x <- as.numeric(x)
		x <- x / 100
		env_combined$clay_field <- x

	### COMBINED: extract FIELD SOIL TEXTURE data to sites
	######################################################

		say('COMBINED: extract FIELD SOIL CHEMISTRY data to sites')
		
		soil_chemistry <- soil_chemistry[-1, ] # blank row
		
		sample_name <- soil_chemistry[ , 1]
		sample_name <- unlist(sample_name)
		sample_name <- gsub(sample_name, pattern = '\\*', replacement = '')
		sample_name <- toupper(sample_name)
		sample_name <- strsplit(sample_name, split = ' ')
		sample_name <- do.call(rbind, sample_name)
		sample_name <- sample_name[ , 1]
		sample_name[nchar(sample_name) == 2] <- paste0(sample_name[nchar(sample_name) == 2], 1)
		soil_chemistry$sample_name <- sample_name

		x <- soil_chemistry$`pH\n`[match(env_combined$location, soil_chemistry$sample_name)]
		x <- as.numeric(x)
		env_combined$ph_field <- x
		
		x <- soil_chemistry$`Total N %`[match(env_combined$location, soil_chemistry$sample_name)]
		x <- as.numeric(x)
		env_combined$nitrogen_field_perc <- x
		
		x <- soil_chemistry$`Total C %`[match(env_combined$location, soil_chemistry$sample_name)]
		env_combined$soc_field_perc <- x
		
	### COMBINED: extract SOILGRIDS SOIL CHEMISTRY data to sites
	############################################################

		say('COMBINED: extract SOILGRIDS SOIL CHEMISTRY data to sites')

		ph <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
		sand <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
		silt <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
		clay <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
		soc <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean_northAmerica.tif'))
		nitrogen <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean.tif'))

		env_combined_sg <- project(env_combined_nad83, ph)
		env_combined_sg_buffer <- buffer(env_combined_sg, 10000)

		ph <- crop(ph, env_combined_sg_buffer)
		sand <- crop(sand, env_combined_sg_buffer)
		silt <- crop(silt, env_combined_sg_buffer)
		clay <- crop(clay, env_combined_sg_buffer)
		soc <- crop(soc, env_combined_sg_buffer)
		nitrogen <- crop(nitrogen, env_combined_sg_buffer)

		soil <- c(ph, sand, silt, clay, soc, nitrogen)
		names(soil) <- c('ph_soilgrids', 'sand_soilgrids', 'silt_soilgrids', 'clay_soilgrids', 'soc_soilgrids_perc', 'nitrogen_soilgrids_perc')

		soil[['ph_soilgrids']] <- soil[['ph_soilgrids']] / 10
		soil[['sand_soilgrids']] <- soil[['sand_soilgrids']] / 1000
		soil[['silt_soilgrids']] <- soil[['silt_soilgrids']] / 1000
		soil[['clay_soilgrids']] <- soil[['clay_soilgrids']] / 1000
		soil[['soc_soilgrids_perc']] <- soil[['soc_soilgrids_perc']] / 1000
		soil[['nitrogen_soilgrids_perc']] <- soil[['nitrogen_soilgrids_perc']] / 10000

		soil_soilgrids_extract <- extract(soil, env_combined_sg, ID = FALSE)
		env_combined <- cbind(env_combined, soil_soilgrids_extract)

	### COMBINED: Andropogon gerardi estimated lambda from species-level SDM
	########################################################################
	
		say('COMBINED: Andropogon gerardi estimated lambda from species-level N-MIXTURE SDM')

		### posterior from SDM
		chains <- readRDS('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_0.95]_[climate_soil]_[priors_ddexp]/sdm_nmixture_chains.rds')

		# subset chain summary to just the lambdas associated with background sites
		summary <- chains$summary$all.chains

		which_lambda <- grepl(rownames(summary), pattern = 'lambda')
		lambda <- summary[which_lambda, ]

		### spatial vector
		ag_vect <- vect('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_0.95]_[climate_soil]_[priors_ddexp]/sdm_nmixture_1961_2020_climate_focus.gpkg')

		# fields <- c('area_km2', 'any_ag_quality1to3', 'num_poaceae_records')
		# ag_vect <- ag_vect[ , fields]

		# ag <- as.data.frame(ag_vect)
		# completes <- complete.cases(as.data.frame(ag_vect))
		# ag_focus <- ag[completes, ]
		# ag_vect_focus <- ag_vect[completes, ]
		
		# # get just counties with data
		# ag_vect_focus$ag_lambda <- lambda[ , 'Mean']
		
		ag_lambda_extract <- extract(ag_vect, env_combined_lambert)
		ag_lambda <- ag_lambda_extract$lambda_mean

		env_combined$ag_lambda_nmixture <- ag_lambda
		
	write.csv(env_combined, paste0(out_dir, '/environment_combined.csv'), row.names = FALSE)

	### COMBINED: study design
	##########################

		study_design_combined <- env_combined[ , c('location', 'plant', 'rhizobiome_or_bulk')]
		write.csv(study_design_combined, paste0(out_dir, './study_design_combined.csv'), row.names = FALSE)

	### DENOUEMENT
	##############

	sink(paste0(out_dir, '/!readme.txt'), split = TRUE)
	say('COLLATED ABUNDANCE, TAXONOMIC, and ENVIRONMENTAL DATA FOR MICROBES ASSOCIATED WITH ANDROPOGON GERARDI')
	say(date(), post = 2)
	say('These files represent:')
	say('* Abundance of microbes associated with Andropogon gerardi roots and in bulk soil samples')
	say('* Taxonomic identities (domain/phylum/class) for the samples')
	say('* Environmental data associated with each site')
	say('* A study design matrix', post = 2)
	say('IMPORTANT: Except for the taxonomic table, the *order* of rows and columns in files needs to coincide. They should be subsetted or re-ordered simultaneously (or not at all, preferably). Specifically, these orders need to be respected:')
	say(' * Rows of the environmental table, the abundances table, and the study design table')
	say(' * Rows of the taxonomic table and columns in the abundance table that have abundances')
	say('This script contains automated checks to ensure proper order is retained between tables and derived abundances match those from the raw data.')
	sink()

# say('#############################################################################################')
# say('### construct raster template and data table of predictors for making spatial predictions ###')
# say('#############################################################################################')
	
# 	### North America
# 	#################

# 	say('extent')
# 	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))

# 	nam <- nam[!(nam$NAME_1 %in% c('Alaska', 'Hawaii'))]
# 	extent <- as.vector(ext(nam))
# 	extent[1] <- -117
# 	extent[3] <- 22
# 	extent[4] <- 55
# 	extent <- ext(extent)
# 	nam <- crop(nam, extent)

# 	# template
# 	climatena <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))

# 	### soil
# 	########

# 	say('SoilGrids')
# 	ph <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
# 	sand <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
# 	silt <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
# 	clay <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
# 	soc <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean_northAmerica.tif'))
# 	nitrogen <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean.tif'))

# 	extent_vect <- as.polygons(extent, crs = crs(nam))
# 	extent_vect_homolosine <- project(extent_vect, ph)
# 	extent_vect_homolosine_buffer <- buffer(extent_vect_homolosine, 500000)
# 	extent_vect_homolosine_buffer <- ext(extent_vect_homolosine_buffer)

# 	ph <- crop(ph, extent_vect_homolosine_buffer)
# 	sand <- crop(sand, extent_vect_homolosine_buffer)
# 	silt <- crop(silt, extent_vect_homolosine_buffer)
# 	clay <- crop(clay, extent_vect_homolosine_buffer)
# 	soc <- crop(soc, extent_vect_homolosine_buffer)
# 	nitrogen <- crop(nitrogen, extent_vect_homolosine_buffer)

# 	nitrogen <- crop(nitrogen, ph)

# 	soil <- c(ph, sand, silt, clay, soc, nitrogen)
# 	names(soil) <- c('ph_soilgrids', 'sand_soilgrids', 'silt_soilgrids', 'clay_soilgrids', 'soc_soilgrids_perc', 'nitrogen_soilgrids_perc')

# 	soil[['ph_soilgrids']] <- soil[['ph_soilgrids']] / 10
# 	soil[['sand_soilgrids']] <- soil[['sand_soilgrids']] / 1000
# 	soil[['silt_soilgrids']] <- soil[['silt_soilgrids']] / 1000
# 	soil[['clay_soilgrids']] <- soil[['clay_soilgrids']] / 1000
# 	soil[['soc_soilgrids_perc']] <- soil[['soc_soilgrids_perc']] / 1000
# 	soil[['nitrogen_soilgrids_perc']] <- soil[['nitrogen_soilgrids_perc']] / 10000

# 	soil <- project(soil, climatena, threads = 4)

# 	### BIOCLIMs & AG lambda: present
# 	#################################

# 	say('ClimateNA - present')

# 	# Get values of ClimateNA from the output from the SDM. This is much faster than repeating the calculation of BIOCLIMs again.
# 	sq_vect_focus <- vect(paste0('./outputs_loretta/sdm_[nmixture]/sdm_nmixture_1961_2020_climate_focus.gpkg'))
# 	sq_vect_complement <- vect(paste0('./outputs_loretta/sdm_[nmixture]/sdm_nmixture_1961_2020_climate_focus_complement.gpkg'))
# 	sq_vect <- rbind(sq_vect_focus, sq_vect_complement)
# 	vars <- sort(unique(c(MASTER_climate_predictor_names, 'lambda_mean')))
# 	if (exists('sq_rasts')) rm(sq_rasts)
# 	for (var in vars) {

# 		sq_rast <- rasterize(sq_vect, climatena, field = var)
# 		names(sq_rast) <- var

# 		if (exists('sq_rasts')) {
# 			sq_rasts <- c(sq_rasts, sq_rast)
# 		} else {
# 			sq_rasts <- sq_rast
# 		}

# 	}
# 	names(sq_rasts)[names(sq_rasts) == 'lambda_mean'] <- 'ag_lambda_nmixture'

# 	# bc <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))
# 	# layers <- c('bio1', 'bio12', MASTER_climate_predictor_names[MASTER_climate_predictor_names %in% names(bc)])
# 	# layers <- unique(layers)
# 	# bc <- bc[[layers]]

# 	# nam_lambert <- project(nam, bc)
# 	# bc <- crop(bc, nam_lambert)

# 	# aridity <- bc$bio1 / (bc$bio12 + 1)
# 	# names(aridity) <- 'aridity'
	
# 	# # AG SDM
# 	# ag_sdm_focus <- vect('./outputs_loretta/sdm_[nmixture]/sdm_nmixture_1961_2020_climate_focus.gpkg')
# 	# ag_sdm_complement <- vect('./outputs_loretta/sdm_[nmixture]/sdm_nmixture_1961_2020_climate_focus_complement.gpkg')
# 	# ag_sdm <- rbind(ag_sdm_focus, ag_sdm_complement)
# 	# ag_sdm <- rasterize(ag_sdm, climatena, field = 'lambda_mean')
# 	# names(ag_sdm) <- 'ag_lambda_nmixture'
# 	# ag_sdm <- crop(ag_sdm, bc)
	
# 	# sq_rasts <- c(bc, aridity, ag_sdm)
# 	# sq_rasts <- sq_rasts[[c(MASTER_climate_predictor_names, 'ag_lambda_nmixture')]]

# 	### BIOCLIMs & AG lambda: future
# 	################################

# 	futs <- c(
# 		'ssp245_2041_2070',
# 		'ssp245_2071_2100',
# 		'ssp370_2041_2070',
# 		'ssp370_2071_2100'
# 	)

# 	clim_futs <- list()
# 	vars <- sort(unique(c(MASTER_climate_predictor_names, 'lambda_mean')))
# 	for (fut in futs) {

# 		say('ClimateNA - ', fut)

# 		# Get future values of ClimateNA from the output from the SDM. This is much faster than repeating the calculation of BIOCLIMs again.
# 		fut_vect <- vect(paste0('./outputs_loretta/sdm_[nmixture]/sdm_nmixture_', fut, '.gpkg'))
# 		if (exists('fut_rasts')) rm(fut_rasts)
# 		for (var in vars) {

# 			fut_rast <- rasterize(fut_vect, climatena, field = var)
# 			names(fut_rast) <- var

# 			if (exists('fut_rasts')) {
# 				fut_rasts <- c(fut_rasts, fut_rast)
# 			} else {
# 				fut_rasts <- fut_rast
# 			}

# 		}
# 		names(fut_rasts)[names(fut_rasts) == 'lambda_mean'] <- 'ag_lambda_nmixture'
# 		clim_futs[[length(clim_futs) + 1]] <- fut_rasts

# 	}
# 	names(clim_futs) <- futs

# 	### aggregate up to make smaller grid cells... we need to do this to obviate memory issues when predict()'ing HMSC
# 	##################################################################################################################
# 	say('combined')

# 	sq_rasts <- crop(sq_rasts, soil)
# 	sq_rasts <- extend(sq_rasts, soil)

# 	soil <- aggregate(soil, 8, 'mean', na.rm = TRUE)

# 	sq_rasts <- aggregate(sq_rasts, 8, 'mean', na.rm = TRUE)
# 	clim_futs <- lapply(clim_futs, aggregate, 8, 'mean', na.rm = TRUE)

# 	env_sq <- c(sq_rasts, soil)
# 	for (i in seq_along(clim_futs)) clim_futs[[i]] <- c(clim_futs[[i]], soil)
	
# 	writeRaster(env_sq, paste0('./outputs_sonny/climatena_1961_2020_soilgrids_aggregated_8x.tif'), overwrite = TRUE)
# 	for (i in seq_along(clim_futs)) {

# 		name <- names(clim_futs)[i]

# 		writeRaster(clim_futs[[i]], paste0('./outputs_sonny/climatena_ensemble_8GCMs_', name, '_soilgrids_aggregated_8x.tif'), overwrite = TRUE)

# 	}

say('DONE!', level = 1, deco = '!')

### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles abundance data for microbes associated with the rhizobiome of Andropogon gerardi and bulk soil samples taken in the plants' vicinities.
###
### source('C:/Ecology/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_02_collate_microbial_data.r')
### source('E:/Adam/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_02_collate_microbial_data.r')
###
### CONTENTS ###
### setup ###
### collate rhizobiome & bulk soil microbial abundances ###
### collapse abundance, phylogeny, and trait data by taxon (ignoring rID tags) ###

#############
### setup ###
#############

rm(list = ls())

# drive <- 'C:/Ecology/'
drive <- 'E:/Adam/'

setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

devtools::load_all(paste0(drive, '/R/airUpThere'))

library(data.table)
library(enmSdmX)
library(readxl)
library(omnibus)

say('######################################################################################')
say('### compile microbial abundance and environmental data into formats needed by HMSC ###')
say('######################################################################################')

	### user-defined
	################
	
		# number of years across which "long-term" variables are calculated
		master_long_term_years <- 10

		# start/end month-day of seasons
		summer_month_day_start <- '06-01'
		summer_month_day_end <- '08-31'

		winter_month_day_start <- '12-01'
		winter_month_day_end <- '02-28'
		
		# pr_dir <- 'D:/ecology/PRISM/working/an81'
		# pr_dir <- 'E:/ecology/PRISM/working/an81'
		pr_dir <- 'F:/ecology/PRISM/working/an81'
		
		climate_years <- (2022 - master_long_term_years + 1):2022 # years across which to calculate long-term variables
	
	# Want:
	# * An "abundances" matrix with each row = sample, each column = taxon abundance.
	# 	* Need to collapse the abundances by Domain/Phylum/Class because some are represented >1 time.
	#	* Need to remove "control" samples (positive and negative controls.)
	# * Three "taxonomy" matrices (rhizobiome, bulk soil, and combined) with each row = one taxon, where there are at least three columns, named "Domain", "Phylum", and "Class". All taxa in "abundances" must be represented. No duplicates allowed in these taxonomy matrices.
	# * An "environment" frame with one row per site/sample and with environmental data, including a flag indicating if a sample was from the bulk or rhizomicrobiome soil.
	#' * A "study design" matrix with factors for each site and each sample within a site.

	### save to
	###########
	
	out_dir <- paste0('./data_from_sonny/collated_for_hmsc')
	dirCreate(out_dir)

	### inputs
	##########
	
	raw_abund_sites_rhizo <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Rhizomicrobiome_Level_7.csv')
	
	raw_abund_sites_bulk <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Bulk_Soil_Level_7_R-Friendly.csv')
	
	soil_chem <- read_xlsx('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Biogeographical_Sampling.xlsx', sheet = 'Soil Chemistry', skip = 1)
	
	soil_chem <- as.data.table(soil_chem)

	### RHIZOBIOME & BULK: remove controls
	######################################
	say('RHIZOBIOME & BULK: remove controls')

	raw_abund_sites_rhizo <- raw_abund_sites_rhizo[index != 'NC_S8']
	cond <- names(raw_abund_sites_bulk) != 'NC-11292023-16s'
	raw_abund_sites_bulk <- raw_abund_sites_bulk[ , ..cond, drop = FALSE]

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
		
		taxon_combined <- paste(c(domain, phylum, class), collapse = '_')

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
	
	taxa_rhizo <- taxa_rhizo[!duplicated(taxa_rhizo)]
	write.csv(taxa_rhizo, paste0(out_dir, '/taxa_rhizobiome.csv'), row.names = FALSE)

	### RHIZOBIOME: rename abundance columns to proper taxon names
	##############################################################
	say('RHIZOBIOME: rename abundance columns to proper taxon names')
	
	cols <- colnames(raw_abund_sites_rhizo)
	taxa <- cols[cols %notin% c('index', 'location', 'plant')]

	for (taxon in taxa) {
		colnames(raw_abund_sites_rhizo)[colnames(raw_abund_sites_rhizo) == taxon] <- taxa_rhizo$taxon[taxa_rhizo$taxon_orig == taxon]
	}

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
		
		raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Rhizomicrobiome_Level_7.csv')
		raw <- raw[index != 'NC_S8']

		# taxon 1
		raw_test_taxon <- 'k__Bacteria;p__Firmicutes;c__Clostridia' # >1 column has this kingdom/phylum/class
		clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class

		these <- grepl(names(raw), pattern = raw_test_taxon)
		these <- raw[ , ..these]
		raw_n <- rowSums(these)

		clean_n <- unlist(collapsed_abund_sites_rhizo[ , ..clean_test_taxon])

		stopifnot(all(raw_n == clean_n))

		# taxon 2
		raw_test_taxon <- 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria' # >1 column has this kingdom/phylum/class
		clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class

		these <- grepl(names(raw), pattern = raw_test_taxon)
		these <- raw[ , ..these]
		raw_n <- rowSums(these)

		clean_n <- unlist(collapsed_abund_sites_rhizo[ , ..clean_test_taxon])

		stopifnot(all(raw_n == clean_n))

	### BULK: clean taxon names
	###########################
	say('BULK: clean taxon names')
	
	taxa_bulk <- raw_abund_sites_bulk[ , c('Kingdom', 'Phylum', 'Class')]	
	colnames(taxa_bulk)[colnames(taxa_bulk) == 'Kingdom'] <- 'Domain'
	colnames(taxa_bulk) <- tolower(colnames(taxa_bulk))
	
	# # remove controls
	# taxa_bulk <- taxa_bulk[domain != 'd__Bacteria' & phylum != '__' & class != '__']	

	taxa_bulk[ , domain_orig := domain]
	taxa_bulk[ , phylum_orig := phylum]
	taxa_bulk[ , class_orig := class]
	
	for (col in c('domain', 'phylum', 'class')) {
	
		taxa_bulk[ , (col) := lapply(.SD, 'trimws'), .SDcols = col]
		
		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = 'd__', replacement = ''), .SDcols = col]
		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = 'p__', replacement = ''), .SDcols = col]
		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = 'c__', replacement = ''), .SDcols = col]
		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = '\\(', replacement = '_'), .SDcols = col]
		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = ')', replacement = ''), .SDcols = col]

		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = '__', replacement = 'unknown'), .SDcols = col]
		taxa_bulk[ , (col) := lapply(.SD, gsub, pattern = '-', replacement = '_'), .SDcols = col]
		
	}
	
	taxa_bulk$taxon <- apply(taxa_bulk[ , c('domain', 'phylum', 'class')], 1, paste, collapse = '_')

	# add "nice" taxon names to bulk
	raw_abund_sites_bulk$domain <- taxa_bulk$domain
	raw_abund_sites_bulk$phylum <- taxa_bulk$phylum
	raw_abund_sites_bulk$class <- taxa_bulk$class
	raw_abund_sites_bulk$taxon <- taxa_bulk$taxon
	
	taxa_bulk <- taxa_bulk[!duplicated(taxa_bulk)]
	
	taxa_bulk <- taxa_bulk[ , c('taxon', 'domain', 'phylum', 'class', 'domain_orig', 'phylum_orig', 'class_orig')]
	write.csv(taxa_bulk, paste0(out_dir, '/taxa_bulk.csv'), row.names = FALSE)
	
	# # ### BULK: sum abundances across rows with the same domain/phylum/class
	# # ######################################################################
	# # say('BULK: sum abundances across rows with the same domain/phylum/class')

	# # sites_bulk <- colnames(raw_abund_sites_bulk)[colnames(raw_abund_sites_bulk) %notin% c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Samples', 'domain', 'phylum', 'class', 'taxon')]

	# # collapsed_abund_sites_bulk <- raw_abund_sites_bulk[ , c('taxon', 'domain', 'phylum', 'class', ..sites_bulk)]
	# # cols <- c('taxon', 'domain', 'phylum', 'class')
	# # collapsed_abund_sites_bulk <- collapsed_abund_sites_bulk[ , lapply(.SD, sum), by = cols]

	# # ### check the collapsing process for some taxa
		
	# # 	raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Bulk_Soil_Level_7_R-Friendly.csv')
	# # 	cond <- names(raw) != 'NC-11292023-16s'
	# # 	raw <- raw[ , ..cond, drop = FALSE]

	# # 	# taxon 1
	# # 	raw_test_domain <- 'd__Bacteria' # >1 row has this kingdom/phylum/class
	# # 	raw_test_phylum <- 'p__Firmicutes' # >1 row has this kingdom/phylum/class
	# # 	raw_test_class <- 'c__Clostridia' # >1 row has this kingdom/phylum/class
	# # 	clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class

	# # 	these <- which(raw$Kingdom == raw_test_domain & raw$Phylum == raw_test_phylum & raw$Class == raw_test_class)
	# # 	these <- raw[these, ..sites_bulk]
	# # 	raw_n <- colSums(these)

	# # 	these <- which(collapsed_abund_sites_bulk$taxon == clean_test_taxon)
	# # 	clean_n <- unlist(collapsed_abund_sites_bulk[these, ..sites_bulk])

	# # 	stopifnot(all(raw_n == clean_n))

	# # 	# taxon 2
	# # 	raw_test_domain <- 'd__Bacteria' # >1 row has this kingdom/phylum/class
	# # 	raw_test_phylum <- 'p__Proteobacteria' # >1 row has this kingdom/phylum/class
	# # 	raw_test_class <- 'c__Gammaproteobacteria' # >1 row has this kingdom/phylum/class
	# # 	clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class

	# # 	these <- which(raw$Kingdom == raw_test_domain & raw$Phylum == raw_test_phylum & raw$Class == raw_test_class)
	# # 	these <- raw[these, ..sites_bulk]
	# # 	raw_n <- colSums(these)

	# # 	these <- which(collapsed_abund_sites_bulk$taxon == clean_test_taxon)
	# # 	clean_n <- unlist(collapsed_abund_sites_bulk[these, ..sites_bulk])

	# # 	stopifnot(all(raw_n == clean_n))

	### COMBINED: taxon list
	########################
	say('COMBINED: taxon list')
	
	# flag for in raw_abund_sites_bulk or not
	tdpc <- c('taxon', 'domain', 'phylum', 'class')
	taxa_combined <- rbind(taxa_bulk[ , ..tdpc], taxa_rhizo[ , ..tdpc]) # order matters here!!!
	taxa_combined <- taxa_combined[!duplicated(taxa_combined)]
	taxa_combined[ , in_bulk := taxa_combined$taxon %in% taxa_bulk$taxon]
	taxa_combined[ , in_rhizobiome := taxa_combined$taxon %in% taxa_rhizo$taxon]

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
	
		raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Rhizomicrobiome_Level_7.csv')
		raw <- raw[index != 'NC_S8']

		# taxon 1
		raw_test_taxon <- 'k__Bacteria;p__Firmicutes;c__Clostridia' # >1 column has this kingdom/phylum/class
		clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class

		these <- grepl(names(raw), pattern = raw_test_taxon)
		these <- raw[ , ..these]
		raw_n <- rowSums(these)

		clean_n <- unlist(abund_rhizo[ , ..clean_test_taxon])

		stopifnot(all(raw_n == clean_n))

		# taxon 2
		raw_test_taxon <- 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria' # >1 column has this kingdom/phylum/class
		clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class

		these <- grepl(names(raw), pattern = raw_test_taxon)
		these <- raw[ , ..these]
		raw_n <- rowSums(these)

		clean_n <- unlist(abund_rhizo[ , ..clean_test_taxon])

		stopifnot(all(raw_n == clean_n))

	write.csv(abund_rhizo, paste0(out_dir, '/abundances_site_by_taxon_rhizobiome.csv'), row.names = FALSE)

	### BULK: collate abundances
	############################

	say('BULK: collate abundances')

	sites_bulk <- colnames(raw_abund_sites_bulk)[colnames(raw_abund_sites_bulk) %notin% c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Samples', 'domain', 'phylum', 'class', 'taxon')]

	sites_bulk_clean <- gsub(sites_bulk, pattern = '-16s', replacement = '')

	sites_bulk_clean <- strsplit(sites_bulk_clean, split = '-')
	sites_bulk_clean <- do.call(rbind, sites_bulk_clean)
	abund_bulk <- data.table(index = sites_bulk, location = sites_bulk_clean[ , 1], sample = sites_bulk_clean[ , 2], plant = sites_bulk_clean[ , 3])

	abund_bulk$rhizobiome_or_bulk <- 'bulk'

	### compile site-by-taxon abundance matrix
	# make one column for each taxon found in either rhizobiome or bulk soil (we need to have columns for taxa only in bulk so we can later stack it onto the rhizobiome abundance matrix)
	taxa <- taxa_combined$taxon
	for (taxon in taxa) {
		abund_bulk[ , DUMMY := 0]
		names(abund_bulk)[ncol(abund_bulk)] <- taxon
	}

	for (i in 1:nrow(abund_bulk)) {
	
		site <- abund_bulk$index[i]

		for (taxon in taxa_bulk$taxon) {
		
			rows <- raw_abund_sites_bulk$taxon == taxon
			n <- raw_abund_sites_bulk[rows, ..site]
			n <- unlist(n)
			n <- sum(n)

			abund_bulk[i, taxon] <- n
		
		}
	
	}

	### check assignment process for some taxa
	
	### check the collapsing process for some taxa
		
		raw <- fread('./data_from_sonny/!ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Bulk_Soil_Level_7_R-Friendly.csv')
		cond <- names(raw) != 'NC-11292023-16s'
		raw <- raw[ , ..cond, drop = FALSE]

		# taxon 1
		raw_test_domain <- 'd__Bacteria' # >1 row has this kingdom/phylum/class
		raw_test_phylum <- 'p__Firmicutes' # >1 row has this kingdom/phylum/class
		raw_test_class <- 'c__Clostridia' # >1 row has this kingdom/phylum/class
		clean_test_taxon <- 'Bacteria_Firmicutes_Clostridia' # >1 column has this kingdom/phylum/class

		these <- which(raw$Kingdom == raw_test_domain & raw$Phylum == raw_test_phylum & raw$Class == raw_test_class)
		these <- raw[these, ..sites_bulk]
		raw_n <- colSums(these)

		these <- which(colnames(abund_bulk) == clean_test_taxon)
		clean_n <- unlist(abund_bulk[ , ..these])

		stopifnot(all(raw_n == clean_n))

		# taxon 2
		raw_test_domain <- 'd__Bacteria' # >1 row has this kingdom/phylum/class
		raw_test_phylum <- 'p__Proteobacteria' # >1 row has this kingdom/phylum/class
		raw_test_class <- 'c__Gammaproteobacteria' # >1 row has this kingdom/phylum/class
		clean_test_taxon <- 'Bacteria_Proteobacteria_Gammaproteobacteria' # >1 column has this kingdom/phylum/class

		these <- which(raw$Kingdom == raw_test_domain & raw$Phylum == raw_test_phylum & raw$Class == raw_test_class)
		these <- raw[these, ..sites_bulk]
		raw_n <- colSums(these)

		these <- which(colnames(abund_bulk) == clean_test_taxon)
		clean_n <- unlist(abund_bulk[ , ..these])

		stopifnot(all(raw_n == clean_n))

	write.csv(abund_bulk, paste0(out_dir, '/abundances_site_by_taxon_bulk.csv'), row.names = FALSE)

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

	### COMBINED: extract climate to sites
	######################################
	say('COMBINED: extract climate to sites')

	
say('DONE!', level = 1, deco = '!')

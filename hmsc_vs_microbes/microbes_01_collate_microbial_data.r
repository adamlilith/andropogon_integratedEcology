### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles abundance data for microbes associated with the rhizobiome of Andropogon gerardi and bulk soil samples taken in the plants' vicinities.
###
### source('E:/Ecology/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_01_collate_microbial_data.r')
### source('E:/Adam/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_01_collate_microbial_data.r')
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

library(data.table)
library(read_xl)
library(omnibus)

# say('###########################################################')
# say('### collate rhizobiome & bulk soil microbial abundances ###')
# say('###########################################################')

	# ### collate rhizobiome & bulk soil microbial abundances
	# #######################################################

	# # We want a matrix with cell values representing abundances
	# # each column representing a taxon
	# # each row a site
	# # column 1: site
	# # column 2: sample ID

	# rhizo <- fread('./data/data_from_sonny/ORIGINALS_do_not_manipulate_only_copy/DATA_Microbiome_ASVs_Rhizomicrobiome_Level_7.csv')
	# bulk <- fread('./data/data_from_sonny/DATA_Microbiome_ASVs_Bulk_Soil_Level_7_R-Friendly.csv')

	# # remove control columns/rows
	# bulk <- bulk[ , `NC-11292023-16s` := NULL]
	# rhizo <- rhizo[index != 'NC_S8']

	# bulk <- bulk[Kingdom != ''] # reads empty rows for some reason

	# microbes <- cbind(
		# data.table(
			# site = rhizo$location,
			# sample = rhizo$plant,
			# type = 'rhizobiome'
		# ),
		# rhizo
	# )
	# colnames(microbes)[colnames(microbes) == 'index'] <- 'site_sample'
	# microbes[ , c('location', 'plant') := NULL]

	# # add bulk samples to rhizobiome
	# bulk_meta_cols <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Samples')
	# sample_cols <- colnames(bulk)[colnames(bulk) %notin% bulk_meta_cols]

	# microbes_meta_cols <- c('site', 'sample', 'type', 'site_sample')
	# microbes_abund_cols <- which(colnames(microbes) %notin% microbes_meta_cols)
	# for (i in seq_along(sample_cols))
	
		# site_sample <- sample_cols[i]
		# site <- substr(site_sample, 1, 3)
		# sample <- substr(site_sample, 5, nchar(site_sample))
	
		# # grab existing row and use it as a template
		# new_row <- microbes[1]
		# new_row$site <- site
		# new_row$sample <- sample
		# new_row$type <- 'bulk'
		# new_row$site_sample <- site_sample
		# new_row[ , (microbes_abund_cols) := 0]

		# # populate cells with abundances	
		# bulk_taxa <- apply(bulk[ , c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')], 1, paste, collapse = ';')
		# microbe_taxa <- colnames(microbes[ , ..microbes_abund_cols])
		# for (j in seq_along(bulk_taxa)) {

			# bulk_taxon <- bulk_taxa[j]
			
			# # this bulk taxon is already in microbe data
			# if (bulk_taxon %in% microbe_taxa) {
			
				# stop()
				# this_col <- which(microbe_taxa == bulk_taxon)
				# new_row[1, this_col] <- bulk[XYZ]

			# # this bulk taxon is not yet in microbe data
			# } else {

				# # add taxon to microbe data
				# microbes[ , DUMMY := 0]
				# colnames(microbes)[ncol(microbes)] <- bulk_taxon

				# # add taxon to new data
				# new_row[1, ..bulk_taxon := bulk[j, ..site_sample]]
				# colnames(new_row)[ncol(new_row)] <- bulk_taxon

				# microbe_taxa <- c(microbe_taxa, bulk_taxon)

			# }

			# microbes <- rbind(microbes, new_row)

		# }

	# fwrite(microbes, './microbes_rhizobiome_bulk.csv', row.names = FALSE)

say('##################################################################################')
say('### collapse abundance, phylogeny, and trait data by taxon (ignoring rID tags) ###')
say('##################################################################################')

	# Task: Collapse the abundances, "trait" data (taxonomic data), and sites by taxon, ignoring the "rID" tags.
	# 	* species_rhz_26JUL2024 (abundances)
	# 	* environment_rhz_26JUL2024
	# 	* rhizo_HMSC_reorder_26JUL2024
	# 	* studyDesign_rhz_26JUL2024
	# 	* traits_rhz_26JUL2024
	
	### general
	###########
	
	out_dir <- paste0('./data/data_from_sonny/compiled_for_modeling_with_hmsc/', Sys.Date(), '_collapsed_by_taxon')
	dirCreate(out_dir)
	
	### add abundances across columns that represent the same taxon
	###############################################################
	
	abundances <- read.csv('./data/data_from_sonny/compiled_for_modeling_with_hmsc/2024_07_26_erica/species_rhz_26JUL2024.csv',  as.is = FALSE, check.names = FALSE)
	taxon_names <- colnames(abundances)
	taxon_names <- taxon_names[taxon_names %notin% c('my.id', 'id')]
	rID_positions <- regexpr(taxon_names, pattern = 'rID_')
	taxon_names <- substr(taxon_names, 1, rID_positions - 1)
	taxon_names <- trimws(taxon_names)
	unique_taxa <- unique(taxon_names)

	collapsed_abundances <- abundances[ , c('my.id', 'id')]
	for (i in seq_along(unique_taxa)) {
	
		taxon <- unique_taxa[i]
		columns_with_taxon <- which(taxon_names == taxon) + 2 # add 2 because of 'my.id' and 'id' columns
		this_abund <- abundances[ , columns_with_taxon, drop = FALSE]
		this_abund <- rowSums(this_abund)
		this_abund <- data.frame(this_abund)
		colnames(this_abund) <- taxon
		
		collapsed_abundances <- cbind(collapsed_abundances, this_abund)
	
	}

	# make nice taxon names
	taxa <- colnames(collapsed_abundances)
	taxa <- strsplit(taxa, split = ' ')

	new_names <- rep(NA_character_, length(taxa))
	for (i in seq_along(taxa)) {

		if (length(taxa[[i]]) == 1) {
			new_names[i] <- taxa[[i]]
		} else {
			this_name <- taxa[[i]]
			this_name <- gsub(this_name, pattern = 'd__', replacement = '')
			this_name <- gsub(this_name, pattern = 'p__', replacement = '')
			this_name <- gsub(this_name, pattern = 'c__', replacement = '')
			this_name <- gsub(this_name, pattern = '__', replacement = 'unknown')
			this_name <- gsub(this_name, pattern = '-', replacement = '_')
			this_name <- paste(this_name, collapse = '_')
			new_names[i] <- this_name
		}

	}

	new_names <- trimws(new_names)
	names(collapsed_abundances) <- new_names
	
	write.csv(collapsed_abundances, paste0(out_dir, '/abundances_', Sys.Date(), '.csv'), row.names = FALSE)
	
	### taxonomy/traits
	###################
	
	taxonomy_traits <- read.csv('./data/data_from_sonny/compiled_for_modeling_with_hmsc/2024_07_26_erica/traits_rhz_26JUL2024.csv')

	# taxon names must match those from "collapsed_abundances"
	taxa <- taxonomy[ , c('Domain', 'Phylum', 'Class')]
	taxa <- apply(taxa, 2, gsub, pattern = 'd__', replacement = '')
	taxa <- apply(taxa, 2, gsub, pattern = 'p__', replacement = '')
	taxa <- apply(taxa, 2, gsub, pattern = 'c__', replacement = '')
	taxa <- apply(taxa, 2, gsub, pattern = '__', replacement = 'unknown')
	taxa <- apply(taxa, 2, gsub, pattern = '-', replacement = '_')
	taxa <- apply(taxa, 2, trimws)
	taxa_together <- apply(taxa, 1, paste, collapse = '_')
	
	taxonomy_traits$taxon <- taxa_together

	taxonomy_traits$Domain <- taxa[ , 'Domain']
	taxonomy_traits$Phylum <- taxa[ , 'Phylum']
	taxonomy_traits$Class <- taxa[ , 'Class']

	# flip "in_bulk" flag in "phylogeny" data to TRUE if at least one taxon occurs in at least one bulk sample
	taxa <- taxonomy_traits$taxon
	for (taxon in taxa) {
	
		index <- which(taxonomy_traits$taxa == taxon)
		in_bulk <- taxonomy_traits$in_bulk[index]
		if (any(in_bulk)) taxonomy_traits$in_bulk[index] <- TRUE
	
	}

	# remove rows with duplicated taxa
	dups <- duplicated(taxa_together)
	taxonomy_traits <- taxonomy_traits[!dups, ]

	taxonomy_traits$Domain <- factor(taxonomy_traits$Domain)
	taxonomy_traits$Phylum <- factor(taxonomy_traits$Phylum)
	taxonomy_traits$Class <- factor(taxonomy_traits$Class)

	rownames(taxonomy_traits) <- taxonomy_traits$taxon

	# check that columns of abundances have same names as taxonomic data
	collapsed_abundances_cols <- colnames(collapsed_abundances)
	collapsed_abundances_cols <- collapsed_abundances_cols[collapsed_abundances_cols %notin% c('my.id', 'id')]
	stopifnot(collapsed_abundances_cols == taxonomy_traits$taxon) # no need for manual inspection

	

	

say('DONE!', level = 1, deco = '!')

### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2024-11
###
### This script constructs a network diagram of microbial taxa that co-occur.
###
### source('C:/Ecology/R/andropogon_integratedEcology/hmsc_vs_microbes/community_network.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/hmsc_vs_microbes/community_network.r')
###
### CONTENTS ###
### setup ###
### integrated SDM-ADMIXTURE model ###

#############
### setup ###
#############

	# using the ChordDiag approach to model rhizo and bulk separately
	# https://r-graph-gallery.com/chord-diagram-interactive.html

	# remove everything in R to ensure reproducibility
	rm(list = ls())

	library(dplyr) # dplyr(!)
	library(omnibus) # helper functions
	library(ggraph) # network plots
	library(igraph) # network plots
	library(omnibus) # network plots
	library(RColorBrewer) # colors

say('#####################################################')
say('### network diagram of co-occurrence between taxa ###')
say('#####################################################')

	### user-specified settings
	###########################

		# remove connections between taxa 
		min_assoc_threshold <- 0.5

		# plot log10 of abundance
		# log10_abund <- TRUE
		log10_abund <- FALSE

	### load data
	#############

	# setwd("/Users/iskander/Desktop/New_R_scripts/AGER_data/data_processed/collated_for_hmsc")
	setwd("C:/Ecology/Research/Andropogon/Andropogon")
	# df <- read.csv("abundances_site_by_taxon_combined.csv", header = TRUE)
	df <- read.csv("./data_from_sonny/collated_for_hmsc_collapsed_to_class/abundances_site_by_taxon_combined.csv", header = TRUE)

	# top left, top right corners of data frame
	corner(df)
	corner(df, 2)

	### data cleaning
	#################

	# remove Eukaryota
	taxa <- names(df)
	bads <- grepl(taxa, pattern = 'Eukaryota')
	df <- df[ , !bads]

	# remove unknowns... NB there are no taxa of the pattern "known_unknown_unknown" beyond those listed here
	unknowns <- c(
		'Unknown_unknown_unknown',
		'Bacteria_unknown_unknown',
		'Archaea_unknown_unknown',
		'Unassigned_unknown_unknown',
		'Eukaryota_unknown_unknown'
	)

	taxa <- names(df)
	bads <- taxa %in% unknowns
	df <- df[ , !bads]

	meta_columns <- c('index', 'location', 'plant', 'sample', 'rhizobiome_or_bulk')


	### cycle through each sample type and make the square matrix

	# list to hold square association matrices (one for rhizobiome, one for bulk)
	assocs <- list()

	sample_types <- c('rhizobiome', 'bulk')
	for (sample_type in sample_types) {

		# get sample of the desire type
		samples <- df[df$rhizobiome_or_bulk == sample_type, ]

		sample_names <- names(samples)
		sample_names <- sample_names[sample_names %notin% meta_columns]
		sample_names <- sub(sample_names, pattern = 'Archaea_', replacement = '')
		sample_names <- sub(sample_names, pattern = 'Bacteria_', replacement = '')
		sample_names <- sub(sample_names, pattern = '_', replacement = ' ')
		names(samples)[(length(meta_columns) + 1):ncol(samples)] <- sample_names

		### remove taxa with 0 total abundance
		taxa_columns <- names(samples) %notin% meta_columns
		total_abund <- colSums(samples[ , taxa_columns])
		keeps <- total_abund > 0
		keeps <- c(
			rep(TRUE, length(meta_columns)),
			keeps
		)
		samples <- samples[ , keeps]

		### make edge matrix with "from/to" for each sub-Kingdom pair
		taxa <- names(samples)
		taxa <- taxa[taxa %notin% meta_columns]
		connect <- expand.grid(
			from = taxa,
			to = taxa
		)

		# remove cases where "from" is same as "to"
		keeps <- connect$from != connect$to
		connect <- connect[keeps, ]

		# remove duplicated cases where "from/to" is "A/B" and another "from/to" is "B/A"
		connect <- connect %>%
			rowwise() %>%
			mutate(pair = list(sort(c(from, to)))) %>%
			ungroup() %>%
			distinct(pair, .keep_all = TRUE) %>%
			select(-pair)

		### calculate pairwise associations between taxa
		################################################

		# within a site, take mean number of plants abundance was >0 for both taxa
		# across sites, take sum of the site-level scores

		# vector of sites
		sites <- unique(samples$location)

		# for each possible pairwise association
		connect$value <- NA
		for (count_assoc in 1:nrow(connect)) {
		
			# names of each taxon
			taxon1 <- connect$from[count_assoc]
			taxon2 <- connect$to[count_assoc]

			taxon1 <- as.character(taxon1)
			taxon2 <- as.character(taxon2)
		
			# vector holds mean association score within sites
			assoc_at_sites <- numeric()
			for (site in sites) {

				this_pair <- samples[samples$location == site, c(taxon1, taxon2)]

				this_pair <- this_pair > 0
				this_pair <- rowSums(this_pair)
				
				# both occur on this plant
				this_pair <- this_pair == 2

				site_assoc <- mean(this_pair)
				assoc_at_sites <- c(assoc_at_sites, site_assoc)


			}

			assoc_across_sites <- sum(assoc_at_sites )

			connect$value[count_assoc] <- assoc_across_sites
		
		}

		# rescale association value to [0, 1]
		connect$value <- connect$value / length(sites)

		# remove associations that are less than this value
		connect <- connect[connect$value >= min_assoc_threshold, ]

		### make edges frame
		####################

		taxa <- names(samples)[(length(meta_columns) + 1):ncol(samples)]
		phyla <- substr(taxa, 1, regexpr(taxa, pattern = ' ') - 1)
		phyla_unique <- unique(phyla)

		edges_origin <- data.frame(
			from = 'origin',
			to = phyla_unique
		)

		edges_taxa <- data.frame(
			from = phyla,
			to = taxa
		)

		# edges_taxa$to <- substr(edges_taxa$to, regexpr(taxa, pattern = ' ') + 1, nchar(edges_taxa$to))

		edges <- rbind(
			edges_origin,
			edges_taxa
		)

		### make vertices frame
		#######################

		vertices_origin <- data.frame(
			name = 'origin',
			group = NA
		)

		vertices_phyla <- data.frame(
			name = phyla_unique,
			group = 'origin'
		)

		vertices_taxa <- data.frame(
			name = taxa,
			group = phyla
		)

		vertices <- rbind(vertices_origin, vertices_phyla, vertices_taxa)
		
		# assign values to vertices (using overall abundance)

		taxa <- vertices$name
		abunds <- samples[ , names(samples) %in% taxa]
		abunds <- colSums(abunds)
		if (log10_abund) abunds <- log10(abunds)

		vertices$value[vertices$name %in% names(abunds)] <- abunds
		vertices$id <- NA
		vertices$angle <- NA
		vertices$hjust <- NA

		#Let's add information concerning the label we are going to add: angle, horizontal adjustment and potential flip
		#calculate the ANGLE of the labels
		myleaves <- which(is.na(match(vertices$name, edges$from)))
		nleaves <- length(myleaves)
		vertices$id[myleaves] <- seq(1:nleaves)
		vertices$angle <- 90 - 360 * vertices$id / nleaves	

		# calculate the alignment of labels: right or left
		# If I am on the left part of the plot, my labels have currently an angle < -90
		vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
		
		# flip angle BY to make them readable
		vertices$angle <- ifelse(vertices$angle < -90, vertices$angle + 180, vertices$angle)

		mygraph <- igraph::graph_from_data_frame(edges, vertices = vertices)

		# The connection object must refer to the ids of the leaves:
		from  <-  match(connect$from, vertices$name)
		to  <-  match(connect$to, vertices$name)

		# Basic usual argument
		assoc <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
			geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.9, aes(color=after_stat(index))) +
			scale_edge_colour_distiller(palette = "RdPu") +

			geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust, color=group), size=2, alpha=1) +

			geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2)) +
			scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
			scale_size_continuous( range = c(0.1,10) ) +

			theme_void() +
			theme(
			legend.position="none",
			plot.margin=unit(c(0,0,0,0),"cm"),
			) +
			expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))

		assocs[[length(assocs) + 1]] <- assoc

	print(assoc)
	

	}






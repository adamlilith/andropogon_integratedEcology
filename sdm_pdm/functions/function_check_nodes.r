check_nodes <- function(model) {

	# check for infinite or NaN likelihoods

	# get all nodes
	all_nodes <- model$getNodeNames()

	# check each stochastic node
	stochastic_nodes <- model$getNodeNames(stochOnly = TRUE)
	say('Checking ', length(stochastic_nodes), ' stochastic nodes...')

	for (node in stochastic_nodes) {
		node_calc <- model$calculate(node)
		if (is.na(node_calc) || is.infinite(node_calc)) {
			say('Problem with node: ', node, ' = ', node_calc)
		}
	}

	# check top-level nodes
	top_nodes <- model$getNodeNames(topOnly = TRUE)
	say('Checking ', length(top_nodes), ' top-level nodes...')

	for (node in top_nodes) {
		node_calc <- model$calculate(node)
		if (is.na(node_calc) || is.infinite(node_calc)) {
			say('Problem with top node: ', node, ' = ', node_calc)
		}
	}

	# check data nodes
	data_nodes <- model$getNodeNames(dataOnly = TRUE)
	say('Checking ', length(data_nodes), ' data nodes...')

	for (node in data_nodes) {
		node_calc <- model$calculate(node)
		if (is.na(node_calc) || is.infinite(node_calc)) {
			say('Problem with data node: ', node, ' = ', node_calc)
		}
	}

}

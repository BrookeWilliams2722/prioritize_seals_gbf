objective.function <- function(w1, c, w2, s){
  y <- (w1*c + w2*s)
  return(y + abs(min(y)) + 1)
}

calculate_spp_delta_er <- function(A0, Ac, speciesweights, z=0.25, delta=0.0001){
  # delta_er <- rep(0, length(A0))
  # spp_delta_er <- ((1 - (Ac / A0)^z) - (1 - ((Ac + (delta * A0)) / A0)^z)) / (delta * A0)
  spp_delta_er <- ((1 - (Ac / A0)^z) - (1 - ((Ac + delta) / A0)^z)) / delta
  spp_delta_er[which(spp_delta_er < 0)] <- 0
  return(spp_delta_er * speciesweights)
}

calculate_pu_delta_er_v2 <- function(pu.df, species.pu.list, spp_delta_er, species.pu.idx.list){
  pu_delta_er <- rep(0, length(pu.df$puid))
  for (i in 1:length(spp_delta_er)){
    pu_delta_er[species.pu.idx.list[[i]]] <- pu_delta_er[species.pu.idx.list[[i]]] + spp_delta_er[i]
  }
  return(pu_delta_er)
}

assign_taxon_weights <- function(groups){
	n <- length(unique(groups))
	# assumes 2 animal groups, the rest plants
	animals <- c("amphibians", "birds", "mammals", "reptiles")
	plants <- unique(groups)[which(!unique(groups) %in% animals)]
	#number of animals
	na <- 4
	#number of plants
	np <- n - na
  # in the current case there are no birds, so just spreading the weights among the animal taxon
	wp <- 0 / np
	wa <- 1 / na
	w <- rep(0, length(groups))
	
	for (i in 1:length(animals)){
		w[which(groups == animals[i])] <- wa / length(which(groups == animals[i]))
	}

	for (i in 1:length(plants)){
		w[which(groups == plants[i])] <- wp / length(which(groups == plants[i]))
	}
	return(w)
}

assign_taxon_weights_test <- function(groups){
  n <- length(unique(groups))
  # assumes 2 animal groups, the rest plants
  animals <- c("bird", "mammal")
  plants <- unique(groups)[which(!unique(groups) %in% animals)]
  #number of animals
  na <- 2
  #number of plants
  np <- n - na
  # in the current case there are no birds, so just spreading the weights among the animal taxon
  wp <- 0 / np
  wa <- 1 / na
  w <- rep(0, length(groups))
  
  for (i in 1:length(animals)){
    w[which(groups == animals[i])] <- wa / length(which(groups == animals[i]))
  }
  
  for (i in 1:length(plants)){
    w[which(groups == plants[i])] <- wp / length(which(groups == plants[i]))
  }
  return(w)
}

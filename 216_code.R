library(remotes)
#remotes::install_github("ropensci/rmangal")
library(rmangal)
help(rmangal)
library(rlang)
library(tidygraph)
library(ggraph)
library(igraph)
library(intergraph)
library(network)
library(ergm)
library(philentropy)
library(rlist)
library(popbio)
library(rstanarm)
library(sna)


mgs <- search_networks("Huizache-Caimanero lagoon")
mgn <- get_collection(mgs)
diet_ig <- as.igraph(mgn)

###########################################################
# DATA PREPROCESSING PART 1: Cleaning up Vertex Attributes#
###########################################################

# delete the irrelevant vertex attributes
diet_ig <- delete_vertex_attr(diet_ig, "node_level")
diet_ig <- delete_vertex_attr(diet_ig, "network_id")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy_id")
diet_ig <- delete_vertex_attr(diet_ig, "created_at")
diet_ig <- delete_vertex_attr(diet_ig, "updated_at")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.id")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.ncbi")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.tsn")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.eol")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.bold")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.gbif")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.col")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.created_at")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy.updated_at")
diet_ig <- delete_vertex_attr(diet_ig, "taxonomy")

# remove the N/A values for taxonomy.name and replace with None 
z <- vertex_attr(diet_ig)["taxonomy.name"]
z <- rapply(z, f=function(x) ifelse(is.na(x), "None", x), how="replace" )
vertex_attr(diet_ig)["taxonomy.name"] <- z
vertex_attr(diet_ig)

# remove the N/A values for taxonomy.rank and replace with None 
z <- vertex_attr(diet_ig)["taxonomy.rank"]
z <- rapply(z, f=function(x) ifelse(is.na(x), "None", x), how="replace" )
vertex_attr(diet_ig)["taxonomy.rank"] <- z
vertex_attr(diet_ig)

###########################################################
# DATA PREPROCESSING PART 2: Cleaning up Edge Attributes###
###########################################################

# delete the irrelevant edge attributes
diet_ig <- delete_edge_attr(diet_ig, "date")
diet_ig <- delete_edge_attr(diet_ig, "direction")
diet_ig <- delete_edge_attr(diet_ig, "method")
diet_ig <- delete_edge_attr(diet_ig, "attr_id")
diet_ig <- delete_edge_attr(diet_ig, "public")
diet_ig <- delete_edge_attr(diet_ig, "network_id")
diet_ig <- delete_edge_attr(diet_ig, "created_at")
diet_ig <- delete_edge_attr(diet_ig, "updated_at")
diet_ig <- delete_edge_attr(diet_ig, "attribute.id")
diet_ig <- delete_edge_attr(diet_ig, "attribute.name")
diet_ig <- delete_edge_attr(diet_ig, "attribute.description")
diet_ig <- delete_edge_attr(diet_ig, "attribute.unit")
diet_ig <- delete_edge_attr(diet_ig, "attribute.created_at")
diet_ig <- delete_edge_attr(diet_ig, "attribute.updated_at")

vertex_attr(diet_ig)["name"] <- vertex_attr(diet_ig)["original_name"]

# using incident edges to determine if the organism acts as a herbivore
# based on this data, create a vertex attribute "herbivore" such that
# 0 = node does not function as a herbivore, 1 = node functions as a herbivore
# as a double vector. Likewise for detritivore and predator.

# For ease of visualization, for each node, we will first create subgraphs of that 
# node and its outgoing edges, such that edges will be labeled with consumer type.
# Note that we only look at outgoing ties as opposed to incoming ties,
# because organism consumer type is determined by the outgoing organism that the 
# node organism consumes (not is consumed by). 

pdf("consumer_plots.pdf")

# create a subgraph of diet_ig only looking at outgoing ties of
# the nodes in my_nodes.
# NOTE: Detritus, Phytoplankton, and Macrophytes do not have any outgoing edges 
# so are not in our nodes list
my_nodes <- c('Scianids', 'Elopids', 'Lutjanids', 'Carangids', 'Centropomids',
              'Ariids', 'Haemulids', 'Pleuronectoids', 'Callinectes',
              'Belonoids', 'Clupeoids', 'Gerreids', 'Poeciliids', 'Gobioids',
              'Mugilids', 'Palaemonids', 'Litopenaeus', 'Bivalves', 'Microcrustaceans',
              'Chanids', 'Polychaetes', 'Gastropods', 'Zooplankton')

for (i in my_nodes) {
  diet_ig_sub <- diet_ig
  out_ties <- incident_edges(diet_ig_sub, i, mode="out")
  # convert list object to type double as specified by subgraph.edges API
  out_ties <- do.call(rbind, lapply(out_ties, as.numeric))
  diet_ig_sub <- subgraph.edges(diet_ig_sub, out_ties, delete.vertices = TRUE)
  # plot graph
  vlab <- vertex_attr(diet_ig_sub, name="name")
  elab <- edge_attr(diet_ig_sub, name="type")
  plot(diet_ig_sub, edge.arrow.size = 0.2, vertex.label = vlab,
       vertex.label.cex = 1,
       edge.label = elab, 
       edge.label.cex = 0.7,
       edge.label.color="darkred"
  )
}

dev.off()


herbivore <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 
               0, 0, 1, 0, 0)

detritivore <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 
                 1, 1, 0, 1, 0, 0)

predator <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              1, 0, 1, 0, 0)

herbivore <- unlist(herbivore)
detritivore <- unlist(detritivore)
predator <- unlist(predator)

diet_ig <- set_vertex_attr(graph=diet_ig, name="herbivore", value=herbivore)
diet_ig <- set_vertex_attr(graph=diet_ig, name="detritivore", value=detritivore)
diet_ig <- set_vertex_attr(graph=diet_ig, name="predator", value=predator)


# we now have our network, with added vertex (nodal) attributes of 
# herbivore/detritivore/predator and plot the diet network of
# each organism's diet by the breakdown of its consumption proportion 

pdf("diet_plots.pdf")

# create a subgraph of diet_ig only looking at outgoing ties of
# the nodes in my_nodes.
# NOTE: Detritus, Phytoplankton, and Macrophytes do not have any outgoing edges 
# so are not in our nodes list
my_nodes <- c('Scianids', 'Elopids', 'Lutjanids', 'Carangids', 'Centropomids',
              'Ariids', 'Haemulids', 'Pleuronectoids', 'Callinectes',
              'Belonoids', 'Clupeoids', 'Gerreids', 'Poeciliids', 'Gobioids',
              'Mugilids', 'Palaemonids', 'Litopenaeus', 'Bivalves', 'Microcrustaceans',
              'Chanids', 'Polychaetes', 'Gastropods', 'Zooplankton')

for (i in my_nodes) {
  diet_ig_sub <- diet_ig
  out_ties <- incident_edges(diet_ig_sub, i, mode="out")
  # convert list object to type double as specified by subgraph.edges API
  out_ties <- do.call(rbind, lapply(out_ties, as.numeric))
  diet_ig_sub <- subgraph.edges(diet_ig_sub, out_ties, delete.vertices = TRUE)
  # plot graph
  vlab <- vertex_attr(diet_ig_sub, name="name")
  elab <- edge_attr(diet_ig_sub, name="value")
  plot(diet_ig_sub, edge.arrow.size = 0.2, vertex.label = vlab,
       vertex.label.cex = 1,
       edge.label = elab, 
       edge.label.cex = 0.7,
       edge.label.color="darkred"
  )
}

dev.off()

# create barplots of the diet consumption per node where height denotes proportion of
# diet consumption

scianids_value <- c(0, 0, 0, 0.005, 0.06, 0, 0.025, 0.01, 0.019, 0.07,
                    0.11, 0.051, 0, 0.003, 0.06, 0.033, 0.035, 0, 
                    0.461, 0, 0.041, 0.005, 0.012, 0, 0, 0)

elopids_value <- c(0, 0, 0.038, 0, 0, 0, 0.108, 0, 0, 0.019, 
                   0, 0, 0, 0.063, 0.262, 0.007, 0.134, 0, 0.244, 
                   0, 0.068, 0.043, 0.014, 0, 0, 0)

lutjanids_value <- c(0, 0, 0, 0, 0, 0, 0.041, 0, 0.028, 0, 
                     0.162, 0.085, 0.018, 0.097, 0.201, 0,
                     0.094, 0.047, 0.094, 0, 0.031, 0.023, 
                     0.079, 0, 0, 0)

carangids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0.018, 0.026, 0.101, 
                     0.231, 0, 0.023, 0.103, 0, 0.103, 0.082, 0.263, 
                     0, 0, 0, 0.039, 0, 0, 0.011)

centropomids_value <- c(0.016, 0, 0, 0.005, 0, 0, 0, 0, 0.019, 0, 0.11, 
                        0.055, 0.002, 0.007, 0.097, 0.022, 0.2, 0, 0.225, 
                        0, 0.11, 0.025, 0.055, 0.042, 0, 0.01)

ariids_value <- c(0, 0, 0, 0, 0.012, 0.025, 0.002, 0.006, 0.074, 
                  0.044, 0.089, 0.029, 0, 0.003, 0.005, 0.011, 
                  0.075, 0.047, 0.386, 0.003, 0.042, 0.041, 
                  0.042, 0, 0, 0.064)

haemulids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.005, 
                     0, 0, 0.103, 0.095, 0.044, 0.454, 0, 
                     0.074, 0.052, 0.02, 0.108, 0, 0.045)

pleuronectoids_value <- c(0, 0, 0, 0, 0, 0, 0, 0.041, 0, 0, 0, 
                          0, 0, 0.02, 0, 0, 0, 0.024, 0.514, 
                          0, 0.185, 0.047, 0.169, 0, 0, 0)

callinectes_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0.026, 0,
                       0.025, 0.003, 0, 0, 0.009, 0.009, 
                       0.016, 0.284, 0.231, 0, 0.079, 0.09, 
                       0.197, 0, 0, 0.031)

belonoids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0.1, 0, 0.002, 0.004, 0, 0.006, 0.007, 0, 0.12, 
                     0, 0.035, 0.154, 0, 0.279, 0, 0.293)

clupeoids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0.214, 0, 0, 0, 0.103, 0.485, 0.1, 0.098)

gerreids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.006,
                    0, 0, 0.007, 0.112, 0.073, 0, 0.105, 0.107, 
                    0.175, 0.188, 0, 0.227)

poeciliids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0.291, 0, 0.282, 0, 0, 0, 0, 0.427)

gobioids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0.157, 0, 0.218, 0.047, 0.324, 0.133, 
                    0, 0.121)

mugilids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.005, 
                    0.24, 0, 0.131, 0.116, 0.409, 0.027, 0, 0.072)

palaemonids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0.008, 0, 0.283, 0, 0.499, 0.105, 0, 0.105)

litopenaeus_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.108, 
                       0, 0.269, 0, 0.478, 0, 0, 0.145)

bivalves_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0.177, 0.292, 0.531, 0)

microcrustaceans_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0.097, 0, 0.173, 0, 0.438, 0.012, 
                            0.074, 0.206)

chanids_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0.042, 0.054, 0.091, 0.027, 0.603, 0.183)

polychaetes_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0.044, 0, 0.01, 0, 0.753, 0, 0, 0.193)

gastropods_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0.438, 0, 0, 0.562)

detritus_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

zooplankton_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0.025, 0, 0, 0, 0.193, 0.162, 0.62, 0)

phytoplankton_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

macrophytes_value <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


diet_ig <- set_vertex_attr(graph=diet_ig, name="scianids_value", value=scianids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="elopids_value", value=elopids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="lutjanids_value", value=lutjanids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="carangids_value", value=carangids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="centropomids_value", value=centropomids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="ariids_value", value=ariids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="haemulids_value", value=haemulids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="pleuronectoids_value", value=pleuronectoids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="callinectes_value", value=callinectes_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="belonoids_value", value=belonoids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="clupeoids_value", value=clupeoids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="gerreids_value", value=gerreids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="poeciliids_value", value=poeciliids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="gobioids_value", value=gobioids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="mugilids_value", value=mugilids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="palaemonids_value", value=palaemonids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="litopenaeus_value", value=litopenaeus_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="bivalves_value", value=bivalves_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="microcrustaceans_value", value=microcrustaceans_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="chanids_value", value=chanids_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="polychaetes_value", value=polychaetes_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="gastropods_value", value=gastropods_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="detritus_value", value=detritus_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="zooplankton_value", value=zooplankton_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="phytoplankton_value", value=phytoplankton_value)
diet_ig <- set_vertex_attr(graph=diet_ig, name="macrophytes_value", value=macrophytes_value)

pdf("diet_histograms.pdf")

edges <- vertex_attr(diet_ig, name="scianids_value")
barplot(edges, main="Diet Breakdown: Scianids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="elopids_value")
barplot(edges, main="Diet Breakdown: Elopids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="lutjanids_value")
barplot(edges, main="Diet Breakdown: Lutjanids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="carangids_value")
barplot(edges, main="Diet Breakdown: Carangids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="centropomids_value")
barplot(edges, main="Diet Breakdown: Centropomids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="ariids_value")
barplot(edges, main="Diet Breakdown: Ariids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="haemulids_value")
barplot(edges, main="Diet Breakdown: Haemulids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="pleuronectoids_value")
barplot(edges, main="Diet Breakdown: Pleuronectoids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="callinectes_value")
barplot(edges, main="Diet Breakdown: Callinectes",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="belonoids_value")
barplot(edges, main="Diet Breakdown: Belonoids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="clupeoids_value")
barplot(edges, main="Diet Breakdown: Clupeoids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="gerreids_value")
barplot(edges, main="Diet Breakdown: Gerreids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="poeciliids_value")
barplot(edges, main="Diet Breakdown: Poeciliids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="gobioids_value")
barplot(edges, main="Diet Breakdown: Gobioids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="mugilids_value")
barplot(edges, main="Diet Breakdown: Mugilids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="palaemonids_value")
barplot(edges, main="Diet Breakdown: Palaemonids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="litopenaeus_value")
barplot(edges, main="Diet Breakdown: Litopenaeus",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="bivalves_value")
barplot(edges, main="Diet Breakdown: Bivalves",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="microcrustaceans_value")
barplot(edges, main="Diet Breakdown: Microcrustaceans",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="chanids_value")
barplot(edges, main="Diet Breakdown: Chanids",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="polychaetes_value")
barplot(edges, main="Diet Breakdown: Polychaetes",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="gastropods_value")
barplot(edges, main="Diet Breakdown: Gastropods",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="detritus_value")
barplot(edges, main="Diet Breakdown: Detritus",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="zooplankton_value")
barplot(edges, main="Diet Breakdown: Zooplankton",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="phytoplankton_value")
barplot(edges, main="Diet Breakdown: Phytoplankton",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

edges <- vertex_attr(diet_ig, name="macrophytes_value")
barplot(edges, main="Diet Breakdown: Macrophytes",
        xlab="Consumed organisms",
        ylab="Proportion of diet consumption", 
        ylim=c(0,1))

dev.off()

# We will compute the KL divergences and eigenvector centralities of the diet probability distributions 
# to measure degree of divergence between the consumption patterns of different organisms. 
# 2 identical distributions (ie. the same organism) will have KL divergence = 0. 
# We find the top 3 KL probability distributions with the smallest KL divergence, and the eigenvector centrality for 
# the graph comprised of the KL_nodes and report the results. 

KL_nodes <- list(scianids_value, elopids_value, lutjanids_value, carangids_value, centropomids_value,
              ariids_value, haemulids_value, pleuronectoids_value, callinectes_value,
              belonoids_value, clupeoids_value, gerreids_value, poeciliids_value, gobioids_value,
              mugilids_value, palaemonids_value, litopenaeus_value, bivalves_value, microcrustaceans_value,
              chanids_value, polychaetes_value, gastropods_value, zooplankton_value)


first <- .Machine$integer.max
second <- .Machine$integer.max
third <- .Machine$integer.max
min_i_first <- scianids_value
min_j_first <- elopids_value
min_i_second <- scianids_value
min_j_second <- elopids_value
min_i_third <- scianids_value
min_j_third <- elopids_value
count <- 0

for (i in KL_nodes){
  for (j in KL_nodes) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      tmp_val <- min(KL(x), KL(y))
      if (tmp_val < first) {
        third <- second
        second <- first
        first <- tmp_val
        min_i_third <- min_i_second
        min_j_third <- min_j_second
        min_i_second <- min_i_first
        min_j_second <- min_j_first
        min_i_first <- i
        min_j_first <- j
      }
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}
cat('1st smallest KL divergence: ', first, '\n')
cat('vector i:', min_i_first, '\n')
cat('vector j:', min_j_first, '\n')

cat('2st smallest KL divergence: ', second, '\n')
cat('vector i:', min_i_second, '\n')
cat('vector j:', min_j_second, '\n')

cat('3st smallest KL divergence: ', third, '\n')
cat('vector i:', min_i_third, '\n')
cat('vector j:', min_j_third, '\n')

##########################################################################

my_nodes <- c('Scianids', 'Elopids', 'Lutjanids', 'Carangids', 'Centropomids',
              'Ariids', 'Haemulids', 'Pleuronectoids', 'Callinectes',
              'Belonoids', 'Clupeoids', 'Gerreids', 'Poeciliids', 'Gobioids',
              'Mugilids', 'Palaemonids', 'Litopenaeus', 'Bivalves', 'Microcrustaceans',
              'Chanids', 'Polychaetes', 'Gastropods', 'Zooplankton')
library(igraph)
first <- .Machine$integer.max
second <- .Machine$integer.max
third <- .Machine$integer.max
min_i_first <- ''
min_i_second <- ''
min_i_third <- ''
count <- 0

evcent_ig <- diet_ig
evcent_ig <- delete_vertices(evcent_ig, "Detritus")
evcent_ig <- delete_vertices(evcent_ig, "Phytoplankton")
evcent_ig <- delete_vertices(evcent_ig, "Macrophytes")

for (i in my_nodes) {
  tmp_val <- evcent(evcent_ig)$vector[i]
  print(tmp_val)
  if (tmp_val < first) {
    third <- second
    second <- first
    first <- tmp_val
    min_i_third <- min_i_second
    min_i_second <- min_i_first
    min_i_first <- i
  }
}

cat('1st smallest eigen centrality: ', first, '\n')
cat('organism:', min_i_first, '\n')

cat('2st smallest eigen centrality: ', second, '\n')
cat('organism:', min_i_second, '\n')

cat('3st smallest eigen centrality: ', third, '\n')
cat('organism:', min_i_third, '\n')

###################################################################################


KL_nodes_herbivores <- list(clupeoids_value, bivalves_value, microcrustaceans_value,
                            chanids_value, zooplankton_value)

KL_nodes_detritivores <- list(scianids_value, elopids_value, lutjanids_value, 
                              carangids_value, centropomids_value, ariids_value, 
                              haemulids_value, pleuronectoids_value, callinectes_value, 
                              clupeoids_value, gerreids_value, gobioids_value, mugilids_value,
                              palaemonids_value, litopenaeus_value, bivalves_value, microcrustaceans_value, 
                              chanids_value, polychaetes_value, gastropods_value, 
                              zooplankton_value)

KL_nodes_predator <- list(scianids_value, elopids_value, lutjanids_value, carangids_value,
                          centropomids_value, ariids_value, haemulids_value, pleuronectoids_value,
                          callinectes_value, belonoids_value, clupeoids_value, gerreids_value,
                          poeciliids_value, gobioids_value, mugilids_value, palaemonids_value,
                          litopenaeus_value, bivalves_value, microcrustaceans_value, chanids_value,
                          polychaetes_value, gastropods_value, zooplankton_value)

# We will now compute KL divergences between different consumer types (herbivore/detritivore/predator) as indicated 
# by the data. We will also randomly assign organisms to one of the 3 consumer types and compute the KL divergences 
# for these randomly assigned organisms. Doing this will enable us to determine if consumer type has a 
# correlation with KL divergence or if the results are simply due to random chance. 

# Use indices from all_KL_nodes to randomly assign herbivore, predator, and detritivore KL nodes.

all_KL_nodes <- c('Scianids', 'Elopids', 'Lutjanids', 'Carangids', 'Centropomids',
               'Ariids', 'Haemulids', 'Pleuronectoids', 'Callinectes',
               'Belonoids', 'Clupeoids', 'Gerreids', 'Poeciliids', 'Gobioids',
               'Mugilids', 'Palaemonids', 'Litopenaeus', 'Bivalves', 'Microcrustaceans',
               'Chanids', 'Polychaetes', 'Gastropods', 'Zooplankton')

herb_index_list <- c()
herb_count <- 5

detritivore_index_list <- c()
detritivore_count <- 21

predator_index_list <- c()
predator_count <- 23

while (herb_count > 0) {
  idx <- floor(runif(1, min=1, max=24))
  while (idx %in% herb_index_list) {
    idx <- floor(runif(1, min=1, max=24))
  }
  herb_index_list <- append(herb_index_list, idx)
  herb_count <- herb_count - 1
}

while (detritivore_count > 0) {
  idx <- floor(runif(1, min=1, max=24))
  while (idx %in% detritivore_index_list) {
    idx <- floor(runif(1, min=1, max=24))
  }
  detritivore_index_list <- append(detritivore_index_list, idx)
  detritivore_count <- detritivore_count - 1
}

while (predator_count > 0) {
  idx <- floor(runif(1, min=1, max=24))
  while (idx %in% predator_index_list) {
    idx <- floor(runif(1, min=1, max=24))
  }
  predator_index_list <- append(predator_index_list, idx)
  predator_count <- predator_count - 1
}

herb_index_list # 11 12 21 19 20
detritivore_index_list # 9 21 12 17 7 13 8 5 6 2 11 15 14 22 19 3 23 1 16 4 10
predator_index_list # 16 23 18 7 9 5 22 12 14 8 2 11 21 13 20 4 17 10 3 15 1 6 19

KL_nodes_herbivores_random <- list(clupeoids_value, gerreids_value, polychaetes_value, microcrustaceans_value, chanids_value)

KL_nodes_detritivores_random <- list(callinectes_value, polychaetes_value, gerreids_value, litopenaeus_value, 
                                     haemulids_value, poeciliids_value, pleuronectoids_value, centropomids_value, 
                                     ariids_value, elopids_value, clupeoids_value, mugilids_value, gobioids_value, 
                                     gastropods_value, microcrustaceans_value, lutjanids_value, zooplankton_value, 
                                     scianids_value, palaemonids_value, carangids_value, belonoids_value)

KL_nodes_predator_random <- list(palaemonids_value, zooplankton_value, bivalves_value, haemulids_value, callinectes_value,
                                 centropomids_value, gastropods_value, gerreids_value, gobioids_value, pleuronectoids_value,
                                 elopids_value, clupeoids_value, polychaetes_value, poeciliids_value, chanids_value,
                                 carangids_value, litopenaeus_value, belonoids_value, lutjanids_value, mugilids_value, 
                                 scianids_value, ariids_value, microcrustaceans_value)



# We will now find the average KL divergence within different herbivores and randomly assigned herbivores.

elems <- list()
count <- 0
for (i in KL_nodes_herbivores){
  for (j in KL_nodes_herbivores) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_herbivores_random){
  for (j in KL_nodes_herbivores_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

####################################################################

# We will now find the average KL divergence between herbivores and predators and between randomly assigned herbivores
# and randomly assigned predators. 

elems <- list()
count <- 0
for (i in KL_nodes_herbivores){
  for (j in KL_nodes_predator) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_herbivores_random){
  for (j in KL_nodes_predator_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

####################################################################

# We will now find the average KL divergence between herbivores and detritivores and between randomly assigned
# herbivores and randomly assigned detritivores.

elems <- list()
count <- 0
for (i in KL_nodes_herbivores){
  for (j in KL_nodes_detritivores) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_herbivores_random){
  for (j in KL_nodes_detritivores_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

####################################################################

# We will now find the average KL divergence within different predators and between randomly assigned predators. 

elems <- list()
count <- 0
for (i in KL_nodes_predator){
  for (j in KL_nodes_predator) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_predator_random){
  for (j in KL_nodes_predator_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

####################################################################

# We will now find the average KL divergence between predators and detritivores and between randomly assigned
# predators and randomly assigned detritivores. 

elems <- list()
count <- 0
for (i in KL_nodes_predator){
  for (j in KL_nodes_detritivores) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_predator_random){
  for (j in KL_nodes_detritivores_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

####################################################################

# We will now find the average KL divergence within different detritivores and within randomly assigned detritivores. 

elems <- list()
count <- 0
for (i in KL_nodes_detritivores){
  for (j in KL_nodes_detritivores) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_detritivores_random){
  for (j in KL_nodes_detritivores_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

###################################################################

# We will now find the average eigenvector centrality scores for the herbivore, detritivore and predator only graphs.

evcent_herbivores <- diet_ig
evcent_herbivores <- delete_vertices(evcent_herbivores, "Scianids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Elopids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Lutjanids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Carangids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Centropomids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Ariids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Haemulids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Pleuronectoids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Callinectes")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Belonoids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Gerreids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Poeciliids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Gobioids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Mugilids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Palaemonids")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Litopenaeus")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Polychaetes")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Gastropods")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Detritus")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Phytoplankton")
evcent_herbivores <- delete_vertices(evcent_herbivores, "Macrophytes")

my_nodes <- c("Clupeoids", "Bivalves", "Microcrustaceans", "Chanids", "Zooplankton")

elems <- list()
count <- 0
for (i in my_nodes){
  elems <- list.append(elems, evcent(evcent_herbivores)$vector[i])
}

cat('Mean eigenvector centrality: ', mean(unlist(elems)), '\n')
cat('Standard Deviation eigenvector centrality: ', sd(unlist(elems)), '\n')

####################################################################

evcent_detritivores <- diet_ig
evcent_detritivores <- delete_vertices(evcent_detritivores, "Belonoids")
evcent_detritivores <- delete_vertices(evcent_detritivores, "Poeciliids")
evcent_detritivores <- delete_vertices(evcent_detritivores, "Detritus")
evcent_detritivores <- delete_vertices(evcent_detritivores, "Phytoplankton")
evcent_detritivores <- delete_vertices(evcent_detritivores, "Macrophytes")

my_nodes <- c("Scianids", "Elopids", "Lutjanids", 
                 "Carangids", "Centropomids", "Ariids", 
                 "Haemulids", "Pleuronectoids", "Callinectes", 
                 "Clupeoids", "Gerreids", "Gobioids", "Mugilids",
                 "Palaemonids", "Litopenaeus", "Bivalves", "Microcrustaceans", 
                 "Chanids", "Polychaetes", "Gastropods", 
                 "Zooplankton")

elems <- list()
count <- 0
for (i in my_nodes){
  elems <- list.append(elems, evcent(evcent_detritivores)$vector[i])
}

cat('Mean eigenvector centrality: ', mean(unlist(elems)), '\n')
cat('Standard Deviation eigenvector centrality: ', sd(unlist(elems)), '\n')

##################################################################

evcent_predator <- diet_ig
evcent_predator <- delete_vertices(evcent_predator, "Detritus")
evcent_predator <- delete_vertices(evcent_predator, "Phytoplankton")
evcent_predator <- delete_vertices(evcent_predator, "Macrophytes")


my_nodes <- c("Scianids", "Elopids", "Lutjanids", "Carangids",
                 "Centropomids", "Ariids", "Haemulids", "Pleuronectoids",
                 "Callinectes", "Belonoids", "Clupeoids", "Gerreids",
                 "Poeciliids", "Gobioids", "Mugilids", "Palaemonids",
                 "Litopenaeus", "Bivalves", "Microcrustaceans", "Chanids",
                 "Polychaetes", "Gastropods", "Zooplankton")

elems <- list()
count <- 0
for (i in my_nodes){
  elems <- list.append(elems, evcent(evcent_predator)$vector[i])
}

cat('Mean eigenvector centrality: ', mean(unlist(elems)), '\n')
cat('Standard Deviation eigenvector centrality: ', sd(unlist(elems)), '\n')


########################################################################

KL_nodes_family <- list(scianids_value, lutjanids_value, carangids_value, 
                        centropomids_value, ariids_value, haemulids_value, pleuronectoids_value,
                        belonoids_value, clupeoids_value, gerreids_value, poeciliids_value, 
                        mugilids_value, palaemonids_value, chanids_value)

KL_nodes_genus <- list(elopids_value, callinectes_value, litopenaeus_value)

KL_nodes_class <- list(bivalves_value, polychaetes_value, gastropods_value)

# We will compute the KL divergences between organisms of the 3 taxonomic types (family/genus/class)
# as indicated by the data. We will also randomly assign organisms to one of the 3 taxonomy types and 
# compute the KL divergences for these randomly assigned organisms. Doing this will enable us to determine 
# if taxonomy type has a correlation with KL divergence or if the results are simply due to random chance. 

# Use indices from all_KL_nodes to randomly assign family, genus, and class KL nodes.

all_KL_nodes <- c('Scianids', 'Elopids', 'Lutjanids', 'Carangids', 'Centropomids',
                  'Ariids', 'Haemulids', 'Pleuronectoids', 'Callinectes',
                  'Belonoids', 'Clupeoids', 'Gerreids', 'Poeciliids', 'Gobioids',
                  'Mugilids', 'Palaemonids', 'Litopenaeus', 'Bivalves', 'Microcrustaceans',
                  'Chanids', 'Polychaetes', 'Gastropods', 'Zooplankton')

family_index_list <- c()
family_count <- 14

genus_index_list <- c()
genus_count <- 3

class_index_list <- c()
class_count <- 3

while (family_count > 0) {
  idx <- floor(runif(1, min=1, max=24))
  while (idx %in% family_index_list) {
    idx <- floor(runif(1, min=1, max=24))
  }
  family_index_list <- append(family_index_list, idx)
  family_count <- family_count - 1
}

while (genus_count > 0) {
  idx <- floor(runif(1, min=1, max=24))
  while (idx %in% genus_index_list) {
    idx <- floor(runif(1, min=1, max=24))
  }
  genus_index_list <- append(genus_index_list, idx)
  genus_count <- genus_count - 1
}

while (class_count > 0) {
  idx <- floor(runif(1, min=1, max=24))
  while (idx %in% class_index_list) {
    idx <- floor(runif(1, min=1, max=24))
  }
  class_index_list <- append(class_index_list, idx)
  class_count <- class_count - 1
}

family_index_list # 4 18 2 8 19 23 21 3 12 20 7 17 6 1
genus_index_list # 21 8 18
class_index_list # 13 7 5

KL_nodes_family_random <- list(carangids_value, bivalves_value, elopids_value, pleuronectoids_value,
                               microcrustaceans_value, zooplankton_value, lutjanids_value,
                               gerreids_value, chanids_value, haemulids_value, litopenaeus_value, 
                               ariids_value, scianids_value)

KL_nodes_genus_random <- list(polychaetes_value, pleuronectoids_value, bivalves_value)

KL_nodes_class_random <- list(poeciliids_value, haemulids_value, centropomids_value)

# We will now find the average KL divergence within different organisms of the same family and within 
# randomly assigned organisms of the same family.

elems <- list()
count <- 0
for (i in KL_nodes_family){
  for (j in KL_nodes_family) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_family_random){
  for (j in KL_nodes_family_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

###########################################################################

# We will now find the average KL divergence within different organisms of the same genus and within randomly 
# assigned organisms of the same genus.

elems <- list()
count <- 0
for (i in KL_nodes_genus){
  for (j in KL_nodes_genus) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_genus_random){
  for (j in KL_nodes_genus_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

##############################################################################

# We will now find the average KL divergence within different organisms of the same class and within randomly assigned
# organisms of the same class.

elems <- list()
count <- 0
for (i in KL_nodes_class){
  for (j in KL_nodes_class) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

### randomly assigned ###
elems <- list()
count <- 0
for (i in KL_nodes_class_random){
  for (j in KL_nodes_class_random) {
    if (!isTRUE(all.equal(i, j))) {
      x <- rbind(i, j)
      y <- rbind(j, i)
      elems <- list.append(elems, min(KL(x), KL(y)))
    }
  }
  cat('Iteration', count, 'done.')
  count <- count + 1
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')
cat('Standard Deviation KL divergence: ', sd(unlist(elems)), '\n')

###################################################################

evcent_family <- diet_ig
evcent_family <- delete_vertices(evcent_family, "Callinectes")
evcent_family <- delete_vertices(evcent_family, "Gobioids")
evcent_family <- delete_vertices(evcent_family, "Litopenaeus")
evcent_family <- delete_vertices(evcent_family, "Bivalves")
evcent_family <- delete_vertices(evcent_family, "Microcrustaceans")
evcent_family <- delete_vertices(evcent_family, "Polychaetes")
evcent_family <- delete_vertices(evcent_family, "Gastropods")
evcent_family <- delete_vertices(evcent_family, "Detritus")
evcent_family <- delete_vertices(evcent_family, "Zooplankton")
evcent_family <- delete_vertices(evcent_family, "Phytoplankton")
evcent_family <- delete_vertices(evcent_family, "Macrophytes")

my_nodes <- c("Scianids", "Lutjanids", "Carangids", 
                "Centropomids", "Ariids", "Haemulids", "Pleuronectoids",
                 "Belonoids", "Clupeoids", "Gerreids", "Poeciliids", 
                 "Mugilids", "Palaemonids", "Chanids")

elems <- list()
count <- 0
for (i in my_nodes){
  elems <- list.append(elems, evcent(evcent_family)$vector[i])
}

cat('Mean eigenvector centrality: ', mean(unlist(elems)), '\n')
cat('Standard Deviation eigenvector centrality: ', sd(unlist(elems)), '\n')

###################################################################

evcent_genus <- diet_ig
evcent_genus <- delete_vertices(evcent_genus, "Scianids")
evcent_genus <- delete_vertices(evcent_genus, "Lutjanids")
evcent_genus <- delete_vertices(evcent_genus, "Carangids")
evcent_genus <- delete_vertices(evcent_genus, "Centropomids")
evcent_genus <- delete_vertices(evcent_genus, "Ariids")
evcent_genus <- delete_vertices(evcent_genus, "Haemulids")
evcent_genus <- delete_vertices(evcent_genus, "Pleuronectoids")
evcent_genus <- delete_vertices(evcent_genus, "Belonoids")
evcent_genus <- delete_vertices(evcent_genus, "Clupeoids")
evcent_genus <- delete_vertices(evcent_genus, "Gerreids")
evcent_genus <- delete_vertices(evcent_genus, "Poeciliids")
evcent_genus <- delete_vertices(evcent_genus, "Gobioids")
evcent_genus <- delete_vertices(evcent_genus, "Mugilids")
evcent_genus <- delete_vertices(evcent_genus, "Palaemonids")
evcent_genus <- delete_vertices(evcent_genus, "Bivalves")
evcent_genus <- delete_vertices(evcent_genus, "Microcrustaceans")
evcent_genus <- delete_vertices(evcent_genus, "Chanids")
evcent_genus <- delete_vertices(evcent_genus, "Polychaetes")
evcent_genus <- delete_vertices(evcent_genus, "Gastropods")
evcent_genus <- delete_vertices(evcent_genus, "Detritus")
evcent_genus <- delete_vertices(evcent_genus, "Zooplankton")
evcent_genus <- delete_vertices(evcent_genus, "Phytoplankton")
evcent_genus <- delete_vertices(evcent_genus, "Macrophytes")

my_nodes <- c("Elopids", "Callinectes", "Litopenaeus")

elems <- list()
count <- 0
for (i in my_nodes){
  elems <- list.append(elems, evcent(evcent_genus)$vector[i])
}

cat('Mean eigenvector centrality: ', mean(unlist(elems)), '\n')
cat('Standard Deviation eigenvector centrality: ', sd(unlist(elems)), '\n')

###################################################################

evcent_class <- diet_ig
evcent_class <- delete_vertices(evcent_class, "Scianids")
evcent_class <- delete_vertices(evcent_class, "Elopids")
evcent_class <- delete_vertices(evcent_class, "Lutjanids")
evcent_class <- delete_vertices(evcent_class, "Carangids")
evcent_class <- delete_vertices(evcent_class, "Centropomids")
evcent_class <- delete_vertices(evcent_class, "Ariids")
evcent_class <- delete_vertices(evcent_class, "Haemulids")
evcent_class <- delete_vertices(evcent_class, "Pleuronectoids")
evcent_class <- delete_vertices(evcent_class, "Callinectes")
evcent_class <- delete_vertices(evcent_class, "Belonoids")
evcent_class <- delete_vertices(evcent_class, "Clupeoids")
evcent_class <- delete_vertices(evcent_class, "Gerreids")
evcent_class <- delete_vertices(evcent_class, "Poeciliids")
evcent_class <- delete_vertices(evcent_class, "Gobioids")
evcent_class <- delete_vertices(evcent_class, "Mugilids")
evcent_class <- delete_vertices(evcent_class, "Palaemonids")
evcent_class <- delete_vertices(evcent_class, "Litopenaeus")
evcent_class <- delete_vertices(evcent_class, "Microcrustaceans")
evcent_class <- delete_vertices(evcent_class, "Chanids")
evcent_class <- delete_vertices(evcent_class, "Detritus")
evcent_class <- delete_vertices(evcent_class, "Zooplankton")
evcent_class <- delete_vertices(evcent_class, "Phytoplankton")
evcent_class <- delete_vertices(evcent_class, "Macrophytes")

my_nodes <- c("Bivalves", "Polychaetes", "Gastropods")

elems <- list()
count <- 0
for (i in my_nodes){
  elems <- list.append(elems, evcent(evcent_class)$vector[i])
}

cat('Mean eigenvector centrality: ', mean(unlist(elems)), '\n')
cat('Standard Deviation eigenvector centrality: ', sd(unlist(elems)), '\n')


######################################################################


#########################
# KL DIVERGENCE RESULTS #
#########################

# The smallest possible KL divergence for all graphs is 0.1073854. 
# Bivalves and Zooplankton have the closest consumption diet proportion compared
# to all other organisms. The 2nd smallest KL divergence for all graphs is 0.174856.
# Litopenaeus and Microcrustaceans have the 2nd closest consumption diet proportion
# compared to all other organisms. The 3rd smallest KL divergence for all graphs
# is 0.268147. Gobioids and Mugilids have the 3rd closest consumption diet proportion
# compared to all other organisms. 

# The mean KL divergence within herbivores is 1.71832 +/- 0.8633437
# The mean KL divergence between herbivores and predators is 5.78527 +/- 3.94843
# The mean KL divergence between herbivores and detritivores is 5.661169 +/- 3.914441
# The mean KL divergence within predators is 4.442644 +/- 3.371918
# The mean KL divergence between predators and detritivores is 4.407891 +/- 3.391976
# The mean KL divergence between detritivores is 4.35195 +/- 3.413382

# The mean KL divergence within randomly assigned herbivores is 2.197416 +/- 0.8321091
# The mean KL divergence between randomly assigned herbivores and randomly assigned predators is 3.971786 +/- 3.09513
# The mean KL divergence between randomly assigned herbivores and randomly assigned detritivores is 4.031168 +/- 3.159512
# The mean KL divergence within randomly assigned predators is 4.442644 +/- 3.371918
# The mean KL divergence between randomly assigned predators and randomly assigned detritivores is 4.158786 +/- 3.119217
# The mean KL divergence within randomly assigned detritivores is 3.802122 +/- 2.743471

# The mean KL divergence within the family taxonomy is 4.046818 +/- 2.907662
# The mean KL divergence within the genus taxonomy is 3.541938 +/- 1.750825
# The mean KL divergence within the class taxonomy is 4.972296 +/- 3.987272
# These results show that the more closely related the organisms are, the
# more similar their diet consumption proportions are. 

# The mean KL divergence within randomly assigned taxonomic family organisms is 5.401385 +/- 3.93155
# The mean KL divergence within randomly assigned taxonomic genus organisms is 6.984435 +/- 4.080252
# The mean KL divergence within randomly assigned taxonomic class organisms is 1.967317 +/- 0.6708289
a <- dnorm(4.3520, 3.4134)
b <- dnorm(3.8021, 2.7435)
x <- rbind(a,b)
y <- rbind(b,a)
max(KL(x), KL(y))

##################################
# EIGENVECTOR CENTRALITY RESULTS #
##################################

# The smallest possible eigenvector centrality is 0.2014359 from Chanids. 
# The 2nd smallest eigenvector centrality is 0.3166119 from Poeciliids.
# The 3rd smallest eigenvector centrality is 0.3892605 from Pleuronectoids.

# The mean eigenvector centrality for herbivores is 0.582667 +/- 0.3467036
# The mean eigenvector centrality for detritivores is 0.6027696 +/- 0.1763133
# The mean eigenvector centrality for predators is 0.5781874 +/- 0.1763821

# The mean eigenvector centrality for the family taxonomy is 0.6293961 +/- 0.2415005
# The mean eigenvector centrality for the genus taxonomy is 0.7489932 +/- 0.2812418
# The mean eigenvector centrality for the class taxonomy is 0.3333333 +/- 0.5773503

##################################################################################

# We will now find the average KL divergence of one organism to all other organisms 
# and store it in avg_KL_divergence

KL_nodes <- list(scianids_value, elopids_value, lutjanids_value, carangids_value, centropomids_value,
                 ariids_value, haemulids_value, pleuronectoids_value, callinectes_value,
                 belonoids_value, clupeoids_value, gerreids_value, poeciliids_value, gobioids_value,
                 mugilids_value, palaemonids_value, litopenaeus_value, bivalves_value, microcrustaceans_value,
                 chanids_value, polychaetes_value, gastropods_value, zooplankton_value)

elems <- list()
i <- zooplankton_value
for (j in KL_nodes) {
  if (!isTRUE(all.equal(i, j))) {
    x <- rbind(i, j)
    y <- rbind(j, i)
    elems <- list.append(elems, min(KL(x), KL(y)))
  }
}

cat('Mean KL divergence: ', mean(unlist(elems)), '\n')


avg_KL_divergence <- c(5.468417, 6.0859, 5.515445, 5.761814, 3.399128, 3.674893, 3.515657,
                       3.876666, 3.535125, 5.15708, 4.460764, 3.263526, 4.458027, 2.74795,
                       2.496422, 3.205786, 3.055368, 8.633731, 2.722405, 6.212586, 3.199651,
                       4.837615, 6.896864)

stan_ig_model <- diet_ig

stan_ig_model <- delete_vertices(stan_ig_model, "Detritus")
stan_ig_model <- delete_vertices(stan_ig_model, "Phytoplankton")
stan_ig_model <- delete_vertices(stan_ig_model, "Macrophytes")
stan_ig_model <- set_vertex_attr(graph=stan_ig_model, name="avg_KL_divergence", value=avg_KL_divergence)

stan_df_model <- as_long_data_frame(stan_ig_model)

stan_herbivore <- vertex_attr(stan_ig_model, name="herbivore")
stan_detritivore <- vertex_attr(stan_ig_model, name="detritivore")
stan_predator <- vertex_attr(stan_ig_model, name="predator")

help(stan_glm)
fit1 <- stan_glm(avg_KL_divergence ~ stan_herbivore + stan_detritivore + stan_predator
                 + stan_herbivore:stan_detritivore + stan_herbivore:stan_predator 
                 + stan_detritivore:stan_predator, data = stan_df_model)
summary(fit1)

# Model Info:
#   function:     stan_glm
# family:       gaussian [identity]
# formula:      avg_KL_divergence ~ stan_herbivore + stan_detritivore + stan_predator + 
#   stan_herbivore:stan_detritivore + stan_herbivore:stan_predator + 
#   stan_detritivore:stan_predator
# algorithm:    sampling
# sample:       4000 (posterior sample size)
# priors:       see help('prior_summary')
# observations: 23
# predictors:   7
# 
# Estimates:
#   mean          sd            10%           50%        
#   (Intercept)                               0.2          16.6         -20.9           0.3
# stan_herbivore                           -0.1           4.1          -5.2          -0.1
# stan_detritivore                          0.0           3.9          -4.9           0.0
# stan_herbivore:stan_detritivore           0.0           3.9          -5.0          -0.1
# stan_herbivore:stan_predator              0.0           4.1          -5.2          -0.1
# stan_detritivore:stan_predator            0.0           3.9          -5.1           0.0
# sigma                           29760659334.8      125360.5 29760501608.7 29760659041.0
# 90%        
# (Intercept)                              21.5
# stan_herbivore                            5.2
# stan_detritivore                          4.9
# stan_herbivore:stan_detritivore           5.1
# stan_herbivore:stan_predator              5.3
# stan_detritivore:stan_predator            5.0
# sigma                           29760817959.3
# 
# Fit Diagnostics:
#   mean          sd            10%           50%           90%        
# mean_PPD   -35682163.8  6285023844.0 -7989016967.5  -143293569.9  8036293178.7
# 
# The mean_ppd is the sample average posterior predictive distribution of the outcome variable (for details see help('summary.stanreg')).
# 
# MCMC diagnostics
# mcse        Rhat        n_eff
# (Intercept)                             0.3         1.0 4085 
# stan_herbivore                          0.1         1.0 3683 
# stan_detritivore                        0.1         1.0 4002 
# stan_herbivore:stan_detritivore         0.1         1.0 3430 
# stan_herbivore:stan_predator            0.1         1.0 2935 
# stan_detritivore:stan_predator          0.1         1.0 3363 
# sigma                                2017.7         1.0 3860 
# mean_PPD                        100383178.8         1.0 3920 
# log-posterior                           0.1         1.0 1489 
# 
# For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).

fit2 <- stan_glm(log(avg_KL_divergence) ~ stan_herbivore + stan_detritivore + stan_predator
                 + stan_herbivore:stan_detritivore + stan_herbivore:stan_predator 
                 + stan_detritivore:stan_predator, data = stan_df_model)
summary(fit2)

# Model Info:
#   function:     stan_glm
# family:       gaussian [identity]
# formula:      log(avg_KL_divergence) ~ stan_herbivore + stan_detritivore + 
#   stan_predator + stan_herbivore:stan_detritivore + stan_herbivore:stan_predator + 
#   stan_detritivore:stan_predator
# algorithm:    sampling
# sample:       4000 (posterior sample size)
# priors:       see help('prior_summary')
# observations: 23
# predictors:   7
# 
# Estimates:
#   mean         sd           10%          50%          90%       
# (Intercept)                              0.1          3.6         -4.5          0.1          4.8
# stan_herbivore                           0.0          0.9         -1.1          0.0          1.1
# stan_detritivore                         0.0          0.9         -1.1          0.0          1.1
# stan_herbivore:stan_detritivore          0.0          0.8         -1.0          0.0          1.1
# stan_herbivore:stan_predator             0.0          0.8         -1.1          0.0          1.1
# stan_detritivore:stan_predator           0.0          0.8         -1.1          0.0          1.0
# sigma                           2877616773.3      18000.1 2877593960.2 2877616812.6 2877640219.4
# 
# Fit Diagnostics:
#   mean         sd           10%          50%          90%       
# mean_PPD    2011507.6  594376832.0 -786534690.2   10043975.5  761327331.9
# 
# The mean_ppd is the sample average posterior predictive distribution of the outcome variable (for details see help('summary.stanreg')).
# 
# MCMC diagnostics
# mcse      Rhat      n_eff
# (Intercept)                           0.1       1.0 3907 
# stan_herbivore                        0.0       1.0 4268 
# stan_detritivore                      0.0       1.0 3783 
# stan_herbivore:stan_detritivore       0.0       1.0 3792 
# stan_herbivore:stan_predator          0.0       1.0 3533 
# stan_detritivore:stan_predator        0.0       1.0 3649 
# sigma                               285.6       1.0 3973 
# mean_PPD                        9603644.6       1.0 3830 
# log-posterior                         0.0       1.0 1770 
# 
# For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).



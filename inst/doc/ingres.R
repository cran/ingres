## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- dependencies, eval=FALSE------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("viper")
#  BiocManager::install("org.Hs.eg.db")
#  BiocManager::install("AnnotationDbi")
#  BiocManager::install("aracne.networks") # optional dependency

## ----setup--------------------------------------------------------------------
library(ingres)

## ----seurat-------------------------------------------------------------------
small_blca_wang


## ----network------------------------------------------------------------------
network

## ----network_conversion-1-----------------------------------------------------
# Load the graphml example file
filename = system.file("extdata", "example_network.graphml", package = "ingres")

#Convert to tidygraph format
graphmlAsTidy(filename)


## ----network-conversion-2-----------------------------------------------------
# Load the GinSim (.zginml) example file
filename = system.file("extdata", "example_ginsim.zginml", package = "ingres")

# Convert to graphml and store it in a temporary file
temp = tempfile()
gml = ginmlToGraphml(filename, dest = temp)
head(gml)

# Convert to tidygraph
graphmlAsTidy(temp)

## ----symbols------------------------------------------------------------------
network_genes

## ----network-genes-template---------------------------------------------------
# store and modify = F just for demonstration.
createNetworkGenesTemplate(network, store = F, modify = F)
# Then, modify as needed.

## ----object-creation----------------------------------------------------------
ing = createIngresObjectFromSeurat(
  seurat.object = small_blca_wang,
  seurat.assay  = "RNA",
  slot          = "data",
  network.genes = network_genes,
  network       = network
)

ing

## ----matrix-idents------------------------------------------------------------
exp = ing@expression
exp[1:2, 1:2]

idents = ing@idents
head(idents)

## ----create-no-seurat---------------------------------------------------------
createIngresObject(
  expression.matrix = exp, 
  idents            = idents,   
  network.genes     = network_genes,
  network           = network)


## ---- include = FALSE---------------------------------------------------------
#needed for following chunks to run, sometimes not readily installed by automatic checks
knitr::opts_chunk$set(
  eval=requireNamespace("aracne.networks")
)

## ----viper--------------------------------------------------------------------
# Using a single regulon for speed
ing = performViper(ing, aracne.networks::regulonblca)
ing

## ----compute-cell-pbn---------------------------------------------------------
ing = computePbnBySingleCell(ing)
ing
head(ing@single.cell.pbn)

head(computePbnBySingleCell(ing, c(-0.5, 0.5))@single.cell.pbn)

## ----compute-cluster-pbn------------------------------------------------------
ing = computePbnByCluster(ing)
ing


## ----produce-network, eval=FALSE----------------------------------------------
#  produceNetworkForCell(ing, "sample1@ACAGCTAAGATCCCGC-1")
#  
#  produceNetworkForCluster(ing, "1")

## ----network-plotting, fig.height=15, fig.width=15, out.width='100%'----------
cellPbnPlot(ing, "sample1@ACAGCTAAGATCCCGC-1")

clusterPbnPlot(ing, "1")

## ----heatmaps, fig.height=5, fig.width=10, out.width='100%'-------------------
cellGenesHeatmap(ing)

clusterGenesHeatmap(ing)


## ----boolnet------------------------------------------------------------------
ing_network = produceNetworkForCluster(ing, "1")

boolnet_network = produceBoolnetNetwork(ing_network)

# Generate a initial state where all nodes have state 1
initial_state = rep(1, length(network))

# Compute one state transition from the given initial state.
transition1 = BoolNet::stateTransition(
  boolnet_network, 
  type  = "probabilistic",
  state = initial_state)

head(transition1)



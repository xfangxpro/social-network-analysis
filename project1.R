# Project 1
# George Fang
# Sharavani Deiva Senapathy

library(igraph)
library(netrw)
library(MASS)
library(ggplot2)

input_file = "/Users/georgef/Documents/grad_school/spring_2016/EE_232E/project1/facebook_combined.txt"


#############################################################################################
#                                         QUESTION 1                                        #
#############################################################################################
# Generate graph based on input file
facebook.graph = read.graph(file = input_file, 
                            format = "edgelist",
                            directed = FALSE)

# Keep original vertex numbers since induced subgraphs lose this information
V(facebook.graph)$id = 1:vcount(facebook.graph)
arp = vcount(facebook.graph)

# Check connectivity
is.conn = is.connected(facebook.graph)

# Check diameter
diam = diameter(facebook.graph)

# Determine the degree distribution
deg.dist = degree.distribution(facebook.graph)

# Plot the degree distribution
graph.degree = degree(facebook.graph)
deg.dist.hist = hist(degree(facebook.graph),
                     main = "Facebook Graph Degree Distribution",
                     breaks = seq(0.0, 
                                  by = 1, 
                                  length.out = max(graph.degree) + 2),
                     xlab = "Number of nodes",
                     ylab = "Degree Frequency")       

plot(deg.dist, 
     main="Facebook Graph Degree Distribution", 
     xlab = "Number of nodes", 
     ylab = "Degree Density",
     type = "o")

# Determine the average degree
mean.degree = mean(graph.degree) 

# Fit in a curve
degree.data.frame = data.frame(x=deg.dist.hist$mids, y=deg.dist.hist$density)

model = nls(y ~ I(exp(1)^(a + b * x)), 
              data = degree.data.frame, 
              start = list(a=0,b=0))

# Plot the fitted curve
fitted.curve = ggplot(degree.data.frame, 
                      aes(x, y)) + 
                      geom_point(size = 2) + 
                      stat_smooth(method = "nls", 
                                  formula = as.formula(model), 
                                  data = degree.data.frame, 
                                  start = list(a=0,b=0), 
                                  size = 1, 
                                  se = FALSE, 
                                  colour = "red") + 
                      xlab("Number of nodes") +
                      ylab("Degree Density") +
                      ggtitle("Best fitted curve")

print(fitted.curve)

# Determine the total mean squared error using residuals
mse = mean(residuals(model)^2)


#############################################################################################
#                                         QUESTION 2                                        #
#############################################################################################
# Determine node 1's neighbors
node.one.neighbors = neighborhood(facebook.graph, 
                                  order = 1,
                                  node = 1)

# Find the induced neighborhood subgraph
facebook.subgraph = induced.subgraph(facebook.graph, 
                        which(((1:vcount(facebook.graph)) %in% node.one.neighbors[[1]])))

# Determine the number of vertices and edges in the induced subgraph
subgraph.vertex.count = vcount(facebook.subgraph)
subgraph.edge.count = ecount(facebook.subgraph)
plot(degree.distribution(facebook.subgraph),
     main = "Personal Network Degree Distribution",
     xlab = "Number of nodes",
     ylab = "Degree Density")

subgraph.degree = degree(facebook.subgraph)
sub.deg.dist.hist = hist(degree(facebook.subgraph),
                         main = "Node 1 Neighborhood Degree Distribution",
                         breaks = seq(0.0, 
                                      by = 1, 
                                      length.out = max(subgraph.degree) + 2),
                         xlab = "Number of nodes", 
                         ylab = "Degree Frequency")   

# Visualize the facebook subgraph (personal network)
vcolor = rep("magenta", vcount(facebook.subgraph))
vcolor[1] = "black"
vsize = rep(2, vcount(facebook.subgraph))
vsize[1] = 4
plot(facebook.subgraph, 
     vertex.color = vcolor, 
     vertex.size = vsize, 
     vertex.label = NA)


#############################################################################################
#                                         QUESTION 3                                        #
#############################################################################################
# Generate a list of nodes that have a neighborhood size > 200
core.nodes = which(neighborhood.size(facebook.graph, 
                                     order = 1, 
                                     nodes = V(facebook.graph)) > 200)

# Determine neighborhood size of each core node (DIAGNOSTIC PURPOSES)
core.neighborhood.size = neighborhood.size(facebook.graph, 
                                           order = 1, 
                                           nodes = core.nodes)

# Determine the number of core nodes
num.core.nodes = length(core.nodes)

# Find the mean degree of core nodes
core.mean.degree = mean(graph.degree[core.nodes])

# Find the induced subgraph for randomly selected core node
chosen.core.index = sample(1:num.core.nodes, 1)

# Try the first node for simplicity sake
chosen.core = core.nodes[1]

# Determine the chosen core's neighbors
chosen.core.neighbors = neighborhood(facebook.graph, 
                                     order = 1,
                                     nodes = chosen.core)

# Find the induced neighborhood subgraph
chosen.core.subgraph = induced.subgraph(facebook.graph, 
                                     which(((1:vcount(facebook.graph)) %in% chosen.core.neighbors[[1]])))

# Use fast greedy community
core.fast.com = fastgreedy.community(chosen.core.subgraph)
plot(core.fast.com,
     chosen.core.subgraph,
     main = "Node 1 Core Communities (Fast Greedy)",
     vertex.size = 4,
     vertex.label = NA)
hist(core.fast.com$membership,
     col="orange",
     main="Community Distribution (Fast greedy)",
     xlab="Community Number",
     ylab="Number of Nodes per Community")

# Use edge betweenness community
core.edge.com = edge.betweenness.community(chosen.core.subgraph)
plot(core.edge.com,
     chosen.core.subgraph,
     main = "Node 1 Core Communities (Edge Betweenness)",
     vertex.size = 4,
     vertex.label = NA)
hist(core.edge.com$membership,
     col="red",
     main="Community Distribution (Edge Betweenness)",
     xlab="Community Number",
     ylab="Number of Nodes per Community")

# Use infomap community
core.info.com = infomap.community(chosen.core.subgraph)
plot(core.info.com,
     chosen.core.subgraph,
     main = "Node 1 Core Communities (Info Map)",
     vertex.size = 4,
     vertex.label = NA)
hist(core.info.com$membership,
     col="blue",
     main="Community Distribution (Info Map)",
     xlab="Community Number",
     ylab="Number of Nodes per Community")


#############################################################################################
#                                         QUESTION 4                                        #
############################################################################################# 
# Remove the core node from the induced subnetwork
coreless.graph = delete.vertices(chosen.core.subgraph, 
                                 v = which(V(chosen.core.subgraph)$id == chosen.core))

# Run community algorithms once again
# Use fast greedy community
coreless.fast.com = fastgreedy.community(coreless.graph)
plot(coreless.fast.com,
     coreless.graph,
     main = "Removed Core Personal Network (Fast Greedy)",
     vertex.size = 4,
     vertex.label = NA)
hist(coreless.fast.com$membership,
     col="orange",
     main="Core-less Community Distribution (Fast greedy)",
     xlab="Community Number",
     ylab="Number of Nodes per Community")

# Use edge betweenness community
coreless.edge.com = edge.betweenness.community(coreless.graph)
plot(coreless.edge.com,
     coreless.graph,
     main = "Removed Core Personal Network (Edge Betweenness)",
     vertex.size = 4,
     vertex.label = NA)
hist(coreless.edge.com$membership,
     col="red",
     main="Core-less Community Distribution (Edge Betweenness)",
     xlab="Community Number",
     ylab="Number of Nodes per Community")

# Use infomap community
coreless.info.com = infomap.community(coreless.graph)
plot(coreless.info.com,
     coreless.graph,
     main = "Removed Core Personal Network (Info Map)",
     vertex.size = 4,
     vertex.label = NA)
hist(coreless.info.com$membership,
     col="blue",
     main="Core-less Community Distribution (Info Map)",
     xlab="Community Number",
     ylab="Number of Nodes per Community")


#############################################################################################
#                                         QUESTION 5                                        #
#############################################################################################
# Function will return the common neighbors
common.neighbors = function(g, u, v) {
  u.neighbors = neighborhood(g,
                             order = 1,
                             nodes = u)[[1]][-1]
  v.neighbors = neighborhood(g,
                             order = 1,
                             nodes = v)[[1]][-1]
  intersection = intersect(u.neighbors,
                           v.neighbors)
  return(intersection)
}

# Embeddedness is number of common neighbors
compute.embeddedness = function(g, u, v) {
  embeddedness = length(common.neighbors(g, u, v))
  return(embeddedness)
}

# Compute the dispersion
compute.dispersion = function(g, u, v) {

  dispersion = 0
  common.u.v = common.neighbors(g, u, v)
  without.u.v = delete.vertices(g, v = c(u, v))

  for (s in common.u.v) {
    for (t in common.u.v) {
      if(s != t &&
         !is.na(match(s, V(without.u.v))) &&
         !is.na(match(t, V(without.u.v))) &&
         !are.connected(without.u.v, s, t) &&
         length(common.neighbors(without.u.v, s, t)) == 0) {

            dispersion = dispersion + 1 
      }
    }
  }
  return(dispersion)
}

# Simplify subgraph creations for a given node (ALTERNATE)
# create.personal.network = function(g, node) {
#   node.neighborhood = neighborhood(facebook.graph, 
#                                    order = 1,
#                                    nodes = node)
#   new.subgraph = induced.subgraph(g, 
#                                   which(((1:vcount(g)) %in% node.neighborhood[[1]])))
#   new.subgraph$name = sort(node.neighborhood[[1]])
#   return(new.subgraph)
# }

# Simplify subgraph creations for a given node
create.personal.network = function(g, node) {
  relevant.nodes = neighborhood(g, 
                                order = 1 ,
                                nodes = node)[[1]]
  nonrelevant.nodes = which( !( (1:vcount(g)) %in% relevant.nodes)  )
  new.subgraph = delete.vertices(g , nonrelevant.nodes)
  new.subgraph$name =  sort(relevant.nodes)
  return(new.subgraph)
}


dispersion.embeddedness.ratio <- function(g,coreNode) {
  high.disp.value = 0
  disp.node = 0
  high.ratio.value = 0
  ratio.node = 0
  high.embed.value = 0
  embed.node = 0
  personal.network.u = create.personal.network(g, coreNode)
  u = which(personal.network.u$name == coreNode)
  
  nodes <- V(personal.network.u)
  for(v in nodes) {
    if(v == u)
      next
    
    disp = compute.dispersion(g,u,v)
    embd = compute.embeddedness(g,u,v)
    
    if (embd > 0)
    {
      rt = disp/embd
      if (rt > high.ratio.value)
      {
        ratio.node = v;
        high.ratio.value=rt;
      }
    }
    
    if (disp > high.disp.value)
    {
      disp.node = v;
      high.disp.value=disp;
    }
    if (embd > high.embed.value)
    {
      embed.node = v
      high.embed.value=embd;
    }
    
  }
  if (disp.node > 0)
  {
    # community detection
    fc = fastgreedy.community(personal.network.u)
    sizes(fc)
    mfc = membership(fc)
    sizeVet = rep(3, length(V(personal.network.u)));
    sizeVet[disp.node] = 8;  
    colEd = rep(8, length(E(personal.network.u)));
    colEd[which(get.edgelist(personal.network.u,name=F)[,1] == disp.node | get.edgelist(personal.network.u,name=F)[,2] == disp.node)] = 3;
    E(personal.network.u)$color = colEd;
    widEd = rep(1, length(E(personal.network.u)));
    widEd[which(get.edgelist(personal.network.u,name=F)[,1] == disp.node | get.edgelist(personal.network.u,name=F)[,2] == disp.node)] = 3;
    dev.new ();
    plot(personal.network.u, vertex.label= NA, vertex.color=mfc,vertex.size=sizeVet, edge.width = widEd,mark.groups = by(seq_along(mfc), mfc, invisible),main="Max dispersion");
  }
  else
  {
    print (paste(c("No high dispersion node", toString(coreNode)), collapse=" "));
  }
  
  
  if (embed.node > 0)
  {
    # community detection
    fc = fastgreedy.community(personal.network.u); sizes(fc)
    mfc = membership(fc)
    sizeVet = rep(3, length(V(personal.network.u)));
    sizeVet[embed.node] = 8;  
    colEd = rep(8, length(E(personal.network.u)));
    colEd[which(get.edgelist(personal.network.u,name=F)[,1] == embed.node | get.edgelist(personal.network.u,name=F)[,2] == embed.node)] = 3;
    E(personal.network.u)$color = colEd;
    widEd = rep(1, length(E(personal.network.u)));
    widEd[which(get.edgelist(personal.network.u,name=F)[,1] == embed.node | get.edgelist(personal.network.u,name=F)[,2] == embed.node)] = 3;
    dev.new ();
    plot(personal.network.u, vertex.label= NA, vertex.color=mfc,vertex.size=sizeVet, edge.width = widEd,mark.groups = by(seq_along(mfc), mfc, invisible),main="Max embeddedness");# ,mark.groups = by(seq_along(mfc), mfc) );
  }
  else
  {
    print (paste(c("No high embeddedness node", toString(coreNode)), collapse=" "));
  }
  
  if (ratio.node > 0)
  {
    
    # community detection
    fc = fastgreedy.community(personal.network.u); sizes(fc)
    mfc = membership(fc)
    
    sizeVet = rep(3, length(V(personal.network.u)));
    sizeVet[ratio.node] = 8;  
    colEd = rep(8, length(E(personal.network.u)));
    colEd[which(get.edgelist(personal.network.u,name=F)[,1] == ratio.node | get.edgelist(personal.network.u,name=F)[,2] == ratio.node)] = 3;
    E(personal.network.u)$color = colEd;
    widEd = rep(1, length(E(personal.network.u)));
    widEd[which(get.edgelist(personal.network.u,name=F)[,1] == ratio.node | get.edgelist(personal.network.u,name=F)[,2] == ratio.node)] = 3;
    dev.new ();
    plot(personal.network.u, vertex.label= NA, vertex.color=mfc,vertex.size=sizeVet, edge.width = widEd,mark.groups = by(seq_along(mfc), mfc, invisible) , main="Max dispersion/embeddedness");
  }
  else
  {
    print (paste(c("No high dispersion node", toString(coreNode)), collapse=" "));
  }
}

# Main functionality
embeddedness.values = list()
dispersion.values = list()
mean.dispersion.values = list()
max.dispersion.values = list()

all.embeddedness.values = c()
all.dispersion.values = c()

counter = 1
 
for (core in core.nodes) {
  core.embeddedness.vector = numeric()
  core.dispersion.vector = numeric()
  this.core.subgraph = create.personal.network(facebook.graph, core)
  # u = which(this.core.subgraph$name == core)
  u = which(V(this.core.subgraph)$id == core)
  this.subgraph.nodes = V(this.core.subgraph)
  
  for (v in this.subgraph.nodes) {
    if (v != u) {
      core.embeddedness = compute.embeddedness(this.core.subgraph, u, v)
      core.dispersion = compute.dispersion(this.core.subgraph, u, v)
      core.embeddedness.vector = append(core.embeddedness.vector, core.embeddedness)
      core.dispersion.vector = append(core.dispersion.vector, core.dispersion)
      all.embeddedness.values = append(all.embeddedness.values, core.embeddedness)
      all.dispersion.values = append(all.dispersion.values, core.dispersion)
      cat("done with node",v,"in core number",counter,"\n")
    }
  }
  embeddedness.values[[counter]] = core.embeddedness.vector
  dispersion.values[[counter]] = core.dispersion.vector
  mean.dispersion.values[[counter]] = mean(core.dispersion.vector)
  max.dispersion.values[[counter]] = max(core.dispersion.vector)
  counter = counter + 1
}

hist (all.embeddedness.values, 
      breaks = seq (-0.5, 
                    by = 1, 
                    length.out = max(all.embeddedness.values) + 2), 
      main ="Embeddedness Distribution")

# hist (all.dispersion.values, 
#       breaks=seq (-0.5, 
#                   by = 1, 
#                   length.out = max(all.dispersion.values) + 2), 
#       main="Dispersion Distribution")

plot(all.dispersion.values,
     xlab = "Node",
     ylab = "Frequency",
     main = "Dispersion Distribution")

dispersion.embeddedness.ratio(facebook.graph, core.nodes[1])
dispersion.embeddedness.ratio(facebook.graph, core.nodes[8])
dispersion.embeddedness.ratio(facebook.graph, core.nodes[12])

# max.embeddedness.values = apply(embeddedness.values, max)


#############################################################################################
#                                         QUESTION 6                                        #
############################################################################################# 
# Initialize feature we will use for comparison
max.avg.degree <- numeric()
max.index.avg.degree <- numeric()
min.avg.degree <- numeric()
min.index.avg.degree <- numeric()

max.cluster.coeff <- numeric()
max.index.cluster.coeff <- numeric()
min.cluster.coeff <- numeric()
min.index.cluster.coeff <- numeric()

max.density <- numeric()
min.density <- numeric()
max.index.density <- numeric()
min.index.density <- numeric()

# For each corenode, run the algorithm
for(i in 1:length(core.nodes))
{
  #Create a Personal network of the core-node
  pn = induced.subgraph(facebook.graph,c(core.nodes[i],neighbors(facebook.graph,core.nodes[i])))
  nodecount = vcount(pn)
  edgecount = ecount(pn)
  cat("\nIndex:",i)
  cat("\nNode count:",nodecount)
  cat("\nEdge count:",edgecount)
  walktrap.com = walktrap.community(pn)
  
  #Find teh sub-communities with >10 members
  subcommunities.greaterthan10 <- numeric()
  for(j in 1:length(walktrap.com))
  {
    cnodes <- length(V(pn)[which(walktrap.com$membership ==j)])
    if(cnodes>10)
    {
      subcommunities.greaterthan10 <-c(subcommunities.greaterthan10,j)
    }
  }
  
  #With the sub-communities of teh PN found, for each subcommunity
  #Find degree,transitivity(golbal), and density
  avg.degree <- numeric()
  cluster.coeff <- numeric()
  density <- numeric()
  for(k in 1:length(subcommunities.greaterthan10))
  {
    #Create a smaller graph for each sub-community
    pn.sub <- induced.subgraph(pn,V(pn)[which(walktrap.com$membership==subcommunities.greaterthan10[k])])
    
    #Get the features
    avg.degree <- c(avg.degree,mean(degree(pn.sub)))
    cluster.coeff <- c(cluster.coeff,transitivity(pn.sub,type="global"))
    density <- c(density,graph.density(pn.sub))
    m <- vector()
    con <- vector()
    # for (s in 0: nrow(walktrap.com$merges)) 
    # {
    #   
    #   memb <- cutat(walktrap.com, steps=s)
    #   m <- c(m, modularity (pn, memb, weights=NULL))
    #   
    #   make_clusters(facebook.graph, memb)
    #   
    #   intra<-0
    #   extra<-0 
    #   for(i in 1:length(E(facebook.graph))) 
    #   {
    #     ifelse(crossing(g2, facebook.graph)[i]==FALSE, intra<-intra+1, extra<-extra+1)
    #   }
    #   con <-c(con, extra/(2*intra+extra)) 
    #   
    # } 
  }
  
  #Calculate density stats
  max.density <- c(max.density,max(density))
  min.density <- c(min.density,min(density))
  max.index.density <- c(max.index.density,subcommunities.greaterthan10[which.max(density)])
  min.index.density <- c(min.index.density,subcommunities.greaterthan10[which.min(density)])
  
  #Calculate average degree stats
  max.avg.degree <- c(max.avg.degree,max(avg.degree))
  min.avg.degree <- c(min.avg.degree,min(avg.degree))
  max.index.avg.degree <- c(max.index.avg.degree,subcommunities.greaterthan10[which.max(avg.degree)])
  min.index.avg.degree <- c(min.index.avg.degree,subcommunities.greaterthan10[which.min(avg.degree)])
  
  #Calculate global transitivity stats
  max.cluster.coeff <-c(max.cluster.coeff,max(cluster.coeff))
  min.cluster.coeff <-c(min.cluster.coeff,min(cluster.coeff))
  max.index.cluster.coeff <-c(max.index.cluster.coeff,subcommunities.greaterthan10[which.max(cluster.coeff)])
  min.index.cluster.coeff <-c(min.index.cluster.coeff,subcommunities.greaterthan10[which.min(cluster.coeff)])
}

#Create data frames
maximum.density <- data.matrix(max.density)
minimum.density <- data.matrix(min.density)
maximum.index.density <- data.matrix(max.index.density)
minimum.index.density <- data.matrix(min.index.density)
colnames(maximum.density) <- "Maximum Density"
colnames(minimum.density) <- "Minimum Density"
colnames(maximum.index.density) <- "Maximum Index Density"
colnames(minimum.index.density) <- "Minimum Index Density"

maximum.avg.degree <- data.matrix(max.avg.degree)
minimum.avg.degree <- data.matrix(min.avg.degree)
maximum.index.avg.degree <- data.matrix(max.index.avg.degree)
minimum.index.avg.degree <- data.matrix(min.index.avg.degree)
colnames(maximum.avg.degree) <- "Maximum Average Degree"
colnames(minimum.avg.degree) <- "Minimum Average Degree"
colnames(maximum.index.avg.degree) <- "Maximum Index Average Degree"
colnames(minimum.index.avg.degree) <- "Minimum Index Average Degree"


maximum.cluster.coeff <- data.matrix(max.cluster.coeff)
minimum.cluster.coeff <- data.matrix(min.cluster.coeff)
maximum.index.cluster.coeff <- data.matrix(max.index.cluster.coeff)
minimum.index.cluster.coeff <- data.matrix(min.index.cluster.coeff)
colnames(maximum.cluster.coeff) <- "Maximum Cluster Co-efficient"
colnames(minimum.cluster.coeff) <- "Minimum Cluster Co-efficient"
colnames(maximum.index.cluster.coeff) <- "Maximum Index Cluster Co-efficient"
colnames(minimum.index.cluster.coeff) <- "Minimum Cluster Co-efficient"


x1 <- data.frame(cbind(minimum.avg.degree,maximum.avg.degree,minimum.index.avg.degree,maximum.index.avg.degree))
write.table(x1,file= "D:/UCLA/Spring 2016/EE232E/Project1/soln6-1.txt",append=FALSE,sep=";")

X2 <- data.frame(cbind(minimum.density,maximum.density,minimum.index.density,maximum.index.density))
write.table(X2,file= "D:/UCLA/Spring 2016/EE232E/Project1/soln6-2.txt",append=FALSE,sep=";")

X3 <- data.frame(cbind(minimum.cluster.coeff,maximum.cluster.coeff,minimum.index.cluster.coeff,maximum.index.cluster.coeff))
write.table(X3,file= "D:/UCLA/Spring 2016/EE232E/Project1/soln6-3.txt",append=FALSE,sep=";")
#Plot graphs

plotdegree = par(mfrow=c(2,2))
hist(max.avg.degree,col = "aquamarine4")
hist(max.index.avg.degree,col = "aquamarine1")
hist(min.avg.degree,col = "aquamarine2")
hist(min.index.avg.degree,col = "aquamarine3")
par(plotdegree)

plotcoeff = par(mfrow=c(2,2))
hist(max.cluster.coeff,col = "darkolivegreen1")
hist(max.index.cluster.coeff,col = "darkolivegreen2")
hist(min.cluster.coeff,col = "darkolivegreen3")
hist(min.index.cluster.coeff,col = "darkolivegreen4")
par(plotcoeff)

plotcoeff = par(mfrow=c(2,2))
hist(max.density,col = "brown1")
hist(min.density,col = "brown2")
hist(max.index.density,col = "brown3")
hist(min.index.density,col = "brown4")
par(plotcoeff)


#############################################################################################
#                                         QUESTION 7                                        #
#############################################################################################
# Path to gplus directory on the computer
gplus_input = "/Users/georgef/Documents/grad_school/spring_2016/EE_232E/project1/gplus/"

# Get unique Ids from the name of teh file

#Method1
#Take all files, find only unique first part, and 
#Use group-capturing to just tke the first part of teh file names
# allfiles <- list.files(gplus_input)
# allIds <-sub("^([^.]*).*", "\\1", allfiles)

#MEthod2
# USe one type of files as an indicator
#Use group-capturing to just tke the first part of teh file names
circleFiles <- list.files(gplus_input, pattern = "*.circles")
circleIds <-sub("^([^.]*).*", "\\1", circleFiles)

egoNodes <- unique(circleIds) 


#Filter Ids to find the list of ppl who have more than two circles
allCircles <- numeric()
for(a in 1:length(circleFiles))
{
  path <- paste(gplus_input,circleFiles[a],sep = "")
  circleFileCon <- file(path, open="r")
  lines <- readLines(circleFileCon)
  allCircles <- append(allCircles,length(lines))
  close(circleFileCon)
}

#Store the egoIds of teh nodes with > 2 circles in the circles vector
circles <- numeric()
circlelen <- numeric()
for(i in 1:length(allCircles))
{
  if(allCircles[i]>2)
  {
    circles <- append(circles,egoNodes[i])
    circlelen <- append(circlelen,allCircles[i])
  }
}

cat("Number of nodes with more than 2 circles is:",length(circles))

#Ego node is not available in the edgelist, so we add the ego node manually
#Running for just one egoId at first

for(i in 1:length(circles))
{
  egoNodeId <- circles[i]
  cat("\nNodeId:",egoNodeId)
  #We need the following two files
  edgelistFile = paste(gplus_input , egoNodeId  , ".edges", sep="")
  circlesFile = paste(gplus_input , egoNodeId , ".circles" , sep="")
  #allcircles contains teh circle count
  g1 = read.graph(edgelistFile, format = "ncol", directed = TRUE)
  #Add ego node to the graph
  g2 = add.vertices(g1, 1, name=egoNodeId)
  
  #Add edge between egoNode and all other nodes
  addEgoEdges = c()
  for (nodeIndex in 1:(vcount(g2)-1)) 
  {
    addEgoEdges = c(addEgoEdges , c(vcount(g2),nodeIndex))
  }
  g2 = add.edges(g2,addEgoEdges)
  
  circleNodes <- list()
  for (x in 1:circlelen[i]) {
    linesplit = strsplit(lines[x],"\t")
    circleNodes[[x]] = linesplit[[1]][-1]
  }
  
  #Get Walktrap community
  walktrap.com = walktrap.community(g2)
  cat("Walktrap\n",egoNodeId," ","Number of communties",":",max(walktrap.com$membership)," ","Number of circles",":",circlelen[i])
  #Write out stats:
  for(j in 1:max(walktrap.com$membership))
  { 
    #Number of communities in a node
    members = c()
    for(k in 1:length( walktrap.com$membership))
    { 
      if(walktrap.com$membership[k]==j)
      {
        members = c(members,(walktrap.com$name[k]))
      }
    }
    
    #Find out the number of members that overlap across the two circles
    member_overlap = c()
    circle_overlap = c()
    for(l in 1:circlelen[i])
    { 
      # no of circles of a node given
      interesction = intersect(members,circleNodes[[l]])
      #Percentage members shared 
      x1 = length(interesction)/length(members)
      #Is it shared with all circles?
      x2 = length(interesction)/length(circleNodes[[l]])
      member_overlap = c(member_overlap, x1)
      circle_overlap= c(circle_overlap, x2)
    }
    
    cat('\nCommunity:',j)
    cat("\nCommunities Overlap:", member_overlap,"\nCircle Overlap:", circle_overlap)
  }
  # plot(walktrap.com,
  #      g2,
  #      main = "Walktrap Community: ",
  #      vertex.size = 5,
  #      edge.arrow.size= 0.2,
  #      vertex.label = NA)
  
  #Get Infomap community
  infomap.com = infomap.community(g2)
  cat("InfoMap\n",egoNodeId," ","Number of communties",":",max(infomap.com$membership)," ","Number of circles",":",circlelen[i])
  #Write out stats:
  for(j in 1:max(infomap.com$membership))
  { 
    #Number of communities in a node
    members = c()
    for(k in 1:length( infomap.com$membership))
    { 
      if(infomap.com$membership[k]==j)
      {
        members = c(members,(infomap.com$name[k]))
      }
    }
    
    #Find out the number of members that overlap across the two circles
    member_overlap = c()
    circle_overlap = c()
    for(l in 1:circlelen[i])
    { 
      # no of circles of a node given
      interesction = intersect(members,circleNodes[[l]])
      #Percentage members shared 
      x1 = length(interesction)/length(members)
      #Is it shared with all circles?
      x2 = length(interesction)/length(circleNodes[[l]])
      member_overlap = c(member_overlap, x1)
      circle_overlap= c(circle_overlap, x2)
    }
    
    cat('\nCommunity:',j)
    cat("\nCommunities Overlap:", member_overlap,"\nCircle Overlap:", circle_overlap)
  }
  plot(infomap.com,
       g2,
       main = "Infomap Community: ",
       vertex.size = 5,
       edge.arrow.size= 0.2,
       vertex.label = NA)
  
}


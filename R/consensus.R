######################################################################
#
# consensus.R
#
# copyright (c) 2001-2002, Hao Wu and Gary A. Churchill, The Jackson Lab.
# written May, 2002
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/maanova package
#
#
######################################################################

consensus <- function(macluster, level=0.8, draw=TRUE)
{

     
  if(class(macluster) != "macluster")
    stop("The input variable is not an object of class macluster.")

  if( (level<0.5) | (level>1) )
    stop("Level must be between 0.5 and 1.")

 # if(macluster$n.perm == 1) {
 #   msg <- "You don't have permutation result in the input object."
 #   msg <- paste(msg, "Set n.perm to be an integer bigger than 1")
 #   msg <- paste(msg, "in macluster function in order to do permutation.")
 #   stop(msg)
 # }

  # local variable
  n.perm <- macluster$n.perm

  # result
  consensus <- NULL
  
  if( macluster$method == "hc")
    #stop("Consensus for HC has not implemented yet.")
    consensus <- consensus.hc(macluster, level, draw)
  else if( macluster$method == "kmean") {
    consensus <- consensus.kmean(macluster, level, draw)
  }


  invisible(consensus)
}




consensus.hc <- function(macluster, level, draw)
{
  # convert the observed hierarchical cluster and store it
  # each clade is represented by three integers:
  # number of nodes, sum of node indices and product of node indices
  obstree <- cluster2num(macluster$cluster.obs$merge)

  # how many nodes in each tree?
  n.node <- length(obstree$tree)
  # make the result to be the observed tree
  result <- obstree
  # total number of nodes in consensus tree
  n.totalnode <- length(result$tree)
  result$count <- rep(1, n.totalnode)
  # loop for the permutation trees
  if(macluster$n.perm > 1) {
    for(i in 2:macluster$n.perm) {
      permtree <- cluster2num(macluster$cluster.perm[[i]])
      # compare the permutation tree with the nodes already in result
      for(j in 1:n.node) {
        # number of leaves in this tree
        n.leaves.tmp <- permtree$num[j,1]
        # sum of leaves numbers
        sum.leaves.tmp <- permtree$num[j,2]
        # find the matches for number of nodes and sum of indices
        # for this perm tree in the result
        idx.match <- which( (result$num[,1]==n.leaves.tmp)&
                           (result$num[,2]==sum.leaves.tmp) )

        flag <- FALSE

        # if there' some match, check if they are really the result we want
        if(length(idx.match) != 0) {
          for(k in 1:length(idx.match)) {
            if(length(setdiff(permtree$tree[[j]], result$tree[[idx.match[k]]])) == 0) {
              # no difference in this two sets. BINGO!
              result$count[idx.match[k]] <- result$count[idx.match[k]] + 1
              flag <- TRUE
              break
            }
          }
        }
        
        # if there's no real match, add this node to the result
        if( (length(idx.match) == 0) | (flag==FALSE) ) {
          n.totalnode <- n.totalnode + 1
          result$tree[[n.totalnode]] <- permtree$tree[[j]]
          result$num <- rbind(result$num, permtree$num[j,])
          result$count <- c(result$count, 1)
        }
      } # finish loop for nodes in this permutation tree
    
    } # finish loop for all trees
  }
  
  # total number of leaves in the tree
  result$nleaves <- n.node+2
  result$count <- result$count/macluster$n.perm

  # the leave names
  result$leaves.name <- macluster$leave.names
#  if(missing(leaves.name)) {
#    if(macluster$what == "gene")
#      idx.tmp <- macluster$idx.gene
#    else
#      idx.tmp <- 1:result$nleaves
#    leaves.name <- paste(macluster$what, idx.tmp, sep="")
#  }
#  result$leaves.name <- leaves.name
  # visualize the result consensus tree

  idx.node <- which(result$count>level)
  # number of inner nodes
  n.nodes <- length(idx.node)

  # make the consensus tree object as list
  ct <- vector("list", n.nodes+1)
  # make the root
  ct[[1]]$depth <- 0
  ct[[1]]$children <- NULL
  ct[[1]]$nleaves <- result$nleaves
  ct[[1]]$occurance <- 1

  if(n.nodes==0) {
    msg <- "There's no inner nodes in the tree."
    msg <- paste(msg, "All leaves are direct children of the root.")
    msg <- paste(msg, "You may need to lower the significant level.")
    stop(msg)
  }
  else { # there are some inner nodes
    # update occurance
    for(i in 1:n.nodes)
      ct[[i+1]]$occurance <- result$count[idx.node[i]]
        
    # take out the nodes with occurance bigger than level
    tree <- result$tree[idx.node]
    
    # convert the selected nodes to a matrix
    # the number of rows is the number of nodes exclude the root
    # the number of columns is the number of total leaves
    binstr <- matrix(rep(0, n.nodes*result$nleaves), n.nodes,
                     result$nleaves)
    # turn on certain positions
    for(i in 1:n.nodes)
      binstr[i, tree[[i]]] <- 1

    # call build tree recursively
    ct <- buildtree(ct, binstr, 0, 1, 1:n.nodes, 1:result$nleaves)
  }

  result$ct <- ct
  result$level <- level
  result$tree <- NULL
  result$count <- NULL
  result$num <- NULL
  class(result) <- "consensus.hc"

  # plot the result
  if(draw)
    plot(result, title=macluster$term)
  
  invisible(result)

}


plot.consensus.hc <- function(x, title, ...)
{
  con <- x
  # visualize the tree
  ###########################
  # calculate the coordinates
  ###########################
  # start from the root
  coord <- NULL
  coord$A <- NULL; coord$B <- NULL;
  coord$xycoord <- NULL; coord$nodecoord <- vector("list", length(con$ct))
  coord$leavesidx <- NULL
  maxdepth <- 0
    
  # find the maximum depth in the tree
  for(i in 1:length(con$ct))
    maxdepth <- max(maxdepth, con$ct[[i]]$depth)
  
  coord <- makeAB(con$ct, coord, 1, 0, maxdepth)
    
  # draw the figure
  if(missing(title)) title <- ""
  title <- paste(title, " - ", con$level*100, "% Consensus tree", sep="")
  plot(coord$A, coord$B, type="n", main=title,
       xlab="", ylab="depth", xaxt="n")
  for(i in 1:dim(coord$A)[2]) 
    lines(coord$A[,i], coord$B[,i])
  # add texts for leaves
  for(i in 1:length(coord$leavesidx))
    text(i-1, 0, con$leaves.name[coord$leavesidx[i]],srt=90, pos=4,
         offset=0.7)
  # add texts for significant levels
  for(i in 1:(length(con$ct)-1)) {
    x <- mean(range(coord$nodecoord[[i+1]][,1]))
    y <- coord$nodecoord[[i+1]][1,2]
    text(x, y, formatC(con$ct[[i+1]]$occurance, format="f", digits=2),
         pos=3)
  }
}



makeAB <- function(ct, coord, treeidx, startx, maxdepth)
{
  # local variable
  A <- coord$A
  B <- coord$B
  xycoord <- coord$xycoord
  nodecoord <- coord$nodecoord
  leavesidx <- coord$leavesidx
  
  subct <- ct[[treeidx]]
  # current x coordinate
  currentx <- startx
  
  # find the direct leaves in children
  directleaves <- subct$children[subct$children<0]
  # find the inner nodes in children 
  innernode <- subct$children[subct$children>0]

  # create the verticle lines for each leaf
  if(length(directleaves) != 0){
    leavesidx <- c(leavesidx, -directleaves)
    for(i in 1:length(directleaves)) {
      x0 <- currentx
      x1 <- currentx
      y0 <- 0
      y1 <- maxdepth - subct$depth + 1
      # append to A and B
      A <- cbind(A, c(x0,x1))
      B <- cbind(B, c(y0,y1))
      # increase startx
      currentx <- currentx + 1
      # update nodecoord for this node
      nodecoord[[treeidx]] <- rbind(nodecoord[[treeidx]], c(x1, y1))
    }
  }
  # if there's no inner nodes, create the horizontal and
  # vertical lines and return
  if(length(innernode) == 0) {
    x0 <- startx
    y0 <- maxdepth - subct$depth + 1
    x1 <- currentx - 1
    y1 <- y0
    A <- cbind(A, c(x0,x1))
    B <- cbind(B, c(y0,y1))
    # create a verticle line to connect the children and parent
    A <- cbind(A, c((x0+x1)/2,(x0+x1)/2))
    B <- cbind(B, c(y1, y1+1))
    # return the midpoint of the horizontal line
    xycoord <-  c((x0+x1)/2, y1+1)
  }
  else {
    # if there are any inner nodes,
    # call makeAB recursively for all inner nodes
    for(i in 1:length(innernode)) {
      if( is.null(A) )
        substartx <- 0
      else
        substartx <- max(A) + 1
      # update input coord object
      coord$A <- A;  coord$B <- B;
      coord$xycoord <- xycoord;
      coord$nodecoord <- nodecoord
      coord$leavesidx <- leavesidx
      # recursively call makeAB
      coord.tmp <- NULL
      coord.tmp <- makeAB(ct, coord, innernode[i], substartx, maxdepth)
      A <- cbind(coord.tmp$A, c(coord.tmp$xycoord[1], coord.tmp$xycoord[1]))
      B <- cbind(coord.tmp$B, c(coord.tmp$xycoord[2]-1, coord.tmp$xycoord[2]))
      nodecoord <- coord.tmp$nodecoord
      nodecoord[[treeidx]] <- rbind(nodecoord[[treeidx]], coord.tmp$xycoord)
      leavesidx <- coord.tmp$leavesidx
    }
    # create a horizontal line to connect all inner nodes and direct leaves 
    x0 <- min(nodecoord[[treeidx]][,1])
    x1 <- max(nodecoord[[treeidx]][,1])
    y0 <- nodecoord[[treeidx]][1,2]
    y1 <- y0
    A <- cbind(A, c(x0,x1))
    B <- cbind(B, c(y0,y1))
    # create a verticle line to connect the children and parent
    # if this node is not root
    if(treeidx != 1) {
      A <- cbind(A, c((x0+x1)/2, (x0+x1)/2))
      B <- cbind(B, c(y1,y1+1))
    }

    # return the midpoint of the horizontal line
    xycoord <- c((x0+x1)/2, y1+1)
  }

  # make the return variable
  result <- NULL
  result$A <- A
  result$B <- B
  result$xycoord <- xycoord
  result$nodecoord <- nodecoord
  result$leavesidx <- leavesidx
  result
}

 


    

buildtree <- function(ct, binstr, depth, parent, idx.node, idx.leave)
{
  # local variable
  n <- dim(binstr)[1]
  l <- dim(binstr)[2]

  # find the missing leaves in binstr and make them the direct
  # children of parent
  missleave <- NULL
  for(i in 1:l) {
    if(all(binstr[,i]==0)) # direct children
      missleave <- c(missleave, i)
  }
  if(!is.null(missleave)) {
    ct[[parent]]$children <- c(ct[[parent]]$children, -idx.leave[missleave])
    binstr <- binstr[,-missleave]
    idx.leave <- idx.leave[-missleave]
  }

  # if there are only one node, we are done
  if(n==1) {
    ct[[idx.node+1]]$children <- -idx.leave
    ct[[idx.node+1]]$depth <- depth + 1
    ct[[idx.node+1]]$nleaves <- length(idx.leave)
    ct[[parent]]$children <- c(ct[[parent]]$children, idx.node+1)
    return(ct)
  }

  # find the node with maximum elements
  nele <- NULL
  for(i in 1:n)
    nele[i] <- sum(binstr[i,]==1)
  nleave <- max(nele)
  maxidx <- which(nele==nleave)[1]

  # update the node object for maxidx
  ct[[idx.node[maxidx]+1]]$depth <- depth + 1
  ct[[idx.node[maxidx]+1]]$nleaves <- nleave
  ct[[parent]]$children <- c(ct[[parent]]$children, idx.node[maxidx]+1)
  
  # use binstr(maxidx) to divide the inner nodes into two groups, e.g.,
  # set1 is the nodes to contain the leaves in binstr(maxidx) but doesn't
  # include maxidx itselt. set2 is the nodes don't contain that.
  leaves <- which(binstr[maxidx,]==1)
  noleaves <- which(binstr[maxidx,]==0)
  set1 <- NULL
  set2 <- NULL
  for(i in 1:dim(binstr)[1]) {
    if(all(binstr[i,leaves]==0))
      set2 <- c(set2,i)
    else
      set1 <- c(set1, i)
  }
  if(length(set1)==1)
    # only one node in set1
    ct[[idx.node[maxidx]+1]]$children <- -idx.leave[leaves]
  
  # exclude maxidx from set1
  set1 <- setdiff(set1, maxidx)

  # create two matrices
  binstr1 <- binstr[set1,leaves]
  binstr2 <- binstr[set2, noleaves]

  # call buildtree recursively
  if(length(binstr1)!=0)
    ct <- buildtree(ct, matrix(binstr1, nrow=length(set1)), ct[[idx.node[maxidx]+1]]$depth,
              idx.node[maxidx]+1, idx.node[set1], idx.leave[leaves])
  if(length(binstr2) != 0)
    ct <- buildtree(ct, matrix(binstr2, nrow=length(set2)), depth, parent, idx.node[set2],
              idx.leave[noleaves])

  ct

}

cluster2num <- function(clust)
{
  # number of leaves, exclude the last one(root)
  n.leaves <- dim(clust)[1] - 1
  result <- NULL
  result$tree <- list(NULL)
  result$num <- matrix(rep(0, n.leaves*2), n.leaves, 2)
  result$tree[[1]] <- -clust[1,]
  result$num[1,] <- c(2, sum(-clust[1,]))
  for(i in 2:n.leaves) {
    tmp1 <- -clust[i,1]
    tmp2 <- -clust[i,2]
    if(tmp1 < 0)
      tmp1 <- result$tree[[-tmp1]]
    if(tmp2 < 0)
      tmp2 <- result$tree[[-tmp2]]
    result$tree[[i]] <- c(tmp1, tmp2)
    result$num[i,] <- c(length(result$tree[[i]]), sum(result$tree[[i]]))
  }

  result
}


#############################################################
# The following functions are using for Kmeans clustering
#############################################################

consensus.kmean <- function(macluster, level, draw)
{
  # local variable
  what <- macluster$what
  n.perm <- macluster$n.perm
  term <- macluster$term
  n.leaves <- length(macluster$cluster.obs)
  n.grp <- macluster$kmean.ngroups
  idx.gene <- macluster$idx.gene
  groups <- matrix(rep(0, n.perm*n.leaves), n.perm, n.leaves)
  
  # fill groups
  groups[1,] <- macluster$cluster.obs
  if(n.perm > 1)
    for(i in 2:n.perm)
      groups[i,] <- macluster$cluster.perm[[i]]

  # initialize output
  result <- list(NULL)
  for(i in 1:(n.grp+1))
    result[[i]] <- numeric(0)
  
  tmp <- rep(0, n.grp)
  # given level, calculate consensus result
  for(i in 1:n.leaves) {
    for(j in 1:n.grp) 
      tmp[j] <- sum(groups[,i]==j)/n.perm

    # find the consensus group for this leaf
    # if level is not 1, use ">" (that's because
    # if level is 0.5, the leaf can be in two groups, each
    # has 50% occurance)
    # if level is 1, use ">="
    if(level != 1)
      idx <- which(tmp>level)
    else 
      idx <- which(tmp>=level)
    
    if(length(idx) != 0) # there's a group for this leaf
      result[[idx]] <- c(result[[idx]], i)
    else
      result[[n.grp+1]] <- c(result[[n.grp+1]], i)
  }
  
  # return variable
  output <- NULL
  output$group <- result
  output$what <- what
  output$data.draw <- macluster$VG
  output$condition.names <- macluster$condition.names
  class(output) <- "consensus.kmean"

  # draw the VG profiles
  if(draw)
    plot(output)
  
  invisible(output)
  
}

plot.consensus.kmean <- function(x, ...)
{
  con <- x
  # local variables
  data.draw <- con$data.draw
  ncondition <- dim(data.draw)[2]
  what <- con$what
  group <- con$group
  n.grp <- length(group) - 1

  # number of figures per plot
  if(n.grp > 8)
    nfig.plot <- 9
  else
    nfig.plot <- n.grp + 1
  
  # determin how many plots will be generated
  # maximum figures on one plot is 9
  # don't forget the not-in-any-group group
  # calculate the number of plots and figures per plot for a pretty layout
  
  if(n.grp > 8)
    n.plot <- ceiling((n.grp+1)/9)
  else
    n.plot <- 1
  
  # determine how many figures horizontally
  if(n.grp<=6)
    n.col <- 2
  else
    n.col <- 3
  # number of cols in figure is n.col, calculate number of rows
  n.row <- ceiling(nfig.plot/n.col)
  
  # save old par
  #old.mar <- par("mar")
  #old.las <- par("las")
  #old.mfrow <- par("mfrow")
  #old.mfcol <- par("mfcol")
  #on.exit(par(las = old.las, mar = old.mar, mfrow = old.mfrow,
  #            mfcol=old.mfcol))

  for(iplot in 1:n.plot) {
    # draw gene profiles
    ylabel <- "Expression level"
    if(what=="gene")
      xlabel <- "Sample"
    else 
      xlabel <- "Gene"

    # starting and ending group in this plot
    istart <- nfig.plot*(iplot-1) + 1
    iend <- nfig.plot * iplot
    if(iend > n.grp)
      iend <-  n.grp + 1
    # number of figures on this plot
    nfig.thisplot <- iend - istart + 1
    
    # generate a new figure for each plot
    get(getOption("device"))()
#    if(.Platform$GUI == "AQUA")
#      quartz()
#    else
#      x11()

    # build a matrix for figure layout
    mat <- matrix(c(1:nfig.thisplot, rep(0, n.col*n.row-nfig.thisplot)),
                  n.row, n.col, byrow=TRUE)
    layout(mat)

    # y-axis limits. It'll be the same for all plots
    ylim <- range(data.draw)
    
    for(i in istart:iend) {
      idx.grp <- group[[i]]
      if(i==(n.grp+1))
        title <- paste("Not in any group \n\n(",length(idx.grp)," ",
                       what, "s)", sep="")
      else
        title <- paste("Expression profile for group ", i,
                       "\n\n (",length(idx.grp)," ", what,"s)", sep="")
      
      # data to be drawn for this group
      data.grp <- data.draw[idx.grp,]
      # create an empty plot
      plot(1:ncondition, ylim=ylim, type="n", xaxt="n",
           xlab=xlabel, ylab=ylabel, main=title)
      
      if(length(idx.grp) == 1)
        # only one gene in this group
        # draw a line
        lines(data.grp, col="blue", type="b")
      else if(length(idx.grp) > 0)
         # draw the lines
        apply(data.grp, 1, function(x) lines(x, col="blue", type="b"))
      
      # draw the axis
      axis(1, at=1:ncondition, labels=con$condition.names)
    }
  }
}

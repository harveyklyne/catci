make_leaf <- function(leaf){
  return(list("root" = leaf,
              "children" = list()))
}
make_parent <- function(children){
  return(list("root" = sort(unlist(lapply(children, function(x){x$root}))),
              "children" = children))
}
make_binary_above <- function(level){
  if(length(level)<2){stop("Level does not contain enough nodes to construct any parents.")}
  out <- lapply(seq(floor(length(level)/2)),
                function(j){make_parent(level[(2*(j-1)+1):(2*j)])})
  if(length(level) %% 2 == 1){
    out <- c(out, level[length(level)])
  }
  return(out)
}

make_binary_tree <- function(d){
  level <- lapply(seq(d), make_leaf)
  while(length(level) > 1){
    level <- make_binary_above(level)
  }
  return(level[[1]])
}

get_siblings <- function(node, tree){
  siblings <- list()
  if (length(setdiff(tree$root, node)) == 0){return(siblings)}
  for (child in tree$children){
    if ((all(node %in% child$root)) & (length(setdiff(child$root, node)) > 0)){
      return(get_siblings(node, child))
    }
    if (any(child$root %in% node)){
      if (!all(child$root %in% node)){stop("Only part of vertex in node.")}
    } else{
      siblings <- c(siblings, list(child))
    }
  }
  return(siblings)
}


get_tree_ind2 <- function(ind1, tree, category){
  if (ind1 >= length(category)) return(integer(0))
  cat1 <- category[[ind1]]
  siblings <- get_siblings(cat1, tree)
  test_inds <- seq(ind1 + 1, length(category))
  return(test_inds[which(sapply(test_inds, function(j){
    cat2 <- category[[j]]
    any(sapply(siblings, function(x){setequal(cat2, x$root)}))
  }))])
}

get_tree_levels <- function(tree, category){
  return(sum(sapply(seq_along(category), function(x){
    length(get_tree_ind2(x, tree, category))
  })))
}

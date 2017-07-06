setwd("~/workspace/rmfpp")
source("subtimeline_maker.r")
source("likelihood_functions.r")
source("partitioned_pcll.r")

library(tidyr)

## Forest modification proposals

## sketch
# proposals occur at:
#   1. leaf nodes
#   2. new trees (stump)
#   3. switch from leaf to distribution

# proposal types:
# a. splits ([event], time)

# replace tree[active,] with
# [1] split
# [2] trueNode
# [3] falseNode
replace_leaf = function(tree, active, internalTrueFalseTree, markReplace=F) {
  if(!exists("replaceLeafCols")) replaceLeafCols = c("condition", "timeFrame", "trueNode","falseNode","distribution",
                                                     "par0","par1","par2","modifiable")
  if(internalTrueFalseTree[[1,"trueNode"]] != 2 || internalTrueFalseTree[[1,"falseNode"]] != 3)
    stop("internalTrueFalseTree must point trueNode to second row and falseNode to third row")
  internalTrueFalseTree[1,c("trueNode","falseNode")] = (dim(tree)[1] + 1:2)
  tree[active,] = internalTrueFalseTree[1,replaceLeafCols] 
  tree = tree %>% bind_rows(internalTrueFalseTree[-1,replaceLeafCols])
  
  if(markReplace) {
    tree = tree %>%
      mutate(replaced=rep(F,nrow(.)))
    tree[c(active,nrow(tree)-1, nrow(tree)),"replaced"] = T
    tree
  } else {
    tree
  }
}

modify_leaf = function(tree, active, leaf) {
  if(!exists("replaceLeafCols")) replaceLeafCols = c("condition", "timeFrame", "trueNode","falseNode","distribution",
                                                     "par0","par1","par2","modifiable")
  tree[active,] = leaf[replaceLeafCols]
  tree
}

# Adds a blank tree to the end of the forest;
#   modifiable: 0 --> fixed, 1 --> modifiable, NA --> raise node (modified through raise (par[3]) only)
add_blank_tree = function(forest, modifiable = 1) {
  forest$tree = forest$tree %>% rbind(list("",NA,NA,NA,1,NA,0,NA,modifiable))
  forest$roots = c(forest$roots, nrow(forest$tree))
  forest
}

# Get a blank tree
get_blank_tree = function(modifiable = 1) {
  tree = data.frame(
    condition=character(),
    timeFrame=numeric(),
    trueNode=numeric(),
    falseNode=numeric(),
    distribution=numeric(),
    par0=character(),
    par1=numeric(),
    par2=numeric(),
    modifiable=numeric(),
    stringsAsFactors = F
  ) %>% tbl_df()
  
  tree[dim(tree)[1]+1,] = c("" ,NA,NA,NA,1 ,NA,0 ,NA, modifiable)
  finish_tree_initialization = function(tree) {
    tree$timeFrame = as.numeric(tree$timeFrame)
    tree$distribution = as.numeric(tree$distribution)
    tree$par0 = as.character(tree$par0)
    tree$par1 = as.numeric(tree$par1)
    tree$par2 = as.numeric(tree$par2)
    tree$trueNode = as.numeric(tree$trueNode)
    tree$falseNode = as.numeric(tree$falseNode)
    tree$modifiable = as.numeric(tree$modifiable)
    tree
  }
  tree = finish_tree_initialization(tree)
  tree
}

# Get a blank forest
get_blank_forest = function(modifiable = 1) {
  forest = list(tree = get_blank_tree(modifiable), roots = 1)
  forest
}

# The raise node adjusts the logmultiplier by par1, which is set by par[3] in non-exp optimization
add_raise_node = function(forest) {
  # if(sum(is.na(forest[["tree"]][,"modifiable"]) &
  #        forest[["tree"]][,"distribution"]==1, na.rm=T)==0) {
    add_blank_tree(forest, modifiable = NA)
  # } else {
  #   warning("forest already has raise node, no raise node added")
  #   forest
  # }
}

# generate proposals
# 1. random generation
# 2. identified-pattern generation
# 3. exhaustive

# generate structure
# identify parameter estimates
# select best
# repeat


# intermediate representation: compute unexpanded forest
# then conduct lesion as necessary

add_internaltruefalsetrees = function(tree, events, timeFrames) {
  n = length(events)
  if(n==0) return(tree)
  tree %>% 
    rbind(list(c(NA           ,"","") %>% rep(n) %>% replace(seq(1,3*n,3),events %>% as.character()),
               c(NA           ,NA,NA) %>% rep(n) %>% replace(seq(1,3*n,3),timeFrames),
               (c(dim(tree)[1],NA,NA) %>% rep(n) ) + 0 + 1:(3*n),
               (c(dim(tree)[1],NA,NA) %>% rep(n) ) + 1 + 1:(3*n),
               c(NA           ,1,1)   %>% rep(n), 
               c(NA           ,NA,NA) %>% rep(n),
               c(NA           ,0,0)   %>% rep(n),
               c(NA           ,NA,NA)  %>% rep(n),
               c(NA           ,1, 1) %>% rep(n)
               ))
    #rbind(list(""   ,        NA,NA,NA,1 ,0 ,NA)) %>%
    #rbind(list(""   ,        NA,NA,NA,1 ,0 ,NA))
}
add_modificationleaves = function(tree, events, distributions, modifiable=1) {
  n = length(distributions)
  chosenEvents = sample(events,size = n)
  tree %>% 
    rbind(list(rep("",n) %>% (function(.) {.[distributions==2]=as.character(chosenEvents[distributions==2]); .})(.),
               rep(NA,n),
               rep(NA,n),
               rep(NA,n),
               distributions,
               chosenEvents, # par0: for distribution == 2, this is the raise value; it is NA for distribution == 1 or NA
               rep(0,n), # par1
               rep(NA,n), # par2
               rep(modifiable, n) # modifiable
               ))
}

reordered_sample = function(ns, dots) {
  do.call("sample", list(x=dots[["t"]], size=ns))
}
noisy_reordered_sample = function(ns, dots) {
  do.call("sample", list(x=dots[["t"]], size=ns)) * 2 * runif(ns)^dots$localityFactor
}

generate_random_proposals = function(tree,
                                     events,
                                     proposals=1,
                                     pars = list(
                                       probSplit=1,
                                       type="default",
                                       localityFactor = 2,
                                       maxLength=100,
                                       d = NA, # distribution function to sample from
                                       dots = NA # parameters to pass to 'd', with nsplits preceding those arguments
                                       )
                                     ) {
  # browser()
  # tis = sample(which(tree[["condition"]]==""), size = proposals, replace=T)
  if(is.null(pars))
    pars = list(
      probSplit=0.5,
      type="default",
      localityFactor = 2,
      maxLength=100,
      d = NA, # distribution function to sample from
      dots = NA # parameters to pass to 'd', with nsplits preceding those arguments
    )
  nsplits = floor(pars$probSplit*proposals)
  nmods = floor((1-pars$probSplit)*proposals)
  if(nsplits + nmods != proposals ) {
    if(proposals*(pars$probSplit-nsplits/proposals) > runif(1))
      nsplits = nsplits + 1
    else
      nmods = nmods + 1
  }

  if(pars$type == "default") {
    proposalTrees = tree[1,] %>%
      add_internaltruefalsetrees(sample(events,nsplits, replace=T),
                                 pars$maxLength^(runif(nsplits)^pars$localityFactor)-1
                                 ) %>%
      .[-1,] %>%
      mutate(proposalId = rep(seq(1,length.out = nsplits), each=3),
             proposalType = 1)
  }
  
  proposalMods = tree[1,] %>%
    add_modificationleaves(events, rep(2,nmods))%>%#, replace = T)) %>%
    .[-1,] %>% 
    mutate(proposalId = seq(3*nsplits+1,length.out = nmods),
           proposalType = 2)
    
  proposalTrees %>% bind_rows(proposalMods)
}

generate_proposal_cuts = function(pars, n=1) {
  pars$maxLength^(runif(n)^pars$localityFactor)-1
}
# given a list of proposals, produce a map object for the reducer to use.
#   the map object will contain all the precompute statistics necessary for the reduce function
#   - for exponential proposal, this is mlss
#   - for loglogistic proposal, this is relativeTimeLLInput
#   proposal should be a list (or df) containing:
#   - $proposalType
#   - $existingTree
#   - $roots
#   - $proposalLocation (in tree)
#   - $proposedTree
evaluate_proposals_map = function(proposal, od, target) {
  od = od %>%
    group_by(id) %>%
    mutate(withProposal = list(
         get_withproposal_object(tree = proposal[["existingForest"]][["tree"]],
                                 proposedTree = proposal[["proposedForest"]][["tree"]],
                                 roots = proposal[["existingForest"]][["roots"]],
                                 proposalLocation = proposal$proposalLocation,
                                 data = pts[[1]],
                                 lbubs = lbubs[[1]])
         )) %>%
    ungroup()
                                        
  # TODO ensure the max and min logRates enter into the MLE (MAP of exp, but REDUCE of ll)
  if(proposal$proposalType == 1 ) { #exp
    #return mlss
    od %>% 
      group_by(id) %>%
      mutate(mapobject = 
           list(get_mlss(proposal[["proposedForest"]][["tree"]], proposal[["proposedForest"]][["roots"]],
                         withProposal[[1]], pts[[1]], lbubs[[1]], target))) %>%
      #TODO caution only takes first of ids in group_by (so problem if nonunique ids)
      ungroup()
  }
}

get_withproposal_object = function(tree, proposedTree, roots, proposalLocation, data, lbubs) {
  # 1. make additional window splits for proposals affecting durations
  # 2. if split only affects one tree, combine factors for other trees
  #    b. compute this for all held out trees
  # 3. per proposal:
  #    filter out unaffected durations
  #    apply proposal change
  
  
  #withProposal[1,"treeIndex"] = 5 # to test having distribution type 2 (log logistic) overwriting
  withProposal = proposedTree %>%
    subtimeline_forest_preexpansion(data,roots,lbubs, 4) %>% 
    do((function(x) if(!("replacedIndex" %in% names(x))) 
      mutate(x, replacedIndex=
               ifelse(treeIndex==proposalLocation,
                      yes = proposalLocation, 
                      no=NA)
             ) else x)(.)) %>%
    mutate(oldIndexInNewTree = (function(ri,ti) { ninsa = !is.na(ri); ti[ninsa] = ri[ninsa]; ti })(replacedIndex,treeIndex)) %>% 
    group_by(treeNumber,treeIndex,replacedIndex,oldIndexInNewTree) %>%
    do((function(x) {ret_subtimeline_tree(tree = tree,
                                         data = data,
                                         active = x$oldIndexInNewTree,
                                         lbubs = x,
                                         leafType = 1)})(.)) %>%
    ungroup() %>%
    select(-oldIndexInNewTree) %>%
    arrange(treeNumber,lb,ub)
    
    #withProposal %>% do((function(x) data.frame(la=mean(x$treeIndex),lo=T))(.)) #combines into tibble if same type/size of tbl
  withProposal
}

get_mlss = function(tree, roots, withProposal, data, lbubs, target) {
  # versus
  #tree %>% subtimeline_forest_preexpansion(data,roots,lbubs, 4) %>% arrange(treeNumber,lb, lb-ub)

  # expanded computation of sum(logRateFactor)
  # collect events: Mp_hat and qp*Tp_hat
  expansion = withProposal %>% arrange(lb) %>% (function(.) fix_timelineCPP(.$lb, .$ub))
  factors = withProposal %>% select(lb,ub,logRateFactor) %>% 
    mutate_active_duration_indices(expansion) %>%
    value_sum_onto_durations(expansion, value="logRateFactor")
  expansionCounts = findInterval((data %>% filter(event==target) %>% .[["time"]] - 1e-14), expansion$lb, left.open=T, rightmost.closed = T) %>%
    vector_of_counts(nrow(expansion))
  # note fudge factor to deal with double precision problems of R's findInterval(...); this limits the precision to 1e-14 (finer will cause uncaught errors!)
  
  ### Proposal type 1: branch exponential calculation
  liCalculation = expansion %>%
    mutate(factors=factors) %>%
    mutate(qptp_hat = exp(factors)*(ub-lb)) %>% # may need to conduct proposal expansion prior to here
    mutate(Mp_hat = expansionCounts) %>% tbl_df() # may need to conduct proposal expansion prior to here

  # convert expansion into summaries in collapsed form (i.e. collapse)
  collapsedCalculation = withProposal %>% mutate_active_duration_indices(expansion) %>%
    mutate(qptp_hat = value_sum_from_durations(.,liCalculation,"qptp_hat"),
           Mp_hat = value_sum_from_durations(.,liCalculation,"Mp_hat")
    )
  
  # filter the replacedIndex intervals, summarise, and get qp', and thus dLL.  
  mlss = collapsedCalculation %>% filter(!is.na(replacedIndex)) %>% # TODO can we move up this filter command to make faster?
    group_by(treeIndex, replacedIndex) %>%
    summarise(M = sum(Mp_hat), QT = sum(qptp_hat)) %>% ungroup() # possible double counting here if overlapping mods, watch out
  
  # maximum likelihood is just those by log qp' and (qp'-1)
  mlss = mlss %>% mutate(log_qp_prime = (function(.) { x = log(.$M)-log(.$QT); x[is.infinite(x)]=-20; x})(.), dll = M*log_qp_prime - expm1(log_qp_prime)*QT)
  mlss 
}

get_ll_mapobject = function(proposedTree, proposalLocation, withProposal, data, target) {
  ### Proposal type 2: log logistic calculation
  # get logmultipliers across active areas in withProposal expansions + data expansions
  # these go into llcll(df = active areas, pc = logmultipliers in active areas)
  # these need to be aligned with respect to lastTime -- couldn't the logmultipliers be different for
  # different active areas? yes ... so then calculation should be over logmultiplier-specific event,
  # i.e. no overlapping intervals. How to compute from this (single df, not {df, pc})?
  activeLogLogistic = withProposal %>% 
    filter(!is.na(replacedIndex)) # i.e. the durations to covers
  if(nrow(activeLogLogistic)==0)
    return(withProposal %>%
             mutate(M=0, lastTime=-Inf) %>%
             select(-replacedIndex,-treeIndex,-treeNumber) %>%
             .[0,])
    
  # expand for the model (multiple trees preexpanded and overlapping) and
  #   sum the logRateFactors *excluding* replaceIndex?
  expansion = withProposal %>% arrange(lb) %>% (function(.) fix_timelineCPP(.$lb, .$ub))
  factors = withProposal %>% 
    filter(is.na(replacedIndex)) %>% # exclude logRateFactors that will be replaced
    select(lb,ub,logRateFactor) %>% 
    mutate_active_duration_indices(expansion) %>%
    value_sum_onto_durations(expansion, value="logRateFactor")
  expansion = expansion %>% mutate(summedLogRateFactor = factors)

  # expand further to account for data events    
  
  expansionOnData = expansion %>%
    insert_events_into_expansion(data %>% filter(event==target) %>% .[["time"]])
  
  expansionLogLogistic = expansionOnData %>%
    mutate(logRateFactor=value_apply_onto_durations(
      expansion %>% mutate_active_duration_indices(expansionOnData), ., "summedLogRateFactor"))
  
  activeLogLogistic = activeLogLogistic %>% 
    mutate_active_duration_indices(expansionLogLogistic)
  
  relativeTimeLLInput = expansionLogLogistic %>%
    # DOTO you can't just throw away overlapping intervals in other trees
    # now filter down to active expanded intervals (active is in withProposal)
    mutate(replacedIndex=value_apply_onto_durations(activeLogLogistic, ., "replacedIndex")) %>%
    filter(replacedIndex > 0) %>%
    mutate(lastTime = (function(lb,eventTimes)
      c(-Inf, eventTimes, Inf)[1+findInterval(lb, eventTimes)]
      )(.$lb,filter(data, data$event==proposedTree[[proposalLocation,"condition"]]) %>% .[["time"]]) ) %>% #parameter to manipulate here # old par0->condition
    mutate(lb = ifelse(is.finite(lb-lastTime), lb-lastTime, lb), 
           ub = ifelse(is.finite(ub-lastTime), ub-lastTime, ub)) %>% 
    #      if "never", do not overwrite times with -Inf, instead let the next function deal with the non-shift
    select(-replacedIndex)

  # browser()
  relativeTimeLLInput
}

# od should contain column of mapobjects
evaluate_proposals_reduce = function(od, proposal) {
  if(proposal$proposalType == 1) {
    reduce_exp_mapobjects(od)
  }
}

reduce_exp_mapobjects = function(od) {
  od %>% select(mapobject) %>% 
    unnest() %>% 
    group_by(treeIndex, replacedIndex) %>% 
    summarise(M = sum(M), QT = sum(QT)) %>%
    ungroup() %>%
    mutate(log_qp_prime = (function(.) { x = log(.$M)-log(.$QT); x[is.infinite(x)]=-100; x})(.),
           dll = M*log_qp_prime - expm1(log_qp_prime)*QT)
}

runTests = F
if(runTests)
  source("tests/forest_modifications_tests.r")

### learn_one_random_proposal.r

proposal_location_generator = function(tree) {
  # browser()
  data.frame(index = 1:nrow(tree), modifiableValue = exp(tree$modifiable)) %>% tbl_df() %>%
    filter(!is.na(modifiableValue)) %>% mutate(modifiableValue = cumsum(modifiableValue)/sum(modifiableValue)) %>%
    mutate(belowIndex = (runif(1) < modifiableValue)*index) %>% .[["belowIndex"]] %>% (function(.) min(.[.>0]))(.)
}

# generates a prop(osal) object
generate_one_proposal_po = function(forest, dats, ...) {
  # proposalLocation = which(forest$tree$modifiable > 0) %>% sample(1)
  proposalLocation = proposal_location_generator(forest$tree)
  dots = list(...)[names(list(...)) %in% names(formals(generate_random_proposals))]
  proposal = do.call("generate_random_proposals",
                     list(forest[["tree"]],
                       dats[[1,"pts"]][["event"]]%>%unique()%>%as.character(),
                       proposals = 1,
                       dots$pars)
  )

  
  # Proposal type 1: branch
  if(proposal$proposalType[1] == 1)
    proposedTree = forest$tree %>% replace_leaf(active = proposalLocation[1],
                                       internalTrueFalseTree = proposal[1:3,], T)
  proposedForest = list(tree = proposedTree, roots = forest$roots)

  # TODO: a default mlp for having no data enter a "true/false" node --> currently unrepresented in mlss 
  
  #   proposal should be a list (or df) containing:
  #   - $proposalType
  #   - $existingTree
  #   - $roots
  #   - $proposalLocation (in tree)
  #   - $proposalTree
  prop = list()
  prop$proposalType = proposal$proposalType[1]
  prop$existingForest = forest
  prop$proposalLocation = proposalLocation
  prop$proposedForest = proposedForest
  # TODO sometimes trueNode or falseNode is not represented because no hits; need to give an MLP
  
  prop
}

learn_pc = function(forest, od, target, iter=10, pf = generate_one_proposal_po, maintainBlankTree=T,...) {
  model = list()
  model$forest = forest
  learn_pc_model(model, od, target, iter, pf, maintainBlankTree, ...)
}
learn_pc_model = function(model, od, target, iter, pf = generate_one_proposal_po, maintainBlankTree,...) {
  if(iter==0) return(model)
  cat(paste0(iter," "))
  
  prop = pf(model$forest, od, ...)
  
  proposalResults = evaluate_proposals_map(prop, od %>% ungroup(), target) %>%
    evaluate_proposals_reduce(prop)
  
  dll = proposalResults$dll %>% sum()
  # browser()
  
  if(is.nan(dll)) {
    # browser()
    warning("proposal caused dll to be NaN")
  } else if(dll > 1) {# AIC criterion (1 par --> 2 pars, net 1)
    # if(dll > 900)
    #   browser()
    # accept proposal
    if(prop$proposalType == 1) { #TODO debug this case
      log_qp = model$forest$tree[proposalResults$replacedIndex,"par1"]
      
      prop$proposedForest$tree = prop$proposedForest$tree %>% select(-replaced)
      prop$proposedForest$tree[proposalResults$treeIndex,"par1"] = proposalResults$log_qp_prime + log_qp
      model$forest = prop$proposedForest
      # model$forest$tree = prop$proposedTree %>% select(-replaced)
      # model$forest$tree[proposalResults$treeIndex, "par1"] = proposalResults$log_qp_prime + log_qp
    }
    
    if(maintainBlankTree) {
      # browser()
      if(sum(
            model[["forest"]][["tree"]][ model[["forest"]][["roots"]], "modifiable" ] > 0 &
            model[["forest"]][["tree"]][ model[["forest"]][["roots"]], "distribution" ]==1,
            na.rm = T
            ) == 0
      ) {
        dots = list(...)[names(list(...)) %in% names(formals(add_blank_tree))]
        if(length(dots)==0)
          model[["forest"]] = model[["forest"]] %>% add_blank_tree()
        else
          model[["forest"]] = model[["forest"]] %>% add_blank_tree(dots$modifiable)
      }
    }
    
    if(with(model,!exists("stats")))
      model$stats = data.frame(iter = iter, dll = dll) %>% tbl_df() %>% mutate(forest = list(model$forest))
    else
      model$stats = model$stats %>% bind_rows(data.frame(iter=iter, dll=dll) %>% tbl_df() %>% mutate(forest=list(model$forest)))
  }
  learn_pc_model(model, od, target, iter-1, pf, maintainBlankTree, ...)

}

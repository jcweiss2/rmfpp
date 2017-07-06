### simulation 1. single event, {periodic, quasiperiodic, quasiperiodic with missingness, scheduled but aperiodic, quasischeduled but aperiodic}

setwd("~/workspace/rmfpp/")
source("forest_modifications.r")
source("learn_one_random_proposal.r")

dats = data.frame( time = rep(seq(0,100,10),10)+rep(1:10,each=11),
                   id = rep(1:10, each=11) %>% as.character(),
                   event = "A") %>% tbl_df() %>%
  group_by(id) %>%
  do(pts = data.frame(.) %>% select(-id))

lbs = data.frame( lb = rep(0,10), ub = rep(110,10), id = 1:10 %>% as.character()) %>%
  group_by(id) %>%
  do(lbubs=data.frame(.) %>% select(-id))

dats = dats %>% left_join(by="id", lbs)

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

tree[dim(tree)[1]+1,] = c("" ,NA,NA,NA,1 ,NA,0 ,NA, 1)
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
forest = list(tree = tree, roots = 1)

set.seed(41); learn_pc(forest, dats, "A", iter=2, maintainBlankTree = F, pars=list(
                                  type="default",
                                  localityFactor = 2,
                                  maxLength=100,
                                  probSplit=1
                                ))


Rprof(tmp <- tempfile(), line.profiling=T)
output = learn_pc(forest,dats, "A",iter=20, maintainBlankTree = F, pars=list(
  type="default",
  localityFactor = 2,
  maxLength=100,
  probSplit=1
  ))
Rprof()
summaryRprof(tmp,lines = "show")

plot_pt_rate_graph = function(dats, model, index=1) {
  require(ggplot2)
  df = dats %>% 
    mutate(sts = list(subtimeline_forest_preexpansion(model[["forest"]][["tree"]],
                                                      pts,
                                                      model[["forest"]][["roots"]],
                                                      lbubs,
                                                      1))) %>% # DOTO shows an improper expansion!
    .[["sts"]] %>% 
    .[[index]] %>% 
    data.frame
  ggplot(df, aes(x=ub, y=logRateFactor)) +
    geom_point(data=dats$pts[[index]], aes(x=time,y=0, color=event)) +
    geom_segment(aes(y=logRateFactor, yend=logRateFactor,x=lb, xend=ub)) + xlab("Time") + theme_get()
}

plot_pt_rate_graph(dats,output,2) # nice result :)
if(exists("od")) plot_pt_rate_graph(od, output, 2)

# conclusions:
#   can we identify the best binary split point once we've selected an event type? that would be an improvement 
#   can we speed up some of the expansion operations (rate-limiting code)
#   can we extend (in a multitude of ways)?

#   can we add roots to make this a multiplicative forest
### subtimeline_maker_tests.r
# source("subtimeline_maker.r")

# basically, wrap everything in do() and compute rowwise along patients (map)
# then, we will aggregate over patients, i.e. we will take the mlp statistic components per patient and aggregate (reduce)

# things to reduce together:
# the mlp statistics M Mhat qt dll
# the timelines for log logistic
#

setwd("~/workspace/rmfpp")
source("subtimeline_maker.r")


# Example data
dat = data.frame(
  id=c("1","1","2","2"),
  event=c("A", "A", "B", "A"),
  time=c(1,2,2.2, 10)
) %>% tbl_df()

# od "object of data":
#   this object contains not only the data but accessory information, such as,
#     precomputes necessary for using extracts from the data
#   it should be static prior to entry into functions
#   requirements on "od" should be specified by functions
od = dat %>% group_by(id) %>% do(pts=data.frame(.) %>% select(-id) %>% tbl_df()) %>% ungroup()

lbs = data.frame(
  id=c(1,2,2,2),
  lb=c(0,1,3,8),
  ub=c(10,2,7,11)
) %>% tbl_df() %>% group_by(id) %>% do(lbubs=data.frame(.) %>% select(-id) %>% tbl_df()) %>% ungroup() %>% select(lbubs)

od = od %>% bind_cols(lbubs = lbs)

# Example tree [condition, value, true node, false node] 
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

tree[dim(tree)[1]+1,] = c("A",3 ,2 ,3 ,NA,NA ,NA,NA,0)
tree[dim(tree)[1]+1,] = c("B",2 ,4 ,5 ,NA,NA ,NA,NA,0)
tree[dim(tree)[1]+1,] = c("B" ,NA,NA,NA,2 ,0,1 ,1,1)
tree[dim(tree)[1]+1,] = c("" ,NA,NA,NA,1 ,NA ,2 ,NA,1)  # dist 1 := exponential
tree[dim(tree)[1]+1,] = c("A" ,NA,NA,NA,2 ,1,3 ,0.1,1) # dist 2 := log-logistic # time since "A", never = 
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

roots = 1 # c(1,...) index list of root nodes
forest = list(tree = tree, roots = roots)

### tests

od = od %>% 
  rowwise() %>% 
  # mutate(pts = (function(.) {y=list(); y$data = .; y})(pts)) %>%
  mutate(st = list(subtimeline_tree_ordered(tree, pts, roots, lbubs)))
subtimeline_tree_ordered(tree, od$pts[[1]], roots, od$lbubs[[1]])
  
product_of_timelines(od[[1,4]], od[[2,4]]) %>% product_of_timelines(od[[2,4]])

partitioned_union_of_durations(od$lbubs[[1]],od$lbubs[[2]]*2)

dat2 = data.frame(a=sort(runif(20)),b=sort(runif(20))) %>% tbl_df() %>% transmute(lb=pmin(a,b),ub=pmax(a,b))
dat2 %>% union_of_durations()

dat2 %>% union_of_durations() %>% intersection_of_durations_df(dat2 %>% union_of_durations())

od %>% mutate(sf = list(subtimeline_forest(tree,pts, roots,lbubs)))
od %>% mutate(sf = list(subtimeline_forest_preexpansion(tree,pts,rep(c(1,2,3),2),lbubs, 3) %>% arrange(lb, lb-ub) %>% data.frame))

od %>% mutate(sf = list(subtimeline_forest_preexpansion(tree,pts,rep(c(1,2,3),2),lbubs, 2) %>% arrange(lb, lb-ub) %>% data.frame))
od %>% mutate(sf = list(subtimeline_forest_preexpansion(tree,pts,rep(c(1,2,3),2),lbubs, 1) %>% arrange(lb, lb-ub) %>% data.frame))
od %>% mutate(sf = list(subtimeline_forest_preexpansion(tree,pts,rep(c(1,2,3),2),lbubs, 4) %>% arrange(lb, lb-ub) %>% data.frame))

# testing value_sum_onto_durations
valued = od$lbubs[[1]] %>% mutate(value=10)
partitions = subtimeline_forest_preexpansion(tree,od$pts[[1]],rep(c(1,2,3),2),od$lbubs[[1]], 3) %>% arrange(lb, lb-ub) %>%
  (function(x) fix_timelineCPP(x[[1]],x[[2]]))
valued %>% bind_rows(valued+2) %>% 
  mutate_active_duration_indices(partitions) %>%
  value_sum_onto_durations(partitions)

od$lbubs[[1]] %>% insert_events_into_expansion(c(2,7))

insert_events_into_expansionCPP(od$lbubs[[1]][[1]],od$lbubs[[1]][[2]], c(1.5,2,2.5,3,3.5))
insert_events_into_expansionCPP(od$lbubs[[2]][[1]],od$lbubs[[2]][[2]], c(1.5,3.5,3.5))
insert_events_into_expansionCPP(od$lbubs[[2]][[1]],od$lbubs[[2]][[2]],  c(0,1.5,3.5,3.5, 12))
insert_events_into_expansionCPP(od$lbubs[[1]][[1]],od$lbubs[[1]][[2]],  c(0,1.5,7,7)) #DOTO fix so M=2@7
insert_events_into_expansionCPP(od$lbubs[[1]][[1]],od$lbubs[[1]][[2]], c(9.5,9.5))


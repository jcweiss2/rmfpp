### subtimeline_maker_tests.r
# source("subtimeline_maker.r")
setwd("~/workspace/rmfpp")
source("subtimeline_maker.r")


# Example data
dat = data.frame(
  event=c("A", "A", "B", "A"),
  time=c(1,2,2.2, 10)
) %>% tbl_df()

# od "object of data":
#   this object contains not only the data but accessory information, such as,
#     precomputes necessary for using extracts from the data
#   it should be static prior to entry into functions
#   requirements on "od" should be specified by functions
od = list()
od$data = dat

lbubs = data.frame(
  lb=c(1,3,8),
  ub=c(2,7,11)
) %>% tbl_df()

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
  stringsAsFactors = F
) %>% tbl_df()

tree[dim(tree)[1]+1,] = c("A",3 ,2 ,3 ,NA,NA ,NA,NA)
tree[dim(tree)[1]+1,] = c("B",2 ,4 ,5 ,NA,NA ,NA,NA)
tree[dim(tree)[1]+1,] = c("B" ,NA,NA,NA,2 ,0,1 ,1)
tree[dim(tree)[1]+1,] = c("" ,NA,NA,NA,1 ,NA ,2 ,NA)  # dist 1 := exponential
tree[dim(tree)[1]+1,] = c("B" ,NA,NA,NA,2 ,0,3 ,0.1) # dist 2 := log-logistic # time since "A", never = 
tree$timeFrame = as.numeric(tree$timeFrame)
tree$distribution = as.numeric(tree$distribution)
tree$par0 = as.character(tree$par0)
tree$par1 = as.numeric(tree$par1)
tree$par2 = as.numeric(tree$par2)
tree$trueNode = as.numeric(tree$trueNode)
tree$falseNode = as.numeric(tree$falseNode)

roots = 1 # c(1,...) index list of root nodes


# tests
product_of_timelinesCPP(c(0,1,2), c(1,2,2.5), c(1,2,3),
                        c(0, 0.5, 1.0), c(0.5,0.6, 2.5), c(10,20,30)
) %>% data.frame() %>% tbl_df() 
partitioned_union_of_durationsCPP(c(0,1,2), c(1,2,2.5), 
                                  c(0, 0.5, 1.0), c(0.5,0.6, 2.5)
)





###### test (random)

st = subtimeline_tree_ordered(tree,od$data,roots,lbubs) 
product_of_timelines(st,st)
product_of_timelines(st,st) %>% product_of_timelines(st)

partitioned_union_of_durations(lbubs,lbubs*2)

dat2 = data.frame(a=sort(runif(20)),b=sort(runif(20))) %>% tbl_df() %>% transmute(lb=pmin(a,b),ub=pmax(a,b))
dat2 %>% union_of_durations()

dat3 = data.frame(a=sort(runif(10)),b=sort(runif(10))) %>% tbl_df() %>% transmute(lb=pmin(a,b),ub=pmax(a,b)) %>% arrange(lb)
dat3 %>% union_of_durations()

dat2 %>% union_of_durations() %>% intersection_of_durations_df(dat3 %>% union_of_durations())
dat3 %>% union_of_durations() %>% intersection_of_durations_df(dat2 %>% union_of_durations())

subtimeline_forest(tree,od$data, roots,lbubs)
subtimeline_forest_preexpansion(tree,od$data,rep(c(1,2,3),2),lbubs, 3) %>% arrange(lb, lb-ub) %>% data.frame
subtimeline_forest_preexpansion(tree,od$data,rep(c(1,2,3),2),lbubs, 2) %>% arrange(lb, lb-ub) %>% data.frame
subtimeline_forest_preexpansion(tree,od$data,rep(c(1,2,3),2),lbubs, 1) %>% arrange(lb, lb-ub) %>% data.frame
subtimeline_forest_preexpansion(tree,od$data,rep(c(1,2,3),2),lbubs, 4) %>% arrange(lb, lb-ub) %>% data.frame

# testing value_sum_onto_durations
valued = lbubs %>% mutate(value=10)
partitions = subtimeline_forest_preexpansion(tree,od$data,rep(c(1,2,3),2),lbubs, 3) %>% arrange(lb, lb-ub) %>%
  (function(x) fix_timelineCPP(x[[1]],x[[2]]))
valued %>% bind_rows(valued+2) %>% 
  mutate_active_duration_indices(partitions) %>%
  value_sum_onto_durations(partitions)

lbubs %>% insert_events_into_expansion(c(2,7))

insert_events_into_expansionCPP(lbubs[[1]],lbubs[[2]], c(1.5,2,2.5,3,3.5))
insert_events_into_expansionCPP(lbubs[[1]],lbubs[[2]], c(1.5,3.5,3.5))
insert_events_into_expansionCPP(lbubs[[1]],lbubs[[2]], c(0,1.5,3.5,3.5, 12))
insert_events_into_expansionCPP(lbubs[[1]],lbubs[[2]], c(0,1.5,7,7))


union_of_durationsCPP(lbubs$lb+10, lbubs$ub+10,lbubs$lb, lbubs$ub)
union_of_durationsCPP(lbubs$lb, lbubs$ub, lbubs$lb+10, lbubs$ub+10)
intersection_of_durationsCPP(lbubs$lb+10, lbubs$ub+10,lbubs$lb, lbubs$ub)
intersection_of_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb+10, lbubs$ub+10)

union_of_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb, lbubs$ub)
union_of_durationsCPP(lbubs$lb, lbubs$ub, lbubs$lb*2, lbubs$ub*2)
intersection_of_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb, lbubs$ub)
intersection_of_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb+0.1, lbubs$ub+0.1)
intersection_of_durationsCPP(lbubs$lb+0.1, lbubs$ub+0.1,lbubs$lb, lbubs$ub)

union_of_duration_columnCPP(c(1,4,5),c(3,4,10))
union_of_durations_df(lbubs)

complement_in_durationsCPP(lbubs$lb+10, lbubs$ub+10,lbubs$lb, lbubs$ub)
complement_in_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb+10, lbubs$ub+10)
complement_in_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb+0.1, lbubs$ub+0.1)
complement_in_durationsCPP(lbubs$lb+0.1, lbubs$ub+0.1,lbubs$lb, lbubs$ub)
complement_in_durationsCPP(lbubs$lb, lbubs$ub,lbubs$lb[0], lbubs$ub[0])
complement_in_durationsCPP(lbubs$lb[0], lbubs$ub[0],lbubs$lb, lbubs$ub)


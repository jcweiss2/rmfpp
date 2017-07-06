library(dplyr)
setwd("~/workspace/pcl")
source("piecewise_constant_converter.r")

# requires 'lb' in first column and 'ub' in second column ordered by lb
union_of_durations = function(durations) {
  durations %>% tbl_df() %>% # TODO slow!
    mutate(largestseen=cummax(ub)) %>%
    mutate(largestseenlagged=lag(largestseen)) %>% 
    filter((largestseenlagged-lb)<0 | is.na(largestseenlagged)) %>%
    mutate(unionub=lead(largestseenlagged,default=durations$ub[dim(durations)[1]])) %>%
    select(lb,ub=unionub)
}

union_of_durations_df = function(durations) {
  union_of_duration_columnCPP(durations[[1]],durations[[2]]) %>%
    (function(.) if(nrow(.)==0) durations[0,1:2] else .)(.) %>%
    tbl_df()
  # union_of_durations(durations)
}

library(Rcpp)
cppFunction(
  '
  DataFrame union_of_duration_columnCPP(NumericVector lb1,
  NumericVector ub1) {
    int n1 = lb1.size();
    int largestUb = DBL_MIN;
    int i1 = 0, oi=0;
    std::vector<double> lb(n1), ub(n1);
    for(int i=0; i<n1; ++i) {
      if(lb1[i]>largestUb) {
        if(oi!=0) // finish previous line
          ub[oi-1] = largestUb;
        lb[oi++] = lb1[i]; // start next line
        largestUb = ub1[i];
      } else
        largestUb = largestUb > ub1[i] ? largestUb : ub1[i];
    }
    if(oi>0)
      ub[oi-1] = largestUb;
    lb.resize(oi);
    ub.resize(oi);
    return DataFrame::create( Named("lb")=lb, Named("ub")=ub );
  // code is untested, returns union (partitioned) of durations with sum of external values provided
  }
  ', includes="#include <cfloat>" 
)

cppFunction(
  '
  DataFrame union_of_durationsCPP(NumericVector lb1,
  NumericVector ub1,
  NumericVector lb2, 
  NumericVector ub2) {
  /** Returns the union of durations of ordered (lb,ub] in O(t+s) **/
    int n1 = lb1.size();
    int n2 = lb2.size();
    int on1 = 0, on2 = 0;
    int i1 = 0, i2 = 0, oi=0;
    std::vector<double> lb(n1+n2), ub(n1+n2);
    if(n1 == 0 || n2 == 0)
      return DataFrame::create();
    if(n1 != ub1.size() || n2 != ub2.size())
      return DataFrame::create();
    for(int i=0; i<n1; ++i)
      if(ub1[i]<lb1[i])
        return DataFrame::create();
    for(int i=0; i<n2; ++i)
      if(ub2[i]<lb2[i])
        return DataFrame::create();
    double v1 = lb1[0], v2 = lb2[0];
    while((i1 <= n1 || on1) && (i2 <= n2 || on2)) {
      if((i1 == n1 && !on1) && (i2 == n1 && !on2))
        break;
      if(v2 < v1 || (v2==v1 && on1)) { // step 2
        if(!on1 && !on2) // going to be on
          lb[oi] = v2;
        if(!on1 && on2) {  // going to be off
          ub[oi] = v2;
          if(ub[oi]!=lb[oi])
            ++oi;
        }
        if(on2 && i2==n2) { // going to be off irrespective of on1
            v2 = DBL_MAX;
            on2 = 0;
            continue;
        }
        v2 = on2? lb2[i2]:ub2[i2]; // advance
        ++on2;
        if(on2 == 2)
          on2 = 0;
        else
          ++i2;
      } else {
        if(!on1 && !on2) // going to be on
          lb[oi] = v1;
        if(on1 && !on2) { // going to be off
          ub[oi] = v1;
          if(ub[oi]!=lb[oi])
            ++oi;
        }
        if(on1 && i1==n1) { // going to be off irrespective of on2
          v1 = DBL_MAX;
          on1 = 0;
          continue;
        }
        v1 = on1? lb1[i1]:ub1[i1]; // advance
        ++on1;
        if(on1 == 2)
          on1 = 0;
        else
          ++i1;
      }
    }
    lb.resize(oi);
    ub.resize(oi);
  
  /*List out(8);  // for debugging
  out[0] = p0;
  out[1] = p1;
  out[2] = product;
  out[3] = poo;
  out[4] = taLong;
  out[5] = saLong;
  return out;
  */
  return DataFrame::create( Named("lb")=lb, Named("ub")=ub );
  // code is untested, returns union (partitioned) of durations with sum of external values provided
  }
  ', includes="#include <cfloat>" 
)


# would be helpful for intersection of durations (complement of union),
#   but this is of type A AND not B (which requires an AND anyways)
# complement_of_durations = function(durations, in_durations) {}

# requires ordered non-overlapping data.frame durations
intersection_of_durations_df = function(durations1,durations2) {
  return(intersection_of_durationsCPP(durations1[[1]],durations1[[2]],durations2[[1]],durations2[[2]])) %>% 
    (function(.) if(nrow(.)==0) durations1[0,] else .)(.) %>%
    tbl_df()
  # return(
  #   intersection_of_durations(
  #     durations1 %>% data.frame(),
  #     durations2 %>% data.frame()
  #   )
  # )
}

intersection_of_durations = function(durations1, durations2) {
  data.frame( # TODO slow!
    lb=c(durations1[,1],
         durations1[,2],
         durations2[,1],
         durations2[,2]
         ),
    markOn=c(rep(1, dim(durations1)[1]),
            rep(0, dim(durations1)[1]),
            rep(1, dim(durations2)[1]),
            rep(0, dim(durations2)[1])
            ),
    d=c(rep(1,dim(durations1)[1]*2),
        rep(2,dim(durations2)[1]*2)
        )
  ) %>% tbl_df() %>%
    arrange(lb) %>%
    mutate(d1ison=(cumsum(markOn==1 & d==1)-cumsum(markOn==0 & d==1))) %>%
    mutate(d2ison=(cumsum(markOn==1 & d==2)-cumsum(markOn==0 & d==2))) %>%
    mutate(intersect=d1ison&d2ison) %>%
    mutate(ub=lead(lb)) %>%
    filter(intersect & lb!=ub) %>%
    select(lb,ub)
}

library(Rcpp)
cppFunction(
  '
  DataFrame intersection_of_durationsCPP(NumericVector lb1,
  NumericVector ub1,
  NumericVector lb2, 
  NumericVector ub2) {
  /** Returns the intersection of durations of ordered (by lb1, lb2) (lb,ub] in O(t+s) **/
  int n1 = lb1.size();
  int n2 = lb2.size();
  int on1 = 0, on2 = 0;
  int i1 = 0, i2 = 0, oi=0;
  std::vector<double> lb(n1+n2), ub(n1+n2);
  if(n1 == 0 || n2 == 0)
    return DataFrame::create();
  if(n1 != ub1.size() || n2 != ub2.size())
    return DataFrame::create();
  for(int i=0; i<n1; ++i)
    if(ub1[i]<lb1[i])
      return DataFrame::create();
  for(int i=0; i<n2; ++i)
    if(ub2[i]<lb2[i])
      return DataFrame::create();
  double v1 = lb1[0], v2 = lb2[0];
  while((i1 < n1 || on1) && (i2 < n2 || on2)) {
      if(v2 < v1) { // step 2
        if(on1 && !on2) // going to be on
          lb[oi] = v2;
        if(on1 && on2) {  // going to be off
          ub[oi] = v2;
          if(ub[oi]!=lb[oi])
            ++oi;
          if(i2==n2)
            break;
        } 
        v2 = on2? lb2[i2]:ub2[i2]; // advance
        ++on2;
        if(on2 == 2)
          on2 = 0;
        else
          ++i2;
        if(i2 == n2 && !on2)
          break;
      } else {
        if(!on1 && on2) // going to be on
          lb[oi] = v1<v2? v1 : v2;
        if(on1 && on2) { // going to be off
          ub[oi] = v1;
          if(ub[oi]!=lb[oi])
            ++oi;
          if(i1==n1)
            break;
        }
        v1 = on1? lb1[i1]:ub1[i1]; // advance
        ++on1;
        if(on1 == 2)
          on1 = 0;
        else
          ++i1;
        if(i1 == n1 && !on1)
          break;
      }
  }
  lb.resize(oi);
  ub.resize(oi);
  
  /*List out(8);  // for debugging
  out[0] = p0;
  out[1] = p1;
  out[2] = product;
  out[3] = poo;
  out[4] = taLong;
  out[5] = saLong;
  return out;
  */
  return DataFrame::create( Named("lb")=lb, Named("ub")=ub );
  // code is untested, returns union (partitioned) of durations with sum of external values provided
  }
  ', includes="#include <cfloat>" 
)


### Combine subtimelines from multiple trees
# partitioned_union_of_durations() with log rate factors summed
cppFunction(
  '
  DataFrame partitioned_union_of_durationsCPP(NumericVector t0,
    NumericVector t1, 
    NumericVector s0,
    NumericVector s1) {
  /** Returns the product of durations of (t0,t1] and (s0,s1] in O(t+s) **/
  int nT = t0.size(), nT2 = nT*2;
  int nS = s0.size(), nS2 = nS*2;
  NumericVector t(nT2), oo(nS2+nT2), s(nS2); //t=[t0_0 t1_0 t0_1 t1_1 ...], oo=[0 1 0 1 0 ...], ...
  for (int i=0, i2=0; i<nT; ++i, i2+=2 ) {
    t[i2] = t0[i];
    t[i2+1] = t1[i];
    oo[i2+1] = 1; // indicates "on prior to t"
  }
  for (int i=0, i2=0; i<nS; ++i, i2+=2) {
    s[i2] = s0[i];
    s[i2+1] = s1[i];
    oo[i2+1] = 1;
  }
  NumericVector product(nS2+nT2), poo(nS2+nT2);
  int pi = 0, ti = 0, si = 0;
  for( ; ti < nT2; ++ti) {
    while( nS!=0 && s[si] < t[ti] ) {
      product[pi] = s[si];
      poo[pi] = ( oo[si] || oo[ti] )? 1 : 0;
      ++pi;
      ++si;
      if( si == nS2 )
        break;
    }
    if(si==nS2) { // finish up ti\'s
      for(; ti<nT2; ++ti) {
        product[pi] = t[ti];
        poo[pi] = oo[ti] ? 1 : 0;
        ++pi;
      }
      break;
    }
  
    // s[si] >= t[ti];
    product[pi] = t[ti];
    poo[pi] = ( oo[si] || oo[ti] )? 1 : 0;
    ++pi;
  }
  if(si!=nS2) { // finish up si\'s
    for( ; si<nS2; ++si) {
      product[pi] = s[si];
      poo[pi] = oo[si] ? 1 : 0; 
      ++pi;
    }
  }
  
  std::vector<double> p0(nS2+nT2), p1(nS2+nT2);
  pi = 0;
  int p01i = 0;
  for( ; pi<nS2+nT2-1; ++pi ) {
    if(poo[pi+1] != 0  &&  product[pi] != product[pi+1]) {
      p0[p01i] = product[pi];
      p1[p01i] = product[pi+1];
      ++p01i;
    }
  }
  p0.resize(p01i);
  p1.resize(p01i);
  
  /*List out(8);  // for debugging
  out[0] = p0;
  out[1] = p1;
  out[2] = product;
  out[3] = poo;
  out[4] = taLong;
  out[5] = saLong;
  return out;
  */
  return DataFrame::create( Named("lb")=p0, Named("ub")=p1 );
  //Named("product")=product, Named("poo")=poo,
  //Named("t01")=t, Named("s01")=s);
  
  // code is untested, returns union (partitioned) of durations with sum of external values provided
  //TODO generalize
  }
  '
)
partitioned_union_of_durations = function(timeline1, timeline2) {
  return(partitioned_union_of_durationsCPP(timeline1[[1]],
                                 timeline1[[2]],
                                 timeline2[[1]],
                                 timeline2[[2]]) %>% data.frame() %>% tbl_df() )
}

cppFunction(
  '
  DataFrame fix_timelineCPP(NumericVector t0, NumericVector t1) {
  /** Resolves overlying durations of (t0,t1] in O(t). Requires lb\'s sorted **/
  int nT = t0.size(), nT2 = nT*2;
  std::vector<double> o0(nT2), o1(nT2);
  double ubs[nT];
  for (int i=0; i<nT; ++i)
    if(t0[i] > t1[i])
      return DataFrame::create();

  int ti=0, oi=0, ubi=0;
  double accepted_past = t0[0];
  double lb = t0[0];
  double ub = t1[0];
  while(1) {
    if( ti+1 < nT ) {
      if(t0[ti+1] >= ub) {  // "row" is in order wrt next row 
        o0[oi] = lb;
        o1[oi] = ub;
        if(lb != ub)
          ++oi;
        std::sort(ubs, ubs+ubi);
        if(ubi>0 && ubs[0] < t0[ti+1]) { // ubs non-empty, fill in
          lb = ub;
          ub = ubs[0];
          ubs[0] = DBL_MAX;
        } else {                // if ubs is empty, jump up
          ++ti;
          lb = t0[ti];

          double lowest = t1[ti];
          int lowesti = -1;
          for(int j=0; j<ubi; ++j)
            if(lowest > ubs[j]) {
              lowest = ubs[j];
              lowesti = j;
            }
          ub = lowest; // should be the min of (ubs,ti[ti])
          if(lowesti != -1) // "remove" ubs[j] from ubs
            ubs[lowesti] = DBL_MAX;
        }
      } else { // "row" is not in order wrt next row
          o0[oi] = lb;
          o1[oi] = t0[ti+1];
          lb = t0[ti+1];
          if(ub > t1[ti+1])
            ubs[ubi++] = ub;
          else
            ubs[ubi++] = t1[ti+1];
          ub = (ub < t1[ti+1])? ub : t1[ti+1];
          if(o0[oi] != o1[oi])
            ++oi;
          ++ti;
      }
    } else if(ti < nT) { // last iteration through t0
      o0[oi] = lb;
      o1[oi] = (ub < t1[ti])? ub : t1[ti];
      ++oi;
      ++ti;
      std::sort(ubs, ubs+ubi);
      for(int j=0; j<ubi; ++j) {
        if(ubs[j] == DBL_MAX)
          break;
        o0[oi] = ub;
        o1[oi] = ubs[j];
        ub = ubs[j];
        if(o0[oi] != o1[oi])
          ++oi;
      }

    } else { // wrap up
      break;
    }
  }
  if(oi>0 && o0[oi-1] == o1[oi-1])
    --oi;
  o0.resize(oi);
  o1.resize(oi);

  return DataFrame::create( Named("lb")=o0, Named("ub")=o1 );

/*
List out(8);  // for debugging
  out[0] = o0;
  out[1] = o1;
  NumericVector ubsv(nT);
  for(int j=0; j<ubi; ++j)
    ubsv[j] = ubs[j];
  out[2] = ubsv;
  return out;
*/

  }
  ', includes="#include <cfloat>" 
) # it works on "bad" = lbubs with lbubs[1,2] = 4 


complement_in_durations_df = function(inDurations, durations) {
  return(complement_in_durationsCPP(inDurations[[1]],inDurations[[2]],durations[[1]],durations[[2]])) %>% 
    (function(.) if(nrow(.)==0) inDurations[0,] else .)(.) %>%
    tbl_df()
  # return(
  #   complement_of_durations(
  #     inDurations %>% data.frame(),
  #     durations %>% data.frame()
  #   )
  # )
}

# maybe simplify this and the other duration functions with Rcpp
complement_of_durations = function(inDurations, durations) {
  data.frame( #TODO slow!
    lb=c(inDurations[,1],
         inDurations[,2],
         durations[,1],
         durations[,2]
    ),
    markOn=c(rep(1, dim(inDurations)[1]),
             rep(0, dim(inDurations)[1]),
             rep(1, dim(durations)[1]),
             rep(0, dim(durations)[1])
    ),
    d=c(rep(1,dim(inDurations)[1]*2),
        rep(2,dim(durations)[1]*2)
    )
  ) %>% tbl_df() %>%
    arrange(lb) %>%
    mutate(idison=(cumsum(markOn==1 & d==1)-cumsum(markOn==0 & d==1))) %>%
    mutate(dison=(cumsum(markOn==1 & d==2)-cumsum(markOn==0 & d==2))) %>%
    mutate(disoff=!dison) %>%
    mutate(intersect=idison&disoff) %>%
    mutate(ub=lead(lb)) %>%
    filter(intersect) %>%
    filter(lb!=ub) %>%
    select(lb,ub)
} 
library(Rcpp)
cppFunction(
'
  DataFrame complement_in_durationsCPP(NumericVector lb1,
    NumericVector ub1,
    NumericVector lb2, 
    NumericVector ub2) {
    /** Returns, within (lb1,ub1], complement of (lb2,ub2] in O(t+s) **/
    int n1 = lb1.size();
    int n2 = lb2.size();
    int on1 = 0, on2 = 0;
    int i1 = 0, i2 = 0, oi=0;
    std::vector<double> lb(n1+n2), ub(n1+n2);
    if(n1 == 0)
      return DataFrame::create();
    if(n1 != ub1.size() || n2 != ub2.size())
      return DataFrame::create();
    for(int i=0; i<n1; ++i)
      if(ub1[i]<lb1[i])
        return DataFrame::create();
    for(int i=0; i<n2; ++i)
      if(ub2[i]<lb2[i])
        return DataFrame::create();
    double v1 = lb1[0], v2;
    if(n2 == 0)
      v2 = DBL_MAX;
    else
      v2 = lb2[0];
    lb[0] = lb1[0];
    while((i1 < n1 || on1) && (i2 <= n2)) {
      if(v2 < v1) { // step 2
        if(on1 && !on2) { // going to be on
          ub[oi] = v2;
          if(ub[oi]!=lb[oi])
            ++oi;
        }
        if(on1 && on2)  // going to be off
          lb[oi] = v2;
        if(on2 && i2==n2) { // going to be off irrespective of on1
          v2 = DBL_MAX;
          on2 = 0;
          continue;
        }
        v2 = on2? lb2[i2]:ub2[i2]; // advance
        ++on2;
        if(on2 == 2)
          on2 = 0;
        else
          ++i2;
        //if(i2 == n2 && !on2)
        //  break;
      } else {
        if(!on1 && !on2) // going to be on
          lb[oi] = v1;
        if(on1 && !on2) { // going to be off and recording
          ub[oi] = v1;
          if(ub[oi]!=lb[oi])
            ++oi;
          if(i1==n1)
            break;
        }
        v1 = on1? lb1[i1]:ub1[i1]; // advance
        ++on1;
        if(on1 == 2)
          on1 = 0;
        else
          ++i1;
        //if(i1 == n1 && !on1)
        //  break;
      }
    }
    lb.resize(oi);
    ub.resize(oi);

/*List out(8);  // for debugging
out[0] = p0;
out[1] = p1;
out[2] = product;
out[3] = poo;
out[4] = taLong;
out[5] = saLong;
return out;
*/
    return DataFrame::create( Named("lb")=lb, Named("ub")=ub );
    // code is untested, returns union (partitioned) of durations with sum of external values provided
}
', includes="#include <cfloat>" 
)



# return types for subtimeline_tree
piecewise_expansion = function(lbubs,leaf) {
  lbubs %>% rowwise() %>% do(piecewise_expansion_row(.,leaf)) %>% ungroup()
}
piecewise_expansion_row = function(lbubs, leaf) {
  if(leaf[["distribution"]]==2) { 
    if(!is.tbl(lbubs)) lbubs = lbubs %>% tbl_df()
    # browser()
    if(is.na(as.numeric(leaf[["par0"]])))
      browser()
    expandedLbubs = piecewise_constant_subtimeline(p=pllog,
                                                   q=qllog,
                                                   scale=leaf[["par1"]],
                                                   shape=leaf[["par2"]],
                                                   t0=lbubs[["mostrecenttime"]],
                                                   lb=lbubs[[1]],
                                                   ub=lbubs[[2]],
                                                   offset = NA,
                                                   as.hazard = T,
                                                   notinf.default = as.numeric(leaf[["par0"]])) # %>%
      #mutate(logRateFactor=logRateFactor+as.numeric(leaf[["par0"]]))
    expandedLbubs = expandedLbubs %>%
      do((function(ex, x) if("replacedIndex" %in% names(x)) ex %>% mutate(replacedIndex=x$replacedIndex) else ex)(.,lbubs)) %>%
    do((function(ex, x) if("treeIndex" %in% names(x)) ex %>% mutate(treeIndex=x$treeIndex) else ex)(.,lbubs)) %>%
      do((function(ex, x) if("oldIndexInNewTree" %in% names(x)) ex %>% mutate(oldIndexInNewTree=x$oldIndexInNewTree) else ex)(.,lbubs)) %>%
      do((function(ex, x) if("treeNumber" %in% names(x)) ex %>% mutate(treeNumber=x$treeNumber) else ex)(.,lbubs)) %>%
      do((function(ex) if(ncol(ex)==3) ex else ex[,c(1,2,4:ncol(ex),3)])(.))
    return(expandedLbubs)
  } else
    stop("Unimplemented distribution in PCL tree")
}
 
# leafType
#   1. expanded
#   2. store distribution parameters
#   3. store distribution parameters with extracted and processed data
#   4. keep treeIndex of tree data structure (and possibly the replacedIndex)
ret_subtimeline_tree = function(tree, data, active, lbubs, leafType) {
  if(length(active) > 1 && (active %>% unique() %>% length()) == 1) active = active[1]
  if(leafType==1) { # expanded
    if(tree[[active,"distribution"]]==1)
      return(lbubs %>% mutate(logRateFactor=tree[[active,"par1"]]))
    else if (tree[[active,"distribution"]]==2) {
      if(nrow(lbubs)==0) return(lbubs %>% mutate(logRateFactor = numeric()))
      lbtemp = data.frame(lb = c(lbubs$lb[1],data$time)) %>% tbl_df() %>%
        mutate(ub = lead(lb, default=lbubs$ub[nrow(lbubs)]),
               mostrecenttime = c(-Inf, lb[-1])) %>%
        (function(x) if(x[[1,"lb"]]       >= x[[1,"ub"]]      ) x[-1,]       else x)(.) %>%
        (function(x) if(x[[nrow(x),"lb"]] >= x[[nrow(x),"ub"]]) x[-nrow(x),] else x)(.)
      if(nrow(lbubs)>1 || lbubs$lb[1] > lbtemp$lb[1] || lbubs$ub[nrow(lbubs)] < lbtemp$ub[nrow(lbtemp)] ) #then intersect it
        lbtemp = lbtemp %>% select(-mostrecenttime) %>%
          intersection_of_durations_df(lbubs) %>%
          rowwise() %>%
          mutate(mostrecenttime=precedingtime(tree[[active,"condition"]], lb, data)) # old "par0"
      # lbtemp = lbubs %>% rowwise() %>% mutate(mostrecenttime=precedingtime(tree[[active,"par0"]], lb, data))
      # TODO not keeping the replacedIndex as required by get_ll_mapobject
      return(piecewise_expansion(lbtemp %>%
                                   mutate(distribution=tree[[active,"distribution"]],
                                          par1=tree[[active,"par1"]],
                                          par2=tree[[active,"par2"]]),
                                 tree[active,]))
    }
  } else if(leafType==2) { # distribution
    return(lbubs %>% mutate(condition = tree[[active,"condition"]],
                            distribution=tree[[active,"distribution"]],
                            par0=tree[[active,"par0"]],
                            par1=tree[[active,"par1"]],
                            par2=tree[[active,"par2"]])
    )
  } else if(leafType==3) { # distribution with log-logistic par0 replaced with most recent time of event-type
    #lbubs %>% rowwise() %>% mutate(mostrecenttime=precedingtime(tree[[active,"par0"]], lb, data))
    lbtemp = lbubs %>% rowwise() %>% mutate(mostrecenttime=precedingtime(tree[[active,"condition"]], lb, data)) #old "par0", (no more going down leaves so no worries) 
    return(lbtemp %>%
             mutate(condition=tree[[active,"condition"]],
               distribution=tree[[active,"distribution"]],
                          par1=tree[[active,"par1"]],
                          par2=tree[[active,"par2"]])
    )
  } else { # indexed
    return(lbubs %>% mutate(treeIndex=active))
    # for a proposal: tree index, modification "suff stats" estimates, not just "countsIn"
  }
}

### Creates the subtimeline for a tree
subtimeline_tree_ordered = function(tree, data, active, lbubs, leafType=1) {
  return(subtimeline_tree(tree, data, active, lbubs, leafType) %>% arrange(lb))
}
subtimeline_tree = function(tree, data, active, lbubs, leafType=1) {
  if(tree[[active,1]]=="" | (!is.na(tree[[active,"distribution"]]) && tree[[active,"distribution"]] == 2))
    return(ret_subtimeline_tree(tree, data, active, lbubs, leafType)) # later relax to quantized distribution;
    # let the data enter? make more complicated data structure? (list[[1]] = data; ...)
  trueBranch = data %>% #TODO slow! consider passing indices instead of actually filtering
    filter(event==tree[[active,1]]) %>% # later relax "==" to "%in%"
    select(-event) %>%
    mutate(end_t=time+tree[[active,2]]) %>%
    rename(lb=time,ub=end_t) %>%
    union_of_durations_df() %>%
    # union_of_durations() %>%
    intersection_of_durations_df(lbubs) # true arm active durations
  falseBranch = lbubs %>% complement_in_durations_df(trueBranch)
  if(leafType==4 && ("replaced" %in% names(tree)) && tree[active,"replaced"]==T) { #TODO slow!
    # possibility of overwriting with last-visited replacement index #i.e. don't sample proposal_locations with replacement
    trueBranch = trueBranch %>% mutate(replacedIndex = active)
    falseBranch = falseBranch %>% mutate(replacedIndex = active)
  }
  return(
    subtimeline_tree(tree, data, tree[[active,"trueNode"]], trueBranch, leafType) %>%
      bind_rows(subtimeline_tree(tree, data, tree[[active,"falseNode"]], falseBranch, leafType)) # %>%
      # arrange(ts) # wait until all parts are back and sort once per tree
  )
}

### Combine subtimelines from multiple trees
# product_of_durations() with log rate factors summed
library(Rcpp)
cppFunction(
  '
DataFrame product_of_timelinesCPP(NumericVector t0, NumericVector t1, NumericVector ta,
                                   NumericVector s0, NumericVector s1, NumericVector sa) {
/** Returns the product of durations of (t0,t1] and (s0,s1] with ta+sa elements in O(t+s) **/
  int nT = t0.size(), nT2 = nT*2;
  int nS = s0.size(), nS2 = nS*2;
  NumericVector t(nT2), oo(nS2+nT2), s(nS2); //t=[t0_0 t1_0 t0_1 t1_1 ...], oo=[1 0 1 0 ...], ...
  NumericVector taLong(nT2), saLong(nS2);
  for (int i=0, i2=0; i<nT; ++i, i2+=2 ) {
    t[i2] = t0[i];
    t[i2+1] = t1[i];
    //taLong[i2] = ta[i];
    taLong[i2+1] = ta[i];
    oo[i2+1] = 1; // indicates "on prior to t"
  }
  for (int i=0, i2=0; i<nS; ++i, i2+=2) {
    s[i2] = s0[i];
    s[i2+1] = s1[i];
    //saLong[i2] = sa[i];
    saLong[i2+1] = sa[i];
    oo[i2+1] = 1;
  }
  NumericVector product(nS2+nT2), poo(nS2+nT2), pxy(nS2+nT2);
  int pi = 0, ti = 0, si = 0;
  for( ; ti < nT2; ++ti) {
    while( nS!=0 && s[si] < t[ti] ) {
      product[pi] = s[si];
      poo[pi] = ( oo[si] || oo[ti] )? 1 : 0;
      pxy[pi] = taLong[ti] + saLong[si];
      ++pi;
      ++si;
      if( si == nS2 )
        break;
    }
    if(si==nS2) { // finish up ti\'s
      for(; ti<nT2; ++ti) {
        product[pi] = t[ti];
        poo[pi] = oo[ti] ? 1 : 0;
        pxy[pi] = taLong[ti];
        ++pi;
      }
      break;
    }

    // s[si] >= t[ti];
    product[pi] = t[ti];
    poo[pi] = ( oo[si] || oo[ti] )? 1 : 0;
    pxy[pi] = taLong[ti] + saLong[si];
    ++pi;
  }
  if(si!=nS2) { // finish up si\'s
    for( ; si<nS2; ++si) {
      product[pi] = s[si];
      poo[pi] = oo[si] ? 1 : 0; 
      pxy[pi] = saLong[si];
      ++pi;
    }
  }

  std::vector<double> p0(nS2+nT2), p1(nS2+nT2), pxyLong(nS2+nT2);
  pi = 0;
  int p01i = 0;
  for( ; pi<nS2+nT2-1; ++pi ) {
    if(poo[pi+1] != 0  &&  product[pi] != product[pi+1]) {
      p0[p01i] = product[pi];
      p1[p01i] = product[pi+1];
      pxyLong[p01i] = pxy[pi+1];
      ++p01i;
    }
  }
  p0.resize(p01i);
  p1.resize(p01i);
  pxyLong.resize(p01i);

  /*List out(8);  // for debugging
  out[0] = p0;
  out[1] = p1;
  out[2] = product;
  out[3] = poo;
  out[4] = taLong;
  out[5] = saLong;
  out[6] = pxy;
  out[7] = pxyLong;
  return out;
  */
  return DataFrame::create( Named("lb")=p0, Named("ub")=p1, Named("logRateFactor")=pxyLong);
                            //Named("product")=product, Named("poo")=poo,
                            //Named("t01")=t, Named("s01")=s);
  
  // code is untested, returns union (partitioned) of durations with sum of external values provided
  //TODO generalize
}
  '
)

# require expanded trees (distributions expanded into approximations)
product_of_timelines = function(timeline1, timeline2) {
  return(product_of_timelinesCPP(timeline1[[1]],
                                 timeline1[[2]],
                                 timeline1[[3]],
                                 timeline2[[1]],
                                 timeline2[[2]],
                                 timeline2[[3]]) %>% data.frame() %>% tbl_df() )
}

product_of_n_timelines = function(timelines, result=NULL) {
  if(length(timelines)==1) return(product_of_timelines(timelines[[1]], result))
  else {
    return(product_of_n_timelines(timelines[-1], timelines[[1]]))
  }
}

#DOTO change to binary split instead of iterative
#TODO change into for loop in C++ or map reduce this
subtimeline_forest_preexpansion = function(tree,
                                           data,
                                           active,
                                           lbubs,
                                           leafType=1,
                                           depth=ceiling(log2(length(active)))-1,
                                           numbering=(0:(length(active)-1))
                                           ) {
  if(leafType==1) return(subtimeline_forest(tree, data, active, lbubs))
  if(length(active)==1) {
    return(subtimeline_tree_ordered(tree, data, active, lbubs, leafType) %>%
             mutate(treeNumber = numbering)
    )
  } else { #TODO use do(...) %>% bind_rows()
    halfway = ceiling(0.5*length(active))
    return(
      subtimeline_forest_preexpansion(tree,
                                      data,
                                      active[1:halfway],
                                      lbubs,
                                      leafType,
                                      depth-1,
                                      numbering[1:halfway]) %>%
        bind_rows(
          subtimeline_forest_preexpansion(tree,
                                          data,
                                          active[-(1:halfway)],
                                          lbubs,
                                          leafType,
                                          depth-1,
                                          numbering[-(1:halfway)]) #%>%
            #mutate(treeNumber=(2^depth)-1+treeNumber)
        )
    )
    
    #return(subtimeline_forest_preexpansion(data, tree, active[-1],lbubs, leafType, depth+1) %>%
    #        bind_rows(subtimeline_tree_ordered(data,tree,active[1],lbubs,leafType) %>%
    #           mutate(treeNumber=depth)
    #        )
    #)
  }
} 

subtimeline_forest = function(tree, data, active, lbubs) { # TODO good for map-reduce over roots
  if(length(active)==1)
    return(subtimeline_tree_ordered(tree, data, active, lbubs))
  else
    return(subtimeline_forest(tree, data, active[-1],lbubs) %>%
                product_of_timelines(
                  subtimeline_tree_ordered(tree, data,active[1],lbubs)
                ) # TODO set limits on subtimeline (max and min accumulated logRateFactors)
    )
}

precedingtime = function(e, t, data) { # TODO use Rcpp for this really simple function
  # browser()
  if(is.na(e)) return(NA)
  if(e=="") return(NA)
  suppressWarnings(data %>% filter(event==e & time<=t) %>% .[["time"]] %>% max())
}

# value_sum_onto_durations takes the [subtimeline] and the full expansion of intervals [partitions] and assigns sum of [value]
# to the full expansion as indicated by the subtimeline
cppFunction(
  '
  NumericVector value_sum_onto_durationsCPP(NumericVector indexlb, NumericVector indexub, NumericVector value, int ontoSize) {
    NumericVector onto(ontoSize);
    int nI = indexlb.size();
    for(int i=0; i<nI; ++i) {
      for(int j=indexlb[i]; j<=indexub[i]; ++j)
        if(j>=0 && j<ontoSize)
          onto[j] += value[i];
    }
    return(onto);
  }
  '
)

# apply value to finer grained partition, assumes no overlap in lbub list otherwise overwriting
cppFunction(
  '
  NumericVector value_apply_onto_durationsCPP(NumericVector indexlb, NumericVector indexub, NumericVector value, int ontoSize) {
    NumericVector onto(ontoSize);
    int nI = indexlb.size();
    for(int i=0; i<nI; ++i) {
      for(int j=indexlb[i]; j<=indexub[i]; ++j)
        if(j>=0 && j<ontoSize)
          onto[j] = value[i]; //potentially overwrites value with last value seen if there are overlaps
    }
    return(onto);
  }
  '
)


cppFunction(
  '
  NumericVector value_sum_from_durationsCPP(NumericVector indexlb, NumericVector indexub, NumericVector value, int fromSize) {
    int nI = indexlb.size();
    NumericVector onto(nI);
    for(int i=0; i<nI; ++i) {
      for(int j=indexlb[i]; j<=indexub[i]; ++j)
        onto[i] += value[j];
    }
    return(onto);
  }
  '
)

# returns 0 indexed indices of partitions that are active in subtimeline intervals (lb,ub]
mutate_active_duration_indices = function(subtimeline, partitions) {
  subtimeline %>%
    mutate(llb = findInterval(lb, partitions$lb, left.open=T), uub = findInterval(ub, partitions$ub, left.open=T)) %>%
    mutate(llb = (function(llb,uub) {doublePrecisionProblem = llb == uub + 1;
                                     llb[doublePrecisionProblem] = uub[doublePrecisionProblem]; llb})(.$llb,.$uub)
           )
}

value_sum_onto_durations = function(subtimeline, partitions, value="value") {
    subtimeline %>%
      (function(s,p,v) value_sum_onto_durationsCPP(s$llb, s$uub, s[[value]], nrow(partitions)))(.,partitions,value)
}

value_apply_onto_durations = function(subtimeline, partitions, value="value") {
  subtimeline %>% 
    (function(s,p,v) value_apply_onto_durationsCPP(s$llb, s$uub, s[[value]], nrow(partitions)))(.,partitions,value)
}
value_sum_from_durations = function(subtimeline, partitions, value="value") {
  subtimeline %>%
    (function(s,p,v) value_sum_from_durationsCPP(s$llb, s$uub, partitions[[value]], nrow(partitions)))(.,partitions,value)
}

## Make a vector with counts based on the indices provided from 1 to maxIndex (1-indexed)
cppFunction(
  '
  NumericVector vector_of_counts(NumericVector indices, int maxIndex) {
    NumericVector out(maxIndex);
    int indexSize = indices.size();
    for(int i=0; i<indexSize; ++i)
      if(indices[i] != 0)
        ++out[indices[i]-1];
    return(out);
  }
  '
)

## Expand subtimeline with target event times and indicator M #TODO if else not consistent!
insert_events_into_expansion = function(subtimeline, targetEventTimes) {
  # hits = match(targetEventTimes, subtimeline$ub) #this requires a sort of each, but already sorted
  # if(hits %>% is.na %>% sum == 0)
  #   subtimeline = subtimeline %>% mutate(M = (function(x,h) {o = rep(0,nrow(x)); o[h] = 1; o })(.,hits))
  insert_events_into_expansionCPP(subtimeline$lb, subtimeline$ub, targetEventTimes) %>% tbl_df()
  # TODO fix so (1) events are not double counted and (2) other columns copied across
  # (careful this only return 3 columns not the whole subtimeline df)
}
cppFunction(
  '
  List insert_events_into_expansionCPP(NumericVector lb, NumericVector ub, NumericVector targets) {
    int nlu = lb.size();
    int nt = targets.size();
    IntegerVector targetSpots(nt); // index of lbub to insert above; ignore if outside of lbub
    IntegerVector targetEquality(nt);
    int nEquality = 0;
    int nDuplicate = 0; //precedence over nEquality
    int nOutside = 0;
    int ti=0;
    for(int i=0; i<nlu; ++i) {
      while(ti<nt && targets[ti] <= lb[i]) {
        targetSpots[ti] = -1;
        ++nOutside;
        ++ti;
      }
      while(ti<nt && targets[ti] <= ub[i]) {
        targetSpots[ti] = i;
        if(ti>0 && targets[ti]==targets[ti-1])
          ++nDuplicate;
        else if(targets[ti] == ub[i]) {
          targetEquality[ti] = 1;
          ++nEquality;
        }
        ++ti;
      }
      if(ti==nt)
        break;
    }
    while(ti<nt) {
      targetSpots[ti] = -1;
      ++nOutside;
      ++ti;
    }

    int nout = nlu + nt - nEquality - nOutside - nDuplicate;
    NumericVector lbout(nout), ubout(nout), indicator(nout);
    int tsi = 0;
    int outi = 0;
    for(int i=0; i<nlu; ++i) {
      while(tsi<nt && targetSpots[tsi]==-1)
        ++tsi;
      if(tsi<nt && targetSpots[tsi]==i && targetEquality[tsi]!=1) {
        // insert target time and break [lb,ub] duration
        lbout[outi] = lb[i];
        ubout[outi] = targets[tsi];
        ++indicator[outi];
        while(tsi+1<nt && targets[tsi+1]==targets[tsi]) {
          ++indicator[outi];
          ++tsi;
        }
        ++outi;
        lbout[outi] = targets[tsi];
        ++tsi;
        while(tsi<nt && targetSpots[tsi]==i && targetEquality[tsi]!=1) {
          // insert [targetTime_0,targetTime_1] duration inside [lb,ub]
          ubout[outi] = targets[tsi];
          ++indicator[outi];
          while(tsi+1<nt && targets[tsi+1]==targets[tsi]) {
            ++indicator[outi];
            ++tsi;
          }
          ++outi;
          lbout[outi] = targets[tsi];
          ++tsi;
        }
        if(tsi<nt && targetSpots[tsi]==i && targetEquality[tsi]==1) {
          ++indicator[outi];
          ++tsi;
          while(tsi<nt && targets[tsi]==targets[tsi-1]) {
            ++indicator[outi];
            ++tsi;
          }
        }
        ubout[outi] = ub[i];
        ++outi;
      } else {
        if(tsi < nt && targetSpots[tsi]==i && targetEquality[tsi]==1) {
          ++indicator[outi];
          ++tsi;
          while(tsi<nt && targets[tsi]==targets[tsi-1]) {
            ++indicator[outi];
            ++tsi;
          }
        }
        lbout[outi] = lb[i];
        ubout[outi] = ub[i];
        ++outi;
      }
    }
    return DataFrame::create( Named("lb")=lbout, Named("ub")=ubout, Named("M")=indicator);
    /*List out(8);  // for debugging
    out[0] = lbout;
    out[1] = ubout;
    out[2] = indicator;
    out[3] = targetEquality;
    out[4] = targetSpots;
    out[5] = outi;
    out[6] = nout;
    out[7] = tsi;
    return out;*/
  
  }
  '
)

# TODO code review to clean up these expansion functions


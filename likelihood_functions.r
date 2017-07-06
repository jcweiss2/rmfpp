library(dplyr)

library(Rcpp)
cppFunction(
  '
NumericVector countIn(NumericVector t0, NumericVector t1, NumericVector events) {
/** Returns the number of events in each (t0,t1] ); runs in O(t+e) **/
  int nT = t0.size();
  int nE = events.size();
  int ti = 0;
  NumericVector counts(nT);
  for(int i=0; i<nE; ++i) {
    while( ti < nT  &&  t1[ti] < events[i] )  ++ti;
    if( ti == nT ) break;
    if( t0[ti] < events[i] )  ++counts[ti];
  }
  return counts;
}
  '
)


# Log likelihood computed for event given by timeline
# Inputs:
#   data as tbl_df with event, time
#   event_name as character
#   subtimeline as tbl_df with ts, tend, and logRateFactor
ll = function(subtimeline, data, event_name) {
  (
  subtimeline %>%
    bind_cols( # workaround for mutate(counts=countIn(...))
      data.frame(
        counts = countIn(
          subtimeline$lb,
          subtimeline$ub,
          (data %>% filter(event==event_name) %>% select(time))[[1]]
        )
      )
    ) %>% mutate(
      #divide up subtimelime into periods of number of events and pass through l^M exp(-lT)
      ll = counts*logRateFactor-exp(logRateFactor)*(ub-lb)
    ) %>% summarise(sum(ll))
  )[[1]]
}

# TODO modify this to create countsIn for proposal subsets.
#      you'll need to document time of event for non-exp dists, e.g. log logistic
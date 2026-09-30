# Source-data and strategy components.
repo_root <- Sys.getenv("TTE_REPO_ROOT", unset = normalizePath(".", mustWork = TRUE))
source(file.path(repo_root, "R", "02_outcome_time_core.R"))

clip_union <- function(x,lo,hi) {
  if(!nrow(x) || hi<=lo)return(data.frame(start_h=numeric(),end_h=numeric()))
  as.data.frame(merge_episodes(pmax(x$start_h,lo),pmin(x$end_h,hi)))
}
intersect_union <- function(x,y) {
  if(!nrow(x)||!nrow(y))return(data.frame(start_h=numeric(),end_h=numeric()))
  z<-lapply(seq_len(nrow(y)),function(j)clip_union(x,y$start_h[j],y$end_h[j]))
  z<-do.call(rbind,z)
  as.data.frame(merge_episodes(z$start_h,z$end_h))
}
union_hours <- function(x)sum(x$end_h-x$start_h)
within_intervals <- function(t,x) {
  ans<-rep(FALSE,length(t))
  for(j in seq_len(nrow(x)))ans<-ans | (t>=x$start_h[j] & t<x$end_h[j])
  ans
}
rass_mode <- function(x) {
  if(!length(x))return(NA_real_)
  v<-sort(unique(x));v[which.max(tabulate(match(x,v)))]
}
primary_qualified <- function(s) {
  if(!nrow(s))return(logical())
  ok <- is.finite(s$rate) & (
    (s$drug_class=='propofol' & s$rateuom=='mcg/kg/min' & s$rate>=5) |
    (s$drug_class=='midazolam' & s$rateuom=='mg/hour' & s$rate>=.5) |
    (s$drug_class=='dexmedetomidine' & s$rateuom=='mcg/kg/hour' & s$rate>=.1))
  ok[is.na(ok)]<-FALSE;ok
}
state_from_rass <- function(rass,qualifying,mv_h) {
  if(mv_h==0)return('paused')
  if(is.na(rass))return('missing')
  if(rass>=-2)return('light')
  if(qualifying)'deep' else 'low_without'
}

build_mv_scoped_nodes <- function(mv,rass,sed,icu_out_h=NA_real_,n_nodes=14L) {
  mv<-as.data.frame(mv);rass<-as.data.frame(rass);sed<-as.data.frame(sed)
  stopifnot(all(c('start_h','end_h') %in% names(mv)),
    all(c('time_h','rass_score') %in% names(rass)),
    all(c('start_h','end_h','drug_class','rate','rateuom') %in% names(sed)))
  stopifnot(all(is.finite(c(mv$start_h,mv$end_h))),all(mv$end_h>mv$start_h))
  # Missing ICU outtime is not converted to an imputed discharge. Supplied MV
  # must be explicitly attributed to the original ICU stay. Unknown discharge
  # remains flagged for release review.
  scope_end<-if(is.na(icu_out_h))n_nodes*12 else min(icu_out_h,n_nodes*12)
  mv<-clip_union(mv,0,scope_end)
  sed<-sed[primary_qualified(sed),,drop=FALSE]
  last_mode<-NA_real_;rows<-list()
  for(k in seq.int(0L,n_nodes-1L)) {
    lo<-12*k;hi<-lo+12
    seg<-clip_union(mv,lo,hi);mh<-union_hours(seg)
    rk<-rass[rass$time_h>=lo & rass$time_h<hi,,drop=FALSE]
    rv<-rk[within_intervals(rk$time_h,seg),,drop=FALSE]
    observed<-rass_mode(rv$rass_score)
    if(!is.na(observed))last_mode<-observed
    # Paused nodes never update the carried exposure summary from awake/off-MV
    # measurements. This follows from the prespecified analysis scope.
    carried<-if(mh>0)last_mode else NA_real_
    sk<-sed[sed$start_h<hi & sed$end_h>lo,,drop=FALSE]
    coverage<-union_hours(intersect_union(sk,seg))
    qfixed<-coverage>=6-1e-6
    qfraction<-mh>0 && coverage>=.5*mh-1e-6
    peaks<-vapply(c('propofol','midazolam','dexmedetomidine'),function(drug) {
      ss<-sk[sk$drug_class==drug,,drop=FALSE]
      if(!nrow(ss))return(NA_real_)
      overlap<-vapply(seq_len(nrow(ss)),function(j)
        union_hours(intersect_union(ss[j,,drop=FALSE],seg))>0,logical(1))
      if(any(overlap))max(ss$rate[overlap]) else NA_real_
    },numeric(1))
    rows[[k+1L]]<-data.frame(node_id=k,node_start_h=lo,node_end_h=hi,
      assessment_h=if(is.na(icu_out_h))hi else min(hi,icu_out_h),
      mv_h=mh,coverage_h=coverage,coverage_fraction=if(mh>0)coverage/mh else NA_real_,
      partial_mv_node=mh>0 & mh<12-1e-7,scope_exit_unknown=is.na(icu_out_h),
      rass_observed=observed,rass_carried=carried,rass_count=nrow(rv),
      off_mv_rass_not_used=nrow(rk)-nrow(rv),
      rass_imputed=mh>0 & nrow(rv)==0 & !is.na(carried),
      state_fixed6=state_from_rass(carried,qfixed,mh),
      state_fraction50=state_from_rass(carried,qfraction,mh),
      max_propofol_rate=peaks[1],max_midazolam_rate=peaks[2],max_dexmedetomidine_rate=peaks[3])
  }
  ans<-do.call(rbind,rows);rownames(ans)<-NULL
  ans$partial_policy_disagreement<-ans$state_fixed6!=ans$state_fraction50
  attr(ans,'clinical_effect_ready')<-FALSE;ans
}

pause_resume_trace <- function(nodes,strategy=c('SLA','DTLS'),transition_h=48,
                              horizon_h=144,endpoint_h=Inf,icu_out_h=NA_real_,
                              partial_policy=NULL,rass_only=FALSE) {
  strategy<-match.arg(strategy)
  if(!is.null(partial_policy))partial_policy<-match.arg(partial_policy,c('fixed6','fraction50'))
  d<-as.data.frame(nodes);d<-d[order(d$node_id),,drop=FALSE]
  stopifnot(all(c('state_fixed6','state_fraction50','assessment_h','mv_h') %in% names(d)))
  counter<-0L;censor_h<-Inf;records<-list();reason<-'fixed_horizon'
  for(i in seq_len(nrow(d))) {
    a<-d[i,]
    if(a$node_start_h>=horizon_h)break
    if(!is.na(icu_out_h) && a$node_start_h>=icu_out_h) {reason<-'index_icu_ended';break}
    t<-min(a$assessment_h,horizon_h)
    if(endpoint_h<=t) {reason<-'terminal_event_first';break}
    if(a$mv_h==0)state<-'paused' else if(rass_only) {
      state<-if(is.na(a$rass_carried))'missing' else if(a$rass_carried>=-2)'light' else 'deep'
    } else if(is.null(partial_policy) && a$partial_policy_disagreement) {
      reason<-'unresolved_partial_coverage_definition';break
    } else state<-if(is.null(partial_policy)||partial_policy=='fixed6')a$state_fixed6 else a$state_fraction50
    expected<-if(strategy=='DTLS' && a$node_start_h<transition_h)'deep' else 'light'
    before<-counter
    # Pause and unresolved RASS are not evidence of adherence. Do not reset the
    # counter or restart the transition/horizon clock. One evaluated GLOBAL
    # node increments at most once, regardless of its number of MV episodes.
    if(!state %in% c('paused','missing'))counter<-if(state==expected)0L else counter+1L
    trigger<-state=='low_without' || counter>1L
    records[[length(records)+1L]]<-data.frame(node_id=a$node_id,assessment_h=t,
      state=state,expected=expected,counter_before=before,counter_after=counter,
      censor_now=trigger,paused=state=='paused',partial_mv_node=a$partial_mv_node)
    if(trigger) {censor_h<-t;reason<-'artificial_censor';break}
  }
  empty<-data.frame(node_id=integer(),assessment_h=numeric(),state=character(),expected=character(),
    counter_before=integer(),counter_after=integer(),censor_now=logical(),paused=logical(),partial_mv_node=logical())
  list(rows=if(length(records))do.call(rbind,records) else empty,censor_h=censor_h,
       stop_reason=reason,partial_policy=partial_policy,clinical_effect_ready=FALSE)
}

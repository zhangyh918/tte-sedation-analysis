# Sensitivity analyses for the adherence horizon, exposure definition,
# RASS aggregation and temporal split.
repo_root <- Sys.getenv("TTE_REPO_ROOT", unset = normalizePath(".", mustWork = TRUE))
source(file.path(repo_root, "R", "05_complete_history_core.R"))

sync_horizon_history <- function(nodes,strategy='SLA',transition_h=48,rass_only=FALSE,
  target='vfd',death_h=Inf,ext_h=Inf,vfd=28,icu_out_h=NA_real_,horizon_h=144) {
  stopifnot(horizon_h %in% seq(12,168,12),target %in% c('vfd','ext'))
  endpoint<-min(672,death_h,if(target=='ext')ext_h else Inf)
  tr<-pause_resume_trace(nodes,strategy,transition_h=transition_h,horizon_h=horizon_h,
    endpoint_h=endpoint,icu_out_h=icu_out_h,partial_policy='fraction50',rass_only=rass_only)
  stopifnot(tr$stop_reason!='unresolved_partial_coverage_definition')
  h<-minimal_event_history(death_h,ext_h,tr$censor_h,target,observed_vfd=vfd,
    censor_tie='endpoint_first',competing_tie='death')
  rr<-tr$rows;rr$model_eligible<-!rr$state %in% c('paused','missing')
  rr$uncensored<-as.integer(!rr$censor_now)
  extra<-setdiff(names(nodes),names(rr))
  rr<-merge(rr,as.data.frame(nodes)[,c('node_id',extra),drop=FALSE],by='node_id',sort=TRUE)
  decisions<-rr[rr$model_eligible,,drop=FALSE];retained<-decisions[decisions$uncensored==1,,drop=FALSE]
  updates<-data.frame(available_h=retained$assessment_h,information_end_h=retained$assessment_h,
    weight=rep(1,nrow(retained)))
  risk<-minimal_weighted_intervals(h,updates);risk$weight<-NULL
  last_h<-if(nrow(retained))max(retained$assessment_h) else 0
  status<-data.frame(strategy=strategy,target=target,censor_h=tr$censor_h,
    endpoint_h=endpoint,finish_h=h$risk_followup_stop,endpoint_observed=h$endpoint_observed,
    endpoint_type=if(is.null(h$endpoint_type))NA_character_ else h$endpoint_type,
    last_retained_decision_h=last_h,stop_reason=tr$stop_reason,n_trace_rows=nrow(rr),
    n_model_decisions=nrow(decisions),n_pauses=sum(rr$paused),n_partial=sum(rr$partial_mv_node),
    scope_exit_unknown=is.na(icu_out_h),observed_vfd=if(target=='vfd')vfd else NA_real_)
  stopifnot(all(rr$assessment_h<endpoint),all(rr$assessment_h<=horizon_h),
    all(decisions$lag_information_end_h<=decisions$node_start_h+1e-7),
    all(risk$stop>risk$start),abs(sum(risk$stop-risk$start)-h$risk_followup_stop)<1e-7,
    !h$endpoint_observed || last_h<endpoint)
  list(trace=rr,decisions=decisions,status=status,risk=risk,clinical_effect_ready=FALSE)
}

sync_build_histories <- function(b,nodes,horizon,scenario) {
  nn<-split(as.data.frame(nodes),nodes$stay_id)
  all<-lapply(seq_len(nrow(b)),function(i){
    rr<-list();j<-0L
    for(s in c('SLA','DTLS'))for(t in c('vfd','ext')) {
      j<-j+1L
      z<-sync_horizon_history(nn[[as.character(b$stay_id[i])]],strategy=s,
        target=t,death_h=b$death_h[i],ext_h=b$ext_h[i],vfd=b$vfd[i],
        icu_out_h=b$icu_out_h[i],horizon_h=horizon)
      for(k in c('trace','decisions','status','risk')) {
        z[[k]]$stay_id<-rep(b$stay_id[i],nrow(z[[k]]))
        z[[k]]$strategy<-rep(s,nrow(z[[k]]))
        z[[k]]$target<-rep(t,nrow(z[[k]]))
      }
      rr[[j]]<-z
    }
    if(i%%1000L==0L){cat(scenario,'histories',i,'/',nrow(b),'\n');flush.console()}
    rr
  })
  flat<-unlist(all,recursive=FALSE)
  out<-list(baseline=data.table::as.data.table(b),
    scenario=data.table::data.table(scenario=scenario,transition=48,rass_only=FALSE),
    horizon_h=horizon,partial_policy='fraction50',clinical_effect_ready=FALSE)
  for(k in c('trace','decisions','status','risk'))
    out[[k]]<-data.table::rbindlist(lapply(flat,`[[`,k),fill=TRUE)
  stopifnot(nrow(out$status)==4*nrow(b),!anyDuplicated(out$baseline$stay_id),
    all(out$decisions$assessment_h<=horizon),all(out$risk$stop>out$risk$start))
  out
}

sync_variant_builder <- function(dose='primary',aggregation='modal') {
  stopifnot(dose %in% c('primary','any_positive','half_primary','double_primary'),
    aggregation %in% c('modal','arithmetic_mean','minimum'))
  env<-new.env(parent=environment(build_mv_scoped_nodes))
  env$primary_qualified<-local({
    selected<-dose
    function(s) {
      if(!nrow(s))return(logical())
      multiplier<-switch(selected,primary=1,half_primary=.5,double_primary=2,any_positive=0)
      above<-function(rate,cut)if(selected=='any_positive')rate>0 else rate>=multiplier*cut
      ok<-is.finite(s$rate)&(
        (s$drug_class=='propofol' & s$rateuom=='mcg/kg/min' & above(s$rate,5)) |
        (s$drug_class=='midazolam' & s$rateuom=='mg/hour' & above(s$rate,.5)) |
        (s$drug_class=='dexmedetomidine' & s$rateuom=='mcg/kg/hour' & above(s$rate,.1)))
      ok[is.na(ok)]<-FALSE;ok
    }
  })
  env$rass_mode<-switch(aggregation,modal=rass_mode,
    arithmetic_mean=function(x){x<-x[!is.na(x)];if(length(x))mean(x) else NA_real_},
    minimum=function(x){x<-x[!is.na(x)];if(length(x))min(x) else NA_real_})
  fn<-build_mv_scoped_nodes;environment(fn)<-env;fn
}

# No horizon is selected from results. These summaries distinguish remaining
# uncensored copies (including terminal outcomes) from still-active decisions.
sync_history_support <- function(h) {
  data.table::rbindlist(lapply(seq(12,h$horizon_h,12),function(hour)
    h$status[,.(initial_copies=.N,
      not_artificially_censored=sum(censor_h>hour),
      uncensored_fraction=mean(censor_h>hour),
      endpoint_already_observed=sum(endpoint_observed & endpoint_h<=hour),
      still_in_followup=sum(finish_h>hour)),by=.(target,strategy)][,hour:=hour]))
}

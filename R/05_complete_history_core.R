# Strategy histories using the prespecified 50% node-coverage rule.
repo_root <- Sys.getenv("TTE_REPO_ROOT", unset = normalizePath(".", mustWork = TRUE))
source(file.path(repo_root, "R", "04_mv_pause_resume_core.R"))
source(file.path(repo_root, "R", "03_event_history_core.R"))

lag_mv_history <- function(nodes,initial_rass=NA_real_,initial_rates=c(0,0,0)) {
  d<-as.data.frame(nodes);d<-d[order(d$node_id),,drop=FALSE]
  stopifnot(length(initial_rates)==3,all(is.finite(initial_rates)),all(initial_rates>=0),
            !anyDuplicated(d$node_id),all(d$node_start_h<d$node_end_h))
  r<-initial_rass;rates<-initial_rates;info_end<-0
  for(i in seq_len(nrow(d))) {
    d$lag_rass[i]<-r;d$lag_propofol[i]<-rates[1]
    d$lag_midazolam[i]<-rates[2];d$lag_dex[i]<-rates[3]
    d$lag_information_end_h[i]<-info_end
    if(d$mv_h[i]>0) {
      if(!is.na(d$rass_carried[i]))r<-d$rass_carried[i]
      rates<-unlist(d[i,c('max_propofol_rate','max_midazolam_rate','max_dexmedetomidine_rate')])
      # These fields mean recorded qualifying infusion, not all sedative burden.
      # No overlapping qualifying record is 0, not an unknown positive dose.
      rates[is.na(rates)]<-0
      info_end<-d$assessment_h[i]
    }
  }
  stopifnot(all(d$lag_information_end_h<=d$node_start_h+1e-7))
  d
}

complete_clone_history <- function(nodes,strategy='SLA',transition_h=48,
                                   rass_only=FALSE,target='vfd',death_h=Inf,
                                   ext_h=Inf,vfd=28,icu_out_h=NA_real_) {
  stopifnot(target %in% c('vfd','ext'),transition_h %in% c(24,48,72))
  endpoint<-min(672,death_h,if(target=='ext')ext_h else Inf)
  tr<-pause_resume_trace(nodes,strategy,transition_h=transition_h,horizon_h=144,
    endpoint_h=endpoint,icu_out_h=icu_out_h,partial_policy='fraction50',rass_only=rass_only)
  stopifnot(tr$stop_reason!='unresolved_partial_coverage_definition')
  h<-minimal_event_history(death_h,ext_h,tr$censor_h,target,observed_vfd=vfd,
    censor_tie='endpoint_first',competing_tie='death')
  rr<-tr$rows
  # All rows include pauses for audit. Only evaluable decisions enter the
  # censoring likelihood. Missing exposure is not fabricated as adherence.
  rr$model_eligible<-!rr$state %in% c('paused','missing')
  rr$uncensored<-as.integer(!rr$censor_now)
  extra<-setdiff(names(nodes),names(rr))
  rr<-merge(rr,as.data.frame(nodes)[,c('node_id',extra),drop=FALSE],by='node_id',sort=TRUE)
  decisions<-rr[rr$model_eligible,,drop=FALSE]
  retained<-decisions[decisions$uncensored==1,,drop=FALSE]
  # Unit weights here are a placeholder for chronology checks, NOT a clinical
  # unweighted analysis. Future fitting must explicitly replace them.
  updates<-data.frame(available_h=retained$assessment_h,
    information_end_h=retained$assessment_h,weight=rep(1,nrow(retained)))
  risk<-minimal_weighted_intervals(h,updates)
  risk$weight<-NULL
  last_h<-if(nrow(retained))max(retained$assessment_h) else 0
  status<-data.frame(strategy=strategy,target=target,censor_h=tr$censor_h,
    endpoint_h=endpoint,finish_h=h$risk_followup_stop,endpoint_observed=h$endpoint_observed,
    endpoint_type=if(is.null(h$endpoint_type))NA_character_ else h$endpoint_type,
    last_retained_decision_h=last_h,stop_reason=tr$stop_reason,
    n_trace_rows=nrow(rr),n_model_decisions=nrow(decisions),
    n_pauses=sum(rr$paused),n_partial=sum(rr$partial_mv_node),
    scope_exit_unknown=is.na(icu_out_h),observed_vfd=if(target=='vfd')vfd else NA_real_)
  stopifnot(all(rr$assessment_h<endpoint),all(rr$assessment_h<=144),
    all(decisions$lag_information_end_h<=decisions$node_start_h+1e-7),
    all(risk$stop>risk$start),abs(sum(risk$stop-risk$start)-h$risk_followup_stop)<1e-7,
    !h$endpoint_observed || last_h<endpoint)
  list(trace=rr,decisions=decisions,status=status,risk=risk,clinical_effect_ready=FALSE)
}

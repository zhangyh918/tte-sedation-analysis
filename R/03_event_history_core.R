# Event-history components used by the analysis pipeline.
# This module reuses the interval builder without invoking its exposure rules.
repo_root <- Sys.getenv("TTE_REPO_ROOT", unset = normalizePath(".", mustWork = TRUE))
source(file.path(repo_root, "R", "01_event_order_core.R"))

minimal_event_history <- function(death_h, successful_ext_h, censor_h,
                                  target=c('ext','vfd'), observed_vfd,
                                  censor_tie=NULL, competing_tie=NULL) {
  target <- match.arg(target)
  times <- c(death_h,successful_ext_h,censor_h)
  stopifnot(length(death_h)==1L,length(successful_ext_h)==1L,length(censor_h)==1L,
            !anyNA(times),all(times>=0),death_h>0,successful_ext_h>0)
  if(target=='ext') {
    end <- min(death_h,successful_ext_h,672)
    type <- if(end==672 && min(death_h,successful_ext_h)>672) 'day28' else
      if(death_h<successful_ext_h) 'death' else 'successful_extubation'
    if(is.finite(death_h) && death_h==successful_ext_h && death_h<=672) {
      if(is.null(competing_tie))stop('Death/extubation tie needs an explicit rule')
      type <- match.arg(competing_tie,c('death','successful_extubation'))
    }
  } else {
    stopifnot(length(observed_vfd)==1L,is.finite(observed_vfd),observed_vfd>=0,observed_vfd<=28)
    if(death_h<=672 && observed_vfd!=0)stop('Recorded death by day 28 must imply VFD=0')
    end <- min(death_h,672)
    type <- if(death_h<=672)'death' else 'day28'
  }
  kept <- end<censor_h
  if(end==censor_h) {
    if(is.null(censor_tie))stop('Endpoint/censor tie needs an explicit rule')
    kept <- match.arg(censor_tie,c('endpoint_first','censor_first'))=='endpoint_first'
  }
  list(endpoint_observed=kept,endpoint_time=if(kept)end else NA_real_,
       endpoint_type=if(kept)type else NA_character_,risk_followup_stop=min(end,censor_h),
       natural_endpoint_h=end,censor_time=censor_h,target=target,
       observed_vfd=if(target=='vfd')observed_vfd else NA_real_,
       decisions=data.frame(decision_time=numeric()))
}

check_predictable_updates <- function(updates) {
  needed <- c('available_h','information_end_h','weight')
  stopifnot(all(needed %in% names(updates)))
  w <- as.data.frame(updates)[,needed,drop=FALSE]
  stopifnot(all(is.finite(as.matrix(w))),all(w$available_h>=0),
            all(w$information_end_h>=0),all(w$weight>0),!anyDuplicated(w$available_h))
  if(any(w$information_end_h>w$available_h))
    stop('Future information cannot be used before it becomes available')
  w[order(w$available_h),,drop=FALSE]
}

minimal_weighted_intervals <- function(history,updates) {
  w <- check_predictable_updates(updates)
  h <- history
  h$decisions <- data.frame(decision_time=w$available_h)
  r <- make_risk_intervals(h,followup_hours=672)
  # The old pure builder marks any terminal endpoint, including day 28.
  # Separate actual events from administrative completion here.
  r$extubation <- r$event*as.integer(isTRUE(h$endpoint_type=='successful_extubation'))
  r$death <- r$event*as.integer(isTRUE(h$endpoint_type=='death'))
  r$event <- NULL
  r$weight <- vapply(r$start,function(t) {
    j <- which(w$available_h<=t)
    if(length(j))w$weight[max(j)] else 1
  },numeric(1))
  attr(r,'clinical_effect_ready') <- FALSE
  r
}

minimal_terminal_row <- function(history,updates) {
  stopifnot(history$target=='vfd')
  w <- check_predictable_updates(updates)
  if(!history$endpoint_observed)
    return(data.frame(endpoint_h=numeric(),vfd=numeric(),weight=numeric()))
  # Uses the left-limit weight: no decision made at/after the outcome is used.
  eligible <- which(w$available_h<history$endpoint_time)
  weight <- if(length(eligible))w$weight[max(eligible)] else 1
  ans <- data.frame(endpoint_h=history$endpoint_time,vfd=history$observed_vfd,weight=weight)
  attr(ans,'clinical_effect_ready') <- FALSE
  ans
}

# A node-level censor trigger is not an exact timestamp. This classifier
# reports only chronology relative to that recorded node; it does not choose
# a start/end censoring convention or evaluate hypothetical new censor rules.
classify_trigger_window <- function(endpoint_h,node_start_h,node_end_h) {
  stopifnot(length(endpoint_h)==length(node_start_h),length(endpoint_h)==length(node_end_h),
            !anyNA(c(endpoint_h,node_start_h,node_end_h)),all(is.finite(endpoint_h)),
            all(endpoint_h>0),all(node_start_h>=0),all(node_end_h>=node_start_h),
            all(is.infinite(node_start_h)==is.infinite(node_end_h)))
  ans <- rep('within_or_at_trigger_node_unresolved',length(endpoint_h))
  ans[endpoint_h<node_start_h] <- 'endpoint_before_trigger_node'
  ans[endpoint_h>node_end_h] <- 'trigger_node_before_endpoint'
  ans[is.infinite(node_start_h)] <- 'no_recorded_artificial_censor'
  ans
}

assert_minimal_contract_complete <- function(contract) {
  required <- c('decision_time_rule','partial_node_rule','post_extubation_strategy',
                'weight_outcome_pair','history_specification','main_cohort_role')
  unresolved <- required[vapply(required,function(k) {
    val <- contract[[k]]
    is.null(val) || length(val)!=1L || !is.character(val) || is.na(val) ||
      !nzchar(trimws(val)) || startsWith(val,'PENDING')
  },logical(1))]
  if(length(unresolved))stop('Effect analysis not released: ',paste(unresolved,collapse=', '))
  invisible(TRUE)
}

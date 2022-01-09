library(survival)
library(ranger)
library(distr6)
library(survivalsvm)
library(mlr3proba)

set.seed(20220109)

## Get data and variables
data = survival::rats
train = sample(nrow(data), nrow(data) * 2/3)
test = setdiff(seq(nrow(data)), train)
test_unique_times = data$time[test]
target = Surv(test_unique_times, data$status[test])
target_train = Surv(data$time[train], data$status[train])

## Train and predict
cox = coxph(Surv(time, status) ~ ., data = data[train,])
p_cox_lp = predict(cox, newdata = data[test,])
p_cox_surv = survfit(cox, newdata = data[test,])
p_ranger = predict(
  ranger(Surv(time, status) ~ ., data = data[train, ]),
  data = data[test, ]
)
p_svm = predict(
  survivalsvm(Surv(time, status) ~ ., data = data[train, ], gamma.mu = 0.1),
  newdata = data[test, ]
)

# Define Antolini's C-index
#  Copied from pycox
#   https://github.com/havakv/pycox/blob/master/pycox/evaluation/concordance.py
is_comparable = function(t_i, t_j, d_i, d_j) 
  ((t_i < t_j) & d_i) | ((t_i == t_j) & (d_i | d_j))

is_concordant = function(s_i, s_j, t_i, t_j, d_i, d_j)
  (s_i < s_j) & is_comparable(t_i, t_j, d_i, d_j)

sum_comparable = function(t, d) {
    count = 0
    for (i in seq_along(t)) {
      for (j in seq_along(t)) {
        if (j != i) {
          count = count + is_comparable(t[i], t[j], d[i], d[j])
        }
      }
    }
    count
}

sum_concordant_disc = function(s, t, d, s_idx) {
  count = 0
  for (i in seq_along(t)) {
    idx = s_idx[i]
    for (j in seq_along(t)) {
      if (j != i) {
        count = count +
          is_concordant(s[idx, i], s[idx, j], t[i], t[j], d[i], d[j])
      }
    }
  }
  count
}

# truth - Surv object corresponding to true survival outcomes
# surv - predicted survival matrix (T x n)
# surv_idx - 'surv_idx[i]' gives index in 'surv' corresponding to
#   the event time of individual 'i'.
antolini = function (truth, surv, surv_idx) {
  durations = truth[, "time"]
  events = truth[, "status"]
  sum_concordant_disc(surv, durations, events, surv_idx) /
    sum_comparable(durations, events)
}

## Calculative Cindex for CPH and SVM
# Harrell
harrell_cph = concordance(target ~ p_cox_lp, reverse = TRUE)$concordance
harrell_rsf = NA
harrell_svm = concordance(target ~ p_svm$predicted[1,])$concordance

# Uno
uno_cph = concordance(target ~ p_cox_lp, reverse = TRUE, timewt = "n/G2")$concordance
uno_rsf = NA
uno_svm = concordance(target ~ p_svm$predicted[1,], timewt = "n/G2")$concordance


## Method 1 - Ensemble mortality - higher value = more deaths = higher risk
ensemble_rsf = concordance(target ~ rowSums(p_ranger$chf), reverse = TRUE)$concordance
ensemble_cph = concordance(target ~ rowSums(-log(t(p_cox_surv$surv))), reverse = TRUE)$concordance
ensemble_svm = NA

## Method 2 - Antolini
antolini_rsf = antolini(
  target, t(p_ranger$survival),
  findInterval(target[, "time"], p_ranger$unique.death.times, all.inside = TRUE)
)
antolini_cph = antolini(
  target, p_cox_surv$surv,
  findInterval(target[, "time"], p_cox_surv$time, all.inside = TRUE)
)
antolini_svm = NA


## Method 3 - Distribution summary (no extrapolation)
cox_surv = t(p_cox_surv$surv)
colnames(cox_surv) = p_cox_surv$time
ranger_surv = p_ranger$survival
colnames(ranger_surv) = p_ranger$unique.death.times

# Higher value = longer expected lifetime = lower risk - Absurd value due to improper distribution
summary_naive_cph = concordance(target ~ distr6::as.Distribution(1 - cox_surv, fun = "cdf")$mean())$concordance
summary_naive_rsf = concordance(target ~ distr6::as.Distribution(1 - ranger_surv, fun = "cdf")$mean())$concordance
summary_naive_svm = NA

## Method 4 - Distribution summary (extrapolation)
cox_surv = cbind(1, cox_surv, 0) # Add probabilities 1 and 0
colnames(cox_surv)[1] = "0"
colnames(cox_surv)[ncol(cox_surv)] = tail(p_cox_surv$time, 1) + 1e-3
ranger_surv = cbind(1, ranger_surv, 0) # Add probabilities 1 and 0
colnames(ranger_surv)[1] = "0"
colnames(ranger_surv)[ncol(ranger_surv)] = tail(p_ranger$unique.death.times, 1) + 1e-3

summary_extr_cph = concordance(target ~ distr6::as.Distribution(1 - cox_surv, fun = "cdf")$mean())$concordance
summary_extr_rsf = concordance(target ~ distr6::as.Distribution(1 - ranger_surv, fun = "cdf")$mean())$concordance
summary_extr_svm = NA

## Method 5 - Single probability comparison
distr_cox = distr6::as.Distribution(1 - cox_surv, fun = "cdf")
distr_rsf = distr6::as.Distribution(1 - ranger_surv, fun = "cdf")

# survival - higher value = higher prob survival = lower risk
cox_prob_concordance = rsf_prob_concordance = numeric(nrow(p_cox_surv$surv))
for (i in seq_along(cox_prob_concordance)) {
  cox_prob_concordance[i] = concordance(target ~ p_cox_surv$surv[i, ])$concordance
}
for (i in seq_along(rsf_prob_concordance)) {
  rsf_prob_concordance[i] = concordance(target ~ p_ranger$survival[, i])$concordance
}

rsf_max = which.max(rsf_prob_concordance)
rsf_min = which.min(rsf_prob_concordance)
rsf_rand = sample(seq_along(rsf_prob_concordance), 1)

prob_cph = cox_prob_concordance[c(rsf_min, rsf_max, rsf_rand)]
prob_rsf = rsf_prob_concordance[c(rsf_min, rsf_max, rsf_rand)]
prob_svm = rep(NA, 3)

stargazer::stargazer(matrix(c(
  harrell_cph, harrell_rsf, harrell_svm,
  uno_cph, uno_rsf, uno_svm,
  antolini_cph, antolini_rsf, antolini_svm,
  prob_cph[1], prob_rsf[1], prob_svm[1],
  prob_cph[2], prob_rsf[2], prob_svm[2],
  prob_cph[3], prob_rsf[3], prob_svm[3],
  summary_naive_cph, summary_naive_rsf, summary_naive_svm,
  summary_extr_cph, summary_extr_rsf, summary_extr_svm,
  ensemble_cph, ensemble_rsf, ensemble_svm
),
ncol = 3, byrow = TRUE,
dimnames = list(c(
  "Harrell", "Uno", "Antolini", "Prob (min)", "Prob (max)",
  "Prob (rand)", "Summary (naive)", "Summary (extr)", "Ensemble"
), c("CPH", "RSF", "SVM"))
), summary = FALSE)

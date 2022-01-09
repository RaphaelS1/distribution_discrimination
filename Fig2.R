library(ggplot2)
library(survival)
library(dplyr)
library(patchwork)

fit <- survfit(Surv(time, status) ~ 1, rats)
dat <- data.frame(S = fit$surv, T = fit$time)

p_orig <- ggplot(dat, aes(x = T, y = S)) +
  geom_line()

dat2 <- rbind(dat, data.frame(T = 105, S = 0))
(mean <- round(as.Distribution(matrix(1 - dat2$S, 1, 57, FALSE, list(NULL, dat2$T)), fun = "cdf")$mean()))
p_drop <- dat2 %>%
  ggplot(aes(x = T, y = S)) +
  geom_line() +
  geom_text(aes(x = x, y = y), data.frame(x = 95, y = 0.6), label = "m=105") +
  annotate("text", x = mean - 5, y = 0.05, label = expression(mu ~ "= 101")) +
  geom_vline(xintercept = 104, lty = 3)

slope <- -diff(range(fit$surv)) / diff(range(fit$time))
intercept <- fit$surv[1] - slope * fit$time[1]
median <- round((0.5 - intercept) / slope)
new_time <- round(unique(-intercept / slope), 2)
dat2 <- rbind(dat, data.frame(T = new_time, S = 0))
(mean <- round(as.Distribution(matrix(1 - dat2$S, 1, 57, FALSE, list(NULL, dat2$T)), fun = "cdf")$mean()))

p_lin <- dat2 %>%
  ggplot(aes(x = T, y = S)) +
  geom_line() +
  geom_text(aes(x = x, y = y), data.frame(x = median + 20, y = 0.6), label = sprintf("m=%s", median)) +
  annotate("text", x = mean - 20, y = 0.05, label = expression(mu ~ "= 386")) +
  geom_vline(xintercept = 104, lty = 3)

jpeg("fig2.jpeg", 10, 13, "cm", res = 300)
p_orig + p_drop + p_lin +
  plot_layout(1, 3) &
  theme_classic() &
  labs(y = "S(T)") &
  ylim(0, 1) &
  geom_hline(yintercept = 0.5, lty = 2)
dev.off()

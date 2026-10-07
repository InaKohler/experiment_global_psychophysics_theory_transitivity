#' ---
#' title: Sample size calculation for magnitude productions
#' author: Dorina Kohler
#' output: html_document
#' ---
#'
#' Sample size refers to number of trials, i.e. number of magnitude productions.
#' The confidence interval for the mean of the adjusted intensities is taken as
#' the precision of the adjustments. The maximum width of the confidence
#' interval is fixed on the basis of  intra-modal JND. SD is taken from earlier
#' data.
#' 1 Block takes approx. 6 minutes
#'
#' ## Set the maximal width to twice the JND in the respective modality
#'
#' 2 * 1.1 dB SPL for loudness
max_width_loud <- 2 * 1.1

#' 2 * 0.7 dB Lambert for brightness ----
max_width_bright <- 2 * 0.7

#' 2 * 2.2 dB displacement for vibration strength ----
max_width_strong <- 2 * 2.2

#' ## SD estimate from old data
#' - brightness (std = 34): 3.8
#' - loudness (std = 76): 7.6
#' - vibration strength (std = 56): 5.2
#' (no data for p = 1, 5.2 is aggregate over p = 2,3,6)

#' ## One simulation for loudness and brightness
# different SD for modalities
stds <- c(76,34,56)
s <- c(3, 4.3, 5.5)
n_sel <- 75
n <- c(45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
cond <- expand.grid(n = n, s = s)
nrep <- 1000

# simulate 95% CI
sim_ci <- function(n, mu, s) {
  x <- rnorm(n, mean = mu, sd = s)
  t.test(x, conf.level = 0.95)$conf.int
}

# get .2 and .8 quantile of width of the 95% CI of simulated magnitude
# productions
width_quant_fun <- function(n, mu = 0, s, nrep = 10000) {
  width <- replicate(nrep, {
    sim_ci(n, mu, s) |> diff()
  })
  setNames(quantile(width, prob = c(.2, .8)), c("qw20", "qw80"))
}

cond$qwidth <- mapply(width_quant_fun,
  n = cond$n, s = cond$s,
  MoreArgs = list(nrep = nrep)
) |> t()

cols <- c("#164467", "#A51E37","#9B7EA6")
# dev.off()
pdf("output/sample_size_calculation_confidence_interval.pdf")
plot(qwidth[, 1] ~ n, cond,
  ylim = c(0, max(5.5, cond$qwidth)),
  type = "n", ylab = ".2 and .8 quantile of simulated CI width",
  log = "x", axes = FALSE
)
axis(1, at = unique(cond$n))
axis(2)
box()
arrows(
  x0 = cond$n, y0 = cond$qwidth[, 1], y1 = cond$qwidth[, 2],
  length = 0.01, code = 3, angle = 90, lty = 1,
  col = rep(cols, each = nrow(cond) / length(s))
)
# match lines to sd values and their modality labels
sd_labels <- c(glue::glue("{s[1]} (luminance)"), 
               glue::glue("{s[2]} (sound pressure)"), 
               glue::glue("{s[3]} (displacement)"))
line_cols  <- cols[1:3]  # first 3 colors correspond to s = c(2, 5, 6)

abline(h = max_width_loud,   lty = 3, col = line_cols[2])
text(49, max_width_loud   + .2, "sound pressure", col = line_cols[2])
abline(h = max_width_bright, lty = 3, col = line_cols[1])
text(47.5, max_width_bright + .2, "luminance",      col = line_cols[1])
abline(h = max_width_strong, lty = 3, col = line_cols[3])
text(48, max_width_strong + .2, "displacement",   col = line_cols[3])

legend("topright",
       legend = sd_labels, col = line_cols, lty = 1,
       title = "sd"
)
dev.off()


#' # Examples
#' ## Sample size for loudness -----------------------------------------------
#' 1. analytic solution
width_l <- max_width_loud
sig_l <- s[2]
(n_l1 <- (2 * qnorm(0.975, 0, 1) * sig_l / width_l)^2 |> ceiling())

#' 2. simulation
nrep <- 2000
ci1 <- replicate(nrep, sim_ci(n = n_l1, mu = stds[2], s = sig_l))
(pow1 <- mean(width_l >= diff(ci1)) |> round(2))

n_l2 <- n_sel
ci2 <- replicate(nrep, sim_ci(n = n_l2, mu = stds[2], s = sig_l))
(pow2 <- mean(width_l >= diff(ci2)) |> round(2))

ci1 <- ci1[, sample(1:nrep, 300)]
ci2 <- ci2[, sample(1:nrep, 300)]
cols <- rep("grey", nrep)
cols[which(diff(ci1) <= width_l)] <- "#A51E37"
par(mfrow = c(1, 2))
matplot(ci1, rbind(1:300, 1:300),
  type = "l", lty = 1,
  xlab = "Coverage area", ylab = "Confidence intervals",
  main = paste("n =", n_l1, ", small CI: ", pow1 * 100, "%"), col = cols,
  xlim = c(stds[2]-3, stds[2]+3)
)
abline(v = stds[2], lty = 1)
cols <- rep("grey", nrep)
cols[which(diff(ci2) <= width_l)] <- "#A51E37"
matplot(ci2, rbind(1:300, 1:300),
  type = "l", lty = 1,
  xlab = "Coverage area", ylab = "Confidence intervals",
  main = paste("n =", n_l2, ", small CI: ", pow2 * 100, "%"), col = cols,
  xlim = c(stds[2]-3, stds[2]+3)
)
abline(v = stds[2], lty = 1)



#' ## Sample size for brightness ----------------------------------------------
width_b <- max_width_bright
sig_b <- s[1]

#' - analytic solution
(n_b1 <- (2 * qnorm(0.975, 0, 1) * sig_b / width_b)^2 |> ceiling())

#' - simulation
nrep <- 2000
ci1 <- replicate(nrep, sim_ci(n_b1, mu = stds[1], s = sig_b))
(pow1 <- mean(width_b >= diff(ci1)))

n_b2 <- n_sel
ci2 <- replicate(nrep, sim_ci(n_b2, mu = stds[1], s = sig_b))
(pow2 <- mean(width_b >= diff(ci2)))

ci1 <- ci1[, sample(1:nrep, 300)]
ci2 <- ci2[, sample(1:nrep, 300)]
cols <- rep("grey", nrep)
cols[which(diff(ci1) <= width_b)] <- "#164467"
par(mfrow = c(1, 2))
matplot(ci1, rbind(1:300, 1:300),
  type = "l", lty = 1,
  xlab = "Coverage area", ylab = "Confidence intervals",
  main = paste("n =", n_b1, ", small CI: ", pow1 * 100, "%"), col = cols,
  xlim =  c(stds[1]-3, stds[1]+3)
)
abline(v = stds[1], lty = 2)
cols <- rep("grey", nrep)
cols[which(diff(ci2) <= width_b)] <- "#164467"
matplot(ci2, rbind(1:300, 1:300),
  type = "l", lty = 1,
  xlab = "Coverage area", ylab = "Confidence intervals",
  main = paste("n =", n_b2, ", small CI: ", pow2 * 100, "%"), col = cols,
  xlim = c(stds[1]-3, stds[1]+3)
)
abline(v = stds[1], lty = 2)


#' ## Sample size for vibration strength -------------------------------------
width_s <- max_width_strong
sig_s <- s[3]

#' - analytic solution
(n_s1 <- (2 * qnorm(0.975, 0, 1) * sig_s / width_s)^2 |> ceiling())

#' - simulation
nrep <- 2000
ci1 <- replicate(nrep, sim_ci(n_s1, mu = stds[3], s = sig_s))
(pow1 <- mean(width_s >= diff(ci1)))

n_s2 <- n_sel
ci2 <- replicate(nrep, sim_ci(n_s2, mu = stds[3], s = sig_s))
(pow2 <- mean(width_s >= diff(ci2)))

ci1 <- ci1[, sample(1:nrep, 300)]
ci2 <- ci2[, sample(1:nrep, 300)]
cols <- rep("grey", nrep)
cols[which(diff(ci1) <= width_s)] <- "#9B7EA6"
par(mfrow = c(1, 2))
matplot(ci1, rbind(1:300, 1:300),
        type = "l", lty = 1,
        xlab = "Coverage area", ylab = "Confidence intervals",
        main = paste("n =", n_b1, ", small CI: ", pow1 * 100, "%"), col = cols,
        xlim = c(stds[3]-4, stds[3]+4)
)
abline(v = stds[3], lty = 2)
cols <- rep("grey", nrep)
cols[which(diff(ci2) <= width_s)] <- "#9B7EA6"
matplot(ci2, rbind(1:300, 1:300),
        type = "l", lty = 1,
        xlab = "Coverage area", ylab = "Confidence intervals",
        main = paste("n =", n_s2, ", small CI: ", pow2 * 100, "%"), col = cols,
        xlim = c(stds[3]-4, stds[3]+4)
)
abline(v = stds[3], lty = 2)


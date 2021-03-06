junk_clust_plot <- function(res, xrange = NULL, sig = NULL,
                            mu = 0, mu_null = 0,
                            junk_mixture = TRUE, null_mixture = TRUE) {
  cbpalette <- c("#000000", "#999999", "#0072B2",
                 "#D55E00", "#CC79A7", "#009E73",
                 "#56B4E9", "#E69F00", "#F0E442")
  n_palette <- length(unique(res$cluster)) - length(cbpalette)
  if (n_palette >= 0) {
    cbpalette <- c(cbpalette,
                   RColorBrewer::brewer.pal(max(3, n_palette + 1),
                  "Spectral")[1:(n_palette + 1)])
  }

  # redefine null cluster as cluster zero
  res$cluster[res$cluster_class == "Null" | res$cluster_class == "Junk"] <- 0
  res$clusters <- as.factor(res$cluster)

  junk_obs <- res$cluster_class == "Junk"
  null_obs <- res$cluster_class == "Null"
  clust_obs <- which(!res$cluster_class == "Null" &
                     !res$cluster_class == "Junk")
  theta <- res$theta
  theta_se <- res$theta_se

  if (is.null(sig)) {
    rng_thet <- range(theta)
    max.disp <- which.max(abs(theta) + 2 * theta_se)
    sig <- (rng_thet[2] - rng_thet[1] + theta_se[max.disp])
  }
  if (is.null(mu)) {
    mu <- 0
  }
  if (sum(null_obs) > 0) {
    sig_null <- max(theta_se[null_obs])
    mu_null <- 0
  } else {
    sig_null <- 1
    mu_null <- 0
  }

  jit_junk <- gen_t(theta[junk_obs], df = 4, mu = mu, sig = sig, log = F)
  jit_null <- stats::dnorm(theta[null_obs], mean = mu_null, sd = sig_null,
                           log = F)
  jit_clusts <- rep(0, length(theta))
  jit_clusts[junk_obs] <- jit_junk * stats::runif(sum(junk_obs), 0, 0.5)
  jit_clusts[null_obs] <- jit_null * stats::runif(sum(null_obs), 0, 0.5)
  if (sum(junk_obs) > 0 | sum(null_obs) > 0) {
    mny <- max(c(jit_clusts[junk_obs], jit_clusts[null_obs]))
  } else {
    mny <- 1
  }
  jit_clusts[clust_obs] <- -stats::runif(length(clust_obs), 0, mny)

  res$jit_clusts <- jit_clusts
  maxx <- max(c(theta + 2.5 * theta_se, theta - 2.5 * theta_se))
  minx <- min(c(theta + 2.5 * theta_se, theta - 2.5 * theta_se))

  maxy <- if (sum(junk_obs | null_obs) > 0) {
    max(jit_clusts[junk_obs | null_obs])
  } else if (sum(junk_obs | null_obs) == 0 & junk_mixture & !null_mixture) {
    max(gen_t(c(minx, maxx), mu = mu, sig = sig, log = FALSE))
  } else if (sum(junk_obs | null_obs) == 0 & !junk_mixture & null_mixture) {
    if (sum(sig_null) > 0) {
      stats::dnorm(0, 0, sig_null)
    } else {
      0
    }
  } else if (sum(junk_obs | null_obs) == 0 & junk_mixture & null_mixture) {
    max(max(gen_t(c(minx, maxx), mu = mu, sig = sig, log = FALSE)),
      if (sum(sig_null) > 0) {
      stats::dnorm(0, 0, sig_null)
      } else {
      0
      }
    )
  }

  miny <- if (length(clust_obs) > 0) {
    min(jit_clusts[clust_obs])
  } else {
    0
  }
  mx_dens <- if (junk_mixture & !null_mixture) {
    gen_t(mu, df = 4, mu = mu, sig = sig, log = FALSE)
  } else if (!junk_mixture & null_mixture) {
    maxy
  } else if (junk_mixture & null_mixture) {
    max(gen_t(mu, mu = mu, sig = sig, log = FALSE),
          stats::dnorm(0, 0, sig_null))
  } else {
    0
  }

  if (maxy < mx_dens) {
    shrink <- maxy / mx_dens
    res$jit_clusts[junk_obs | null_obs] <- (shrink *
                                              jit_clusts[junk_obs | null_obs])
    minx <- 1.25 * minx
  } else {
    shrink <- 1
  }

  if (!is.null(xrange)) {
    maxx <- xrange[2]
    minx <- xrange[1]
  }

  annotations <- data.frame(
    xpos = c(minx, minx),
    ypos = c(-miny, miny / 1.25),
    annotateText = c("Junk/null\nestimates", "Clustered\nestimates"),
    hjustvar = c(0, 0),
    vjustvar = c(1, 0)
  )


p <- ggplot2::ggplot(data = res, ggplot2::aes(theta, jit_clusts))
p <- p + ggplot2::geom_point(ggplot2::aes(colour = cluster_class), size = 1.55)
p <- p + ggplot2::scale_color_manual(values = cbpalette)
p <- p + ggplot2::geom_errorbarh(
     ggplot2::aes(xmin = theta - 1.96 * theta_se,
                   xmax = theta + 1.96 * theta_se,
                   color = cluster_class), linetype = "solid")
p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                        color = cbpalette[1], lwd = 0.5) 
    if (junk_mixture) {
p <- p + ggplot2::geom_line(stat = "function", fun = gen_t_scale,
        args = with(res, c(df = 4, mu = mu, sig = sig, log = FALSE,
                           scale = shrink)),
        ggplot2::aes(colour = "Junk density")
      )
    }
p <- p + ggplot2::ylim(miny, maxy)
  clusts <- unique(res$cluster_class)
  clusts <- clusts[order(clusts)]
  ind <- which(!clusts %in% "Null" & !clusts %in% "Junk")
  if (sum(ind) > 0) {
    clust_label <- c("Junk cluster")
  } else {
    clust_label <- NULL
  }
  if (sum(ind) > 0) {
    for (i in seq_len(length(ind))) {
      mn <- res$cluster_mean[res$cluster == as.numeric(clusts[ind[i]])][1]
      clust_label <- c(clust_label, paste0("cluster_", ind[i]))
      col <- cbpalette[i]
      p <- p + ggplot2::geom_segment(x = mn, y = miny, xend = mn, yend = 0,
                                     color = col, linetype = "dotted",
                                     lwd = 0.2)
    }
  }
  p <- p + ggplot2::geom_text(data = annotations,
                            ggplot2::aes(x = xpos, y = ypos, hjust = hjustvar,
                                           vjust = vjustvar,
                                           label = annotateText), size = 3) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Cluster",
                                                  labels = clust_label)) +
    ggplot2::xlim(minx, maxx)
  p <- p +
    ggplot2::ylab("") + ggplot2::xlab("Two-stage ratio estimate") +
    ggplot2::ggtitle("Null and/or junk observations (top)\nand cluster
                     partitioned ratio estimates (bottom)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(), ## <- this line
      axis.text.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(colour = "white"),
      plot.background = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  if (null_mixture & sum(null_obs) > 0) {
    p <- p + ggplot2::geom_line(stat = "function",
      fun = scaled_dnorm,
      args = with(res, c(mean = mu_null, sd = sig_null, log = FALSE,
                         scale = shrink)),
      ggplot2::aes(colour = "Null density")) +
      ggplot2::ylab("") + ggplot2::xlab("Two-stage ratio estimate") +
      ggplot2::ggtitle("Null and/or junk observations (top)\nand cluster
                       partitioned ratio estimates (bottom)") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(), ## <- this line
        axis.text.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(colour = "white"),
        plot.background = ggplot2::element_blank(),
        legend.position = "bottom"
      )
  }

  return(p)
}

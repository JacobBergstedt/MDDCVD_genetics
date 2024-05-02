

forestplot_mod <- function (df, name = name, estimate = estimate, se = se, pvalue = NULL, 
                            colour = NULL, shape = NULL, logodds = FALSE, psignif = 0.05, 
                            ci = 0.95, ...) 
{
  stopifnot(is.data.frame(df))
  stopifnot(is.logical(logodds))
  name <- enquo(name)
  estimate <- enquo(estimate)
  se <- enquo(se)
  pvalue <- enquo(pvalue)
  colour <- enquo(colour)
  shape <- enquo(shape)
  args <- list(...)
  const <- stats::qnorm(1 - (1 - ci)/2)
  df <- df %>% dplyr::mutate(`:=`(!!name, factor(!!name, levels = !!name %>% 
                                                   unique() %>% rev(), ordered = TRUE)), .xmin = !!estimate - 
                               const * !!se, .xmax = !!estimate + const * !!se, .filled = TRUE, 
                             .label = sprintf("%.2f", !!estimate))
  if (logodds) {
    df <- df %>% mutate(.xmin = exp(.data$.xmin), .xmax = exp(.data$.xmax), 
                        `:=`(!!estimate, exp(!!estimate)))
  }
  if (!quo_is_null(pvalue)) {
    df <- df %>% dplyr::mutate(.filled = !!pvalue < !!psignif)
  }
  g <- ggplot2::ggplot(df, aes(x = !!estimate, y = !!name))
  if (logodds) {
    if ("xtickbreaks" %in% names(args)) {
      g <- g + scale_x_continuous(trans = "log10", breaks = args$xtickbreaks)
    }
    else {
      g <- g + scale_x_continuous(trans = "log10", breaks = scales::log_breaks(n = 7))
    }
  }
  g <- g + theme_forest() + scale_colour_ng_d() + scale_fill_ng_d() + 
    geom_stripes() + geom_vline(xintercept = ifelse(test = logodds, 
                                                    yes = 1, no = 0), linetype = "solid", size = 0.4, colour = "black")
  g <- g + geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                    colour = !!colour, shape = !!shape, filled = .data$.filled), 
                       position = ggstance::position_dodgev(height = 0.5)) + 
    ggplot2::scale_shape_manual(values = c(21L, 22L, 23L, 
                                           24L, 25L)) + guides(colour = guide_legend(reverse = TRUE), 
                                                               shape = guide_legend(reverse = TRUE))
  if ("title" %in% names(args)) {
    g <- g + labs(title = args$title)
  }
  if ("subtitle" %in% names(args)) {
    g <- g + labs(subtitle = args$subtitle)
  }
  if ("caption" %in% names(args)) {
    g <- g + labs(caption = args$caption)
  }
  if ("xlab" %in% names(args)) {
    g <- g + labs(x = args$xlab)
  }
  if (!"ylab" %in% names(args)) {
    args$ylab <- ""
  }
  g <- g + labs(y = args$ylab)
  if ("xlim" %in% names(args)) {
    g <- g + coord_cartesian(xlim = args$xlim)
  }
  if ("ylim" %in% names(args)) {
    g <- g + ylim(args$ylim)
  }
  if ("xtickbreaks" %in% names(args) & !logodds) {
    g <- g + scale_x_continuous(breaks = args$xtickbreaks)
  }
  g
}
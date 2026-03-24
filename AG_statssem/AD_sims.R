dd <- readRDS(here::here("data", "trait3_subtree.rds"))
if (!is.data.frame(dd)) dd <- as.data.frame(t(dd))
dd$ldiff <- dd$RTMB_loglik - dd$orig_loglik
tt <- t.test(dd$ldiff)
qq <- quantile(dd$ldiff, c(0.25, 0.75))

ddtb <- dd |>
  transmute(
    seed, ntaxa, ntrait, model,
    ldiff,
    opt_RTMB = RTMB_opt.time,
    opt_orig = orig_opt.time,
    tot_RTMB = RTMB_tot.time,
    tot_orig = orig_tot.time
  ) |>
  pivot_longer(
    cols = c(opt_RTMB, opt_orig, tot_RTMB, tot_orig),
    names_to      = c("type", "method"),
    names_pattern = "(opt|tot)_(RTMB|orig)",
    values_to     = "value"
  ) |>
  mutate(
    type   = factor(type,   levels = c("opt", "tot")),
    method = factor(method, levels = c("orig", "RTMB")),

    good = case_when(
      is.na(ldiff) ~ NA_character_,
      abs(ldiff) < 1e-6 ~ "OK",  
      ldiff < 0 ~ "RTMB",  
      ldiff > 0 ~ "orig"
    ),
    good = factor(good, levels = c("OK", "orig", "RTMB"))
  )

ggplot(ddtb, aes(x = ntaxa, y = value, colour = method)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ type, scales = "free_y") +
  geom_smooth(
    method = "lm",  
    se = FALSE,
    size = 0.8
  ) +
  labs(
    x = "Number of taxa",
    y = "time (seconds)",
    colour   = "method"
  ) +
    scale_colour_discrete_qualitative() +
  theme_minimal()

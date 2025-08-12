library(ggplot2)
library(dplyr)

#####################################
#### Experimental data analysis #####
#####################################

## Fluorescence reading data 

data <- data.frame(
  dataset = rep(c("eco3R_tr3","eco3R_tr4","eco3R_tr6","eco3R_tr1",
                  "eco3R_tr2","eco3R_tr7","eco3R_P","eca1_P","eca1_tr1"), each = 6),
  group = rep(rep(c("Blank","Test"), each = 3), times = 9),
  value = c(
    c(45,66,39), c(5698,5258,9622),
    c(23,23,21), c(2071,2399,3235),
    c(22,27,18), c(1596,5276,6122),
    c(32,23,27), c(1659,2895,1726),
    c(37,34,31), c(2078,1970,3429),
    c(22,21,17), c(1714,2652,5067),
    c(19,16,17), c(517,1528,1894),
    c(37,24,27), c(4531,3313,1142),
    c(41,51,48), c(2385,5633,3635)
  )
)

# Summary stats for each aptamer variant and group
summary_stats <- data %>%
  group_by(dataset, group) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )

# P-values for each aptamer dataset
p_values <- data %>%
  group_by(dataset) %>%
  summarise(
    p_val = t.test(value ~ group, var.equal = FALSE)$p.value,
    max_y = max(value)
  )

# Plot
ggplot(data, aes(x = group, y = value)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  geom_point(data = summary_stats, 
             aes(x = group, y = mean), 
             color = "red", size = 3, inherit.aes = FALSE) +
  geom_errorbar(data = summary_stats,
                aes(x = group, ymin = mean - sd, ymax = mean + sd),
                width = 0.15, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ dataset, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Blank vs Test Fluorescence Readings of aptamer variants",
       y = "Fluorescence", x = "") +
  geom_text(data = p_values,
            aes(x = 1.5, y = max_y * 0.95,
                label = paste0("p = ", signif(p_val, 3))),
            inherit.aes = FALSE, size = 4)

library('ggplot2')
library('dplyr')
library('tidyr')
library('data.table')
library('latex2exp')
library('scales')
library('cowplot')
library('RColorBrewer')
library('showtext')

# Add Avenir font, comment out if the font file is not on the machine
font_add("Avenir", "Avenir.ttc")

# Set working directory
setwd('~/decode/figures/')

#############################
# Set up general aesthetics #
#############################
text_size = 8

presentation <-  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,
                                  family = "Avenir",
                                  # face = "bold",
                                  size = text_size),
        axis.title = element_text(family = "Avenir",
                                  # face = "bold",
                                  size = text_size),
        axis.text = element_text(family = "Avenir",
                                 # face = "bold",
                                 size = text_size - 2),
        legend.text = element_text(family = "Avenir",
                                   # face = "bold",
                                   size = text_size),
        legend.title = element_text(family = "Avenir",
                                    # face = "bold",
                                    size = text_size),
        strip.text = element_text(family = "Avenir",
                                  # face = "bold",
                                  size = text_size),
        strip.text.x = element_text(vjust = 0))

cowplot_config <- theme(legend.position='top',
                        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
                        axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
                        axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0)))

############
# Figure 2 #
############

rosetta_target_data = data.frame(method = c('DeCoDe', 'SwiftLib'), n_covered  = c(191, 173))

rosetta_target_data$method = factor(rosetta_target_data$method, levels = c('SwiftLib', 'DeCoDe'))

# Plot panel A
panel_A <- ggplot(rosetta_target_data) +
  geom_bar(aes(x=method, y=n_covered, fill=method), stat='identity') +
  scale_fill_manual(values = c('DeCoDe'='#E41A1C','SwiftLib'='#377EB8')) +
  geom_text(aes(x=method, y=100, label=n_covered),
            family = "Avenir", size=2, color='white') +
  labs(x='Method',
       y='Targets covered',
       fill='Method:') +
  coord_flip() +
  presentation + 
  cowplot_config +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

rosetta_data <- fread('../results/sl_comparison/kmers.csv') %>% 
  filter(task == '1xbi') %>% 
  mutate(SLiT_frac = SLiT / target_kmers,
         DCiT_frac = DCiT / target_kmers)

rosetta_data_melt <- rosetta_data %>% 
  select(lib_limit, sublibs, k, SLiT_frac, DCiT_frac) %>% 
  gather('fraction', 'value', SLiT_frac, DCiT_frac) %>% 
  mutate(fraction = as.factor(fraction))

### Scatter

rosetta_data <- fread('../results/sl_comparison/kmers.csv') %>% 
  filter(task == '1xbi') %>% 
  mutate(SLiT_frac_unique = SLiT / target_kmers,
         DCiT_frac_unique = DCiT / target_kmers,
         SLiT_frac_all = SLiT_all / target_kmers_all,
         DCiT_frac_all = DCiT_all / target_kmers_all)

rosetta_kmer_unique <- rosetta_data %>% 
  select(k, SLiT_frac_unique, DCiT_frac_unique) %>% 
  mutate(measure = 'unique')
  
rosetta_kmer_all <- rosetta_data %>% 
  select(k, SLiT_frac_all, DCiT_frac_all) %>% 
  mutate(measure = 'all')

panel_B <- ggplot(rosetta_kmer_all, aes(x=SLiT_frac_all, y = DCiT_frac_all, color=as.factor(k))) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) + 
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  labs(x="SwiftLib",
       y="DeCoDe",
       color = "k:") + 
  ggtitle("Fraction of total\nk-mers covered") +
  scale_y_continuous(limits = c(0.895, 1), breaks=seq(0.9, 1, by = .05)) +
  scale_x_continuous(limits = c(0.895, 1), breaks=seq(0.9, 1, by = .05)) +
  presentation + 
  cowplot_config +
  coord_fixed() +
  theme(legend.position='none',
        axis.text.x = element_text(angle = 90, vjust = 0.5))

panel_C <- ggplot(rosetta_kmer_unique, aes(x=SLiT_frac_unique, y = DCiT_frac_unique, color=as.factor(k))) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) + 
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  labs(x="SwiftLib",
       y="DeCoDe",
       color = "k:") +
  ggtitle("Fraction of unique\nk-mers covered") +
  scale_y_continuous(limits = c(0.895, 1), breaks=seq(0.9, 1, by = .05)) +
  scale_x_continuous(limits = c(0.895, 1), breaks=seq(0.9, 1, by = .05)) +
  presentation + 
  cowplot_config +
  coord_fixed() +
  theme(legend.position='bottom',
        legend.margin=margin(t = 0, b = 0, unit='mm'),
        legend.text.align = 0.5,
        axis.text.x = element_text(angle = 90, vjust = 0.5))


# Assemble grid
subgrid_A <- plot_grid(panel_A,
                       labels = c(""),
                       scale=0.9)

subgrid_BC <- plot_grid(panel_B,
                     panel_C + theme(legend.position="none"),
                     ncol = 2,
                     scale = 1,
                     labels = c("B", "C"),
                     label_fontfamily = "Avenir",
                     hjust = 0,
                     vjust = 1.6)

legend <- get_legend(
    panel_C +
    theme(legend.key.width=unit(0.5, 'cm'),
          legend.text.align=0.5,
          legend.spacing.x = unit(0.2, 'cm'))
)

grid <- plot_grid(subgrid_A,
                  subgrid_BC,
                  legend,
                  labels = c("A", "", ""),
                  nrow = 3,
                  rel_heights = c(1, 1.5, 0.1),
                  scale = 1,
                  label_fontfamily = "Avenir",
                  hjust = 0,
                  vjust = 1)

# Save the plot
ggsave('figures/figure_2.png', plot = grid, device = 'png',
       scale = 1, width = 86, height = 75, units = 'mm',
       dpi = 350)

###########


levels(rosetta_data_melt$fraction) <- c('DeCoDe', 'SwiftLib')

ylab <- expression(bold(paste("Fraction of target\n", bold(italic("k")), "-mers covered")))

panel_B <- ggplot(rosetta_data_melt, aes(x=as.factor(k), y = value, fill=fraction)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(x=as.factor(k), y=0.5, label=round(value, 3)),
            family = "Avenir", size=2, angle=90, position=position_dodge2(width=0.9),
            hjust=0.5, color='white') +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, 1) + 
  labs(x='k',
       y="Fraction of unique\ntarget k-mers covered",
       fill='Method:') + 
  presentation + 
  theme(legend.position='bottom',
        legend.margin=margin(t = 0, b = 0, unit='mm'),
        legend.spacing.x = unit(1, 'lines'),
        legend.text.align = 0.5)

# Assemble grid
grid <- plot_grid(panel_A,
                  panel_B,
                  labels = "AUTO",
                  nrow = 2,
                  rel_heights = c(2, 3),
                  scale = .95,
                  label_fontfamily = "Avenir",
                  hjust = 0,
                  vjust = 1)

# Save the plot
ggsave('figures/figure_2.png', plot = grid, device = 'png',
       scale = 1, width = 86, height = 90, units = 'mm',
       dpi = 350)

############
# Figure 3 #
############

# Read in the data
data <- fread('../results/sl_comparison/results.csv')

# Refactor sublibrary count
data$sublibs <- factor(data$sublibs)

levels(data$sublibs) <- c('1 sublibrary', '2 sublibraries')

# Plot panel A
panel_A <- ggplot(data %>% filter(sublibs == '1 sublibrary')) +
  geom_bar(aes(x=as.factor(lib_limit), y=n_covered, fill=method), position='dodge', stat='identity') +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  geom_text(aes(x=as.factor(lib_limit), y=n_covered, label=n_covered), family = "Avenir",
            position=position_dodge2(width=0.9), size=2, vjust=-0.25) +
  ylim(0, 70) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='Targets covered',
       fill='Method:') +
  presentation + 
  cowplot_config +
  theme(axis.text.x = element_text(family = "Avenir", size = 6),
        panel.grid.major.x = element_blank())

# Plot panel B
panel_B <- ggplot(data %>% filter(sublibs == '2 sublibraries')) +
  geom_bar(aes(x=as.factor(lib_limit), y=n_covered, fill=method), position='dodge', stat='identity') +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  geom_text(aes(x=as.factor(lib_limit), y=n_covered, label=n_covered), family = "Avenir",
            position=position_dodge2(width=0.9), size=2, vjust=-0.25) +
  ylim(0, 70) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='Targets covered',
       fill='Method:') +
  presentation + 
  cowplot_config +
  theme(axis.text.x = element_text(family = "Avenir", size = 6),
        panel.grid.major.x = element_blank())

# Get kmer data
kmer_data <- fread('../results/sl_comparison/kmers.csv') %>% 
  filter(task == 'gfp') %>% 
  mutate(SwiftLib = SLiT / target_kmers,
         DeCoDe = DCiT / target_kmers,
         sublibs = as.factor(sublibs),
         lib_limit = as.factor(lib_limit),
         k = as.factor(k)) %>% 
  select(lib_limit, sublibs, k, SwiftLib, DeCoDe) %>% 
  gather('method', 'fraction', SwiftLib, DeCoDe)

levels(kmer_data$sublibs) <- c('1 sublibrary', '2 sublibraries')
levels(kmer_data$k) <- c('k=2', 'k=3', 'k=4')


###### Scatter

kmer_data <- fread('../results/sl_comparison/kmers.csv') %>% 
  filter(task == 'gfp') %>% 
  mutate(SwiftLib_unique = SLiT / target_kmers,
         DeCoDe_unique = DCiT / target_kmers,
         SwiftLib_all = SLiT_all / target_kmers_all,
         DeCoDe_all = DCiT_all / target_kmers_all,
         sublibs = as.factor(sublibs),
         lib_limit = as.factor(lib_limit),
         k = as.factor(k)) %>% 
  select(lib_limit, sublibs, k, SwiftLib_unique, DeCoDe_unique, SwiftLib_all, DeCoDe_all)

levels(kmer_data$k) <- c('k=2', 'k=3', 'k=4')

panel_C <- ggplot(kmer_data, aes(x=SwiftLib_all, y = DeCoDe_all, color = lib_limit, shape = sublibs)) +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2) +
  geom_point() +
  scale_color_brewer(palette = 'Set2') +
  facet_grid(.~k) +
  labs(x='SwiftLib',
       y='DeCoDe',
       color='Library size limit:') +
  ggtitle("Fraction of total k-mers covered") +
  presentation +
  cowplot_config +
  coord_fixed(xlim = c(0.75, 1), ylim = c(0.75, 1)) +
  scale_x_continuous(breaks = seq(.7, 1, .1)) +
  scale_y_continuous(breaks = seq(.7, 1, .1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

panel_D <- ggplot(kmer_data, aes(x=SwiftLib_unique, y = DeCoDe_unique, color = lib_limit, shape = sublibs)) +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2) +
  geom_point() +
  scale_color_brewer(palette = 'Set2',
                     labels=c('100000' = parse(text = TeX('$10^5$')),
                              '1000000' = parse(text = TeX('$10^6$')),
                              '10000000' = parse(text = TeX('$10^7$')),
                              '100000000' = parse(text = TeX('$10^8$')),
                              '1000000000' = parse(text = TeX('$10^9$')))) +
  facet_grid(.~k) +
  labs(x='SwiftLib',
       y='DeCoDe',
       color='Library size limit:',
       shape="Sublibrary count: ") +
  ggtitle("Fraction of unique k-mers covered") +
  presentation +
  cowplot_config +
  coord_fixed(xlim = c(0.2, 0.55), ylim = c(0.2, 0.55)) +
  theme(legend.box = "vertical",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

legend1 <- get_legend(
  panel_A +
    theme(legend.key.width=unit(0.5, 'cm'),
          legend.text.align=0.5,
          legend.spacing.x = unit(0.5, 'cm'))
)

legend2 <- get_legend(
  panel_D +
    theme(legend.key.width=unit(0.05, 'cm'),
          legend.text.align=0,
          legend.spacing.x = unit(0.5, 'cm'),
          # legend.spacing.y = unit(.0001, 'cm'),
          legend.margin = unit(0.05, 'cm'),
          legend.text = element_text(
            margin = margin(l = -10, r = 0, unit = "pt"))))

# Plot the grid
grid <- plot_grid(panel_A + theme(legend.position="none"),
                  panel_C + theme(legend.position="none"),
                  panel_B + theme(legend.position="none"),
                  panel_D + theme(legend.position="none"),
                  legend1,
                  legend2,
                  labels = c("A", "C", "B", "D", "", ""),
                  ncol = 2,
                  scale = 0.95,
                  rel_widths = c(1, 1),
                  rel_heights = c(1, 1, .3),
                  label_fontfamily = "Avenir",
                  hjust = 0,
                  vjust = 1)

# Save the plot
ggsave('figures/figure_3.png', plot = grid, device = 'png',
       scale = 1, width = 178, height = 100, units = 'mm',
       dpi = 350)
###################

# Plot panel C
kmer_1sublib <- kmer_data %>% 
  filter(sublibs == '1 sublibrary')

panel_B <- ggplot(kmer_1sublib, aes(x=lib_limit, y = fraction, fill=method)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = 'Set1') +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  facet_grid(sublibs~k) +
  labs(x='Total library size limit',
       y='Fraction of unique\ntarget k-mers covered',
       fill='Method:') +
  # ylim(0, 0.6) + 
  presentation +
  cowplot_config +
  theme(axis.text.x = element_text(family = "Avenir",
        size = 6))

# Plot panel D
kmer_2sublib <- kmer_data %>% 
  filter(sublibs == '2 sublibraries')

panel_D <- ggplot(kmer_2sublib, aes(x=lib_limit, y = fraction, fill=method)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = 'Set1') +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  facet_grid(sublibs~k) +
  labs(x='Total library size limit',
       y='Fraction of unique\ntarget k-mers covered',
       fill='Method:') +
  # ylim(0, 0.6) + 
  presentation +
  cowplot_config +
  theme(axis.text.x = element_text(family = "Avenir",
                                 size = 6))

legend <- get_legend(
  panel_A +
    theme(legend.key.width=unit(0.5, 'cm'),
          legend.text.align=0.5,
          legend.spacing.x = unit(0.5, 'cm'))
)

# Plot the grid
grid <- plot_grid(panel_A + theme(legend.position="none"),
                  panel_B + theme(legend.position="none"),
                  panel_C + theme(legend.position="none"),
                  panel_D + theme(legend.position="none"),
                  labels = "AUTO",
                  ncol = 2,
                  rel_widths = c(2, 3),
                  label_fontfamily = "Avenir",
                  hjust = 0,
                  vjust = 1)

grid <- plot_grid(grid, legend,
                  ncol = 1,
                  rel_heights = c(3, 0.15),
                  label_fontfamily = "Avenir",
                  hjust = 0,
                  vjust = 1)

# Save the plot
ggsave('figures/figure_3.png', plot = grid, device = 'png',
       scale = 1, width = 178, height = 100, units = 'mm',
       dpi = 350)

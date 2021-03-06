library('ggplot2')
library('dplyr')
library('tidyr')
library('data.table')
library('latex2exp')
library('scales')
library('cowplot')

text_size = 8

# Global aesthetics
presentation <-  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,
                                  family = "Avenir",
                                  face = "bold",
                                  size = text_size),
        axis.title = element_text(family = "Avenir",
          face = "bold",
          size = text_size),
        axis.text = element_text(family = "Avenir",
          # face = "bold",
          size = text_size),
        legend.text = element_text(family = "Avenir",
          # face = "bold",
          size = text_size),
        legend.title = element_text(family = "Avenir",
          face = "bold",
          size = text_size),
        strip.text = element_text(family = "Avenir",
          # face = "bold",
          size = text_size),
        strip.text.x = element_text(vjust = 0))

cowplot_config <- theme(panel.grid.major.x = element_blank(),
                        legend.position='top',
                        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# Read in the data
setwd('decode/figures/figure2/')

data <- fread('../../results/sl_comparison/results.csv')

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
  cowplot_config

# Plot panel C
panel_C <- ggplot(data %>% filter(sublibs == '2 sublibraries')) +
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
  cowplot_config

# Get kmer data
kmer_data <- fread('../../results/sl_comparison/kmers.csv') %>% 
  filter(task == 'gfp') %>% 
  mutate(SwiftLib = SLiT / target_kmers,
         DeCoDe = DCiT / target_kmers,
         sublibs = as.factor(sublibs),
         lib_limit = as.factor(lib_limit)) %>% 
  select(lib_limit, sublibs, k, SwiftLib, DeCoDe) %>% 
  gather('method', 'fraction', SwiftLib, DeCoDe)

levels(kmer_data$sublibs) <- c('1 sublibrary', '2 sublibraries')


# Construct labeller
levels(kmer_data$lib_limit) <- c(
  '10^5',
  '10^6',
  '10^7',
  '10^8',
  '10^9'
)

custom_labeller = labeller(lib_limit = as_labeller(kmer_data$lib_limit, label_parsed))

# Plot panel C
kmer_1sublib <- kmer_data %>% 
  filter(sublibs == '1 sublibrary')

panel_B <- ggplot(kmer_1sublib, aes(x=as.factor(k), y = fraction, fill=method)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = 'Set1') +
  geom_text(aes(x=as.factor(k), y=fraction + .01, label=round(fraction, 3)),
            family = "Avenir", size=2, angle=90, position=position_dodge2(width=0.9),
            hjust=0) +
  facet_grid(sublibs~lib_limit, labeller = custom_labeller) +
  labs(x='k-mer size',
       y='Fraction of target\nk-mers covered',
       fill='Method:') +
  ylim(0, 0.7) + 
  presentation +
  cowplot_config

# Plot panel D
kmer_2sublib <- kmer_data %>% 
  filter(sublibs == '2 sublibraries')

panel_D <- ggplot(kmer_2sublib, aes(x=as.factor(k), y = fraction, fill=method)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = 'Set1') +
  geom_text(aes(x=as.factor(k), y=fraction + .01, label=round(fraction, 3)),
            family = "Avenir", size=2, angle=90, position=position_dodge2(width=0.9),
            hjust=0) +
  facet_grid(sublibs~lib_limit, labeller = custom_labeller) +
  labs(x='k-mer size',
       y='Fraction of target\nk-mers covered',
       fill='Method:') +
  ylim(0, 0.7) + 
  presentation +
  cowplot_config

# Handle legend
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
                  rel_widths = c(2, 3))

grid <- plot_grid(grid, legend, ncol = 1, rel_heights = c(3, 0.15))

grid

# Save the plot
ggsave('figure_2.png', plot = last_plot(), device = 'png',
       scale = 1, width = 178, height = 100, units = 'mm',
       dpi = 350)

####################
#### SI Figures ####
####################

ggplot(data) +
  geom_bar(aes(x=as.factor(lib_limit), y=total_lib_size, fill=method), position='dodge', stat='identity') +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='Total library size',
       fill='Method:') +
  presentation + 
  theme(panel.grid.major.x = element_blank(),
        legend.position='top')

ggsave('Fig_2_SI_1.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)

# Runtime
ggplot(data) +
  geom_bar(aes(x=as.factor(lib_limit), y=time, fill=method), position='dodge', stat='identity') +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept=172800, color='black', linetype=2, lwd=1) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='Runtime (s)',
       fill='Method:') +
  presentation + 
  theme(panel.grid.major.x = element_blank(),
      legend.position='top')

ggsave('Fig_2_SI_2.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)

# Memory usage
ggplot(data) +
  geom_bar(aes(x=as.factor(lib_limit), y=max_mem / 1000000, fill=method), position='dodge', stat='identity') +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='Maximum memory usage (GB)',
       fill='Method:') +
  presentation + 
  theme(panel.grid.major.x = element_blank(),
        legend.position='top')

ggsave('Fig_2_SI_3.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)


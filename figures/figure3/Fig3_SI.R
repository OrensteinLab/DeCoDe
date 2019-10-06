library('ggplot2')
library('latex2exp')
library('data.table')
library('scales')
library('dplyr')
library('tidyr')

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

# Read in the data
setwd('decode/figures/figure3/')

data <- fread('../../results/multi_sublib/results.csv')

# Refactor sublibrary count
data$sublibs <- factor(data$sublibs)

# Log data
log_data <- fread('../../results/multi_sublib/log.csv')

####################
#### SI Figures ####
####################

# Coverage
ggplot(data) +
  geom_bar(aes(x=sublibs, y=n_covered), stat='identity') +
  geom_text(aes(x=sublibs, y=n_covered, label=n_covered), family = "Avenir", vjust=-0.25) +
  geom_hline(yintercept=131, color='darkgreen', linetype=2, lwd=1) +
  labs(x='Number of sublibraries',
       y='Number of target sequences covered') +
  presentation

ggsave('Fig_3_SI_1.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)

ggplot(data) +
  geom_bar(aes(x=sublibs, y=n_covered), stat='identity', fill=c('#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F')) +
  geom_hline(yintercept=131, color='black', linetype=2, lwd=0.5) +
  labs(x='Sublibraries',
       y='Targets covered') +
  presentation +
  theme(panel.border = element_blank(),
        # axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.5, color='darkgray'),
        # panel.grid.minor.y = element_line( size=.5, color='gray'),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave('Fig_3_SI_1.svg', plot = last_plot(),
       scale = 1, width = 45, height = 25, units = 'mm',
       dpi = 300)

# Total library size
ggplot(data) +
  geom_bar(aes(x=sublibs, y=total_lib_size), stat='identity') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept=10^7, color='red', linetype=2, lwd=1) +
  labs(x='Number of sublibraries',
       y='Total library size') +
  presentation

ggsave('Fig_3_SI_2.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)

# Runtime
ggplot(data) +
  geom_bar(aes(x=sublibs, y=time), stat='identity') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept=172800, color='red', linetype=2, lwd=1) +
  labs(x='Number of sublibraries',
       y='Runtime (s)') +
  presentation

ggsave('Fig_3_SI_3.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)

# Maximum memory usage
ggplot(data) +
  geom_bar(aes(x=sublibs, y=max_mem / 1000000), stat='identity') +
  labs(x='Number of sublibraries',
       y='Maximum memory usage (GB)') +
  presentation

ggsave('Fig_3_SI_4.png', plot = last_plot(), device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 300)

# Log traces
log_data_reformed <- log_data %>%
  select(sublibs, time, solution, bound) %>% 
  group_by(sublibs, time) %>% 
  gather('var', 'val', solution, bound)

ggplot(log_data_reformed) +
  geom_line(aes(x=time / 60, y=val, color=var)) +
  facet_wrap(.~sublibs, nrow=2, scales='free_x') +
  scale_color_brewer(palette = "Set1") +
  labs(x='Time (m)',
       y='ILP solution size',
       color='Solution\nelement') +
  presentation +
  theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5))

ggsave('Fig_3_SI_5.png', plot = last_plot(), device = 'png',
       scale = 1, width = 160, height = 100, units = 'mm',
       dpi = 300)


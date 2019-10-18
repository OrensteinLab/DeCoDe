library('ggplot2')
library('dplyr')
library('tidyr')
library('data.table')
library('latex2exp')
library('scales')
library('cowplot')

# Set working directory
setwd('~/decode/figures/')

#############################
# Set up general aesthetics #
#############################
text_size = 8

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

##############
# SI Figures #
##############

###############
# 1 through 3 #
###############

# Read in the data
data <- fread('../results/sl_comparison/results.csv')

# Refactor sublibrary count
data$sublibs <- factor(data$sublibs)

levels(data$sublibs) <- c('1 sublibrary', '2 sublibraries')

fig_s1 <- ggplot(data) +
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

ggsave('figures/figure_s1.png', plot = fig_s1, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

# Runtime
fig_s2 <- ggplot(data) +
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

ggsave('figures/figure_s2.png', plot = fig_s2, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

# Memory usage
fig_s3 <- ggplot(data) +
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

ggsave('figures/figure_s3.png', plot = fig_s3, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

###############
# 5 through 9 #
###############

data <- fread('../results/multi_sublib/results.csv')

# Refactor sublibrary count
data$sublibs <- factor(data$sublibs)

# Log data
log_data <- fread('../results/multi_sublib/log.csv')

# Coverage
fig_s5 <- ggplot(data) +
  geom_bar(aes(x=sublibs, y=n_covered), stat='identity') +
  geom_text(aes(x=sublibs, y=n_covered, label=n_covered), family = "Avenir", vjust=-0.25) +
  geom_hline(yintercept=131, color='darkgreen', linetype=2, lwd=1) +
  labs(x='Number of sublibraries',
       y='Number of target sequences covered') +
  presentation

ggsave('figures/figure_s5.png', plot = fig_s5, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

# Total library size
fig_s6 <- ggplot(data) +
  geom_bar(aes(x=sublibs, y=total_lib_size), stat='identity') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept=10^7, color='red', linetype=2, lwd=1) +
  labs(x='Number of sublibraries',
       y='Total library size') +
  presentation

ggsave('figures/figure_s6.png', plot = fig_s6, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

# Runtime
fig_s7 <- ggplot(data) +
  geom_bar(aes(x=sublibs, y=time), stat='identity') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept=172800, color='red', linetype=2, lwd=1) +
  labs(x='Number of sublibraries',
       y='Runtime (s)') +
  presentation

ggsave('figures/figure_s7.png', plot = fig_s7, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

# Maximum memory usage
fig_s8 <- ggplot(data) +
  geom_bar(aes(x=sublibs, y=max_mem / 1000000), stat='identity') +
  labs(x='Number of sublibraries',
       y='Maximum memory usage (GB)') +
  presentation

ggsave('figures/figure_s8.png', plot = fig_s8, device = 'png',
       scale = 1, width = 125, height = 125, units = 'mm',
       dpi = 350)

# Log traces
log_data_reformed <- log_data %>%
  select(sublibs, time, solution, bound) %>% 
  group_by(sublibs, time) %>% 
  gather('var', 'val', solution, bound)

fig_s9 <- ggplot(log_data_reformed) +
  geom_line(aes(x=time / 60, y=val, color=var)) +
  facet_wrap(.~sublibs, nrow=2, scales='free_x') +
  scale_color_brewer(palette = "Set1") +
  labs(x='Time (m)',
       y='ILP solution size',
       color='Solution\nelement') +
  presentation +
  theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5))

ggsave('figures/figure_s9.png', plot = fig_s9, device = 'png',
       scale = 1, width = 160, height = 100, units = 'mm',
       dpi = 350)

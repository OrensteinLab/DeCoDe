library('ggplot2')
library('dplyr')
library('data.table')
library('latex2exp')
library('scales')

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
setwd('decode/figures/figure2/')

data <- fread('../../results/sl_comparison/results.csv')

# Refactor sublibrary count
data$sublibs <- factor(data$sublibs)

levels(data$sublibs) <- c('1 sublibrary', '2 sublibraries')

# Plot Figure 2
ggplot(data) +
  geom_bar(aes(x=as.factor(lib_limit), y=n_covered, fill=method), position='dodge', stat='identity') +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='Number of target sequences covered',
       fill='Method:') +
  presentation + 
  theme(panel.grid.major.x = element_blank(),
        legend.position='top')

# Save the plot
ggsave('figure_2.png', plot = last_plot(), device = 'png',
       scale = 1, width = 84, height = 125, units = 'mm',
       dpi = 300)

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


# Subset data
decode_data <- data %>% filter(method == 'DeCoDe')

# Runtime
ggplot(decode_data) +
  geom_bar(aes(x=as.factor(lib_limit), y=time), stat='identity') +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='DeCoDe runtime (s)') +
  presentation

# Memory usage
ggplot(decode_data) +
  geom_bar(aes(x=as.factor(lib_limit), y=max_mem / 1000000), stat='identity') +
  scale_x_discrete(labels=c('100000' = parse(text = TeX('$10^5$')),
                            '1000000' = parse(text = TeX('$10^6$')),
                            '10000000' = parse(text = TeX('$10^7$')),
                            '100000000' = parse(text = TeX('$10^8$')),
                            '1000000000' = parse(text = TeX('$10^9$')))) +
  facet_grid(sublibs~.) +
  labs(x='Total library size limit',
       y='DeCoDe maximum memory usage (GB)') +
  presentation


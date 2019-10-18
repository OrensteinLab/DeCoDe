library('ggplot2')
library('latex2exp')
library('data.table')
library('scales')
library('dplyr')
library('tidyr')

text_size = 8

# Global aesthetics
presentation <- theme_bw() +
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
        strip.text.x = element_text(vjust = 0),
        panel.grid.major.x = element_blank())

setwd('~/decode')
gfp_data <- fread('results/sl_comparison/kmers.csv') %>% 
  filter(task == 'gfp') %>% 
  mutate(SLiT_frac = SLiT / target_kmers,
         DCiT_frac = DCiT / target_kmers)

ggplot(gfp_data, aes(x=SLiT_frac, y=DCiT_frac, color=as.factor(k))) +
  geom_point() +
  facet_grid(sublibs~lib_limit) +
  geom_abline(slope=1, intercept=0, color = 'red', linetype=2) + 
  xlim(0, 1) +
  ylim(0, 1) + 
  coord_fixed() +
  presentation

gfp_data_melt <- gfp_data %>% 
  select(lib_limit, sublibs, k, SLiT_frac, DCiT_frac) %>% 
  gather('fraction', 'value', SLiT_frac, DCiT_frac)

ggplot(gfp_data_melt, aes(x=as.factor(k), y = value, fill=fraction)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(sublibs~lib_limit) +
  ylim(0, 1) + 
  presentation

rosetta_data <- fread('results/sl_comparison/kmers.csv') %>% 
  filter(task == '1xbi') %>% 
  mutate(SLiT_frac = SLiT / target_kmers,
         DCiT_frac = DCiT / target_kmers)
  
rosetta_data_melt <- rosetta_data %>% 
  select(lib_limit, sublibs, k, SLiT_frac, DCiT_frac) %>% 
  gather('fraction', 'value', SLiT_frac, DCiT_frac) %>% 
  mutate(fraction = as.factor(fraction))

levels(rosetta_data_melt$fraction) <- c('DeCoDe', 'SwiftLib')

ggplot(rosetta_data_melt, aes(x=as.factor(k), y = value, fill=fraction)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(x=as.factor(k), y=0.5, label=round(value, 3)),
            family = "Avenir", size=2, angle=90, position=position_dodge2(width=0.9),
            hjust=0.5, color='white') +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, 1) + 
  labs(x='k-mer size',
       y='Fraction of target\nk-mers covered',
       fill='Method:') + 
  presentation + 
  theme(legend.position='bottom',
        legend.margin=margin(t = 0, b = 0, unit='mm'))

# Save the plot
ggsave('1xbi_kmers.png', plot = last_plot(), device = 'png',
       scale = 1, width = 86, height = 55, units = 'mm',
       dpi = 350)
  
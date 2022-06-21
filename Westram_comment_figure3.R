## load libraries
library(ggplot2)
library(cowplot)
library(viridis)
theme_set(theme_cowplot())

## model parameters
q = 0 # favored allele in B
p = 1-q # unfavored allele in B
m = 0.001 # migration rate
s = 0.25 # power of selection, must be (s >> m)
r = seq(0.01,1,0.01) # recombination rate per generation
v = c(0,0.2,0.8,1) # rate of transmisison of epigenetic locus

## dataframe to store calulcated values
df = data.frame(r = rep(r,4),
                v = rep(v,each=100),
                RI = NA)

## calculating RI for different r and v values
for(i in 1:length(df[,1])) {
  rv=df$r[i]+(1-df$r[i])*(1-df$v[i]) # functional recombination rate given transmission rate
  me = m*(q+(p*(1-s)*rv)/(1-(1-s)*(1-rv))) # effective migration rate
  ri = 1-me/m # reproductive isolation
  df$RI[i] = ri
}

df$Transmission = paste0("V=",df$v)

## plotting output
p1 = ggplot(data = df, aes(x = r, y = RI, color = Transmission)) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  geom_line(size = 3) +
  theme(text = element_text(size = 20)) +
  ylim(c(0,1))

ggsave("figure3.png",p1,width = 6,height = 6)

library(tidyverse)
data  <-  tibble(x = runif(10000), y = runif(10000))
ggplot(data)+geom_massogram(aes(x=x,mass=y))

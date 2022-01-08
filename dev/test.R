library(tidyverse)
data  <-  tibble(x = runif(10000), y = runif(10000))
ggplot(data)+geom_massogram(aes(x=x,mass=y))



dbrk_surf <- function(n=100, dist=access)
{
  function(d) {
    if(d%in%names(dist))
      quantile(x = dist[[d]], probs=0:n/n, na.rm=TRUE)
    else
      NA
  }
}
access <- arrow::read_parquet("dev/access.parquet")
dffh <- arrow::read_parquet("dev/dffh.parquet")
ggplot(dffh |> drop_na(to50k, dstoth_m, dstoth_p))+
  geom_massogram(aes(x=to50k, mass = dstoth_m), cuts = dbrk_surf(100), fill="red", alpha=0.3, lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE)+
  geom_massogram(aes(x=to50k, mass = dstoth_p), cuts = dbrk_surf(100), fill="green", alpha=0.3, lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE)

library(raster)
library(sf)
library(patchwork)
library(tidyverse)
library(scales)
library(quantogram)
library(glue)
library(vroom)
library(colorspace)

source("dev/f.map utils.r")
source("dev/f.iso2time.r")

localdata <- "G:\\Mon Drive\\DVFdata\\Rda\\"

## fabrique les breaks pour distance, résolution 200m ------------------------------
csvsources <- "{localdata}/csv sources" %>% glue
SVG <- "svg"
c200  <- qs::qread(str_c(localdata, "c200.rda"))
iris <- qs::qread(str_c(localdata, "iris15.rda"))

iris <- iris %>%
  mutate(x=st_coordinates(st_centroid(iris))[,1],
         y=st_coordinates(st_centroid(iris))[,2],
         emp=EMP09,
         pop=P15_POP)
bs <- iris %>%
  filter(UU2010=="00851") %>%
  ungroup() %>%
  as_tibble() %>%
  summarise(mx = mean(x, na.rm=TRUE),
            my = mean(y, na.rm=TRUE),
            mwx=weighted.mean(x,w=pop, na.rm=TRUE),
            mwy=weighted.mean(y,w=pop, na.rm=TRUE),
            mwx_e=weighted.mean(x,w=emp, na.rm=TRUE),
            mwy_e=weighted.mean(y,w=emp, na.rm=TRUE))

bs <- bind_rows(
  tibble(mean="unweighted", x = bs$mx, y=bs$my),
  tibble(mean="weighted pop", x = bs$mwx, y=bs$mwy),
  tibble(mean="weighted emp", x = bs$mwx_e, y=bs$mwy_e)) %>%
  st_as_sf(coords=c("x", "y"), crs=3035)

paris <- st_centroid(st_union(iris %>% filter(DEP=="75"))) %>% st_sf() %>% mutate(mean="75")

bs <- bind_rows(bs, paris)

nd_line <- function(xy, alpha)
{ xy[,2]>=tan(alpha)*(xy[,1]-3760095)+2889557}

uu851 <- iris %>% filter(UU2010=="00851") %>% st_union()
c200851 <- c200 %>% filter(st_within(., uu851, sparse=FALSE))
rm(c200)

tR5emp09 <- qs::qread(str_c(localdata,"tr_r5_2020.rda"))$EMP09 %>%
  iso2time(c(100,1000)*1000) %>%
  r2dt(resolution=200) %>%
  as_tibble() %>%
  rename(t2emp09_100k = to100k,
         t2emp09_1M = to1M)

tR5pop <- qs::qread(str_c(localdata,"popc200_2020.rda"))$pop %>%
  iso2time(c(100,1000)*1000) %>%
  r2dt(resolution=200) %>%
  as_tibble() %>%
  rename(t2pop_100k = to100k,
         t2pop_1M = to1M)

tR5emp17 <- qs::qread(str_c(localdata,"tr_r5emp17_2020.rda"))$emp17 %>%
  iso2time(c(100,1000)*1000) %>%
  r2dt(resolution=200) %>%
  as_tibble() %>%
  rename(t2emp17_100k = to100k,
         t2emp17_1M = to1M)

c851 <- c200851 %>%
  st_centroid() %>%
  mutate(Adultes=Ind_18_24+Ind_25_39+Ind_40_54+Ind_55_64+Ind_65_79+Ind_80p,
         dep = str_sub(Depcom, 1,2)) %>%
  dplyr::select(idINS200, Ind, Men, Adultes, dep, Depcom)

c851 <- c851 %>%
  left_join(tR5emp09, by="idINS200") %>%
  left_join(tR5pop, by="idINS200") %>%
  left_join(tR5emp17, by="idINS200")

dd <- (st_distance(c851, bs) %>% units::drop_units())/1000
xy <- st_coordinates(c851)
c851 <- c851 %>% mutate(dist2b=dd[,1]/max(dd[,1], na.rm=TRUE)*100,
                        dist2wbp=dd[,2]/max(dd[,2], na.rm=TRUE)*100,
                        dist2wbe=dd[,3]/max(dd[,3], na.rm=TRUE)*100,
                        dist275 = dd[,4]/max(dd[,4], na.rm=TRUE)*100,
                        NS0 = if_else(nd_line(xy, -0.12), 1,-1),
                        NS1 = if_else(nd_line(xy, -0.62), 1,-1),
                        EO = if_else(nd_line(xy, pi/2), 1,-1),
                        dist275_NS0 = dist275*NS0,
                        dist275_NS1 = dist275*NS1,
                        dist275_EO = dist275*EO,
                        t2emp17_100k_NS0 = t2emp17_100k * NS0,
                        t2emp17_100k_NS1 = t2emp17_100k * NS1,
                        t2emp17_100k_EO = t2emp17_100k * EO) %>%
  st_drop_geometry()

# tmap_mode("plot")
# ns0 <- tm_shape(c851 %>% select(idINS200, NS0) %>% dt2r(resolution=200))+tm_raster(palette=c("grey", "black"), legend.show = FALSE )+tm_layout(frame=FALSE)
# ns1 <- tm_shape(c851 %>% select(idINS200, NS1) %>% dt2r(resolution=200))+tm_raster(palette=c("grey", "black"), legend.show = FALSE )+tm_layout(frame=FALSE)
# eo <- tm_shape(c851 %>% select(idINS200, EO) %>% dt2r(resolution=200))+tm_raster(palette=c("grey", "black"), legend.show = FALSE )+tm_layout(frame=FALSE)
#
# tmap_save(tm=ns0, "ns0.svg")
# tmap_save(tm=ns1, "ns1.svg")

## fonctions de break -------------------------------
dbrk_surf <- function(n=100, dist=c851)
{
  function(d) {
    if(d%in%names(dist))
      quantile(x = dist[[d]], probs=0:n/n, na.rm=TRUE)
    else
      NA
  }
}

dbrk_ind <- function(n=100, dist=c851)
{
  function(d) {
    if(d%in%names(dist))
      weighted_quantile_x(x = dist[[d]], w=dist$Ind, probs=0:n/n, na.rm=TRUE)
    else
      NA
  }
}

cut_pivot <- function(data, ddata, distances, n=100)
{
  names(distances) <- distances
  newx <- map_dfc(distances, ~{
    cuts <- dbrk_surf(n=n, dist=ddata)(.x)
    labels_x <- (tail(cuts, -1)+head(cuts,-1))/2
    x_cutted <- findInterval(data[[.x]], cuts, all.inside=TRUE)
    labels_x[x_cutted]
  })
  out <- data %>% select(-all_of(distances)) %>% bind_cols(newx)
  out %>% pivot_longer(all_of(distances), names_to="metric", values_to = "distance")
}

trans_cuts <- function(cuts)
{
  ncuts <- names(cuts)
  icuts <-str_extract(ncuts, pattern = "[:digit:]*\\.?[:digit:]*") %>% as.numeric()
  scales::trans_new(
    "qtrans",
    approxfun(y=icuts, x=cuts),
    approxfun(x=icuts, y=cuts),
    breaks = extended_breaks(),
    minor_breaks = regular_minor_breaks(),
    format = scales::format_format(),
    domain = c(min(cuts), max(cuts)))
}

# difference entre distances ----------------

gd <- ggplot(c851)+
  geom_point(aes(y=dist275, x=t2emp17_100k, size=Ind, col=dep), alpha=0.05)+
  xlab("Time to 100k (r5, transit)")+
  ylab("Distance to Notre Dame")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
gd2 <- ggplot(c851)+
  geom_point(aes(y=dist275, x=t2emp17_1M, size=Ind, col=dep), alpha=0.05)+
  xlab("Time to 1M (r5, transit)")+
  ylab("Distance to Notre Dame")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
gg <- gd + gd2+plot_layout(guides = 'collect')
ggsave(plot=gg, "{SVG}/distances comp.jpg" %>% glue, width=25, height = 17, units="cm")
library(ggcorrplot)
ggcorrplot(cor(c851 %>% select(dist275, t2emp17_100k, t2emp17_1M), use="complete.obs"), lab = TRUE)

round(cor(c851 %>% select(dist275, t2emp17_100k, t2emp17_1M), use="complete.obs"), 2)
# prix et valeur  -----------------------------------

dv3f <- qs::qread(str_c(localdata, "dv3fv411.rda"))
fn_l <- map(set_names(10:18, str_c("dpp", 10:18)), function(x) (substitute( ~ .x/.data[[xx]]-1, list(xx=str_c("lvm_20",x)))))
data <- cut_pivot(dv3f %>% filter(inseenotaires) %>% select(dist275,t2emp17_100k, lvm=log_prix, valeurfonc, date_ma),
                  c851, c("dist275", "t2emp17_100k")) %>%
  group_by(date_ma,distance, metric) %>%
  drop_na(lvm, valeurfonc) %>%
  summarize(lvm=median(lvm), valeurfonc=median(valeurfonc), svf = sum(valeurfonc), nt = n()) %>%
  pivot_wider(id_cols = c(distance, metric), names_from = date_ma, values_from = c(lvm, valeurfonc, svf))
datalvm <- data %>%
  dplyr::filter(metric=="dist275") %>%
  dplyr::select(distance, lvm_2010:lvm_2019) %>%
  dplyr::mutate(across(.cols=lvm_2018:lvm_2010, ~exp((.data[["lvm_2019"]]-.x)/(19-as.numeric(str_sub(cur_column(), -2, -1))))-1, .names="d{.col}")) %>%
  tidyr::pivot_longer(cols=starts_with("dlvm"))
gp <- ggplot(datalvm)+geom_line(aes(x=distance, y=value, col=name, group=name), show.legend=FALSE)+
  scale_y_continuous(labels=scales::label_percent(1))+
  coord_trans(x=trans_cuts(dbrk_surf(100)("dist275")))+
  facet_wrap(~name, labeller=as_labeller(set_names( str_c("2019-", 2010:2018), str_c("dlvm_", 2010:2018))))+
  ylab("Annualized growth rate of median price per m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")
datalvm2 <- data %>%
  filter(metric=="t2emp17_100k") %>%
  select(distance, lvm_2010:lvm_2019) %>%
  mutate(across(.cols=lvm_2018:lvm_2010, ~exp((.data[["lvm_2019"]]-.x)/(19-as.numeric(str_sub(cur_column(), -2, -1))))-1, .names="d{.col}")) %>%
  pivot_longer(cols=starts_with("dlvm"))
gp2 <- ggplot(datalvm2)+geom_line(aes(x=distance, y=value, col=name, group=name), show.legend=FALSE)+
  scale_y_continuous(labels=percent)+
  coord_trans(x=trans_cuts(dbrk_surf(100)("t2emp17_100k")))+
  facet_wrap(~name, labeller=as_labeller(set_names( str_c("2019-", 2010:2018), str_c("dlvm_", 2010:2018))))+
  ylab("Annualized growth rate of median price per m²")+
  xlab("Time to 100k (r5, transit)")

gp <- ggplot(datalvm %>% filter(name=="dlvm_2010"))+geom_line(aes(x=distance, y=value), col="darkorange1", show.legend=FALSE)+
  scale_y_continuous(labels=label_percent(1), limits=c(-0.02,0.05))+
  coord_trans(x=trans_cuts(dbrk_surf(100)("dist275")))+
  facet_wrap(~name, labeller=as_labeller(set_names( str_c("2019-", 2010:2018), str_c("dlvm_", 2010:2018))))+
  ylab("Annualized growth rate of median price per m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")
gp2 <- ggplot(datalvm2 %>% filter(name=="dlvm_2010"))+geom_line(aes(x=distance, y=value), col="dodgerblue3", show.legend=FALSE)+
  scale_y_continuous(labels=label_percent(1), limits=c(-0.02,0.05))+
  coord_trans(x=trans_cuts(dbrk_surf(100)("t2emp17_100k")))+
  facet_wrap(~name, labeller=as_labeller(set_names( str_c("2019-", 2010:2018), str_c("dlvm_", 2010:2018))))+
  ylab("Annualized growth rate of median price per m²")+
  xlab("Time to 100k (r5, transit)")
gp2level <- ggplot(datalvm2 %>% filter(name=="dlvm_2018"))+
  geom_line(aes(x=distance, y=exp(lvm_2019)/quantile(exp(lvm_2019),0.05)), col="dodgerblue3", show.legend=FALSE)+
  geom_line(aes(x=distance, y=exp(lvm_2010)/quantile(exp(lvm_2010),0.05)), col="darkolivegreen3", show.legend=FALSE)+
  scale_y_continuous()+
  annotate("text", x=21, y=2.33, label="2019", col="dodgerblue3", size=3, hjust=0)+
  annotate("text", x=21, y=1.33, label="2010", col="darkolivegreen3", size=3, hjust=0)+
  coord_trans(x=trans_cuts(dbrk_surf(100)("t2emp17_100k")))+
  ylab("Price ratio (5% lowest)")+
  xlab("Time to 100k (r5, transit)")+
  theme_minimal()
graph2svg(gp2level, "Price ratio 2019 et 2010", height = 12, width=25)

gggg <- gp+ggtitle("Euclidean")+gp2+ggtitle("Job Accessibility")
ggsave(plot=gggg, ".svg/comp distance dprix 2019 10.svg", width=25, height = 17, units="cm")

ggsave(plot=gp, ".svg/dist275 dprix par an.svg", width=25, height = 17, units="cm")
ggsave(plot=gp2, ".svg/t2emp17_100k dprix par an.svg", width=25, height = 17, units="cm")

gp <- ggplot(dv3f)+geom_quantogram(aes(x=dist275,  mass=exp(lvm), col=date_ma), cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE)+
  scale_y_continuous(labels=f2si2)+
  ylab("price per m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  facet_wrap(~date_ma)
gp2 <- ggplot(dv3f)+geom_quantogram(aes(x=t2emp17_100k,  mass=exp(lvm), col=date_ma), cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE)+
  scale_y_continuous(labels=f2si2)+
  ylab("price per m²")+
  xlab("Time to 100k (r5, transit)")+
  facet_wrap(~date_ma)

ggsave(plot=gp, ".svg/dist275 prix par an.svg", width=25, height = 17, units="cm")
ggsave(plot=gp2, ".svg/emp17_100k prix par an.svg", width=25, height = 17, units="cm")

gp <- ggplot(dv3f)+geom_quantogram(aes(x=dist275,  mass=valeurfonc, col=date_ma), cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE)+
  scale_y_continuous(labels=f2si2)+
  ylab("value of transaction")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  facet_wrap(~date_ma)
ggsave(plot=gp, ".svg/dist275 valeur foncière par an.svg", width=25, height = 17, units="cm")

ggplot(dv3f)+geom_quantogram(aes(x=ttr5dvf_emp09_50k,  mass=lvm, col=date_ma), lines=FALSE, cuts = dbrk_surf(100), probs=c(0.5, 0.75))+facet_wrap(~date_ma)
ggplot(dv3f)+geom_quantogram(aes(x=tcardvf_emp09_500k,  mass=lvm, col=date_ma), cuts = dbrk_surf(100), probs=c(0.5, 0.75))+facet_wrap(~date_ma)
ggplot(dv3f)+geom_quantogram(aes(x=ttr5dvf_emp09_1M,  mass=lvm, col=date_ma), cuts = dbrk_surf(100), probs=c(0.5, 0.75))+facet_wrap(~date_ma)

ggplot(dv3f %>% filter(date_ma==2019))+
  geom_quantogram(aes(x=ttr5dvf_emp09_50k,  measure=(lvm)), cuts = dbrk_surf(100), lines=FALSE)

ggplot(dv3f %>% filter(date_ma==2019))+
  geom_quantogram(aes(x=dist275,  mass=(lvm)), lines=FALSE, col="red", bins=50)+
  geom_quantogram(aes(x=dist275,  mass=(lvm)), lines=FALSE, col="blue", cuts=dbrk_ind(50))+
  geom_quantogram(aes(x=dist275,  mass=(lvm)), lines=FALSE, col="green", cuts=dbrk_surf(50))

ggplot(dv3f %>% filter(date_ma==2019))+
  geom_quantogram(aes(x=ttr5dvf_emp09_50k,  mass=(lvm)), lines=FALSE, col="red", bins=50)+
  geom_quantogram(aes(x=ttr5dvf_emp09_50k,  mass=(lvm)), lines=FALSE, col="blue", cuts=dbrk_ind(50))+
  geom_quantogram(aes(x=ttr5dvf_emp09_50k,  mass=(lvm)), lines=FALSE, col="green", cuts=dbrk_surf(50))

distances <- c("ttr_to50k_emp09", "ttr_to1M_emp09", "dist275")
pdv3f <- cut_pivot(dv3f %>% st_drop_geometry() %>% select(all_of(c(distances,"log_prix", "valeurfonc", "surface", "date_ma"))), c851, distances, n=50)
gprix <- ggplot(pdv3f)+geom_quantogram(aes(x=distance,  mass=log_prix, col=metric), lines=TRUE)+
  facet_grid(row=vars(metric), cols=vars(date_ma), scales="free_x")
gvaleur <- ggplot(pdv3f)+geom_quantogram(aes(x=distance,  mass=valeurfonc, col=metric), lines=TRUE)+
  scale_y_continuous(labels=f2si2)+
  facet_grid(row=vars(metric), cols=vars(date_ma), scales="free_x")
gsurface <- ggplot(pdv3f)+geom_quantogram(aes(x=distance,  mass=surface, col=metric), lines=FALSE)+
  scale_y_continuous(labels=f2si2)+
  facet_grid(row=vars(metric), cols=vars(date_ma), scales="free_x")

ggsave("prix versus distance.svg", g, width=24, height=17, units="cm")

ggplot(pdv3f)+geom_quantogram(aes(x=distance, mass=lvm, col=date_ma), lines=TRUE, bars=FALSE)+
  scale_color_discrete_diverging()+
  facet_wrap(vars(metric), scales="free_x")

ggplot(pdv3f %>% filter(metric=="ttr5dvf_emp09_1M"))+
  geom_quantogram(aes(x=distance, mass=lvm, col=date_ma, fill=date_ma), lines=TRUE, bars=FALSE)+
  scale_color_discrete_diverging()

# fichiers fonciers, surfaces et destinations

# locaux 2019 ---------------------
locaux.ff2019 <- vroom("{localdata}/locaux.ff2019.csv" %>% glue ,
                       col_types = cols_only(idbat="c", idlocal="c", idcom="c",
                                             dteloc="i", dteloctxt="c", ccodep="c", slocal="d",
                                             cconlc="c", dnatlc="c", typeact="c", actvac="c", fburx="c", cconac="c", cconactxt="c",
                                             stoth="d", stotd="d", stotdsueic="d", sprincp="d", habitat="c", ccthp="c",
                                             ssecp="d", ssecncp="d", sparkp="d", sparkncp="d", X="d", Y="d"), n_max=Inf)
locaux.ff2019 <- locaux.ff2019 %>%
  drop_na(X,Y)

locaux.ff2019 <- locaux.ff2019 %>%
  mutate(idINS50=idINS3035(X,Y, resolution=50),
         idINS200=idINS3035(X,Y, resolution=200))

activite.19 <- locaux.ff2019 %>%
  select(actvac, idINS50, sprincp, typeact, dnatlc, cconac, cconactxt, fburx, idcom, habitat, ccthp, X, Y) %>%
  filter(sprincp>0) %>%
  mutate(act = is.na(actvac),
         cconac2 = str_sub(cconac,1,2),
         nonaf=is.na(cconac),
         dep = str_sub(idcom,1,2),
         typeact3 =str_sub(typeact, 1,3)) %>%
  st_as_sf(coords=c("X", "Y"), crs=3035)
rm(locaux.ff2019)

# locaux 2018 -----------------
locaux.ff2018 <- vroom("{csvsources}/locaux.ff2018.csv" %>% glue ,
                       col_types = cols_only(idbat="c", idlocal="c", idcom="c",
                                             dteloc="i", dteloctxt="c", ccodep="c", slocal="d",
                                             cconlc="c", dnatlc="c", typeact="c", actvac="c",
                                             loghvac="c", fburx="c", cconac="c", cconactxt="c", logh="c",
                                             proba_rprs = "c", hlmsem="c", dnatlc = "c", loghlls = "c",
                                             stoth="d", stotd="d", stotdsueic="d", sprincp="d", habitat="c", ccthp="c",
                                             ssecp="d", ssecncp="d", sparkp="d", sparkncp="d", X="d", Y="d", st_geomloc = "c"), n_max=Inf)

XY <- locaux.ff2018$st_geomloc |> str_split(pattern = '\\|', n=2) |> unlist() |> as.numeric() |> matrix(ncol=2, byrow=TRUE)
XY <- sf_project(from = st_crs(2154), to = st_crs(3035), pts = XY)
locaux.ff2018 <- locaux.ff2018 %>%
  mutate(X=XY[,1], Y=XY[,2]) |>
  drop_na(X,Y) |>
  mutate(idINS50=idINS3035(X,Y, resolution=50),
         idINS200=idINS3035(X,Y, resolution=200))

locaux.18.m <- cbind(locaux.ff2018$X, locaux.ff2018$Y)
dd <- rdist::cdist(X = locaux.18.m, bs %>% st_coordinates()) / 1000
tr_r5_2020 <- qs::qread(str_c(localdata,"tr_r5_2020_dvf.rda"))[, .(idINS50=idINS3035(x,y,resolution=50), EMP09, temps)] [,.(EMP09=mean(EMP09)), by=c("temps","idINS50")]
tr_r5_2020 <- dcast(tr_r5_2020, idINS50~temps, value.var="EMP09")
setnames(tr_r5_2020, c(str_c(1:120)), str_c("iso",1:120, "m"))
trr5emp09 <- (qs::qread(str_c(localdata, "tr_r5_2020.rda"))$EMP09 %>% iso2time(c(100,1000)*1000)) %>% r2dt(resolution=200)
setnames(trr5emp09, c("to100k", "to1M"), c("t2emp09_100k", "t2emp09_1M"))
trr5emp17 <- qs::qread(str_c(localdata, "tr_r5emp17_2020.rda"))$emp17 %>% iso2time(c(100,1000)*1000) %>% r2dt(resolution = 200)
setnames(trr5emp17, c("to100k", "to1M"), c("t2emp17_100k", "t2emp17_1M"))
trr5pop <- qs::qread(str_c(localdata,"popc200_2020.rda"))$pop %>% iso2time(c(100,1000)*1000) %>% r2dt(resolution = 200)
setnames(trr5pop, c("to100k", "to1M"), c("t2pop_100k", "t2pop_1M"))

locaux.18 <- locaux.ff2018 %>%
  left_join(trr5emp17 %>% as_tibble(), by="idINS200") %>%
  left_join(trr5emp09 %>% as_tibble(), by="idINS200") %>%
  left_join(trr5pop %>% as_tibble(), by="idINS200") %>%
  left_join(c851 %>% as_tibble() %>% select(-t2emp09_100k,-t2emp09_1M,-t2pop_100k,-t2pop_1M,-t2emp17_100k,-t2emp17_1M), by="idINS200")

locaux.18 <- locaux.18 %>%
  mutate(
    actnv = is.na(actvac),
    lognv = is.na(loghvac),
    cconac2 = str_sub(cconac,1,2),
    nonaf=is.na(cconac),
    dep = str_sub(idcom,1,2),
    typeact3 =str_sub(typeact, 1,3))

logcat <- c("Owner Occ.", "Sec. Home", "Tenant", "Low rent",
            "Vacant LR", "Vacant H.",
            "Others", "Others T.", "Outbuilding", "Business")
actcat <- c("Housing", "Outbuilding", "Health&Education",
            "Services", "Industrial",
            "Offices", "Vacant", "Others")
scat <- c("Housing", "Vacant H.", "Outbuilding", "Business",
          "Vacant B.", "Others")

locaux.18 <- locaux.18 %>%
  mutate(
    hab    = stoth>0,
    dep    = stoth==0&sprincp==0&stotd>0,
    act    = sprincp>0,
    surface = case_when(
      hab ~ stoth,
      dep ~ stotd,
      act ~ sprincp),
    hlm    = hab&loghlls%in%c("OUI", "OUI PROBABLE"),
    ccthp  = if_else(is.na(ccthp), "NA", ccthp),
    hlmsem  = if_else(is.na(hlmsem), "NA", hlmsem),
    prop   = ccthp=="P",
    loc    = ccthp%in%c("L", "G", "F"),
    locaut = ccthp%in%c("B", "R", "U", "X")) %>%
  mutate(
    type_s = case_when(
      hab&!lognv ~ "Vacant H.",
      hab&lognv ~ "Housing",
      act&!actnv ~ "Vacant B.",
      act&actnv ~ "Business",
      dep ~ "Outbuilding",
      TRUE ~ "Others"
    ),
    type_a = case_when(
      hab ~ "Housing",
      dep ~ "Outbuilding",
      !actnv ~ "Vacant",
      typeact3%in%c("HOT", "MAG", "SPE") ~ "Services",
      typeact3%in%c("ATE", "IND", "DEP") ~ "Industrial",
      typeact3%in%c("CLI", "ENS") ~ "Health&Education",
      typeact3%in%c("BUR") ~ "Offices",
      TRUE ~ "Others"),
    type_h = case_when(
      hab&prop&proba_rprs!="RS" ~ "Owner Occ.",
      hab&!hlm&loc ~ "Tenant",
      hab&!hlm&locaut ~ "Others T.",
      hab&prop&proba_rprs=="Sec. Home" ~ "Sec. Home",
      hlm ~ "Low rent",
      hab&!lognv&!hlmsem%in%c(5,6) ~ "Vacant",
      hab&!lognv&hlmsem%in%c(5,6) ~ "Vacant LR",
      hab ~ "Others",
      dep ~ "Outbuilding",
      TRUE ~ "Business")) %>%
  mutate(
    type_s = factor(
      type_s,
      levels = scat),
    type_a = factor(
      type_a,
      levels = actcat),
    type_h = factor(
      type_h,
      levels = logcat))

m851 <- locaux.18 %>%
  filter(hab) %>%
  group_by(idINS200) %>%
  summarize(surface=sum(surface)/0.04)
mm851 <- c851 %>%
  mutate(dh = Men/0.04) %>%
  left_join(m851, by="idINS200")

gd <- ggplot(mm851)+
  geom_quantogram(aes(x=dist275, y=after_stat(median), mass=dh),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, col="blue")+
  annotate("text", x=12.5, y=15000, label="Household density (per km²)", col="blue", size=3, hjust=0)+
  geom_quantogram(aes(x=dist275, y=after_stat(median), mass=surface/100),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, col="green")+
  annotate("text", x=20, y=800, label="m² of housing density (per km²)", col="green", size=3, hjust=1)+
  scale_y_continuous(labels=f2si2, sec.axis = sec_axis(~.*100, labels=f2si2))+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  ylab("Density")+
  theme_minimal(base_size=9)

gd2 <- ggplot(mm851)+
  geom_quantogram(aes(x=dist275, y=after_stat(median), mass=surface/dh),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, col="orange")+
  scale_y_continuous(labels=f2si2)+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  ylab("Average size of housing per household")+
  theme_minimal(base_size=9)

dd <- gd/gd2
graph2svg(dd, "dist275 household density")

ggd <- ggplot(mm851)+
  geom_quantogram(aes(x=t2emp17_100k, y=after_stat(median), mass=dh),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, col="blue")+
  annotate("text", x=20, y=15000, label="Household density (per km²)", col="blue", size=3, hjust=0)+
  geom_quantogram(aes(x=t2emp17_100k, y=after_stat(median), mass=surface/100),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, col="green")+
  annotate("text", x=20, y=800, label="m² of housing density (per km²)", col="green", size=3, hjust=1)+
  scale_y_continuous(labels=f2si2, sec.axis = sec_axis(~.*100, labels=f2si2))+
  xlab("Time to 100k (r5, transit)")+
  ylab("Density")+
  theme_minimal(base_size=9)

ggd2 <- ggplot(mm851)+
  geom_quantogram(aes(x=t2emp17_100k, y=after_stat(median), mass=surface/dh),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, col="orange")+
  scale_y_continuous(labels=f2si2)+
  xlab("Time to 100k (r5, transit)")+
  ylab("Average size of housing per household")+
  theme_minimal(base_size=9)

ddd <- ggd/ggd2
ggsave(plot=ddd, ".svg/t2emp17_100k household density.svg", width=25, height = 17, units="cm")

logcol <- sequential_hcl(8,"Blues 2")
names(logcol) <- logcat[1:8]
logcol <- c(logcol, "Outbuilding"="burlywood1", "Business"=sequential_hcl(6,"Greens 3")[[3]])
actcol <- sequential_hcl(6,"Greens 3")
names(actcol) <- actcat[3:8]
actcol <- c("Housing"=sequential_hcl(8,"Blues 2")[[4]], "Outbuilding"="burlywood1", actcol)
scol <- c("Housing"=sequential_hcl(8,"Blues 2")[[1]], "Vacant H."=sequential_hcl(8,"Blues 2")[[6]], "Outbuilding"="burlywood1",
          "Business"=sequential_hcl(6,"Greens 3")[[1]], "Vacant B."=sequential_hcl(6,"Greens 3")[[5]],
          "Others"=sequential_hcl(6,"Greens 3")[[6]])

# locaux.18.long <- locaux.18 %>%
#   cut_pivot(c851, distances=c("dist275", "t2emp17_100k", "t2emp17_1M", "t2emp09_100k", "t2emp09_1M"), n=100)

g_f <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_s, col=type_s),
                 alpha=0.5, position="fill", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=scol,aesthetics = c("fill", "color"))+
  ylab("Share of built m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=percent)+
  theme_minimal(base_size=9)
g_s <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_s, col=type_s),
                 alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=scol,aesthetics = c("fill", "color"))+
  ylab("Sum of built m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
gg <- g_f+g_s+plot_layout(guides="collect")
ggsave(plot=gg, ".svg/dist275 structure simple.svg", width=25, height = 17, units="cm")

g_f <- ggplot(locaux.18)+
  geom_massogram(aes(x=t2emp17_100k, mass=surface, fill=type_s, col=type_s),
                 alpha=0.5, position="fill", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=scol,aesthetics = c("fill", "color"))+
  ylab("Share of built m²")+
  xlab("Time to 100k (r5, transit)")+
  scale_y_continuous(labels=percent)+
  theme_minimal(base_size=9)
g_s <- ggplot(locaux.18)+
  geom_massogram(aes(x=t2emp17_100k, mass=surface, fill=type_s, col=type_s),
                 alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=scol,aesthetics = c("fill", "color"))+
  ylab("Sum of built m²")+
  xlab("Time to 100k (r5, transit)")+
  scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
gg <- g_f+g_s+plot_layout(guides="collect")
ggsave(plot=gg, ".svg/t2emp17_100k structure simple.svg", width=25, height = 17, units="cm")

loc200 <- locaux.18 %>%
  group_by(idINS200) %>%
  summarize(Housing = sum(surface*(type_s%in%c("Housing", "Vacant H.", "Outbuilding")), na.rm=TRUE),
            Activity = sum(surface*(!type_s%in%c("Housing", "Vacant H.", "Outbuilding")), na.rm=TRUE)) %>%
  filter(idINS200 %in% c851$idINS200) %>%
  mutate(x = idINS2point(idINS200)[, 1],
         y = idINS2point(idINS200)[, 2])

iris <- load_DVF("iris15") %>% filter(UU2010=="00851") %>% st_union()
mg <- tm_shape(iris)+tm_fill(col="grey80")+
  tm_shape((loc200 %>% dt2r(resolution=200))/0.04/1000)+tm_raster(style="kmeans", n=10, title="m² density ('000/km²)")

mg <- ggplot(iris)+
  geom_sf(fill="grey97", color=NA)+
  geom_Raster(dt2r(loc200), aes(fill=value), long=TRUE, style="kmeans", k=10)+
  scale_fill_discrete_sequential(palette="Terrain 2", na.value=NA, rev=TRUE)+
  geom_sf(data=iris,fill="transparent", col="black", size=0.1)+
  theme_void(base_size = 9)+
  labs(fill="m² density\n1000m²/km²\nkmeans")+
  theme(legend.title=element_text(size=7), legend.text=element_text(size=7))+
  facet_wrap(~variable)
ggsave(plot=mg, ".svg/dist275 map housing activity.svg", width=25, height = 17, units="cm")

mg2 <- ggplot(loc200 %>% mutate(Housing = if_else(Housing<=quantile(Housing, 0.999), Housing, quantile(Housing, 0.999))))+
  geom_tile(aes(x=x, y=y, fill=Housing))+
  labs(fill="density")+
  scale_fill_continuous_sequential(palette="Terrain 2", na.value=0, rev=FALSE)+
  coord_fixed()
rayshader::plot_gg(mg2, raytrace = TRUE, width = 3.5, multicore = TRUE, windowsize = c(1800, 1800))
rayshader::render_snapshot(clear = TRUE)

g_f <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_a, col=type_a),
                 alpha=0.5, position="fill", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=actcol,aesthetics = c("fill", "color"), name="Usage")+
  ylab("Share of built m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=percent)+
  theme_minimal(base_size=9)

g_f_r5 <- ggplot(locaux.18)+
  geom_massogram(aes(x=t2emp17_100k, mass=surface, fill=type_a, col=type_a),
                 alpha=0.5, position="fill", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=actcol,aesthetics = c("fill", "color"), name="Usage")+
  ylab("Share of built m²")+
  xlab("Time to 100k (r5, transit)")+
  scale_y_continuous(labels=percent)+
  theme_minimal(base_size=9)
ggsave(plot=g_f+ggtitle("Euclidean Distance")+g_f_r5+ggtitle("Job Accessibility")+plot_layout(guides = "collect"), ".svg/comp distance structure act.svg", width=25, height = 17, units="cm")

g_s <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_a, col=type_a),
                 alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=actcol,aesthetics = c("fill", "color"))+
  ylab("Sum of built m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
gg <- g_f+g_s+plot_layout(guides="collect")
ggsave(plot=gg, ".svg/dist275 structure actvitivité.svg", width=25, height = 17, units="cm")

g_f <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_h, col=type_h),
                 alpha=0.5, position="fill", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=logcol,aesthetics = c("fill", "color"), name="Usage")+
  ylab("Share of built m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=percent)+
  theme_minimal(base_size=9)
g_f_r5 <- ggplot(locaux.18)+
  geom_massogram(aes(x=t2emp17_100k, mass=surface, fill=type_h, col=type_h),
                 alpha=0.5, position="fill", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=logcol,aesthetics = c("fill", "color"), name="Usage")+
  ylab("Share of built m²")+
  xlab("Time to 100k (r5, transit)")+
  scale_y_continuous(labels=percent)+
  theme_minimal(base_size=9)
g_s <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_h, col=type_h),
                 alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=logcol,aesthetics = c("fill", "color"), name="Usage")+
  ylab("Sum of built m²")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
gg <- g_f+g_s+plot_layout(guides="collect")
ggg <- g_f+ggtitle("Euclidean Distance")+g_f_r5+ggtitle("Job Accessibility")+plot_layout(guides="collect")
ggsave(plot=gg, ".svg/dist275 structure logement.svg", width=25, height = 17, units="cm")
ggsave(plot=ggg, ".svg/comp distance structure logement.svg", width=25, height = 17, units="cm")

g_s <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_s, col=type_s),
                 alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=FALSE)+
  scale_fill_manual(values=scol,aesthetics = c("fill", "color"))+
  ylab("Sum of built m²")+
  xlab("Distance to Notre Dame (% of max)")+
  scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
g_s_t <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275, mass=surface, fill=type_s, col=type_s),
                 alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_fill_manual(values=scol,aesthetics = c("fill", "color"))+
  ylab("Sum of built m²")+
  xlab("centiles of ground surface per distance to Notre Dame")+
  scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
gg <- g_s+g_s_t+plot_layout(guides="collect")
ggsave(plot=gg, ".svg/dist275 mystery.svg", width=25, height = 8.5, units="cm")

## symetry mystery  ----------

gd275ns0 <- ggplot(locaux.18 )+
  geom_massity(aes(x=dist275_NS0, mass=surface, fill=type_a, col=type_a),
               alpha=0.5, position="stack", adjust=2)+
    scale_y_continuous(labels=f2si2)+scale_color_discrete_sequential("Terrain 2", aesthetics = c("fill", "color"))+
  theme_minimal(base_size=9)

gd275ns0 <- ggplot(locaux.18)+
  geom_massogram(aes(x=dist275_NS0, mass=slocal, fill=type_s, col=type_s), alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_y_continuous(labels=f2si2, limits=c(0,10e+6))+
  scale_color_manual(values=scol, aesthetics = c("fill", "color"), name="Usage")+
  scale_x_continuous(breaks=-4:4*20)+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  ylab("Sum of built m²")+
  theme_minimal(base_size=9)

gemp17ns0 <- ggplot(locaux.18)+
  geom_massogram(aes(x=t2emp17_100k_NS0, mass=slocal, fill=type_s, col=type_s), alpha=0.5, position="stack", cuts=dbrk_surf(100), trans=TRUE)+
  scale_y_continuous(labels=f2si2, limits=c(0,10e+6))+
  scale_color_manual(values=scol, aesthetics = c("fill", "color"), name="Usage")+
  scale_x_continuous(breaks=-4:4*20)+
  xlab("Time to 100k (r5, transit)")+
  ylab("Sum of built m²")+
  theme_minimal(base_size=9)

ggg <- gd275ns0+ggtitle("Euclidean")+gemp17ns0+ggtitle("Job Accessibility")+plot_layout(guides="collect")
ggsave(plot=ggg, ".svg/comp distance NS1.svg", width=25, height = 17, units="cm")

gd275ns1 <- ggplot(locaux.18 )+
  geom_massity(aes(x=dist275_NS1, mass=surface, fill=type_a, col=type_a), alpha=0.5, position="stack", adjust=2)+
    scale_y_continuous(labels=f2si2, limits=c(0,10e+6))+
  theme_minimal(base_size=9)
gd275ns0_f <- ggplot(locaux.18 )+
    geom_massity(aes(x=dist275_NS0, mass=surface, fill=type_a, col=type_a), alpha=0.5, position="fill", adjust=2)+
    scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)
gd275ns1_f <- ggplot(locaux.18 )+
    geom_massity(aes(x=dist275_NS1, mass=surface, fill=type_a, col=type_a), alpha=0.5, position="fill", adjust=2)+
    scale_y_continuous(labels=f2si2)+
  theme_minimal(base_size=9)

ggsave(plot=gd275ns0, ".svg/dist275 NS0 graph.svg", width=25, height = 17, units="cm")
ggsave(plot=gd275ns1, ".svg/dist275 NS1 graph.svg", width=25, height = 17, units="cm")
ggsave(plot=gd275ns0_f, ".svg/dist275 NS0 f graph.svg", width=25, height = 17, units="cm")
ggsave(plot=gd275ns1_f, ".svg/dist275 NS1 f graph.svg", width=25, height = 17, units="cm")

## test data ---------------

test_data <- tibble(X=round(runif(1000000,min = 3697551 , max= 3844185)), Y=round(runif(1000000,min = 2806553  , max= 2935780)), mass = round(exp(rnorm(1000000, mean=5, sd=0.1))))
test_data <- test_data %>% mutate(idINS200 = idINS3035(X,Y)) %>% left_join(c851, by="idINS200") %>% drop_na(dist275)
trans <- ggplot(test_data )+
  geom_massogram(aes(x=dist275, mass=mass), alpha=0.5, cuts=dbrk_surf(100), trans=TRUE, lines=FALSE, col="blue", show.legend=FALSE)+
  theme_minimal(base_size=9)+
  scale_y_continuous(limits =c(0, 3.5e+5), labels=f2si2)+
  ylab("Sum of mass")+
  xlab("Distance to Notre Dame (% max)")
notrans <- ggplot(test_data )+
  geom_massogram(aes(x=dist275, mass=mass), alpha=0.5, cuts=dbrk_surf(100), trans=FALSE, lines=FALSE, col="blue", show.legend=FALSE)+
  theme_minimal(base_size=9)+
  scale_y_continuous(limits =c(0,3.5e+5), labels=f2si2)+
  ylab("Sum of mass")+
  xlab("Centiles of ground surface per distance to Notre Dame")

gg <- notrans+trans
ggsave(plot=gg, ".svg/dist275 simple mystery.svg", width=21, height = 7.5, units="cm")

# parcelles 2018 -----------------------

parcelles.ff2018 <- vroom("{csvsources}/parcelles.ff2018.csv" %>% glue,
                          col_types = cols_only(idpar="c", idparref="c",
                                                ctpdl="c", gpdl="c", pdlmp="c",
                                                nlogh="d", nlocal="d",
                                                jannatmax="d", jannatmin="d",
                                                nloghloue="d",nloghpp="d", nloghvac2a="d",
                                                nloghlm="d",ssuf="d", stoth="d", slocal="d", X="d", Y="d"))

setDT(parcelles.ff2018)

par_ref <- parcelles.ff2018[pdlmp=="R"]

parcelles_copro <- merge(parcelles.ff2018, par_ref, by.x="idparref", by.y="idpar", suffixes=c("", ".mp"), all.x=TRUE)
parcelles_copro[, part_par:=fifelse(!is.na(idparref)&!is.na(ssuf.mp)&ssuf.mp>0, ssuf/ssuf.mp, 1)]
parcelles_copro[, `:=`(nlogh = ifelse(is.na(idparref), nlogh, nlogh.mp*part_par),
                       nlocal = ifelse(is.na(idparref), nlocal, nlocal.mp*part_par),
                       nloghloue = ifelse(is.na(idparref), nloghloue, nloghloue.mp*part_par),
                       nloghpp = ifelse(is.na(idparref), nloghpp, nloghpp.mp*part_par),
                       nloghvac2a = ifelse(is.na(idparref), nloghvac2a, nloghvac2a.mp*part_par),
                       nloghlm = ifelse(is.na(idparref), nloghlm, nloghlm.mp*part_par),
                       jannatmin = ifelse(is.na(idparref), jannatmin, jannatmin.mp),
                       jannatmax = ifelse(is.na(idparref), jannatmax, jannatmax.mp),
                       ssuf = ifelse(is.na(idparref), ssuf, ssuf.mp*part_par),
                       stoth = ifelse(is.na(idparref), stoth, stoth.mp*part_par),
                       slocal = ifelse(is.na(idparref), slocal, slocal.mp*part_par)) ]

parcelles_copro[, `:=`(hlm = ifelse(nlogh>0, nloghlm/nlogh, NA_real_),
                       vacant =ifelse(nlogh>0, nloghvac2a/nlogh, NA_real_),
                       prop =ifelse(nlogh>0, 1-nloghloue/nlogh, NA_real_),
                       ageb_min = ifelse(jannatmin>0, jannatmin, NA_real_),
                       ageb_max = ifelse(jannatmax>0, jannatmax, NA_real_),
                       hsurl = stoth/ssuf,
                       ssurl = slocal/ssuf)]

parcelles_baties <- parcelles_copro %>%
  as_tibble() %>%
  drop_na(X,Y) %>%
  filter(ssuf>0, slocal>0) %>%
  mutate(idINS200 = idINS3035(sf_project(pts = cbind(.data$X, .data$Y), from=st_crs(2154), to=st_crs(3035)),
                               resolution=200)) %>%
  left_join(mm851, by="idINS200")

gp <- ggplot(parcelles_baties)+
  geom_quantogram(aes(x=dist275,  mass=hsurl), cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, color="blue")+
  geom_quantogram(aes(x=dist275,  mass=ssurl), cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, color="green")+
  scale_y_continuous(labels=f2si2)+
  annotate("text", x=11, y=3.5, label="All built/Land", col="green", size=2, hjust=0)+
  annotate("text", x=10, y=0.5, label="Housing/Land", col="blue", size=2, hjust=1)+
  ylab("H/L")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  theme_minimal(base_size=9)
ggsave(plot=gp, ".svg/dist275 H sur L.svg", width=25, height = 17, units="cm")

idpars <- dv3f$l_idpar %>% str_extract_all("(?<=\\{|\\,)[:alnum:]*(?=\\}|\\,)", simplify=TRUE)
ii <- map(1:ncol(idpars), ~match(idpars[,.x], parcelles_baties$idpar))
ii <- do.call(cbind, ii)
dv3f <- dv3f %>%
  mutate(
    hsurl = rowMeans(apply(ii, 2, function(x) parcelles_baties$hsurl[x]), na.rm=TRUE),
    ssurl = rowMeans(apply(ii, 2, function(x) parcelles_baties$ssurl[x]), na.rm=TRUE))

gp <- ggplot(dv3f)+
  geom_quantogram(aes(x=dist275,  mass=hsurl),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, color="blue")+
  annotate("text", x=10, y=0.5, label="Housing/Land", col="blue", size=2, hjust=1)+
  geom_quantogram(aes(x=dist275,  mass=ssurl),
                  cuts = dbrk_surf(100), lines=FALSE, probs=c(0.5), show.legend = FALSE, trans=TRUE, color="green")+
  annotate("text", x=13, y=3.5, label="All built/Land", col="green", size=2, hjust=0)+
  scale_y_continuous(labels=f2si2)+
  ylab("H/L")+
  xlab("Centiles of ground surface per distance to Notre Dame")+
  theme_minimal(base_size=9)
ggsave(plot=gp, ".svg/dist275 H sur L dv3f.svg", width=25, height = 17, units="cm")

library(ggpmisc)
gg <- ggplot(dv3f %>% filter(inseenotaires, hsurl>0, abs(hsurl/ssurl-1)<0.3), aes(x=log(ssurl), y=lvm, col=date_ma))+
  geom_hex(show.legend=FALSE, binwidth=c(0.1, 0.05), col=NA)+
  scale_fill_continuous_sequential(palette="OrRd", alpha=0.75, rev=TRUE, begin=0.2)+
  geom_smooth(method="lm", show.legend=FALSE, size=0.5)+
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")),
               label.x = "right", label.y = 0.97,
               eq.with.lhs = "italic(log(p))~`=`~",
               eq.x.rhs = "~italic(log(k))",
               formula = y~x, parse = TRUE, size = 2) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right", label.y = 0.9,
               formula = y~x, parse = TRUE, size = 2)+
  facet_wrap(~date_ma)+
  theme_minimal(base_size=9)
ggsave(plot=gg, ".svg/log(hl) log(p).svg", width=25, height = 17, units="cm")

data <- dv3f %>%
  filter(inseenotaires, hsurl>0, abs(hsurl/ssurl-1)<0.3) %>%
  mutate(new =if_else(is.na(ancbat), "Old" , if_else(ancbat>as.numeric(levels(date_ma)[date_ma])-10, "New", "Old")))

gg <- ggplot(data, aes(x=log(ssurl), y=lvm))+
  geom_hex(show.legend=FALSE, binwidth=c(0.2, 0.2/3), col=NA)+
  scale_fill_continuous_sequential(palette="OrRd", alpha=0.75, rev=TRUE, begin=0.2)+
  geom_smooth(method="lm", show.legend=FALSE, size=0.5)+
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")),
               label.x = "right", label.y = 0.97,
               eq.with.lhs = "italic(log(p))~`=`~",
               eq.x.rhs = "~italic(log(k))",
               formula = y~x, parse = TRUE, size = 2) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right", label.y = 0.9,
               formula = y~x, parse = TRUE, size = 2)+
  facet_wrap(~new)+
  theme_minimal(base_size=9)+
  coord_equal(3)
ggsave(plot=gg, ".svg/log(hl) log(p) Old vs New.svg", width=25, height = 12, units="cm")

# CORINE ----------------------

uu851 <- load_DVF("uu851")
coco <- load_DVF("kfoot_corine_50c")
coco <- iso2time(coco, c(0.01))
cocomap <- uu851$fdc+tm_shape(coco)+tm_raster(style="cont", palette=colorspace::sequential_hcl(palette="Terrain 2",n=20))+uu851$hdc+uu851$riv
tmap_save(tm=cocomap, filename="{DVFdata}/presentation/vv/corine.svg" %>% glue, width=25, height = 17, units="cm")

coco200 <- left_join(r2dt(coco, resolution=200), c851, by="idINS200") %>% as_tibble()

coco <- lload_DVF("kfoot_corine_50")
coco <- coco/cellStats(coco, max)
coco <- iso2time(coco, c(0.1,0.2,0.3,0.4,0.5,0.9))
tm_shape(coco)+tm_raster()
ggplot(coco200)+geom_massity(aes(x=to10m,mass=Ind,y=after_stat(density),col=dep ))
ggplot(coco200)+geom_massogram(aes(x=to10m,mass=Ind,y=after_stat(cumsum)/after_stat(totalmass),col=dep), position="identity")

# déciles ------------------------------------

deciles <- fread("{DVFdata}/rd9-CASD/kernel 2018.2/kernel osrm deciles 2018.csv" %>% glue) %>%
  as_tibble() %>%
  filter(iso=="iso1m") %>%
  mutate(
    d10_b = mean(d10/poids, na.rm=TRUE),
    d1_b = mean(d1/poids, na.rm=TRUE),
    d1_5 = d1+d2+d3+d4+d5,
    d1_5_b = mean(d1_5/poids, na.rm=TRUE),
    or10 = ifelse(d10>0, log(d10/(poids-d10))-log(d10_b/(1-d10_b)), NA),
         or10 = ifelse(d10==0, min(or10, na.rm=TRUE), or10),
         or1 = ifelse(d1>0, log(d1/(poids-d1))-log(d1_b/(1-d1_b)), NA),
         or1 = ifelse(d1==0, min(or1, na.rm=TRUE), or1),
         or5 = ifelse(d1_5>0, log(d1_5/(poids-d1_5))- log(d1_5_b/(1-d1_5_b)), NA),
         or5 = ifelse(d1_5==0, min(or5, na.rm=TRUE), or5))

deciles_r <- deciles %>%
  select(or1, or5, or10, idINS200) %>%
  dt2r(resolution=200)

uu851 <- load_DVF("uu851")
IDFOR10 <- uu851.fdc+tm_shape(deciles_r) + tm_raster(n=9, showNA=FALSE, title="Odd ratio")+
  tm_layout(panel.labels = c("First Decile (poorest)", "First five deciles", "Last Decile (richest)")) +uu851.hdc
graph2svg(IDFOR10)

deciles <- fread("{DVFdata}/rd9-CASD/kernel 2018.2/kernel osrm deciles 2018.csv" %>% glue) %>%
  as_tibble() %>%
  filter(iso=="iso1m") %>%
  left_join(c851, by="idINS200") %>%
  mutate(across(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10), ~.x/poids*Ind))

deciles_l <- deciles %>%
  pivot_longer(cols=d1:d10, names_to = "deciles", values_to = "individus") %>%
  mutate(deciles = factor(deciles,levels = str_c("d", 1:10)))

gdec_f <- ggplot(deciles_l )+
  geom_quantogram(aes(x=t2emp17_100k, mass=individus, col=deciles, fill=deciles),
                  cuts=dbrk_surf(100), trans=TRUE, se=FALSE,lines=FALSE, alpha=.5, position="fill")+
  theme_minimal()+
  xlab("Time to 100k (r5, transit)")+
  ylab("Share of individuals")

gdec_facets <- ggplot(deciles_l )+
  geom_quantogram(aes(x=t2emp17_100k, mass=individus, col=deciles, fill=deciles),
                  cuts=dbrk_surf(100), trans=TRUE, se=FALSE,lines=FALSE, alpha=.5)+
  theme_minimal()+
  facet_wrap(~deciles)+
  xlab("Time to 100k (r5, transit)")+
  ylab("Density of individuals")

gdec <- ggplot(deciles_l)+
  geom_quantogram(aes(x=t2emp17_100k, mass=individus, col=deciles, fill=deciles),
                  cuts=dbrk_surf(100), trans=TRUE, se=FALSE,lines=FALSE, alpha=.5, position="stack")+
  theme_minimal()+
  xlab("Time to 100k (r5, transit)")+
  ylab("Density of individuals")

gmen <- ggplot(deciles)+
  geom_quantogram(aes(x=t2emp17_100k, mass=Adultes/Men),
                  cuts=dbrk_surf(100), trans=TRUE, se=FALSE,lines=FALSE, alpha=.5, col="purple")+
  geom_quantogram(aes(x=t2emp17_100k, mass=(Ind-Adultes)/Men),
                  cuts=dbrk_surf(100), trans=TRUE, se=FALSE,lines=FALSE, alpha=.5, col="salmon")+
  annotate("text", x=48, y=1.85, label="Adults per household", col="purple", size=3, hjust=1)+
  annotate("text", x=48, y=0.6, label="Children per household", col="salmon", size=3, hjust=1)+
  theme_minimal()+
  xlab("Time to 100k (r5, transit)")+
  ylab("Size of household")


g <- gdec_f+gdec
graph2svg(g, "deciles IdF")
graph2svg(gdec_facets, "deciles IdF en facets")
graph2svg(gmen, "h size IdF")

deciles %>% summarize(across(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10), sum))

##### BRIERS-LOUW ET AL. 2023 #####

#Single-session spatial capture-recapture models: 2021 survey only 

##### STEP 1: PREPARE DATA #####

### PREP WORKSPACE
install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("rgdal", dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("secr", dependencies = TRUE, INSTALL_opts = '--no-lock')

require(tidyverse)
require(devtools)
require(rgdal)
require(secr)

setwd(".")

#READ IN DATA

#Read in traps. Coordinates are in UTM
traps= read.traps ("traps2021.txt", detector="count", covnames = c("HumanIndex","HabitatFineScale","PreyIndex"), binary.usage = F)

#Read in captures
captures <- read.table("captures2021.txt",
                       header = F, 
                       col.names = c("Session", "ID", "Occasion", "Detector", "Sex"))

captures$Sex <- as.factor(captures$Sex) #convert to factor: Female = 1; male = 2; unknown = NA

### MAKE CAPTURE HISTORY OBJECT
#combine captures and traps into capthist
c11.zd <- make.capthist(captures, 
                        traps,noccasions = 1,
                        covnames = c("Sex"))

summary(c11.zd, terse = TRUE)

c11.zd

### MAKE MASKS
#make masks with different sized buffers

mask1= make.mask(traps, buffer = 10000, 
            spacing = 1000, 
            type ="trapbuffer", poly = NULL)

mask2= make.mask(traps, buffer = 12000, 
                 spacing = 1000, 
                 type ="trapbuffer", poly = NULL)

mask3= make.mask(traps, buffer = 15000, 
                 spacing = 1000, 
                 type ="trapbuffer", poly = NULL)

mask4= make.mask(traps, buffer = 20000, 
                 spacing = 1000, 
                 type ="trapbuffer", poly = NULL)

##### STEP 2: FIT PRELIMINARY MODEL #####

###FIT MODEL
#Vary buffer size. Adjust number of cores (argument "ncores") as needed.
session_mod10 <- secr.fit(capthist = c11.zd, 
                        model = list(D ~ 1, 
                                     lambda0 ~ 1, 
                                     sigma ~ 1), 
                        mask = mask1, 
                        detectfn = "HHN", #Hazard half-normal
                        binomN = 0, #poisson process
                        details = list(fastproximity = FALSE), 
                        method = "Nelder-Mead", 
                        trace = F, 
                        ncores = 4, 
                        control = list(maxit = 19999))

session_mod10$fit$convergence
print(session_mod10)

session_mod12 <- secr.fit(capthist = c11.zd, 
                          model = list(D ~ 1, 
                                       lambda0 ~ 1, 
                                       sigma ~ 1), 
                          mask = mask2, 
                          detectfn = "HHN", #Hazard half-normal
                          binomN = 0, #poisson process
                          details = list(fastproximity = FALSE), 
                          method = "Nelder-Mead", 
                          trace = F, 
                          ncores = 4, 
                          control = list(maxit = 19999))

session_mod12$fit$convergence
print(session_mod12)

session_mod15 <- secr.fit(capthist = c11.zd, 
                          model = list(D ~ 1, 
                                       lambda0 ~ 1, 
                                       sigma ~ 1), 
                          mask = mask3, 
                          detectfn = "HHN", #Hazard half-normal
                          binomN = 0, #poisson process
                          details = list(fastproximity = FALSE), 
                          method = "Nelder-Mead", 
                          trace = F, 
                          ncores = 4, 
                          control = list(maxit = 19999))

session_mod15$fit$convergence
print(session_mod15)

session_mod20 <- secr.fit(capthist = c11.zd, 
                          model = list(D ~ 1, 
                                       lambda0 ~ 1, 
                                       sigma ~ 1), 
                          mask = mask4, 
                          detectfn = "HHN", #Hazard half-normal
                          binomN = 0, #poisson process
                          details = list(fastproximity = FALSE), 
                          method = "Nelder-Mead", 
                          trace = F, 
                          ncores = 4, 
                          control = list(maxit = 19999))

session_mod20$fit$convergence
print(session_mod20)

#Mixture models
session_modM1 <- secr.fit(capthist = c11.zd, 
                         model = list(D ~ 1,
                                      lambda0 ~ h2, 
                                      sigma ~ h2), 
                         mask = mask1, 
                         detectfn = "HHN", #Hazard half-normal
                         binomN = 0, #poisson process
                         details = list(fastproximity = FALSE), 
                         method = "Nelder-Mead", 
                         CL=FALSE,
                         trace = F, 
                         ncores = 4, 
                         control = list(maxit = 19999),
                         hcov=)

session_modM1$fit$convergence
print(session_modM1)

session_modM2 <- secr.fit(capthist = c11.zd, 
                          model = list(D ~ 1,
                                       lambda0 ~ h2, 
                                       sigma ~ h2), 
                          mask = mask2, 
                          detectfn = "HHN", #Hazard half-normal
                          binomN = 0, #poisson process
                          details = list(fastproximity = FALSE), 
                          method = "Nelder-Mead", 
                          CL=FALSE,
                          trace = F, 
                          ncores = 4, 
                          control = list(maxit = 19999),
                          hcov=)

session_modM2$fit$convergence
print(session_modM2)

session_modM3 <- secr.fit(capthist = c11.zd, 
                          model = list(D ~ 1,
                                       lambda0 ~ h2, 
                                       sigma ~ h2), 
                          mask = mask3, 
                          detectfn = "HHN", #Hazard half-normal
                          binomN = 0, #poisson process
                          details = list(fastproximity = FALSE), 
                          method = "Nelder-Mead", 
                          CL=FALSE,
                          trace = F, 
                          ncores = 4, 
                          control = list(maxit = 19999),
                          hcov=)

session_modM3$fit$convergence
print(session_modM3)

session_modM4 <- secr.fit(capthist = c11.zd, 
                          model = list(D ~ 1,
                                       lambda0 ~ h2, 
                                       sigma ~ h2), 
                          mask = mask4, 
                          detectfn = "HHN", #Hazard half-normal
                          binomN = 0, #poisson process
                          details = list(fastproximity = FALSE), 
                          method = "Nelder-Mead", 
                          CL=FALSE,
                          trace = F, 
                          ncores = 4, 
                          control = list(maxit = 19999),
                          hcov=)

session_modM4$fit$convergence
print(session_modM4)


##### STEP 3: EXTRACT COVARIATE VALUES #####

###EXTRACT RASTER VALUES
#Load spatial packages
library(raster)
library(sf)
library(nngeo)
library(sp)
library(ggmap)

###Read in covariate layers
#Load each separately -  for correlations or dropping types
#CRS: WGS 84 UTM zone 36S

#Land cover tifs
grass <- raster("./secr - shapefiles/Land Cover/GrasslandUTM.tif")
plot(grass)
grass@crs

other <- raster("./secr - shapefiles/Land Cover/OtherUTM.tif")
other@crs

shrubs <- raster("./secr - shapefiles/Land Cover/ShrublandUTM.tif")
shrubs@crs

trees <- raster("./secr - shapefiles/Land Cover/Trees_coverUTM.tif")
trees@crs

#Community distance
community <- st_read("./secr - shapefiles/Community 2021/CommunityUTM.shp")
st_crs(community)

###extract covariate values
msk.covs <- list() #empty list to populate with spatial points of mask with covariate values

#convert mask (15km) to spatial points
#Tibbles: convert into dataframe before input into secr
msk_pts <- mask3 %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs =32736, remove = F)

#sample covariates
msk.covs <- msk_pts %>% mutate(Grass = raster::extract(grass, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Trees = raster::extract(trees, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Shrubs = raster::extract(shrubs, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Other = raster::extract(other, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Community = simplify(st_nn(msk_pts, community, k=1, returnDist = T)[[2]]),
                               Comm_log = log(replace(Community, Community < 1, 1)))
msk.covs %>%
write_csv("./outputs/msk.covs.csv")

###scale covariates
#data frame with continuous predictors scaled to have mean 0 and std dev 1
mask.covs_st <- 
  msk.covs %>% 
  mutate_at(vars(c("Grass", "Trees", "Shrubs", "Other", "Comm_log")), function(x) scale(x)[,1])

###check correlations
cor_df <- mask.covs_st %>% 
  dplyr::select(Grass, Trees, Shrubs, Other, Comm_log)
cor_df <- st_drop_geometry(cor_df)
cor_tbl <- round(cor(cor_df), digits = 2)
cor_tbl2 <- cor(cor_df, use="pairwise.complete.obs")
write.csv(cor_tbl2, "covariate_correlations.csv")

###add covariates to mask object
covariates(mask3) <- mask.covs_st%>%
  as.data.frame()
summary(covariates(mask3))

##### STEP 4: MODEL FITTING #####

### SPECIFY MODELS
#specify competing models

#null
secr1=secr.fit(c11.zd, model=list(D~1, lambda0~1, sigma~1),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, method = "Nelder-Mead", 
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999))

#sex specific
secr2=secr.fit(c11.zd, model=list(D~1, lambda0~h2, sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE,method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific & human activity
secr3a=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+HumanIndex), sigma~h2),
               mask=mask3,detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific & prey activity
secrP1=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+PreyIndex), sigma~h2),
                mask=mask3,detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
                method = "Nelder-Mead",
                details = list(fastproximity = FALSE), 
                control = list(maxit = 19999), hcov="Sex")
#sex specific & fine scale habitat
secr4a=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+HabitatFineScale), sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific, human activity & fine scale habitat
secr5a=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+HabitatFineScale+HumanIndex), sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific: community
secr6a=secr.fit(c11.zd, model=list(D~Comm_log, lambda0~h2, sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific: tree cover
secr7a=secr.fit(c11.zd, model=list(D~Trees, lambda0~h2, sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific: tree cover + community
secr8a=secr.fit(c11.zd, model=list(D~Trees + Comm_log, lambda0~h2, sigma~h2),
                mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
                method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")
#sex specific: tree cover + community
secr9a=secr.fit(c11.zd, model=list(D~Trees*Comm_log, lambda0~h2, sigma~h2),
                mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr2,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov="Sex")


### EVALUATE RESULTS
secrmods <- secrlist(secr2, secr3a, secrP1, secr4a, secr5a, secr6a, secr7a, secr8a, secr9a)
str(secrmods)
print(secrmods)

#check convergence
print(secrmods$fit$convergence)

#check compatibility
AICcompatible(secrmods)

#check support
AIC(secrmods, criterion = "AICc", dmax = 7)

#save AIC table
AIC_tblL <- as_tibble(AIC(secrmods, criterion = "AICc", dmax = 7), rownames = "Hyp") %>% 
  mutate(model = str_sub(model, 1, -46),
         dAICc = round(dAICc, 2), Weight = round(AICcwt, 2)) %>% 
  dplyr::select(-detectfn, -logLik, -AICc) %>%
  write_csv("./outputs/AICc_tableSS.csv")

saveRDS(secrmods, file = "./outputs/secrmodelSS..rds")
save(secrmods, file = "./outputs/secrmodelsSS.RData")

#hcov
secr10=secr.fit(c11.zd, model=list(D~1, lambda0~h2, sigma~h2),
               mask=mask3, detectfn = "HHN",binomN = 0, CL=FALSE, trace=FALSE, method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov=)
#hcov & human activity
secr11a=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+HumanIndex), sigma~h2),
               mask=mask3, detectfn = "HHN",binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov=)
#hcov & prey activity
secrP2=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+PreyIndex), sigma~h2),
                 mask=mask3, detectfn = "HHN",binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
                 method = "Nelder-Mead",
                 details = list(fastproximity = FALSE), 
                 control = list(maxit = 19999), hcov=)
#hcov & fine scale habitat
secr12a=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+HabitatFineScale), sigma~h2),
                mask=mask3, detectfn = "HHN",binomN = 0, CL=FALSE, trace=FALSE, start = secr10, 
                method = "Nelder-Mead",
                details = list(fastproximity = FALSE), 
                control = list(maxit = 19999), hcov=)
#hcov, human activity & fine scale habitat
secr13a=secr.fit(c11.zd, model=list(D~1, lambda0~(h2+HabitatFineScale+HumanIndex), sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov=)
#hcov: community
secr14a=secr.fit(c11.zd, model=list(D~Comm_log, lambda0~h2, sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov=)
#hcov: tree cover
secr15a=secr.fit(c11.zd, model=list(D~Trees, lambda0~h2, sigma~h2),
               mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = FALSE), 
               control = list(maxit = 19999), hcov=)
#hcov: tree cover + community
secr16a=secr.fit(c11.zd, model=list(D~Trees + Comm_log, lambda0~h2, sigma~h2),
                mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
                method = "Nelder-Mead",
                details = list(fastproximity = FALSE), 
                control = list(maxit = 19999), hcov=)
#hcov: tree cover + community
secr17a=secr.fit(c11.zd, model=list(D~Trees*Comm_log, lambda0~h2, sigma~h2),
                mask=mask3, detectfn = "HHN", binomN = 0, CL=FALSE, trace=FALSE, start = secr10,
                method = "Nelder-Mead",
                details = list(fastproximity = FALSE), 
                control = list(maxit = 19999), hcov=)

### EVALUATE RESULTS
secrmods2 <- secrlist(secr10, secr11a, secr12a, secr13a, secr14a, secr15a, secr16a, secr17a)
str(secrmods2)
print(secrmods2)

#check convergence
print(secrmods2$fit$convergence)

#check compatibility
AICcompatible(secrmods2)

#check support
AIC(secrmods2, criterion = "AICc", dmax = 7)

#save AIC table
AIC_tbl2 <- as_tibble(AIC(secrmods2, criterion = "AICc", dmax = 7), rownames = "Hyp") %>% 
  mutate(model = str_sub(model, 1, -46),
         dAICc = round(dAICc, 2), Weight = round(AICcwt, 2)) %>% 
  dplyr::select(-detectfn, -logLik, -AICc) %>%
  write_csv("./outputs/AICc_tableNS.csv")

saveRDS(secrmods2, file = "./outputs/secrmodelsNS..rds")
save(secrmods2, file = "./outputs/secrmodelsNS.RData")

#Goodness of fit tests of top models
secrSS <- secr.test (secr4a, fit=TRUE)
secrSS

secrNS <- secr.test (secr12a, fit=TRUE)
secrNS


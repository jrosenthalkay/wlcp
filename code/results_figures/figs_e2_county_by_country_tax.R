## experiment : global warming 
## figures for The Winners and Losers of Climate Policies 
# jordan rosenthal-kay

# libraries 
library(haven)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(patchwork)

# read in data 
data <- read_csv(file.path(path,'welfare_carbontax_unilateral_i_l_50.csv'))

data <- data %>%
  mutate(across(-iso3, ~ .x * 100))

data <- data %>%
  mutate(dW_winsorized = case_when(
    dW <= quantile(dW, 0.02, na.rm = TRUE) ~ quantile(dW, 0.02, na.rm = TRUE),
    dW >= quantile(dW, 0.98, na.rm = TRUE) ~ quantile(dW, 0.98, na.rm = TRUE),
    TRUE ~ dW
  ))

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Join the data with the map
world_data <- world %>%
  left_join(data, by = c("iso_a3_eh" = "iso3"))

# Plot the map with Robinson projection
p1 <- ggplot() +
  geom_sf(data = sf::st_graticule(ndiscr = 10000, margin = 0.05), 
          color = "grey70", size = 0.05, alpha = 0.5) +
  geom_sf(data = world_data, aes(fill = dW_winsorized), color = "grey90", size = 0.1) +
  scale_fill_gradient2(
    low = "firebrick3", 
    high = "dodgerblue3",
    midpoint=0,
    limits=c(-0.75,0.75),
    breaks=c(-0.6,-0.3,0,0.3,0.6),
    name = expression(Delta * "Welfare"),
    na.value = "grey85"
  ) +
  coord_sf(crs = st_crs("+proj=robin"), # Robinson projection
           datum = NA) +
  theme_minimal() +
  labs(
    title = "Each nation's own percent change in welfare",
    subtitle = "from a unilateral $50/ton global carbon tax"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.box = "vertical",
    legend.margin = margin(t = 10),
    legend.box.margin = margin(0, 0, 0, 0)
  )
p1

savename <- 'welfare_carbontax_50usd_Wi_map.png'
ggsave(file.path(savepath,savename), plot = p1, width = 8, height = 6, dpi = 300)

# now world 

# read in data 
data2 <- read_csv(file.path(path,'welfare_carbontax_unilateral_w_l_50.csv'))

data2 <- data2 %>%
  mutate(across(-iso3, ~ .x * 100))

data2 <- data2 %>%
  mutate(dW_winsorized = case_when(
    dW <= quantile(dW, 0.02, na.rm = TRUE) ~ quantile(dW, 0.02, na.rm = TRUE),
    dW >= quantile(dW, 0.98, na.rm = TRUE) ~ quantile(dW, 0.98, na.rm = TRUE),
    TRUE ~ dW
  ))

data$dW_winsorized[data$dW_winsorized < -0.75] <- -0.75
data$dW_winsorized[data$dW_winsorized > 0.75] <- 0.75

# Join the data with the map
world_data <- world %>%
  left_join(data2, by = c("iso_a3_eh" = "iso3"))

# Plot the map with Robinson projection
p2 <- ggplot() +
  geom_sf(data = sf::st_graticule(ndiscr = 10000, margin = 0.05), 
          color = "grey70", size = 0.05, alpha = 0.5) +
  geom_sf(data = world_data, aes(fill = dW_winsorized), color = "grey90", size = 0.1) +
  scale_fill_gradient2(
    low = "firebrick3", 
    high = "dodgerblue3",
    midpoint=0,
    name = expression(Delta * "Welfare"),
    na.value = "grey85",
    limits=c(-0.75,0.75)
  ) +
  coord_sf(crs = st_crs("+proj=robin"), # Robinson projection
           datum = NA) +
  theme_minimal() +
  labs(
    title = "The global welfare gain from unilateral carbon taxation",
    subtitle = "from a unilateral $50/ton carbon tax"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.box = "vertical",
    legend.margin = margin(t = 10),
    legend.box.margin = margin(0, 0, 0, 0)
  )
p2

savename <- 'welfare_carbontax_50usd_Wu_map.png'
ggsave(file.path(savepath,savename), plot = p2, width = 8, height = 6, dpi = 300)

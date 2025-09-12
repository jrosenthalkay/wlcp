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
data <- read_csv(file.path(path,'welfare_carbontaxCHN_l_50.csv'))

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
    name = expression(Delta * "Welfare"),
    na.value = "grey85"
  ) +
  coord_sf(crs = st_crs("+proj=robin"), # Robinson projection
           datum = NA) +
  theme_minimal() +
  labs(
    title = "Percent change in welfare",
    subtitle = "from a $50/ton carbon tax in China"
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

savename <- 'welfare_carbontaxCHN_50usd_map.png'
ggsave(file.path(savepath,savename), plot = p1, width = 8, height = 6, dpi = 300)

# Filter for the specified countries
selected_countries <- c("USA", "CHN", "DEU","GBR", "JPN", "RUS", "IND", "IDN" ,"SAU","BRA","NGA","CAN","WLD")
filtered_data <- data %>%
  filter(iso3 %in% selected_countries)

# Reshape data to long format for stacking
long_data <- filtered_data %>%
  pivot_longer(
    cols = c(dW_D, dW_p, dW_e, dW_π),
    names_to = "component",
    values_to = "value"
  )

# Set the order of the components for stacking
long_data$component <- factor(
  long_data$component,
  levels = c("dW_D", "dW_p", "dW_π", "dW_e")  # Order by typical magnitude
)

# Create a nice label lookup for the components
component_labels <- c(
  "dW_D" = "Direct productivity",
  "dW_p" = "Trade",
  "dW_e" = "Energy cost",
  "dW_π" = "Energy rent"
)

# Plot a stacked bar chart
p2 <- ggplot(long_data, aes(x = iso3, y = value, fill = component)) +
  geom_col(position = "stack") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c("dW_D" = "#0072B2", "dW_p" = "#009E73", "dW_e" = "#CC79A7", "dW_π" = "#E69F00"),
    labels = component_labels,
    name = ""
  ) +
  labs(
    title = "Decomposition of welfare changes",
    subtitle = "from a $50/ton carbon tax in China",
    x = "Country",
    y = expression(Delta * "Welfare Components")
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
  ) +
  guides(fill = guide_legend(nrow = 1))

savename <- 'welfare_carbontaxCHN_50usd_bar.png'
ggsave(file.path(savepath,savename), plot = p2, width = 8, height = 6, dpi = 300)

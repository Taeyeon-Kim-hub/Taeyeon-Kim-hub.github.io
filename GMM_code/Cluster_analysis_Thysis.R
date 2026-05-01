# Cluster analysis
# Sample year check (20 years above)

# Packages
library(tidyverse)
library(mclust)
library(sf)
library(rnaturalearth)
library(gridExtra)
library(grid)
library(cowplot) 

# Data
raw_data <- read.csv("D:\\Project\\Thysis\\Data\\Data_suhyup.csv") %>% drop_na()
coord_data <- read.csv("D:\\Project\\Thysis\\Data\\Gridcell_lat_log.csv")

# Code
sf_use_s2(TRUE)
bbox <- st_bbox(c(xmin = 123, ymin = 28, xmax = 134, ymax = 38.5), crs = 4326)
world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid()
world_clip <- st_crop(world, bbox)

map_list <- list()
legend_only <- NULL

# Checking sea block data (Length 20 ~ 25)
for (y_limit in 20:25) {
  current_df <- raw_data %>%
    group_by(Gridcell_large) %>%
    filter(n_distinct(Year) >= y_limit) %>% 
    mutate(Effort = Tonnage * Fishing_day) %>%
    summarise(
      Total_Yield = sum(Yield_kg, na.rm = TRUE),
      Total_Effort = sum(Effort, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Weighted_CPUE = Total_Yield / Total_Effort) %>%
    inner_join(coord_data, by = c("Gridcell_large" = "Large_grid"))
  
  gmm_input <- current_df %>% select(Weighted_CPUE, longitude, latitude) %>% scale()
  gmm_best <- mclustBIC(gmm_input)
  mod <- Mclust(gmm_input, x = gmm_best)
  
  current_df$Original_Cluster <- mod$classification
  cluster_rank <- current_df %>%
    group_by(Original_Cluster) %>%
    summarise(Mean_CPUE = mean(Weighted_CPUE)) %>%
    arrange(Mean_CPUE) %>% 
    mutate(New_Label = row_number()) 
  
  current_df <- current_df %>%
    left_join(cluster_rank, by = "Original_Cluster") %>%
    mutate(Cluster = factor(New_Label)) 
  
  p <- ggplot() +
    geom_sf(data = world_clip, fill = "grey90", color = "grey50", linewidth = 0.3) +
    
    geom_tile(data = current_df, 
              aes(x = longitude, y = latitude, fill = Cluster), 
              alpha = 0.8, color = "black", linewidth = 0.1) +
    
    scale_fill_brewer(palette = "Set1", name = "Cluster Group") + 
    
    coord_sf(xlim = c(123, 134), ylim = c(28, 38.5), expand = FALSE) +
    
    theme_minimal() +
    theme(
      legend.position = "none",

      axis.title = element_blank(),
      axis.text = element_text(size = 11, color = "black"), 
      axis.ticks = element_line(color = "black"),
      
      panel.grid = element_line(color = "grey90", linetype = "dotted"),
      panel.background = element_rect(fill = "aliceblue", color = NA),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
    )
  map_list[[paste0("Y", y_limit)]] <- p
}
grid.arrange(grobs = map_list, ncol = 3, nrow = 2)

## Select data (above 20 years, 80% limit)
target_limit <- 20
df_final <- raw_data %>%
  group_by(Gridcell_large) %>%
  filter(n_distinct(Year) >= target_limit) %>% 
  mutate(Effort = Tonnage * Fishing_day) %>%
  summarise(
    Total_Yield = sum(Yield_kg, na.rm = TRUE),
    Total_Effort = sum(Effort, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(Weighted_CPUE = Total_Yield / Total_Effort) %>%
  inner_join(coord_data, by = c("Gridcell_large" = "Large_grid"))

gmm_input_final <- df_final %>% select(Weighted_CPUE, longitude, latitude) %>% scale()
gmm_best_final <- mclustBIC(gmm_input_final)
mod_final <- Mclust(gmm_input_final, x = gmm_best_final)
mod_final$modelName
#-------------------------------------------------------------------------------------#
# ─────────────────────────────────────────────────────────────────────────────
# [논문용] GMM Model Selection Table 추출 (AIC, BIC, ICL 포함)
# ─────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(mclust)

bic_matrix <- as.matrix(gmm_best_final)
bic_df <- as.data.frame(as.table(bic_matrix)) %>%
  rename(Clusters = Var1, Model = Var2, mclust_BIC = Freq) %>%
  drop_na() %>%
  arrange(desc(mclust_BIC))

top_5_models <- head(bic_df, 20)
table_list <- list()

for(i in 1:nrow(top_5_models)) {
  G_val <- as.numeric(as.character(top_5_models$Clusters[i]))
  model_name <- as.character(top_5_models$Model[i])
  temp_mod <- Mclust(gmm_input_final, G = G_val, modelNames = model_name)
  logLik_val <- temp_mod$loglik
  df_val <- temp_mod$df       # 파라미터 개수
  n_val <- temp_mod$n         # 샘플 사이즈
  std_AIC <- -2 * logLik_val + 2 * df_val
  std_BIC <- -2 * logLik_val + df_val * log(n_val)
  mclust_icl_val <- icl(temp_mod)
  std_ICL <- -1 * as.numeric(mclust_icl_val)
  table_list[[i]] <- data.frame(
    Rank = i,
    Model_ID = paste0("M", i),
    Model = model_name,
    Clusters = G_val,
    Parameters = df_val,
    LogLikelihood = round(logLik_val, 2),
    AIC = round(std_AIC, 2),
    BIC = round(std_BIC, 2),
    ICL = round(std_ICL, 2)
  )
}

model_selection_table <- bind_rows(table_list)
print(model_selection_table)
write.csv(model_selection_table, "GMM_Model_Selection_Table_Full.csv", row.names = FALSE)
#------------------------------------------------------------------------------#

df_final$Original_Cluster <- mod_final$classification
cluster_rank_final <- df_final %>%
  group_by(Original_Cluster) %>%
  summarise(Mean_CPUE = mean(Weighted_CPUE)) %>%
  arrange(Mean_CPUE) %>% 
  mutate(New_Label = row_number()) 

df_final <- df_final %>%
  left_join(cluster_rank_final, by = "Original_Cluster") %>%
  mutate(Cluster = factor(New_Label)) 

summary(mod_final)
print(mod_final$parameters$pro)

# Checking cluster covariance and correlation
sigmas <- mod_final$parameters$variance$sigma
num_clusters <- mod_final$G
var_names <- c("Weighted_CPUE", "Longitude", "Latitude") 

for (k in 1:num_clusters) {
  cat(paste0("\n[ Cluster ", k, " ]\n"))
  if (is.matrix(sigmas)) {
    cov_mat <- sigmas
  } else {
    cov_mat <- sigmas[, , k]
  }
  cor_mat <- cov2cor(cov_mat)
  rownames(cor_mat) <- var_names
  colnames(cor_mat) <- var_names
  print(round(cor_mat, 3))
}

# Checking each cluster's mean (Long,Lat,CPUE)
real_params <- df_final %>%
  group_by(Cluster) %>%
  summarise(
    Count = n(),                            
    Prob_Weight = n() / nrow(df_final),      
    Mean_CPUE = mean(Weighted_CPUE),        
    Mean_Lon = mean(longitude),              
    Mean_Lat = mean(latitude),               
    SD_CPUE = sd(Weighted_CPUE)              
  ) 
print(real_params)

# East (blue), West (Red), South (green)
my_colors <- c("1" = "#4575B4",   
               "2" = "#74C476",   
               "3" = "#D73027")   

ggplot() +
  geom_sf(data = world_clip, fill = "grey80", color = "grey30", linewidth = 0.3) +
  geom_tile(data = df_final, 
            aes(x = longitude, y = latitude, fill = Cluster), 
            alpha = 1,             
            color = "black",       
            linewidth = 0.05) +    
  
  scale_fill_manual(
    values = my_colors, 
    name = "Cluster group", 
    labels = c("Cluster 1", 
               "Cluster 2", 
               "Cluster 3")
  ) + 
  
  coord_sf(xlim = c(123, 134), ylim = c(28, 38.5), expand = FALSE) +
  labs(
    #title = "Spatial Stratification",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = "white", color = "grey80"), 
    legend.key.size = unit(0.8, "cm"),
    
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey80"),
    
    panel.background = element_rect(fill = "#F0F8FF", color = NA), 
    panel.grid = element_line(color = "grey80", linetype = "dotted")
  ) +
  
  guides(fill = guide_legend(override.aes = list(color = "black", linewidth = 0.2)))

#------------------------------------------------------------------------------#
# CPUE comparison
grid_map <- df_final %>% 
  select(Gridcell_large, Cluster) %>% 
  distinct()

yearly_data <- raw_data %>%
  inner_join(grid_map, by = "Gridcell_large")

cluster_trend <- yearly_data %>%
  mutate(Effort = Tonnage * Fishing_day) %>%  
  group_by(Year, Cluster) %>% 
  summarise(
    Sum_Yield_kg = sum(Yield_kg, na.rm = TRUE), 
    Sum_Effort = sum(Effort, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Sum_Yield_Ton = Sum_Yield_kg / 1000,        
    Cluster_CPUE = Sum_Yield_Ton / Sum_Effort   
  ) %>%
  arrange(Cluster, Year)

print(head(cluster_trend))

c1 <- cluster_trend %>% filter(Cluster == "1") # Blue
c2 <- cluster_trend %>% filter(Cluster == "2") # Green
c3 <- cluster_trend %>% filter(Cluster == "3") # Red

y_max <- max(cluster_trend$Cluster_CPUE) * 1.15
x_range <- range(cluster_trend$Year)

plot(x = x_range, y = c(0, y_max), 
     type = "n", 
     xlab = "Year", 
     ylab = "CPUE (MT / Effort)",
     main = "Annual CPUE trend by cluster",
     cex.main = 1.5,
     cex.lab = 1.3,
     las = 1
     ) 

#grid(col = "lightgray", lty = "dotted")

lines(c1$Year, c1$Cluster_CPUE, col = "#4575B4", lwd = 3, type = "o", pch = 19)

lines(c2$Year, c2$Cluster_CPUE, col = "#74C476", lwd = 3, type = "o", pch = 19)

lines(c3$Year, c3$Cluster_CPUE, col = "#D73027", lwd = 3, type = "o", pch = 19)

legend("topleft", 
       legend = c("Cluster 3 ", "Cluster 2 ", "Cluster 1 "),
       col = c("#D73027", "#74C476", "#4575B4"),
       lty = 1, lwd = 3, pch = 19,
       bty = "n", 
       cex = 1.2)


#-----------------------------------------------------------------------------#
# Create file (Regional CPUE and Yield)
# ─────────────────────────────────────────────────────────────────────────────
# 1. Cluster mapping
# ─────────────────────────────────────────────────────────────────────────────
grid_map <- df_final %>% 
  select(Gridcell_large, Cluster) %>% 
  distinct() %>% 
  mutate(Region = case_when(
    Cluster == 1 ~ "East",   # Blue
    Cluster == 2 ~ "South",  # Green
    Cluster == 3 ~ "West"    # Red
  ))

yearly_data <- raw_data %>%
  inner_join(grid_map, by = "Gridcell_large")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Yearly, by grid cell mapping
# ─────────────────────────────────────────────────────────────────────────────
region_trend <- yearly_data %>%
  mutate(Effort = Tonnage * Fishing_day) %>% 
  group_by(Year, Region) %>%
  summarise(
    Sum_Yield_kg = sum(Yield_kg, na.rm = TRUE), 
    Sum_Effort = sum(Effort, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Sum_Yield_Ton = Sum_Yield_kg / 1000,          
    Region_CPUE = Sum_Yield_Ton / Sum_Effort      
  )
total_trend <- yearly_data %>%
  mutate(Effort = Tonnage * Fishing_day) %>% 
  group_by(Year) %>%  
  summarise(
    Total_Yield_kg = sum(Yield_kg, na.rm = TRUE), 
    Total_Effort = sum(Effort, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Total_Yield_Ton = Total_Yield_kg / 1000,
    Total_CPUE = Total_Yield_Ton / Total_Effort 
  )
# ─────────────────────────────────────────────────────────────────────────────
# 3. CSV
# ─────────────────────────────────────────────────────────────────────────────
df_yield_save <- region_trend %>%
  select(Year, Region, Sum_Yield_Ton) %>%
  pivot_wider(names_from = Region, values_from = Sum_Yield_Ton) %>%
  arrange(Year) %>%
  select(Year, East, South, West) 

df_cpue_save <- region_trend %>%
  select(Year, Region, Region_CPUE) %>%
  pivot_wider(names_from = Region, values_from = Region_CPUE) %>%
  left_join(select(total_trend, Year, Total_CPUE), by = "Year") %>%
  rename(Total = Total_CPUE) %>% 
  arrange(Year) %>%
  select(Year, East, South, West, Total) 

write.csv(df_yield_save, "Yearly_Yield_by_Region.csv", row.names = FALSE)
write.csv(df_cpue_save, "Yearly_CPUE_by_Region.csv", row.names = FALSE)

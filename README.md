# SARS-CoV-2_CADR
This R script is designed to generate figures for the manuscript: Overview of SARS-CoV-2 Genomic Surveillance in Central America and the Dominican Republic from February 2020 to January 2023: The Impact of PAHO and COMISCA's Collaborative Efforts.

## Load required libraries
```R
library(ggplot2)
library(RColorBrewer)
library(lubridate)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library(sf)
library(raster)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggnewscale)
library(colorspace)
library(tibble)
library(tidyr)
library(conflicted)

# Prioritize dplyr functions
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("group_by", "dplyr")

# Define a combined color palette
num_required_colors <- 150 # Generate more colors than potentially needed
base_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928") # A qualitative palette like Paired from ColorBrewer

# Create a color generating function
color_generator <- colorRampPalette(base_colors)

# Generate the custom color palette
custom_colors <- color_generator(num_required_colors)
```
## Create Figure 1. Timeline of processes implemented to strengthen genomic surveillance of SARS-CoV-2 in Central America and the Dominican Republic.

```R
# Update this path to their desired location
read_path <- "/content/path/centralamerica.tsv"

# Define the activities and dates
activities <- c(
  'Assessment of local capacity needs',
  'Acquisition of sequencer INCIENSA',
  'Acquisition of sequencer ICGES',
  'Acquisition of laboratory supplies',
  'SE-COMISCA World Courier account opening',
  'Shipment El Salvador-Costa Rica',
  'Shipment Guatemala-Costa Rica',
  'Shipment Dominican Republic-Panama',
  'Shipment Guatemala-Panama',
  'Shipment Honduras-Panama',
  'Shipment El Salvador-Panama',
  'Shipment Belize-Panama',
  'Acquisition bioinformatic software',
  'Acquisition of computer equipment',
  'Workshop on Bioinformatic Analysis')

dates <- as.Date(c(
  '2021-01-01', '2021-03-22', '2021-03-22',
  '2021-03-01', '2023-04-01', '2021-08-23',
  '2021-09-20', '2021-11-24', '2021-06-01',
  '2021-06-01', '2021-07-01', '2021-07-01',
  '2022-06-15', '2022-02-05', '2022-08-08'))
# Create a data frame for plotting
dataset <- tibble(
  Activity = activities,
  Start = dates,
  End = c(
    '2021-01-31', '2021-05-21', '2021-07-24',
    '2023-02-01', '2023-08-18', '2022-08-22',
    '2022-08-22', '2022-04-30', '2022-03-22',
    '2022-09-29', '2021-08-29', '2023-01-03',
    '2022-09-23', '2022-06-23', '2022-08-12'))

# Convert to tibble and create sorting date column
dataset <- dataset %>%
  mutate(Start = as.Date(Start, "%Y-%m-%d"),
         End = as.Date(End, "%Y-%m-%d"),
         Activity = factor(Activity, levels = rev(Activity))) %>%
  arrange(desc(Start)) %>%
  filter(Start >= as.Date("2020-02-01"), End <= as.Date("2023-01-31")) # Apply date filtering

# The bimonthly grouping was causing issues with geom_segment.
# Reverting to using the original Start and End dates for plotting,
# but keeping the bimonthly breaks on the x-axis.

# Generate colors from custom_colors with an interval for activities
num_activities <- length(unique(dataset$Activity))
color_indices_activities <- floor(seq(1, length(custom_colors), length.out = num_activities))
plot_colors_activities <- custom_colors[color_indices_activities]


# Plotting Gantt chart
timeline<-ggplot(dataset, aes(x = Start, xend = End, y = Activity, yend = Activity, color = Activity)) + # Use original Start and End dates
  geom_segment(size = 10) +
  scale_color_manual(values = plot_colors_activities, guide = FALSE) + # Use plot_colors_activities manual scale
  labs(x = "Date", y = NULL) +
  theme_minimal(base_family = "Arial") + # Set base font family to Arial
  theme(panel.grid.major.x = element_line(color = "#3c3c3c", size = 0.2)) +
  theme(panel.grid.minor.x = element_line(color = "darkgray", size = 0.1)) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold", vjust=-1)) +
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12, hjust = 1)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") # Keep bimonthly breaks on x-axis

timeline
# Extract the directory path from read_path
add_path <- dirname(read_path)

# Create the directory if it doesn't exist
if (!dir.exists(add_path)) {
  dir.create(add_path, recursive = TRUE)
}

ggsave(paste0(add_path, "/Figure 1.tiff"), width=16.1, height=5.74)
ggsave(paste0(add_path, "/Figure 1.png"), width=16.1, height=5.74)
```

## Create Figure 2.  Individual SARS-CoV-2 sequences obtained by country in Central America and Do-minican Republic from February 2020 to January 2023. The top 14 most prevalent lineages were individually labeled, and all the remaining lineages were labeled as other.
```R
# Load data
data0 <- read_tsv(paste0("/content/path/centralamerica.tsv"))

# Descending country count
df_count<- data0 %>%
  count(country) %>%
  arrange(desc(n))

# Renaming data 0
df_central_america<-data0

# Convert 'date' to Date data type using lubridate::parse_date_time with multiple orders
df_central_america$date <- lubridate::parse_date_time(df_central_america$date, orders = c("mdy", "dmy"), quiet = TRUE)

# Apply date filtering
df_central_america <- df_central_america %>%
  filter(date >= as.Date("2020-02-01"), date <= as.Date("2023-01-31"))

# Grouping data by day, 2 weeks, and 1 month (keeping for potential future use or if needed elsewhere)
df_central_america$days<-as.Date(cut(df_central_america$date,breaks = "day",start.on.monday = FALSE))
df_central_america$date2<-as.Date(cut(df_central_america$date,breaks = "2 weeks",start.on.monday = FALSE))
df_central_america$date3<-as.Date(cut(df_central_america$date,breaks = "1 month",start.on.monday = FALSE))


# Create date3 column with format YYYY-MMM
df_central_america$date3 <- format(df_central_america$date, "%Y-%b")


# Ascending country count
df_count <- df_central_america %>% count(country)
names(df_count)[names(df_count) == "country"] <- "country"
names(df_count)[names(df_count) == "n"] <- "Count"


# Add Country Count on the left side of the table
df_central_america = df_central_america %>%
  left_join(df_count, by = c("country" = "country"))

# Define lineages to be individually labeled (top 14 + Other)
lineages <- c('BA.1.1','XBB.1','AY.113','A.2.5','B.1','AY.100','BA.2',
              'BA.2.9','BE.1','BA.5.2.23','BA.5.2','BQ.1.1','BA.4.6','BA.2.12.1','Other')

# Create a new column grouping lineages into 'Other' or specific lineage and convert to factor
df_central_america <- df_central_america %>%
  mutate(pango_lineage_grouped = ifelse(pango_lineage %in% lineages, pango_lineage, 'Other')) %>%
  mutate(pango_lineage_grouped = factor(pango_lineage_grouped, levels = lineages))

# Generate colors from custom_colors with an interval for lineages
num_lineages <- length(unique(df_central_america$pango_lineage_grouped))
color_indices <- floor(seq(1, length(custom_colors), length.out = num_lineages))
plot_colors <- custom_colors[color_indices]

# Plot
panelopas <- ggplot(data=df_central_america) +
  theme_classic(base_family = "Arial")+ # Set base font family to Arial
  geom_segment(aes(x=as.Date("2020-02-01"), y=reorder(country,Count), xend=as.Date("2023-02-01"), yend=country, group=country), colour="grey80", size=5) +

  geom_point(aes(fill=pango_lineage_grouped, x=days, y=reorder(country,Count)),
             position = position_jitter(height=0.2), # Removed width=0.2
             shape=21, stroke=0.05, col='grey70', size=3, na.rm = TRUE) + # Added na.rm = TRUE

  ylab('')+ xlab('month')+ ggtitle('') +

  scale_fill_manual(values = plot_colors, name = 'Pangolin Lineages') + # Use plot_colors manual scale
  theme(legend.position="bottom") +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold", vjust=-0.5)) +
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold")) +
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(panel.grid.major.x = element_line(color = "#3c3c3c", size = 0.2)) +
  theme(panel.grid.minor.x = element_line(color = "darkgray", size = 0.08)) +
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month", date_minor_breaks = "1 month", limits = as.Date(c("2020-02-01", "2023-03-01")))+ # Explicitly set x-axis limits
  guides(fill = guide_legend(override.aes = list(size=5), nrow = 2))+
  xlab('Sample Collection Date')


# Extract the directory path from read_path (assuming read_path is defined earlier)
add_path <- dirname("/content/path/centralamerica.tsv") # Define add_path locally if read_path is not globally available

ggsave(paste0(add_path, "/Figure 2.tiff"), panelopas, width=18.8, height=5.74)
ggsave(paste0(add_path, "/Figure 2.png"), panelopas, width=18.8, height=5.74)
```
## Create Figure 3. Relative percentages of SARS-CoV-2 lineages circulating in Central America and Dominican Republic from February 2020 to January 2023. Lineages with a frequency exceding 120 occurrences (n>120) per month were selected.
```R
# --- 1. Load and prepare data ---
data0 <- read_tsv(paste0("/content/path/centralamerica.tsv"))

# Print unique values and summary of the date column
print("Unique values of date column before conversion:")
print(unique(data0$date))
print("Summary of date column before conversion:")
print(summary(data0$date))

# Correct the date format string using parse_date_time with multiple orders
# Adding a check for NA dates after conversion
data0$date <- lubridate::parse_date_time(data0$date, orders = c("mdy", "dmy"), quiet = TRUE)
if (any(is.na(data0$date))) {
  warning("Some dates could not be parsed and were set to NA.")
}

# Print data0 after date conversion
print("data0 after date conversion:")
print(head(data0))


df <- data0 %>%
  filter(date >= as.Date("2020-02-01"), date <= as.Date("2023-01-31")) %>% # Explicitly convert filter dates
  filter(pango_lineage != "None", pango_lineage != "NA", pango_lineage != "Unassigned") %>%
  filter(!is.na(date)) %>% # Filter out rows with NA dates after parsing
  mutate(month = format(date, "%Y-%b")) # Keep original month format for reference if needed

# Create a bimonthly grouping variable
df$bimonth <- floor_date(df$date, unit = "month")
# Adjust bimonth to start on odd months (Jan, Mar, May, etc.)
df$bimonth <- df$bimonth - months((month(df$bimonth) - 1) %% 2)
df$bimonth <- as.Date(df$bimonth)


# Print df after initial filtering and bimonth creation
print("df after initial filtering and bimonth creation:")
print(head(df))


# --- 2. Summarize data for total and country plots ---

# Summarize for total plot
df_total_summary <- df %>%
  group_by(pango_lineage, bimonth) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pango_lineage) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(pango_lineage_adjusted = ifelse(total <= 35, "Others", pango_lineage)) %>%
  select(-total) %>%
  group_by(bimonth, pango_lineage_adjusted) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  complete(bimonth, pango_lineage_adjusted, fill = list(n = 0)) %>%
  group_by(bimonth) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup()


# Summarize for country plots
df_country_summary <- df %>%
  group_by(country, pango_lineage, bimonth) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(country, bimonth, pango_lineage) %>%
  summarise(n = sum(n), .groups = "drop") %>% # Ensure correct grouping before calculating total per lineage for adjustment
  group_by(pango_lineage) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(pango_lineage_adjusted = ifelse(total <= 35, "Others", pango_lineage)) %>%
  select(-total) %>%
  group_by(country, bimonth, pango_lineage_adjusted) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  complete(country, bimonth, pango_lineage_adjusted, fill = list(n = 0)) %>%
  group_by(country, bimonth) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup()

# Filter out countries with no data before creating plots
countries_with_data <- df_country_summary %>%
  group_by(country) %>%
  summarise(total_n = sum(n)) %>%
  filter(total_n > 0) %>%
  pull(country)

df_country_summary_filtered <- df_country_summary %>%
  filter(country %in% countries_with_data)


# Print summarized dataframes
print("df_total_summary head:")
print(head(df_total_summary))
print("df_country_summary_filtered head:")
print(head(df_country_summary_filtered))


# Define color palette based on unique lineages across all final data (total and country)
all_lineages_final <- unique(c(df_total_summary$pango_lineage_adjusted, df_country_summary_filtered$pango_lineage_adjusted))

# Ensure the levels for total and country plots match the generated color names
df_total_summary$pango_lineage_adjusted <- factor(df_total_summary$pango_lineage_adjusted, levels = all_lineages_final)
df_country_summary_filtered$pango_lineage_adjusted <- factor(df_country_summary_filtered$pango_lineage_adjusted, levels = all_lineages_final)

# Generate colors from custom_colors with an interval that covers all unique lineages
num_lineages_total <- length(all_lineages_final)
color_indices <- floor(seq(1, length(custom_colors), length.out = num_lineages_total))
plot_colors <- custom_colors[color_indices]
names(plot_colors) <- all_lineages_final # Name the colors with the lineage names for correct mapping


# --- 3. Create plots ---

# Create country plots
country_plots <- df_country_summary_filtered %>%
  split(.$country) %>%
  map(~{
    ggplot(.x, aes(x = bimonth, y = percentage, fill = pango_lineage_adjusted)) +
      geom_area() +
      scale_fill_manual(values = plot_colors, name = "Pangolin Lineages", guide = FALSE) + # Use the shared plot_colors manual scale and set guide = FALSE
      scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.5)) +
      scale_x_date(date_labels = "%b\n%Y", date_breaks = "2 months") +
      labs(x = NULL, y = NULL, title = unique(.x$country)) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.8),
        panel.grid.minor = element_blank(),
        legend.position = "none" # Redundant with guide=FALSE, but good practice
      )
  })

# Create the total plot
Q_total <- df_total_summary %>%
  ggplot(aes(x = bimonth, y = percentage, fill = pango_lineage_adjusted)) +
  geom_area() +
  theme_bw(base_family = "Arial") +
  scale_fill_manual(values = plot_colors, name = "Pangolin Lineages") + # Use the shared plot_colors manual scale and set name
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, size = 10)) +
  theme(axis.title.x = element_text(color = "black", size = 15, face = "bold", vjust = -1.5)) +
  theme(axis.title.y = element_text(color = "black", size = 15, face = "bold", vjust = 1.5)) +
  theme(axis.text.y = element_text(color = "black", size = 10)) +
  xlab("Bimonthly Period") +
  ylab("Percentage") +
  labs(fill = "Pangolin Lineages") + # Ensure consistent legend title
  scale_y_continuous(breaks = seq(0, 1, by = 0.10), labels = scales::percent(seq(0, 1, by = 0.10))) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months")


# --- 4. Combine and save plots ---

# Combine into a composite figure
if (length(country_plots) > 0) {
  combined_plot <- (Q_total / wrap_plots(country_plots, ncol = 4)) +
    plot_layout(heights = c(2, 2), guides = "collect") # Use guides = "collect" here

  display_plot <- combined_plot & theme(legend.position = "bottom") # Apply legend position to combined plot
  display_plot <- display_plot + plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")") # Apply annotation separately
  # Removed explicit collect_guides call as plot_layout with guides="collect" should handle it

} else {
  display_plot <- Q_total + plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")") & theme(legend.position = "bottom")
  print("No country plots generated due to no data.")
}


# Show
print(display_plot)

# Save
if (exists("display_plot")) {
    ggsave(filename = paste0(add_path, "/Figure_3_total_plus_countries.tiff"),
           plot = display_plot, width = 16, height = 12, dpi = 300)
    ggsave(filename = paste0(add_path, "/Figure_3_total_plus_countries.png"),
           plot = display_plot, width = 16, height = 12, dpi = 300)
}
```
## Create Figure 5.

```R
# Read data from CSV
# Make sure the CSV file path is correct
df <- read_csv("/content/path/Combined_data.csv")

# Separate the single column into multiple columns
df <- df %>%
  separate(`Country;Month;Parameter;Value`, into = c("Country", "Month", "Parameter", "Value"), sep = ";") %>%
  mutate(Value = as.numeric(Value)) # Convert Value to numeric

# Check columns and first rows to understand structure
glimpse(df)
head(df)

# Pivot to have each parameter in a column
# Aggregate before pivoting to avoid duplicate issues in pivot_wider
df_aggregated <- df %>%
  group_by(Country, Month, Parameter) %>%
  summarise(Value = sum(Value, na.rm = TRUE)) %>%
  ungroup()

df_wide <- df_aggregated %>%
  pivot_wider(names_from = Parameter, values_from = Value, values_fill = 0)


# Convert Month to date and create numeric time variable
# Repair month format by adding day and then parse as date
df_wide <- df_wide %>%
  mutate(Month = paste0(Month, "-01")) %>% # Add day 01
  mutate(Month = ymd(Month), # Convert to Date
         time_num = as.integer(format(Month, "%Y")) * 12 + as.integer(format(Month, "%m")),
         time_num = time_num - min(time_num, na.rm = TRUE) + 1)

# Identify lineage columns (adjust pattern if necessary)
# This pattern may need adjustment depending on how lineages are named in your data
linajes_cols <- grep("^B\\.|^BA\\.|^AY\\.|^XBB|^Delta|^Omicron", colnames(df_wide), value = TRUE)


# Calculate lineage proportion for each row
df_wide <- df_wide %>%
  rowwise() %>%
  mutate(total_linajes = sum(c_across(all_of(linajes_cols)), na.rm = TRUE), # Added na.rm=TRUE
         across(all_of(linajes_cols), ~ ifelse(total_linajes == 0, 0, . / total_linajes), .names = "prop_{col}")) %>%
  ungroup()

# Show first rows of the prepared dataframe
head(df_wide)
df_wide <- df_wide %>%
  rowwise() %>%
  mutate(
    total_linajes = sum(c_across(all_of(linajes_cols)), na.rm = TRUE),
    across(all_of(linajes_cols), ~ ifelse(total_linajes == 0, 0, . / total_linajes), .names = "prop_{col}")
  ) %>%
  ungroup()
df_wide$total_linajes <- rowSums(df_wide[, linajes_cols], na.rm = TRUE)
df_wide <- df_wide %>%
  mutate(across(all_of(linajes_cols),
                ~ ifelse(total_linajes == 0, 0, . / total_linajes),
                .names = "prop_{col}"))
head(df_wide)
df_wide <- df_wide %>%
  mutate(
    Month = as.character(Month),                
    Month = ym(Month),                          
    time_num = as.integer(format(Month, "%Y")) * 12 + as.integer(format(Month, "%m"))
  ) %>%
  mutate(time_num = time_num - min(time_num, na.rm = TRUE) + 1)
df_wide <- df_aggregated %>%
glimpse(df_aggregated)
unique(df_aggregated$Month)
# Calculate total counts per lineage and filter those with more than 35 total cases
linajes_cols <- grep("^B\\.|^BA\\.|^AY\\.|^XBB|^Delta|^Omicron", colnames(df_wide), value = TRUE)

linajes_totales <- df_wide %>%
  dplyr::select(all_of(linajes_cols)) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "linaje", values_to = "total") %>%
  filter(total > 35) %>%
  arrange(desc(total))

linajes_top <- linajes_totales$linaje

# Prepare long-format dataframe with selected lineages and sum monthly cases
df_long_top <- df_wide %>%
  dplyr::select(Month, all_of(linajes_top)) %>%
  pivot_longer(cols = all_of(linajes_top), names_to = "linaje", values_to = "casos") %>%
  group_by(Month, linaje) %>%
  summarise(total_casos = sum(casos, na.rm = TRUE), .groups = "drop")

# Calculate monthly proportions per lineage
df_prop <- df_long_top %>%
  group_by(Month) %>%
  mutate(prop = total_casos / sum(total_casos)) %>%
  ungroup()

# Get total deaths and cases per month
df_muertes <- df_wide %>%
  group_by(Month) %>%
  summarise(total_deaths = sum(total_deaths, na.rm = TRUE), .groups = "drop")

df_casos <- df_wide %>%
  group_by(Month) %>%
  summarise(total_cases = sum(case_number, na.rm = TRUE), .groups = "drop")

# Merge proportions with deaths and cases
df_plot <- left_join(df_prop, df_muertes, by = "Month") %>%
           left_join(df_casos, by = "Month")
write_csv(df_plot, "/content/path/df_plot_completo.csv")


# Define custom color palette - Use the global custom_colors
# Generate colors from custom_colors with an interval that covers all unique lineages in df_plot
num_lineages_plot <- length(unique(df_plot$linaje))
color_indices_plot <- floor(seq(1, length(custom_colors), length.out = num_lineages_plot))
plot_colors_fig5 <- custom_colors[color_indices_plot]
names(plot_colors_fig5) <- unique(df_plot$linaje) # Name colors with lineage names

# Generate plots
grafico_muertes <- ggplot(df_plot, aes(x = Month)) +
  geom_area(aes(y = prop, fill = linaje), position = "stack", alpha = 0.7, color = "black", size = 0.1) +
  geom_line(aes(y = total_deaths / max(total_deaths, na.rm = TRUE), group = 1),
            color = "black", size = 1.2, linetype = "dashed") +
  scale_y_continuous(
    name = "Lineage proportion",
    sec.axis = sec_axis(~ . * max(df_muertes$total_deaths, na.rm = TRUE), name = "Total deaths")
  ) +
  scale_fill_manual(values = plot_colors_fig5, name = "Lineage", guide = FALSE) + # Use plot_colors_fig5 and set name, add guide = FALSE
  labs(title = "Evolution of lineage proportion with death curve",
       x = "Year", 
       fill = "Lineage") + # Ensure consistent legend title
  theme_minimal(base_family = "Arial") + # Set base font family to Arial
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grafico_casos <- ggplot(df_plot, aes(x = Month)) +
  geom_area(aes(y = prop, fill = linaje), position = "stack", alpha = 0.7, color = "black", size = 0.1) +
  geom_line(aes(y = total_cases / max(total_cases, na.rm = TRUE), group = 1),
            color = "black", size = 1.2, linetype = "dashed") +
  scale_y_continuous(
    name = "Lineage proportion",
    sec.axis = sec_axis(~ . * max(df_casos$total_cases, na.rm = TRUE),
                        name = "Total cases", labels = scales::comma)
  ) +
  scale_fill_manual(values = plot_colors_fig5, name = "Lineage") + # Use plot_colors_fig5 and set name
  labs(title = "Evolution of lineage proportion with case curve",
       x = "Year",
       fill = "Lineage") + # Ensure consistent legend title
  theme_minimal(base_family = "Arial") + # Set base font family to Arial
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") # Keep legend for this plot


plot_final = (grafico_muertes / grafico_casos) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") # Collect guides


print(plot_final)

# Save
ggsave(filename = "/content/path/Plot_final.tiff",
       plot = plot_final,
       width = 12, height = 10, dpi = 300, units = "in") # Increased height to accommodate legend
ggsave(filename = "/content/path/Plot_final.png",
       plot = plot_final,
       width = 12, height = 10, dpi = 300, units = "in") # Increased height to accommodate legend
```
## Create Complementary Figure 1. Individual SARS-CoV-2 sequences obtained by country in Central America and Do-minican Republic from February 2020 to January 2023. (30, 50, 100)

```R
# Install and load required packages
install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
install_if_missing("tidyverse")
install_if_missing("patchwork")
install_if_missing("RColorBrewer")
install_if_missing("viridis") # Install viridis package
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(lubridate) # Load lubridate for dmy()
library(viridis) # Load viridis library


# Data upload
add_path <- "/content/path"  # Update this path to your correct location
df_central_america <- read_tsv(paste0("/content/path/centralamerica.tsv"))

# Convert dates using dmy() before filtering and filter out NA dates
df_central_america$date <- lubridate::dmy(df_central_america$date)
df_central_america <- df_central_america %>%
  filter(!is.na(date)) %>% # Filter out rows with NA dates
  filter(date >= as.Date("2020-02-01"), date <= as.Date("2023-01-31")) %>% # Apply date filtering
  mutate(days = date) # Directly assign the date column to days


# Create a bimonthly grouping variable
df_central_america$bimonth <- floor_date(df_central_america$date, unit = "month")
month(df_central_america$bimonth) <- ifelse(month(df_central_america$bimonth) %% 2 == 0, month(df_central_america$bimonth) - 1, month(df_central_america$bimonth))
df_central_america$bimonth <- as.Date(df_central_america$bimonth)


# Function to create graph by number of lineages
make_lineage_plot <- function(data, top_n = 14) {
  # Calculate quantity by country
  df_count <- data %>% count(country) %>% rename(Count = n)
  data <- data %>% left_join(df_count, by = "country")

  # Top linajes
  lineage_freq <- data %>% count(pango_lineage, sort = TRUE)
  top_lineages <- lineage_freq$pango_lineage[1:min(top_n, nrow(lineage_freq))]

  data <- data %>%
    mutate(pango_lineage_grouped = if_else(pango_lineage %in% top_lineages, pango_lineage, "Other"))

  # Generate colors from custom_colors with an interval
  num_lineages <- length(unique(data$pango_lineage_grouped))
  color_indices <- floor(seq(1, length(custom_colors), length.out = num_lineages))
  plot_colors <- custom_colors[color_indices]

  # Ensure that the levels for the fill aesthetic match the selected colors
  data$pango_lineage_grouped <- factor(data$pango_lineage_grouped, levels = unique(data$pango_lineage_grouped))

  # Create graphic
  ggplot(data) +
    theme_classic(base_family = "Arial") + # Set base font family to Arial
    geom_segment(aes(x = as.Date("2020-02-01"), y = reorder(country, Count),
                     xend = as.Date("2023-02-01"), yend = country, group = country),
                 colour = "grey80", size = 5) +
    geom_point(aes(fill = pango_lineage_grouped, x = bimonth, y = reorder(country, Count)), # Use bimonth for x-axis
               position = position_jitter(width = 0.2, height = 0.2),
               shape = 21, stroke = 0.1, col = 'black', size = 3, alpha = 0.7) +
    ylab('') + xlab('Sample Collection Date') +
    scale_fill_manual(values = plot_colors, name = 'Pangolin Lineages') + # Use plot_colors manual scale
    scale_x_date(date_labels = "%b\n%Y", date_breaks = "2 months", date_minor_breaks = "1 month", limits = as.Date(c("2020-02-01", "2023-03-01")))+ # Explicitly set x-axis limits
    coord_cartesian(clip = "off") +  # ← Prevent labels from being cut off
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      axis.title.x = element_text(color = "black", size = 15, face = "bold", vjust = -0.5),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
      axis.title.y = element_text(color = "black", size = 15, face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      panel.grid.major.x = element_line(color = "#3c3c3c", size = 0.2),
      panel.grid.minor.x = element_line(color = "darkgray", size = 0.08),
      plot.margin = margin(t = 10, r = 90, b = 70, l = 10)  # ← Increased right margin
    ) +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 4, byrow = TRUE)) +
    ggtitle(paste0("Top ", top_n, " Pangolin Lineages"))
}

# Create the three selected graphics
p30 <- make_lineage_plot(df_central_america, 30) + theme(legend.position = "bottom")
p50 <- make_lineage_plot(df_central_america, 50) + theme(legend.position = "bottom")
p100 <- make_lineage_plot(df_central_america, 100) + theme(legend.position = "bottom")

# Organize vertically
mosaic_plot <- p30 / p50 / p100 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Show in console
print(mosaic_plot)

# Save final image
ggsave(filename = paste0(add_path, "/Supplement_of_Figure_2.tiff"),
       plot = mosaic_plot,
       width = 20,
       height = 26,
       units = "in",
       dpi = 300
)
```


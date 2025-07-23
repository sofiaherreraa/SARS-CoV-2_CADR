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

# Define the desired color palette
base_colors <- c("lightgrey", "#117864", "#1d6996", "#e15759", "#edad08", "#9c755f", "#5f4690")

# Create a color generating function
color_generator <- colorRampPalette(base_colors)
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

# Generate colors from the defined color_generator for activities
num_activities <- length(unique(dataset$Activity))
plot_colors_activities <- color_generator(num_activities)
names(plot_colors_activities) <- unique(dataset$Activity) # Name colors with activity names


# Plotting Gantt chart
timeline<-ggplot(dataset, aes(x = Start, xend = End, y = Activity, yend = Activity, color = Activity)) + # Use original Start and End dates
  geom_segment(linewidth = 10) + # Changed size to linewidth
  scale_color_manual(values = plot_colors_activities, guide = "none") + # Use plot_colors_activities manual scale, changed guide = FALSE to "none"
  labs(x = "Date", y = NULL) +
  theme_minimal(base_family = "Arial") + # Set base font family to Arial
  theme(panel.grid.major.x = element_line(color = "#3c3c3c", linewidth = 0.2)) + # Changed size to linewidth
  theme(panel.grid.minor.x = element_line(color = "darkgray", linewidth = 0.1)) + # Changed size to linewidth
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

## Create Figure 2.  Individual SARS-CoV-2 sequences obtained by country in Central America and Dominican Republic from February 2020 to January 2023. The top 14 most prevalent lineages were individually labeled, and all the remaining lineages were labeled as other.
```R
# Load libraries
library(tidyverse)
library(lubridate) # Load lubridate for dmy()

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

# Generate colors from the defined color_generator for lineages
num_lineages <- length(unique(df_central_america$pango_lineage_grouped))
plot_colors <- color_generator(num_lineages)
names(plot_colors) <- unique(df_central_america$pango_lineage_grouped) # Name colors with lineage names


# Plot
panelopas <- ggplot(data=df_central_america) +
  theme_classic(base_family = "Arial")+ # Set base font family to Arial
  geom_segment(aes(x=as.Date("2020-02-01"), y=reorder(country,Count), xend=as.Date("2023-02-01"), yend=country, group=country), colour="grey80", linewidth=5) + # Changed size to linewidth

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
  theme(panel.grid.major.x = element_line(color = "#3c3c3c", linewidth = 0.2)) + # Changed size to linewidth
  theme(panel.grid.minor.x = element_line(color = "darkgray", linewidth = 0.08)) + # Changed size to linewidth
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month", date_minor_breaks = "1 month", limits = as.Date(c("2020-02-01", "2023-03-01")))+ # Explicitly set x-axis limits
  guides(fill = guide_legend(override.aes = list(size=5), nrow = 2))+
  xlab('Sample Collection Date')


# Extract the directory path from read_path (assuming read_path is defined earlier)
add_path <- dirname("/content/path/centralamerica.tsv") # Define add_path locally if read_path is not globally available

ggsave(paste0(add_path, "/Figure 2.tiff"), panelopas, width=18.8, height=5.74)
ggsave(paste0(add_path, "/Figure 2.png"), panelopas, width=18.8, height=5.74)
```
## Create Figure 3. Relative percentages of SARS-CoV-2 lineages circulating in Central America and the Dominican Republic from February 2020 to January 2023. 
```R

add_path <- "/content/path"
input_path <- file.path(add_path, "centralamerica.tsv")

# Load and clean data
df <- read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, format = "%d/%m/%y")) %>%
  filter(!is.na(date)) %>%
  filter(date > as.Date("2020-01-01") & date < as.Date("2023-02-01")) %>%
  filter(!pango_lineage %in% c("None", "NA", "Unassigned", NA)) %>%
  filter(!is.na(country))

# Global general table
p <- df %>%
  mutate(month = format(date, "%Y-%b")) %>%
  filter(!is.na(month)) %>%
  group_by(pango_lineage, month) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pango_lineage_adj = ifelse(n <= 35, "Others", as.character(pango_lineage)))

others <- p %>%
  filter(pango_lineage_adj == "Others") %>%
  group_by(month) %>%
  summarise(n = sum(n), pango_lineage_adj = "Others", .groups = "drop")

p_final <- bind_rows(
  p %>% filter(pango_lineage_adj != "Others"),
  others
) %>%
  filter(!is.na(pango_lineage_adj), !is.na(month))

p_final <- p_final %>%
  group_by(month) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup()

# Prepare variables to order and factor months
year_month_split <- strsplit(p_final$month, "-")
years <- sapply(year_month_split, `[`, 1)
months <- sapply(year_month_split, `[`, 2)

month_levels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

p_final <- p_final %>%
  mutate(
    year = as.numeric(years),
    month_f = factor(months, levels = month_levels),
    year_month = paste(year, month_f, sep = "-")
  ) %>%
  filter(!is.na(month_f)) %>%
  arrange(year, month_f)

p_final$year_month <- factor(p_final$year_month, levels = unique(p_final$year_month))

# Valid lineages without "Others" and sorted
linajes_validos <- unique(p_final$pango_lineage_adj)
linajes_validos <- linajes_validos[linajes_validos != "Others"]
linajes_validos <- sort(linajes_validos)

niveles_lineages <- c("Others", linajes_validos)

# Manual color palette
base_colors <- c("#117864", "#1d6996", "#e15759", "#edad08", "#9c755f", "#5f4690")

paletas_por_color <- lapply(base_colors, function(col) {
  colorRampPalette(c(lighten(col, amount = 0.6), col, darken(col, amount = 0.3)))(10)
})

paleta_60 <- unlist(paletas_por_color)

num_linajes <- length(linajes_validos)
if(num_linajes > length(paleta_60)) {
  warning("There are more lineages than available colors, some colors will be repeated")
  paleta_asignada <- rep(paleta_60, length.out = num_linajes)
} else {
  paleta_asignada <- paleta_60[1:num_linajes]
}

palette_colors <- c(Others = "lightgrey", setNames(paleta_asignada, linajes_validos))

# Ordered factor for p_final
p_final$pango_lineage_adj <- factor(p_final$pango_lineage_adj, levels = niveles_lineages)

# --- New country-month calculation with strict filtering ---

# Get valid months and lineages from general to filter country-month data
meses_validos <- unique(p_final$month)
linajes_validos_full <- unique(p_final$pango_lineage_adj)

# Count by country, month, and lineage (without grouping "Others" yet)
df_counts <- df %>%
  mutate(month = format(date, "%Y-%b")) %>%
  filter(month %in% meses_validos) %>%   # Only valid months from general
  group_by(country, month, pango_lineage) %>%
  summarise(n = n(), .groups = "drop")

# Create auxiliary table with valid lineages and months and their general classification
linajes_mes_general <- p %>%
  select(pango_lineage, month, pango_lineage_adj) %>%
  distinct()

# Join to assign general classification to each lineage-country-month
df_counts_joined <- df_counts %>%
  inner_join(linajes_mes_general, by = c("pango_lineage", "month")) %>%
  rename(pango_lineage_adj = pango_lineage_adj)

# Sum "Others" by country and month
others_country_month <- df_counts_joined %>%
  filter(pango_lineage_adj == "Others") %>%
  group_by(country, month) %>%
  summarise(n = sum(n), pango_lineage_adj = "Others", .groups = "drop")

# Lineages different from Others
linajes_country_month <- df_counts_joined %>%
  filter(pango_lineage_adj != "Others")

# Combine all
df_final_country <- bind_rows(linajes_country_month, others_country_month)

# Complete with all possible combinations to show months without data in each country
all_countries <- unique(df_final_country$country)
all_months <- unique(p_final$month)
all_lineages <- niveles_lineages

df_final_country_complete <- df_final_country %>%
  complete(
    country = all_countries,
    month = all_months,
    pango_lineage_adj = all_lineages,
    fill = list(n = 0)
  ) %>%
  group_by(country, month) %>%
  mutate(total_n = sum(n)) %>%
  ungroup() %>%
  mutate(proportion = if_else(total_n == 0, 0, n / total_n)) %>%
  select(-total_n)

# Factorize for plots
df_final_country_complete <- df_final_country_complete %>%
  mutate(
    month_f = factor(month, levels = unique(p_final$month)),
    pango_lineage_adj = factor(pango_lineage_adj, levels = niveles_lineages),
    country = factor(country, levels = sort(unique(country)))
  )

# Create labels like a), b), c) for facets
labels <- data.frame(
  country = levels(df_final_country_complete$country),
  tag = paste0(letters[seq_along(levels(df_final_country_complete$country))], ")")
)

labeller_with_tag <- function(variable, value) {
  tags <- labels$tag[match(value, labels$country)]
  paste0(tags, " ", value)
}

# General plot without legend
g_general <- ggplot(p_final, aes(x = year_month, y = percentage, fill = pango_lineage_adj)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = palette_colors) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(t=10, r=10, b=10, l=10)
  ) +
  labs(
    x = "Months",
    y = "Percentage",
    title = "a) General"
  )

# Country plot with legend and labels a), b), c), ...
country_levels <- sort(unique(df_final_country_complete$country))
country_labels <- paste0(letters[2:(length(country_levels)+1)], ") ", country_levels)
names(country_labels) <- country_levels

g_country <- ggplot(df_final_country_complete, aes(x = month_f, y = proportion, fill = pango_lineage_adj)) +
  geom_col() +
  facet_wrap(~country, ncol = 4, scales = "free_x",
             labeller = labeller(country = country_labels)) +
  scale_fill_manual(values = palette_colors, name = "Lineages") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.6, "cm"),
    plot.margin = margin(t=10, r=10, b=10, l=10)
  ) +
  guides(fill = guide_legend(ncol = 8, byrow = TRUE)) +
  labs(
    x = "Months",
    y = "Relative Frequency"
  )

# Combine and save
final_plot <- g_general / g_country + patchwork::plot_layout(heights = c(2, 3))

ggsave(filename = file.path(add_path, "Figure_3.tiff"), plot = final_plot, width = 16, height = 14, dpi = 300)
ggsave(filename = file.path(add_path, "Figure_3.png"), plot = final_plot, width = 16, height = 14, dpi = 300)
print(final_plot)
```
## Create Figure 5. Temporal dynamics of SARS-CoV-2 lineage proportions with COVID-19 statistical curves in Central America and the Dominican Republic.

```R
# Read data from CSV
# Make sure the CSV file path is correct
df <- read_csv("/content/path/Combined_data.csv")

df <- df %>%
  separate(`Country;Month;Parameter;Value`, into = c("Country", "Month", "Parameter", "Value"), sep = ";") %>%
  mutate(Value = as.numeric(Value))

df_aggregated <- df %>%
  group_by(Country, Month, Parameter) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

df_wide <- df_aggregated %>%
  pivot_wider(names_from = Parameter, values_from = Value, values_fill = 0) %>%
  mutate(Month = ymd(paste0(Month, "-01")))

linajes_cols <- grep("^B\\.|^BA\\.|^AY\\.|^XBB|^Delta|^Omicron", colnames(df_wide), value = TRUE)

df_wide <- df_wide %>%
  rowwise() %>%
  mutate(
    total_linajes = sum(c_across(all_of(linajes_cols)), na.rm = TRUE),
    across(all_of(linajes_cols), ~ ifelse(total_linajes == 0, 0, . / total_linajes), .names = "prop_{col}")
  ) %>%
  ungroup()

linajes_totales <- df_wide %>%
  summarise(across(all_of(linajes_cols), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "linaje", values_to = "total") %>%
  filter(total > 35)

linajes_top <- linajes_totales$linaje

df_long_top <- df_wide %>%
  select(Month, all_of(linajes_top)) %>%
  pivot_longer(cols = all_of(linajes_top), names_to = "linaje", values_to = "casos") %>%
  group_by(Month, linaje) %>%
  summarise(total_casos = sum(casos, na.rm = TRUE), .groups = "drop")

df_prop <- df_long_top %>%
  group_by(Month) %>%
  mutate(prop = total_casos / sum(total_casos)) %>%
  ungroup()

df_muertes <- df_wide %>%
  group_by(Month) %>%
  summarise(total_deaths = sum(total_deaths, na.rm = TRUE), .groups = "drop")

df_casos <- df_wide %>%
  group_by(Month) %>%
  summarise(total_cases = sum(case_number, na.rm = TRUE), .groups = "drop")

df_plot <- df_prop %>%
  left_join(df_muertes, by = "Month") %>%
  left_join(df_casos, by = "Month")

# Colors
base_colors <- c("lightgrey", "#117864", "#1d6996", "#e15759", "#edad08", "#9c755f", "#5f4690")
color_generator <- colorRampPalette(base_colors)
num_lineages_plot <- length(unique(df_plot$linaje))
plot_colors_fig5 <- color_generator(num_lineages_plot)
names(plot_colors_fig5) <- unique(df_plot$linaje)

# Total deaths
grafico_muertes <- ggplot(df_plot, aes(x = Month)) +
  geom_area(aes(y = prop, fill = linaje), position = "stack", alpha = 0.7, color = "black", size = 0.1) +
  geom_line(aes(y = total_deaths / max(total_deaths, na.rm = TRUE), group = 1),
            color = "black", size = 1.2, linetype = "dashed") +
  scale_y_continuous(
    name = "Lineage proportion",
    sec.axis = sec_axis(~ . * max(df_muertes$total_deaths, na.rm = TRUE), name = "Total deaths")
  ) +
  scale_fill_manual(values = plot_colors_fig5, name = "Lineage") +
  labs(title = "Evolution of lineage proportion with death curve", x = "Year") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Total cases
grafico_casos <- ggplot(df_plot, aes(x = Month)) +
  geom_area(aes(y = prop, fill = linaje), position = "stack", alpha = 0.7, color = "black", size = 0.1) +
  geom_line(aes(y = total_cases / max(total_cases, na.rm = TRUE), group = 1),
            color = "black", size = 1.2, linetype = "dashed") +
  scale_y_continuous(
    name = "Lineage proportion",
    sec.axis = sec_axis(~ . * max(df_casos$total_cases, na.rm = TRUE), name = "Total cases", labels = comma)
  ) +
  scale_fill_manual(values = plot_colors_fig5, name = "Lineage") +
  labs(title = "Evolution of lineage proportion with case curve", x = "Year") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Use patchwork for combind legends
final_plot <- grafico_muertes / grafico_casos + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')

print(final_plot)

# Save
ggsave("/cloud/project/Figure_5.tiff", plot = final_plot, width = 12, height = 10.5, dpi = 300, units = "in")
ggsave("/cloud/project/Figure_5.png", plot = final_plot, width = 12, height = 10.5, dpi = 300, units = "in")


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

  # Generate colors from the defined color_generator for lineages
  num_lineages <- length(unique(data$pango_lineage_grouped))
  plot_colors <- color_generator(num_lineages)
  names(plot_colors) <- unique(data$pango_lineage_grouped) # Name colors with lineage names

  # Ensure that the levels for the fill aesthetic match the selected colors
  data$pango_lineage_grouped <- factor(data$pango_lineage_grouped, levels = unique(data$pango_lineage_grouped))

  # Create graphic
  ggplot(data) +
    theme_classic(base_family = "Arial") + # Set base font family to Arial
    geom_segment(aes(x = as.Date("2020-02-01"), y = reorder(country, Count),
                     xend = as.Date("2023-02-01"), yend = country, group = country),
                 colour = "grey80", linewidth = 5) + # Changed size to linewidth
    geom_point(aes(fill = pango_lineage_grouped, x = days, y = reorder(country, Count)), # Use days for x-axis
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
      panel.grid.major.x = element_line(color = "#3c3c3c", linewidth = 0.2), # Changed size to linewidth
      panel.grid.minor.x = element_line(color = "darkgray", linewidth = 0.08), # Changed size to linewidth
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
ggsave(filename = paste0(add_path, "/Supplement_of_Figure_2.png"),
       plot = mosaic_plot,
       width = 20,
       height = 26,
       units = "in",
       dpi = 300
)
```


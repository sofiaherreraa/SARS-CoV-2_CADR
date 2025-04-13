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
```
## Create Figure 1. Timeline of processes implemented to strengthen genomic surveillance of SARS-CoV-2 in Central America and the Dominican Republic.

```R
# Update this path to their desired location
add_path = "/Users/yourpath"

#Define the activities and dates
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
  '2021-03-01', '2021-04-01', '2021-08-23',
  '2021-09-20', '2021-11-24', '2021-06-01', 
  '2021-06-01', '2021-07-01', '2021-07-01',
  '2022-06-15', '2022-02-05', '2022-08-08'))
# Create a data frame for plotting
dataset <- tibble(
  Activity = activities,
  Start = dates,
  End = c(
    '2021-01-31', '2021-05-21', '2021-07-24',
    '2023-02-01', '2021-08-18', '2022-08-22',
    '2022-08-22', '2022-04-30', '2022-03-22', 
    '2022-09-29', '2021-08-29', '2023-01-03',
    '2022-09-23', '2022-06-23', '2022-08-12'))

# Convert to tibble and create sorting date column
dataset <- dataset %>%
  mutate(Start = as.Date(Start, "%Y-%m-%d"),
         End = as.Date(End, "%Y-%m-%d"),
         Activity = factor(Activity, levels = rev(Activity))) %>%
  arrange(desc(Start))

# Plotting Gantt chart
timeline<-ggplot(dataset, aes(x = Start, xend = End + 2, y = Activity, yend = Activity, color = Activity)) +
  geom_segment(size = 10) +
  scale_color_manual(values = rep(c( "#666666", "#e7896f","#cc503e","#f28e2c","#59a14f",
                                     "#117864", 'cadetblue3',"#1d6996", "skyblue2","lightsteelblue1", 
                                     "darkslateblue","#6f4070","#b993ac","#994e95", "#d37295", "#e9a8b3"), 2), guide = FALSE) +
  labs(x = "Date", y = NULL) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month", date_minor_breaks = "1 month")+
  theme(panel.grid.major.x = element_line(color = "#3c3c3c", size = 0.2)) +
  theme(panel.grid.minor.x = element_line(color = "darkgray", size = 0.1)) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold", vjust=-1)) +
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12, hjust = 1))

timeline
ggsave(paste0(add_path, "/Figure 1.tiff"), width=16.1, height=5.74)
```

## Create Figure 2.  Individual SARS-CoV-2 sequences obtained by country in Central America and Do-minican Republic from February 2020 to January 2023. The top 14 most prevalent lineages were individually labeled, and all the remaining lineages were labeled as other.
```R
#Load data for Figure 2 and 3
data0 <- read_tsv(paste0(add_path, "/centralamerica.tsv"))

#Descending country count 
names(data0)
df_count<- data0 %>% 
  count(country) %>%
  arrange(desc(n))
df_count

#Renaming data 0
df_central_america<-data0

#Convert 'date' to Date data type
df_central_america$date<-as.Date(df_central_america$date)

#Grouping data by date 1 month
df_central_america$days<-as.Date(cut(df_central_america$date,breaks = "day",start.on.monday = FALSE))
df_central_america$date2<-as.Date(cut(df_central_america$date,breaks = "2 weeks",start.on.monday = FALSE))
df_central_america$date3<-as.Date(cut(df_central_america$date,breaks = "1 month",start.on.monday = FALSE))

# Create date3 column with format YYYY-MMM
df_central_america$date3 <- format(df_central_america$date, "%Y-%b")

#Determine time frame 2020-01-01 to 2023-02-01
df_central_america = df_central_america%>%
  filter(date> "2020-01-01" , date <"2023-02-01")

#Ascending country count 
df_count <- df_central_america %>% count(country)
names(df_count)[names(df_count) == "country"] <- "country"
names(df_count)[names(df_count) == "n"] <- "Count"
df_count[order(df_count$Count),]

#Add Country Count on the left side of the table
df_central_america = df_central_america %>% 
  left_join(df_count, by = c("country" = "country"))

#Creating data set with Pangolin Lineage Count - Ascending
pangolin_count<-as.data.frame(table(df_central_america$pango_lineage))
pangolin_count<-pangolin_count[order(pangolin_count$Freq), ]
panelopas<-ggplot(data=df_central_america)  + theme_classic()+
  geom_segment(aes(x=as.Date("2020-02-01"), y=reorder(country,Count), xend=as.Date("2023-02-01"), yend=country, group=country), colour="grey80", size=5) +
  geom_point(aes(x=days, fill='Other',y=reorder(country,Count)),position = position_jitter(width=0.2, height=0.2), shape=21,stroke=0.05, col='grey70', size=3)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.1.1'),aes(fill='BA.1.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='XBB.1'),aes(fill='XBB.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='AY.113'),aes(fill='AY.113',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='A.2.5'),aes(fill='A.2.5',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='B.1'),aes(fill='B.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='AY.100'),aes(fill='AY.100',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.2'),aes(fill='BA.2',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.2.9'),aes(fill='BA.2.9',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BE.1'),aes(fill='BE.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.5.2.23'),aes(fill='BA.5.2.23',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.5.2'),aes(fill='BA.5.2',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BQ.1.1'),aes(fill='BQ.1.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.4.6'),aes(fill='BA.4.6',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  geom_point(data=subset(df_central_america, pango_lineage=='BA.2.12.1'),aes(fill='BA.2.12.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=3, alpha=0.5)+
  ylab('')+ xlab('month')+ ggtitle('') +
  scale_fill_manual(values=c('grey','peachpuff3','#edad08','cadetblue3',"#1d6996","#994e95","lightsteelblue1","#e7896f", "#e9a8b3","lightyellow1","#117864","#cc503c","#59a14f","thistle2","#5f4690"), name='Pangolin Lineages')+theme(legend.position="bottom") +
  #scale_fill_manual(values=c("darkgrey","azure2","lightyellow1", "darkseagreen1", "aquamarine1","turquoise2",'cadetblue3',"lightsteelblue1","steelblue","#1d6996","darkslateblue", "mediumpurple4","#994e95", "thistle2", "pink"), name='Pangolin Lineages')+theme(legend.position="bottom") +
  #scale_fill_manual(values=c( "#666666", "#b993ac","#cc503e","skyblue2","#59a14f","#117864", "steelblue", "#edad08","darkblue","darkslateblue", "#6f4070","#e7896f", "thistle2","lightgrey","#d37295"), name='Pangolin Lineages')+theme(legend.position="bottom") +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold", vjust=-0.5)) +
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold")) +
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(panel.grid.major.x = element_line(color = "#3c3c3c", size = 0.2)) +
  theme(panel.grid.minor.x = element_line(color = "darkgray", size = 0.08)) +
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month", date_minor_breaks = "1 month")+
  guides(fill = guide_legend(override.aes = list(size=5), nrow = 2))+
  xlab('Sample Collection Date')

panelopas
ggsave(paste0(add_path, "/Figure 2.tiff"), width=18.8, height=5.74)

```
## Create Figure 3. Relative percentages of SARS-CoV-2 lineages circulating in Central America and Dominican Republic from February 2020 to January 2023. Lineages with a frequency exceding 120 occurrences (n>120) per month were selected.
```R
total<- df_central_america%>%
  filter(pango_lineage !="None", pango_lineage !="NA", pango_lineage !="Unassigned") %>%
  filter(!is.na(Date))%>%
  unique()%>%
  summarise(total = n())

#Removing duplicates from pango_lineages
p <- df_central_america %>%
  group_by(pango_lineage) %>%
  filter(pango_lineage != "None", pango_lineage != "NA", pango_lineage !="Unassigned" ) %>%
  filter(!is.na(Date)) %>%
  unique() %>%
  
  #Adding month column with format YY-MM, grouping and filter by date - lineage
  mutate(month = format(date, "%Y-%b")) %>%
  group_by(pango_lineage, month) %>%
  unique() %>%
  
  
  # Filter out unwanted month values
  filter(month != "None", month != "NA") %>%
  summarise(n = n())

  
 #Identify lineages with counts <=35 and replace them with "Others"
  p <- p %>%
  mutate(pango_lineage_adjusted = ifelse(n <=35, "Others", as.character(pango_lineage)))

 # Filter for "Others" lineages
 others <- p %>%
  filter(pango_lineage_adjusted == "Others") %>%
  
  # Group by month
  group_by(month) %>%
  
  # Summarize to get a single row per month for "Others"
  summarise(n = sum(n),
            pango_lineage_adjusted = "Others")

# Combine with the rest of the data
p <- bind_rows(p %>% filter(pango_lineage_adjusted != "Others"), others)
 summarise(n = n())

  # Calculate percentage
 p <- p %>%
   group_by(month) %>%
   mutate(percentage = n / sum(n))
  
# Convert month to character
p$month <- as.character(p$month)

# Extract year and month separately
year_month <- strsplit(p$month, "-")
years <- sapply(year_month, `[`, 1)
months <- sapply(year_month, `[`, 2)

# Define the order of months
month_order <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Create a factor with correct order of levels for months
p$month <- factor(months, levels = month_order)

# Convert years to numeric
p$year <- as.numeric(years)

# Sort by year and month
p <- p[order(p$year, p$month), ]

# Combine year and month into a single column
p$year_month <- paste(p$year, p$month, sep = "-")

# Create a factor with correct order of levels for year_month
p$year_month <- factor(p$year_month, levels = unique(p$year_month))


#Degradation color palette for Figure 3
# Define the number of colors
num_colors <- 55
# Generate color scale from green to blue to purple
colors <- (colorRampPalette(c("lightgrey", "#666666","#bab0ab",'peachpuff3',"#9c755f", "#e7896f",
                            "#cc503e","#f28e2c","#edad08", "#0f8554","#59a14f",
                             "#117864", 'cadetblue3',"#1d6996","steelblue", "skyblue2","lightsteelblue1","darkblue", "darkslateblue", "#5f4690","#6f4070",
                             "#94346e", "#b993ac", "thistle2","#994e95", "#d37295",
                             "#e9a8b3",  "#ff9da7","#e15759"))(num_colors))

# Ensure pango_lineage_adjusted is a factor with desired order
p$pango_lineage_adjusted <- factor(p$pango_lineage_adjusted, levels = unique(p$pango_lineage_adjusted))
Q<-p%>%
  ggplot(aes(x=year_month, y=percentage,fill=pango_lineage_adjusted, label=scales::percent(abs(percentage)))) +
  geom_col() +
  theme_bw() + 
  geom_text(alpha=0)+
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1,size = 10 )) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold", vjust=-1.5)) +
  theme(axis.title.y = element_text(color="black", size=15, face="bold", vjust=1.5)) +
  theme(axis.text.y = element_text(color="black", size=10))+
  xlab("Months") +
  ylab("Percentage") +
  labs(fill = paste("Lineages"))

# Adjust y-axis breaks and labels
Q <- Q + scale_y_continuous(breaks = seq(0, 1, by = 0.10), labels = scales::percent(seq(0, 1, by = 0.10)))

# Adjust the label for "Lineages n"
Q <- Q +
  annotate("text", x = Inf, y = Inf, label = paste("n = (", total$total, ")", sep = ""), 
           hjust = 1, vjust = 1.5, size = 4, fontface = "bold")
Q
ggsave(paste0(add_path, "/Figure 3.tiff"), width=12, height=6)
```


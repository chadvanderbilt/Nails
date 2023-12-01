#!/opt/R/4.1.3/bin/Rscript


## Load required packages
# Data manipulation
library(dplyr)        

# Data tidying
library(tidyr)        

# Data visualization
library(ggplot2)      

# String manipulation
library(stringr)      

# Variable handling
library(varhandle)      

# Additional ggplot features
library(ggforce)       

# Excel file reading
library(readxl)        

# Publication-ready ggplot extensions
library(ggpubr)        

# Model fitting and summary
library(gmodels)       

# Data manipulation with tibbles
library(tibble)

# Factor handling
library(forcats)

# Grid and table manipulation
library(grid)
library(gtable)
library(gridExtra)

# PNG image reading
library(png)

## Set working directory
setwd("/home/User/Nails/data")
plot_dir <- sub(getwd(), pattern = "data", replacement = "plots")

## Figure 1b - Comparing Fragment Size Between Extraction Methods

## Read in data
Fragment_size <- read.csv("Fragment_size.csv")

## Process data
Fragment_size_processed <- Fragment_size %>%
  separate(case, sep = " ", c("Mnumber", "Method", NA)) %>%
  mutate(Method = ifelse(Method == "BB", "Method 2", "Method 1")) %>%
  rename(Size = size)
level_order <- c("Method 1", "Method 2")
pdf(paste0(plot_dir,"/FIG1B.pdf"), width = 8, height = 12)
## Create the plot
ggplot(Fragment_size_processed, aes(x = factor(Method, level = level_order,) , y = Size)) +
  geom_boxplot(fill =c("#b2b2f2", "#9cf1c2")) +
  labs(y = "Fragment Size (bp)", x = "Extraction Method") +  # Change y-axis label
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20), # Adjusted size for better visibility
    axis.title = element_text(size = 22, face = "bold")
  ) +
  ## Add statistical comparison
  stat_compare_means(
    method = "t.test", 
    comparisons = list(c("Method 1", "Method 2")),
    label = "p.signif",
    paired = TRUE,
    inherit.aes = TRUE
  )
dev.off()


## Figure 1d - Comparing DNA Yield Between Extraction Methods


extraction <- read.csv("extraction.csv", header = T)

pdf(paste0(plot_dir,"/FIG1D.pdf"), width = 9, height = 10)
extraction %>%  mutate(method = ifelse(test_type == "old", "Method 1", "Method 2")) %>% 
  ggplot(aes(x = factor(method, level = level_order), y = yield)) +
  geom_boxplot(fill =c("#b2b2f2", "#9cf1c2")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20), # Adjusted size for better visibility
    axis.title = element_text(size = 22, face = "bold"), 
    panel.border = element_blank()
  ) +
  ## Add statistical comparison
  stat_compare_means(
    method = "t.test", 
    comparisons = list(c("Method 1", "Method 2")),
    label = "p.signif",
    paired = FALSE, 
    inherit.aes = TRUE, 
    label.size = 6
  ) +
  labs(y = "DNA Extraction Yield (ng)", x = "Extraction Method") +  # Change y-axis label  
  facet_zoom(ylim = c(0, 1300), zoom.size = 0.2, horizontal = T) 
dev.off()


## Figure 1e - Comparing DNA Fragment Size Between Extraction Methods by Quartiles

quartiles <- read.csv("fragment_quartiles.csv", header = T)
level_order_2 <- c("Method 1", "Method 2")
pdf(paste0(plot_dir,"/FIG1E.pdf"), width = 9, height = 10)
quartiles %>%
  gather(-sampleid, -method, key = stat, value = size) %>%
  mutate(stat = str_replace_all(stat, "\\.", " ")) %>%
  mutate(stat = str_replace(stat, "e ", "e: ")) %>%
  mutate(method = ifelse(method == "Bead Blaster", "Method 2", "Method 1")) %>%
  ggplot(aes(x = factor(method, level = level_order_2), y = size)) +
  geom_boxplot(fill =c("#b2b2f2", "#9cf1c2", "#b2b2f2", "#9cf1c2","#b2b2f2", "#9cf1c2")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20), # Adjusted size for better visibility
    axis.title = element_text(size = 22, face = "bold"), 
    strip.text = element_text(size = 12, face = "bold") #facet label size
  ) +
  facet_wrap(~stat) +
  labs(y = "Fragment Size (bp)", x = "Extraction Method") +
  ## Add statistical comparison
  stat_compare_means(
    method = "t.test", 
    comparisons = list(c("Method 1", "Method 2")),
    label = "p.signif",
    paired = FALSE, 
    inherit.aes = TRUE, 
    label.size = 6
  )
dev.off()

## Figure 1f - Comparing DNA Fragment Size Patterns Between Extraction Methods

## Read in data
fragment_size_density <- 
read.table("val_run.txt", sep = "\t",  header = T) %>%
rename(method1 = bc1269, beadblaster_sheared = bc1279, beadblaster_nonsheared = bc1280, bonemarrow = bc1283) %>% select(-contains("bc")) %>%
gather(-insert_size, key = sample_type, value = proportion) %>%
filter(insert_size != "Peak") %>% 
mutate(insert_size = as.numeric(insert_size))


## Create the plot

pdf(paste0(plot_dir,"/FIG1F.pdf"), width = 18, height = 9)

fragment_size_density %>%
ggplot(aes(x=insert_size, color = sample_type, y = proportion)) +
  geom_line(size =2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20), # Adjusted size for better visibility
    axis.title = element_text(size = 30, face = "bold"),
    legend.position = c(0.7, 0.8), # Adjust the position of the legend
    legend.background = element_rect(fill = "white", color = NA), # Make legend background white
    legend.box.background = element_rect(fill = "white", color = NA), # Make legend box background white
    legend.key = element_rect(fill = "white", color = NA), # Make legend key background white
    legend.title = element_text(size = 40), # Adjust legend title size
    legend.text = element_text(size = 32, color = "black", face = "italic"), # Adjust legend text size and appearance
    legend.margin = margin(0, 0, 0, 0) # Adjust legend margin
  ) +
  labs(y = "Proportion of Fragments", x = "Fragment Size (bp)") +
  scale_color_discrete(
    name = "Sample Type / Extraction Method",
    type = c(bonemarrow="black", beadblaster_nonsheared ="#b2e0f2", beadblaster_sheared="#9cf1c2", method1="#b2b2f2"),
    breaks = c("beadblaster_nonsheared", "beadblaster_sheared", "method1", "bonemarrow"),
    labels = c(beadblaster_nonsheared = "Method 2: Non-Sheared Nail DNA", beadblaster_sheared= "Method 2: Sheared Nail DNA", method1="Method 1: Sheared Nail DNA", bonemarrow="Bone Marrow DNA")
  )+
  scale_x_continuous(breaks = seq(0, 500, 100)) +
  scale_y_continuous(breaks = seq(0, 0.1, 0.01))  
dev.off()




## Figure 2a/2b - Comparing mutations detected in Nail by tumor type



Fig2a_b_Patient_level <- read.table("Fig2a_b_Patient_level.csv", sep = ",", header = T)

Fig2a_b_Mutation_level <- read.table("Fig2a_b_Mutation_level.csv", sep = ",", header = T)

# using fig2c_d, make a stacked bar chart where Patients with and without mutation in nail stack up to 100%
pdf(paste0(plot_dir,"/FIG2A_B.pdf"), width = 40, height = 20)



desired_order <- c("ALAL","AML","BPDCN","HDCN","MCD","MDS","MDS/MPN","MDS Workup","MLNER","MPN","MPN Workup","BLL","HL","LATL","MBN","MTNN","PCM","TLL","Myeloid","Lymphoid","Grand Total")

B_Cell_Other <- c("BPDCN", "MLNER", "HL") #combine these and regenerate plot. 

tmp <- Fig2a_b_Patient_level %>%
  mutate(OncoTree_Code = str_replace_all(OncoTree_Code, pattern = "_", " ")) %>%
  mutate(totals = ifelse(str_detect(OncoTree_Code, "Myeloid|Lymphoid|Grand Total") == TRUE, "Total", "Individual Tumor Types")) %>%
  mutate(OncoTree_Code = fct_relevel(OncoTree_Code, desired_order)) %>%
  group_by(OncoTree_Code) %>%
  mutate(x_label = paste(OncoTree_Code, " (n =", min(tally), "/", sum(tally), ")")) %>%
  arrange(factor(OncoTree_Code, levels = rev(desired_order))) %>% ungroup()
new_order <- unique(tmp$x_label)
patient_level <- tmp %>%
  mutate(x_label = fct_relevel(x_label, rev(new_order))) %>%
  ggplot(aes(x = x_label, y = tally, fill = patient_with_mutations)) +
  geom_bar(position = "fill", stat = "identity", alpha = 0.9) +
  scale_fill_manual(name = "Patients with Mutations in Nail",
                    values = c("false" = "#428ed5", "true" = "#f3af50"),
                    labels = c("false" = "No Mutations", "true" = "With Mutations")) +
  scale_alpha_continuous(guide = 'none') +
  theme_minimal() +
  labs(
    y = "Percentage",
    x = "OncoTree Code"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust= 1, hjust = 1, size = 30), 
        axis.text.y = element_text(size = 30),
        legend.position = "bottom",
        plot.title = element_text(size = 35, hjust = 0),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20), 
        axis.title = element_text(size = 30), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  scale_x_discrete(
  breaks = rev(new_order),
  labels = rev(new_order),
  expand = c(0.1, 0.1, 0.1, 0.3)  # Adjust the expansion as needed
) +
  coord_flip()

tmp_mutations <- Fig2a_b_Mutation_level %>%
  mutate(OncoTree_Code = str_replace_all(OncoTree_Code, pattern = "_", " ")) %>%
  mutate(totals = ifelse(str_detect(OncoTree_Code, "Myeloid|Lymphoid|Grand Total") == TRUE, "Total", "Individual Tumor Types")) %>%
  mutate(OncoTree_Code = fct_relevel(OncoTree_Code, desired_order)) %>%
  group_by(OncoTree_Code) %>%
  mutate(x_label = paste("(n =", min(tally), "/", sum(tally), ")")) %>%
  arrange(factor(OncoTree_Code, levels = rev(desired_order))) %>% ungroup()
new_order <- unique(tmp_mutations$x_label)
mutation_level <- tmp_mutations %>%
  mutate(x_label = fct_relevel(x_label, rev(new_order))) %>%
  ggplot(aes(x = x_label, y = tally, fill = patient_with_mutations)) +
  geom_bar(position = "fill", stat = "identity", alpha = 0.9) +
  scale_fill_manual(name = "Mutations in Nail",
                    values = c("false" = "#428ed5", "true" = "#f3af50"),
                    labels = c("false" = "No Mutations", "true" = "With Mutations")) +
  scale_alpha_continuous(guide = 'none') +
  theme_minimal() +
  labs(
    y = "Percentage",
    x = "OncoTree Code"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust= 1, hjust = 1, size = 30), 
        axis.text.y = element_text(size = 30),
        legend.position = "bottom",
        plot.title = element_text(size = 35, hjust = 0),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20), 
        axis.title.x = element_text(size = 24), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  scale_x_discrete(
  breaks = rev(new_order),
  labels = rev(new_order),
  expand = c(0.1, 0.1, 0.1, 0.3)  # Adjust the expansion as needed
) +
  coord_flip()

combined_plot <- grid.arrange(patient_level, mutation_level, ncol = 2)
combined_plot
# Print the combined plot
combined_plot

dev.off()


## Figure 2c - Comparing mutations VAF detected in Nail by tumor type



Fig2C<- read.table("Fig2C.csv", sep = ",", header = T)

Fig2C_gather <- Fig2C %>% select(1,2, 6,7) %>%
                  gather(key="mut_source", value = "tally", -VAF_categories)

desired_order <- rev(c("90-100","80-89","70-79","60-69","50-59","40-49","30-39","20-29","10-19","5-9","2-4","1-2"))



pdf(paste0(plot_dir,"/FIG2C.pdf"), width = 20, height = 15)
Fig2C_gather %>% 
          mutate(VAF_categories = fct_relevel(VAF_categories, desired_order)) %>%
          ggplot() +
          geom_bar(stat = "identity",
                    position = "dodge", 
                    width = 1,
                    aes(x= VAF_categories, 
                    y = tally,
                    fill = mut_source)) +
          # scale_fill_manual(values = c("#f3af50", "blue","green","#428ed5")) +
          scale_fill_manual(values = c( Nail_lymphoid_mutation_tally="#eaea8f",Nail_myeloid_mutation_tally="#c9829a",Tumor_mutations_tally ="#428ed5"),
                            labels = c(Nail_lymphoid_mutation_tally = "Nail (Lymphoid Tumors)", Nail_myeloid_mutation_tally = "Nail (Myeloid Tumor)", Tumor_mutations_tally = "Tumor")) +
          coord_flip() +
          theme_minimal()+
  labs(y = "Number of Mutations",
       fill = "Tissue Type",
       x = "Variant Allele Frequency Percentage Bin") +
  theme(axis.text.x = element_text(angle = 90, vjust= 1, hjust = 1, size = 30), 
        axis.text.y = element_text(size = 30),
        legend.position = "bottom",
        # legend.title = element_text()
        plot.title = element_text(size = 35, hjust = 0.5),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        axis.title = element_text(size = 32),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
dev.off()

## figure 2d - Comparing mutations detected in Nail by tumor class

all_cv3_nail_variants_13Dec2021 <- read.delim("supplemental_data_mutations.csv", sep = ",", header = T)
Tumor_Type <- read.csv("Tumor_Type.csv")

pdf(paste0(plot_dir,"/FIG2D.pdf"), width = 20, height = 10)
VAF_myeloid_lymphoid <- all_cv3_nail_variants_13Dec2021 %>% 
  select(OncoTree_Code, Gene, VF, Nail_AltFreq, Called_nail_mutation) %>% 
  filter(Called_nail_mutation == TRUE) %>%
  mutate(Tumor_over_Nail = VF/Nail_AltFreq) %>%
  left_join(Tumor_Type, relationship = "many-to-many") %>% 
  filter(is.na(M_L)==F) %>%
  # group_by(M_L) %>%
  # summarise(mean_tumor_nail = mean(Tumor_over_Nail), mean_nail = mean(Nail_AltFreq), mean_tumor = mean(VF))
  ggplot(aes(fill = M_L, x = Tumor_over_Nail)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c(L="#eaea8f",M="#c9829a" ), 
                    labels = c(M = "Nail (Myeloid Tumors)", L = "Nail (Lymphoid Tumor)")) +
  labs(x = "Ratio of VAF (Tumor over Nail)", y = "Density", fill = "Tissue Type") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 24),
        legend.position = "bottom",
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))
VAF_myeloid_lymphoid
dev.off()

## Figure 3a/3b - Comparing VAF in Nail by to in Tumor. 



Fig3A_B_myeloid <- 
  read_excel("Fig3A_B.xlsx", 
             sheet = "FIGURE 2F ALL MYELOID") %>%
             rownames_to_column(var = "case") %>%
             rename(tumor_vaf = 2, nail_vaf =3) %>%
             arrange(desc(tumor_vaf))





Fig3A_B_myeloid2 <- 
  read_excel("Fig3A_B.xlsx", 
             sheet = "FIGURE 2F ALL MYELOID") %>%
             rownames_to_column(var = "case") %>%
             rename(tumor_vaf = 2, nail_vaf =3) %>%
             arrange(desc(nail_vaf))

desired_order <- Fig3A_B_myeloid$case



Fig3A_B_lymphoid <- 
  read_excel("Fig3A_B.xlsx", 
             sheet = "FIGURE 2G LYMPHOID") %>%
             rownames_to_column(var = "case")%>%
             rename(tumor_vaf = 2, nail_vaf =3) %>%
             arrange(desc(tumor_vaf))






Fig3A_B_lymphoid2 <- 
  read_excel("Fig3A_B.xlsx", 
             sheet = "FIGURE 2G LYMPHOID") %>%
             rownames_to_column(var = "case")%>%
             rename(tumor_vaf = 2, nail_vaf =3) %>%
             arrange(desc(nail_vaf))


list_frames_inner <- c("Fig3A_B_myeloid2", "Fig3A_B_lymphoid2")
list_frames_outer <- c("Fig3A_B_myeloid", "Fig3A_B_lymphoid")

for (i in 1:length(list_frames_outer)) {
   

desired_order <-  get(list_frames_outer[[i]])$case
plot <- get(list_frames_outer[i]) %>% 
          mutate(case = fct_relevel(case, desired_order)) %>%
          ggplot() +
          geom_bar(stat = "identity",width = 1, aes(x= case, 
                                      y = tumor_vaf, 
                                      fill = "VAF in Tumor", 
                                      alpha = 0.4)) +
          geom_bar(stat = "identity",width = 1, aes(x= case, 
                                      y = nail_vaf, 
                                      fill = "VAF in Nails"))+
          scale_fill_manual(values = c("#b66a00", "#b3bfcb")) +
          # coord_flip() +
          scale_alpha_continuous(guide = 'none') +
          theme_minimal()+
  labs(y = "VAF",
      fill = "Tissue Type",
      x = "Mutation") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20),
        legend.position = "bottom",
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20), 
        axis.title = element_text(size = 24),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())





  ggsave(plot = plot, filename = paste0(plot_dir,"/figure_new", list_frames_outer[i], ".pdf"), device = "pdf", width = 20, height = 10, units = "cm")
  plot <- plot + #remove the legend
  theme(legend.position="none")
  ggsave(plot = plot, filename = paste0(plot_dir,"/figure_new", list_frames_outer[i], ".png"), device = "png", width = 20, height = 10, units = "cm", bg = "transparent")
  ggsave(plot = plot, filename = paste0(plot_dir,"/figure_new_nolegend", list_frames_outer[i], ".pdf"), device = "pdf", width = 20, height = 10, units = "cm")
  # Load the vectorized image using readPNG (for example, you may need to use a different function depending on your vectorized image format)
img <- readPNG(paste0(plot_dir,"/figure_new", list_frames_outer[i], ".png"))

# Extract red, green, blue, and alpha channels
red_channel <- img[, , 1]
green_channel <- img[, , 2]
blue_channel <- img[, , 3]
alpha_channel <- img[, , 4]

# Enhance contrast (example: multiply pixel values by 2)
contrast_factor <- 1
red_channel_contrast <- pmax(0, pmin(1, red_channel * contrast_factor))
green_channel_contrast <- pmax(0, pmin(1, green_channel * contrast_factor))
blue_channel_contrast <- pmax(0, pmin(1, blue_channel * contrast_factor))

# Create a 3D array for the new image
img_contrast <- array(0, dim = c(nrow(img), ncol(img), 4))
img_contrast[, , 1] <- red_channel_contrast
img_contrast[, , 2] <- green_channel_contrast
img_contrast[, , 3] <- blue_channel_contrast
img_contrast[, , 4] <- alpha_channel

desired_order <-  get(list_frames_inner[[i]])$case
plot <- get(list_frames_outer[i]) %>% 
          mutate(case = fct_relevel(case, desired_order)) %>%
          ggplot() +
          geom_bar(stat = "identity",width = 1, aes(x= case, 
                                      y = tumor_vaf, 
                                      fill = "VAF in Tumor", 
                                      alpha = 0.6)) +
          geom_bar(stat = "identity",width = 1, aes(x= case, 
                                      y = nail_vaf, 
                                      fill = "VAF in Nails"))+
          scale_fill_manual(values = c("#f3af50", "#428ed5")) +
          # coord_flip() +
          scale_alpha_continuous(guide = 'none') +
          theme_bw()+
  labs(y = "VAF",
      fill = "Tissue Type",
      x = "Mutation") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 40),
        legend.position = "bottom",
        legend.title = element_text(size = 48),
        legend.text = element_text(size = 40), 
        axis.title = element_text(size = 48),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(plot = plot, filename = paste0(plot_dir,"/figure_new", list_frames_inner[i], ".pdf"), device = "pdf", width = 20, height = 10, units = "cm")

# Combine the ggplot and the vectorized image
# combined_plot <- ggplot_gtable(ggplot_build(plot))
# combined_plot <- gtable_add_grob(combined_plot, 
#                                  grobs = list(rasterGrob(image, width = 0.1, height = 0.1)),
#                                  t = 1, l = 4, b = 1, r = 5)
combined_plot <- plot +
  annotation_custom(
    rasterGrob(img_contrast, width = unit(0.7, "npc"), height = unit(0.7, "npc"))
  )
ggsave(plot = combined_plot, filename = paste0(plot_dir, "/FIG3_", list_frames_inner[i], ".pdf"), device = "pdf", width = 80, height = 40, units = "cm")

}




## Figures 3b/3c/3e/3f - Comparing mutations detected in Nail by tumor type


#generate bines
VF_bins <- all_cv3_nail_variants_13Dec2021 %>% 
  select(OncoTree_Code, Gene, VF, Nail_AltFreq) %>% 
  mutate(VF_bins = ifelse(VF >= 0.9, "90-100", 
                          ifelse(VF >= 0.8, "80-89", 
                                 ifelse(VF >= 0.7, "70-79", 
                                        ifelse(VF >= 0.6, "60-69", 
                                               ifelse(VF >= 0.5, "50-59", 
                                                      ifelse(VF >= 0.4, "40-49", 
                                                             ifelse(VF >= 0.3, "30-39", 
                                                                    ifelse(VF >= 0.2, "20-29", 
                                                                           ifelse(VF >= 0.1, "10-19", 
                                                                                  ifelse(VF >= 0.05, "5-9", 
                                                                                         ifelse(VF >= 0.02, "2-4", 
                                                                                                ifelse(VF >= 0.01, "1-2", NA))))))))))))) %>%
  left_join(Tumor_Type, relationship = "many-to-many") 

# make bar chart with flipped axis and facet by M_L and put labels on top of bars saying (n = x)
pdf(paste0(plot_dir,"/FIG3B_C_E_F.pdf"), width = 18, height = 15)
M_L <- c("L", "M")
desired_order <- rev(c("90-100","80-89","70-79","60-69","50-59","40-49","30-39","20-29","10-19","5-9","2-4","1-2"))

fig3b_d <- VF_bins %>% 
mutate(VF_bins = fct_relevel(VF_bins, desired_order)) %>%
filter(is.na(M_L)==F, is.na(VF_bins)==F) %>%
group_by(VF_bins, M_L) %>%
summarise(counts = n()) %>%
mutate(M_L = ifelse(M_L=="M", "1M", M_L)) %>%
ggplot(aes(x = VF_bins, y = counts)) +
geom_bar(stat = "identity", fill = "#428ed5", alpha = 0.6) +
geom_text(aes(label = paste0("(n=",counts, ")"), hjust = -0.5)) +
ylim(c(0, 1500)) +
facet_wrap(~M_L, scales = "free_x", strip.position = "left", ncol=1) +
coord_flip() + 
theme_minimal() +
labs(x = "Variant Allele Frequency Percentage Bin", y = "Number of Mutations (Tumor)") +
theme(axis.text.x = element_text(vjust= 1, hjust = 1, size = 20, face = "bold"), 
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 24), 
      strip.text = element_blank(),
      text = element_text(size = 30, face = "bold"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()) 



fig3c_f <- VF_bins %>% 
mutate(VF_bins = fct_relevel(VF_bins, desired_order)) %>%
filter(is.na(M_L)==F, is.na(VF_bins)==F) %>%
mutate(nail_mut = ifelse(Nail_AltFreq >= 0.01, TRUE, FALSE)) %>%
mutate(M_L = ifelse(M_L=="M", "1M", M_L)) %>%
group_by(VF_bins, M_L, nail_mut) %>%
summarise(counts = n()) %>% ungroup() %>%
ggplot(aes(x = VF_bins, y = counts, fill = nail_mut, alpha = nail_mut)) +
geom_bar(stat = "identity", position = "fill") +
scale_fill_manual(values = c("#f3af50", "#428ed5")) +
scale_alpha_manual(values = c(0.99, 0.6)) +
facet_wrap(~M_L, scales = "free_x", strip.position = "left", ncol = 1) +
coord_flip() + 
theme_minimal() +
labs(x = "", y = "Fraction of Mutations Present in Nail") +
theme(axis.text.x = element_text(vjust= 1, hjust = 1, size = 20, face = "bold"), 
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_blank(),
      axis.title = element_text(size = 24), 
      legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()) 
#attach the two plots
fig3b_d
fig3c_f
combined_plot <- grid.arrange(fig3b_d, fig3c_f, ncol = 2)
combined_plot

dev.off()


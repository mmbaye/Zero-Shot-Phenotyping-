# ════════════════════════════════════════════════════════════
#  ANALYSE BIOMASSE + ACP DES TRAITS SAM
#  Deux analyses complémentaires pour le papier
# ════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(FactoMineR)   
library(factoextra)   

select <- dplyr::select
filter <- dplyr::filter

# Palette
PAL <- c(
  NGANDA="#D85A30", PAYENNE="#854F0B", NGUINTH="#1D9E75",
  GOLOBE="#085041", ICRS14001="#534AB7", FAOUROU="#3C3489",
  ICSV745="#791F1F", S35="#5F5E5A"
)

#  Charger le dataset mergé 
df <- read_csv("dataset_SAM_LIDAR_Yield_merged.csv",
               show_col_types = FALSE) %>%
  mutate(genotype = factor(genotype, levels = names(PAL)))

# ════════════════════════════════════════════════════════════
#  PARTIE 1  ANALYSE BIOMASSE
# ════════════════════════════════════════════════════════════


# Statistiques biomasse
df %>%
  select(genotype, dry_biom_mean, fresh_biom_mean,
         dry_pan_mean, seed_yield_mean) %>%
  arrange(desc(dry_biom_mean)) %>%
  print(n=Inf)

#  Fig A : Biomasse sèche par génotype 
fig_biom <- df %>%
  filter(!is.na(dry_biom_mean)) %>%
  ggplot(aes(x = fct_reorder(genotype, dry_biom_mean, .desc=TRUE),
             y = dry_biom_mean, fill = genotype)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_text(aes(label=round(dry_biom_mean,2)),
            vjust=-0.4, size=3.5, fontface="bold") +
  scale_fill_manual(values=PAL, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(title="A  Dry biomass by genotype",
       x=NULL, y="Dry biomass (t/ha)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.x=element_blank(),
        plot.title=element_text(face="bold"))

# Fig B : Seed yield par génotype 
fig_yield <- df %>%
  filter(!is.na(seed_yield_mean)) %>%
  ggplot(aes(x = fct_reorder(genotype, seed_yield_mean, .desc=TRUE),
             y = seed_yield_mean, fill = genotype)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_text(aes(label=round(seed_yield_mean,3)),
            vjust=-0.4, size=3.5, fontface="bold") +
  scale_fill_manual(values=PAL, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(title="B  Seed yield by genotype",
       subtitle="⚠️ Affected by grain abortion at end of cycle",
       x=NULL, y="Seed yield (t/ha)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.x=element_blank(),
        plot.title=element_text(face="bold"),
        plot.subtitle=element_text(colour="#E24B4A", size=9))

# Corrélation SAM leaf area vs biomasse 
cat("\nCorrélation SAM traits × Biomasse :\n")
traits_sam <- c("sam_pan_L","sam_pan_A",
                "sam_leaf_L_mean","sam_leaf_l_mean",
                "sam_leaf_A_mean","sam_angle_mean")

map_dfr(traits_sam, function(ts) {
  x <- df[[ts]]; y <- df$dry_biom_mean
  v <- !is.na(x) & !is.na(y)
  if(sum(v)<3) return(NULL)
  ct <- cor.test(x[v], y[v])
  tibble(SAM_trait=ts, r=round(ct$estimate,3),
         p=round(ct$p.value,4),
         signif=case_when(ct$p.value<0.05~"*",
                          ct$p.value<0.10~".",TRUE~"ns"))
}) %>% arrange(p) %>% print()

#  Scatter SAM leaf area vs dry biomass 
fig_scatter_biom <- df %>%
  filter(!is.na(sam_leaf_A_mean), !is.na(dry_biom_mean)) %>%
  ggplot(aes(x=sam_leaf_A_mean, y=dry_biom_mean,
             colour=genotype, label=genotype)) +
  geom_smooth(method="lm", colour="grey50", fill="grey80",
              alpha=0.2, linewidth=1, se=TRUE,
              inherit.aes=FALSE,
              aes(x=sam_leaf_A_mean, y=dry_biom_mean)) +
  geom_point(size=4, alpha=0.9) +
  geom_text_repel(size=3.5, fontface="bold",
                  box.padding=0.4, point.padding=0.3) +
  scale_colour_manual(values=PAL, guide="none") +
  labs(title="C  SAM leaf area vs. dry biomass",
       x="SAM mean leaf area (cm²)",
       y="Dry biomass (t/ha)") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"))

# Scatter SAM pan_L vs seed yield 
ct_pan_yield <- cor.test(df$sam_pan_L, df$seed_yield_mean)
fig_scatter_yield <- df %>%
  filter(!is.na(sam_pan_L), !is.na(seed_yield_mean)) %>%
  ggplot(aes(x=sam_pan_L, y=seed_yield_mean,
             colour=genotype, label=genotype)) +
  geom_smooth(method="lm", colour="#E24B4A", fill="#E24B4A",
              alpha=0.15, linewidth=1, se=TRUE,
              inherit.aes=FALSE,
              aes(x=sam_pan_L, y=seed_yield_mean)) +
  geom_point(size=4, alpha=0.9) +
  geom_text_repel(size=3.5, fontface="bold",
                  box.padding=0.4, point.padding=0.3) +
  annotate("label", x=Inf, y=Inf,
           label=sprintf("r = %.3f\np = %.3f",
                         ct_pan_yield$estimate,
                         ct_pan_yield$p.value),
           hjust=1.1, vjust=1.3, size=3.5,
           family="mono", fill="white", colour="grey30") +
  scale_colour_manual(values=PAL, guide="none") +
  labs(title="D  SAM panicle length vs. seed yield",
       subtitle="Grain abortion noted at end of cycle",
       x="SAM-estimated panicle length (cm)",
       y="Seed yield (t/ha)") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"),
        plot.subtitle=element_text(colour="#E24B4A", size=9))

# Assembler Fig7
fig7 <- (fig_biom | fig_yield) / (fig_scatter_biom | fig_scatter_yield) +
  plot_annotation(
    title="Agronomic performance and SAM trait associations",
    theme=theme(plot.title=element_text(face="bold", size=13))
  )

ggsave("Fig7_biomasse_yield_SAM.png", fig7,
       width=13, height=10, dpi=300, bg="white")
ggsave("Fig7_biomasse_yield_SAM.pdf", fig7,
       width=13, height=10, bg="white")
cat("\n✅ Fig7 sauvegardée\n")

# ════════════════════════════════════════════════════════════
#  PARTIE 2 — ACP DES TRAITS SAM
# ════════════════════════════════════════════════════════════


# Préparer la matrice pour l'ACP
# Utiliser les moyennes par génotype
df_pca <- df %>%
  filter(!is.na(sam_pan_L)) %>%
  select(genotype,
         `Panicle L`  = sam_pan_L,
         `Leaf L`     = sam_leaf_L_mean,
         `Leaf l`     = sam_leaf_l_mean,
         `Leaf A`     = sam_leaf_A_mean,
         `Angle`      = sam_angle_mean) %>%
  column_to_rownames("genotype")

# NA
df_pca <- df_pca[complete.cases(df_pca), ]

cat(sprintf("ACP sur %d génotypes × %d traits SAM\n\n",
            nrow(df_pca), ncol(df_pca)))

# Calculer l'ACP (centrée-réduite)
res_pca <- PCA(df_pca, scale.unit=TRUE, graph=FALSE)

# Variance expliquée
eig <- get_eigenvalue(res_pca)
cat("Variance expliquée par composante :\n")
print(round(eig[1:4,], 2))

# Coordonnées variables (loadings)

var_coords <- get_pca_var(res_pca)
print(round(var_coords$contrib, 2))

# Corrélation PC1 avec yield et biomasse
ind_coords <- get_pca_ind(res_pca)
df_pc <- ind_coords$coord %>%
  as_tibble(rownames="genotype") %>%
  rename(PC1=Dim.1, PC2=Dim.2, PC3=Dim.3) %>%
  left_join(df %>% select(genotype, seed_yield_mean,
                           dry_biom_mean), by="genotype")

#Corrélation PC1 × agronomie 
for(var in c("seed_yield_mean","dry_biom_mean")) {
  v <- !is.na(df_pc[[var]])
  if(sum(v)>=3) {
    ct <- cor.test(df_pc$PC1[v], df_pc[[var]][v])
    cat(sprintf("  PC1 × %-20s : r=%.3f  p=%.4f\n",
                var, ct$estimate, ct$p.value))
  }
}

#Corrélation PC2 × agronomie
for(var in c("seed_yield_mean","dry_biom_mean")) {
  v <- !is.na(df_pc[[var]])
  if(sum(v)>=3) {
    ct <- cor.test(df_pc$PC2[v], df_pc[[var]][v])
    cat(sprintf("  PC2 × %-20s : r=%.3f  p=%.4f\n",
                var, ct$estimate, ct$p.value))
  }
}

# Fig 8A : Cercle de corrélations 
pct_var <- round(eig$variance.percent[1:2], 1)

fig_pca_var <- fviz_pca_var(
  res_pca,
  col.var    = "contrib",
  gradient.cols = c("#5F5E5A","#534AB7","#D85A30"),
  repel      = TRUE,
  labelsize  = 4
) +
  labs(
    title    = "A  PCA variable loadings",
    subtitle = sprintf("PC1 = %.1f%%  |  PC2 = %.1f%% of variance",
                       pct_var[1], pct_var[2]),
    colour   = "Contribution (%)"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"))

# Fig 8B : Biplot génotypes 
fig_pca_ind <- fviz_pca_biplot(
  res_pca,
  repel         = TRUE,
  col.ind       = df_pca %>% rownames(),
  palette       = PAL[rownames(df_pca)],
  col.var       = "grey40",
  label         = "all",
  pointsize     = 4,
  labelsize     = 4,
  arrowsize     = 0.8,
  alpha.ind     = 0.9
) +
  labs(
    title    = "B  PCA biplot — genotypes × SAM traits",
    subtitle = sprintf("PC1 = %.1f%%  |  PC2 = %.1f%% of variance",
                       pct_var[1], pct_var[2])
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"),
        legend.position="none")

# Fig 8C : Scree plot 
fig_scree <- fviz_eig(
  res_pca,
  addlabels = TRUE,
  barfill   = "#185FA5",
  barcolor  = "white",
  linecolor = "#D85A30"
) +
  labs(title="C  Scree plot", x="Principal component",
       y="% variance explained") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"))

#  Fig 8D : PC1 vs dry biomass 
fig_pc1_biom <- df_pc %>%
  filter(!is.na(dry_biom_mean)) %>%
  ggplot(aes(x=PC1, y=dry_biom_mean,
             colour=genotype, label=genotype)) +
  geom_smooth(method="lm", colour="grey50",
              fill="grey80", alpha=0.2, se=TRUE,
              inherit.aes=FALSE,
              aes(x=PC1, y=dry_biom_mean)) +
  geom_point(size=4, alpha=0.9) +
  geom_text_repel(size=3.5, fontface="bold",
                  box.padding=0.4) +
  scale_colour_manual(values=PAL, guide="none") +
  labs(title="D  PC1 vs. dry biomass",
       x=sprintf("PC1 (%.1f%% variance)", pct_var[1]),
       y="Dry biomass (t/ha)") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"))

# Assembler Fig8
fig8 <- (fig_pca_var | fig_pca_ind) / (fig_scree | fig_pc1_biom) +
  plot_annotation(
    title = "Principal Component Analysis of SAM-extracted architectural traits",
    theme = theme(plot.title=element_text(face="bold", size=13))
  )

ggsave("Fig8_ACP_SAM_traits.png", fig8,
       width=14, height=11, dpi=300, bg="white")

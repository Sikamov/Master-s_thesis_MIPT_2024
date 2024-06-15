# Загрузка необходимых библиотек
library(immunarch)
library(ggplot2)
library(FactoMineR)
library(lattice)
library(RColorBrewer)
library(pheatmap)

# Загрузка данных
DATA = repLoad("/home/sikamov/TCR/res/")
DATA$meta <- DATA$meta[-c(12,63,64),]
names(DATA$data)[names(DATA$data) == DATA$meta$Sample] <- DATA$meta$Name
names(DATA$meta)[names(DATA$meta) == "Sample"] <- "Library"
names(DATA$meta)[names(DATA$meta) == "Name"] <- "Sample"

# Фильтрация данных
group <- c("HC","UC","CD")
point <- c("K06","K03")
data <- repFilter(DATA, .method = "by.meta", .query = list(Group = include(group), Point = exclude(point)))$data
data_a <- repFilter(DATA, .method = "by.meta", .query = list(Chain = include("α"),Group = include(group), Point = exclude(point)))$data
data_b <- repFilter(DATA, .method = "by.meta", .query = list(Chain = include("β"),Group = include(group), Point = exclude(point)))$data
meta_a <- repFilter(DATA, .method = "by.meta", .query = list(Chain = include("α"),Group = include(group), Point = exclude(point)))$meta
meta_b <- repFilter(DATA, .method = "by.meta", .query = list(Chain = include("β"),Group = include(group), Point = exclude(point)))$meta

#Разнообразие репертуаров
repExplore(data_b, "volume") %>%  vis() +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Количество уникальных клонотипов") +
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(30))

repExplore(data_b, "clones") %>%  vis() +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Количество клонотипов") +
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(30))

repExplore(data_b, "volume") %>%  vis(.by = c("Patient"), .meta = DATA$meta, .test = F) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Количество уникальных клонотипов") +
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(9))

repExplore(data_b, "clones") %>%  vis(.by = c("Patient"), .meta = DATA$meta, .test = F) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Количество клонотипов") +
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(9))

repExplore(data_b, "volume") %>%  vis(.by = c("Group"), .meta = DATA$meta, .test = T) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Количество уникальных клонотипов") +
  scale_fill_manual(values = c("firebrick3", "#238B37", "dodgerblue3")) +  
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman"))

repExplore(data_b, "clones") %>%  vis(.by = c("Group"), .meta = DATA$meta, .test = T) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Количество клонотипов") +
  scale_fill_manual(values = c("firebrick3", "#238B37", "dodgerblue3")) +  
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman"))

#Метрики разнообразия
repDiversity(data_b, "hill", .norm = TRUE) %>%
  vis(.by = c("Group"), .meta = DATA$meta, .test = TRUE) +
  labs(title = NULL, subtitle = NULL, y = "Число Хилла", x = "Параметр α") +
  scale_color_manual("Группа",
    values = c("firebrick3", "#238B37", "dodgerblue3"),
    labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК")) +  
  theme(text = element_text(size = 18, family = "Times New Roman"), 
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)) 

repDiversity(data_b, "gini.simp", .norm = TRUE)  %>%  vis(.by = c("Group"), .meta = DATA$meta, .test = T) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Индекс Джини-Симпсона") +
  scale_fill_manual(values = c("firebrick3", "#238B37", "dodgerblue3")) +  
  scale_x_discrete(labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК"))  + 
  theme(text = element_text(size = 18, family = "Times New Roman"))

#Гомеостаз
repClonality(data_b, "homeo", .clone.types = c(Низкая = 1e-05, Средняя = 0.001, Высокая = 0.01, `Очень высокая` = 1)) %>% 
  vis(.by = c("Group"), .meta = DATA$meta, .test = F)  +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = "Частота встречаемости в репертуаре") +
  scale_fill_manual("Группа",
                     values = c("firebrick3", "#238B37", "dodgerblue3"),
                     labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК")) +  
  theme(text = element_text(size = 18, family = "Times New Roman"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14)) 

repClonality(data_b, "homeo", .clone.types = c(Низкая = 1e-05, Средняя = 0.001, Высокая = 0.01, `Очень высокая` = 1)) %>% vis()  +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  guides(fill=guide_legend(title="Частота встречаемости в репертуаре")) +
  theme(text = element_text(size = 18, family = "Times New Roman"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "top",
        legend) 

#Топ
repClonality(data_b, "top", .head = c(10, 100, 1000, 3000)) %>% 
  vis(.by = c("Group"), .meta = DATA$meta, .test = F)  +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = "Ранг в репертуаре") +
  scale_fill_manual("Группа",
                    values = c("firebrick3", "#238B37", "dodgerblue3"),
                    labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК")) +  
  theme(text = element_text(size = 18, family = "Times New Roman"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "top") 

repClonality(data_b, "top", .head = c(10, 100, 1000, 3000)) %>% vis()  +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  guides(fill=guide_legend(title="Ранг в репертуаре")) +
  theme(text = element_text(size = 18, family = "Times New Roman"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "top")

#Теловая диаграмма с индексами Мориситы-Хорна
ov <- repOverlap(data_a, "morisita")

annotations <- data.frame(Group = meta_a$Group)
rownames(annotations) <- meta_a$Sample

pheatmap(
  ov,
  annotation_col = annotations,
  annotation_colors = list(Group = c("CD" = "firebrick3", "HC" = "#238B37", "UC" = "dodgerblue3")),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  #border_color = NA,
  #border = TRUE,
  main = "",
  fontsize = 14,
  fontsize_row = 14,
  fontsize_col = 14,
  annotation_legend = TRUE,
  colorRampPalette(brewer.pal(11, "Spectral") )(100),  
  fontfamily = "serif",
)

#PCA
ov[is.na(ov)] <- 1
res_pca <- PCA(scale(ov), ncp = 5, graph = FALSE)
pca <- data.frame(res_pca$ind$coord, Sample = row.names(res_pca$ind$coord))
pca$Group <- DATA$meta$Group[match(pca$Sample, DATA$meta$Sample)]
pca$Cluster <- kmeans(pca[, 1:3], centers = 3)$cluster
v <- round(res_pca$eig[, 1] / sum(res_pca$eig[, 1]) * 100, 1)

hulls <- do.call(rbind, lapply(split(pca, pca$Cluster), function(df) df[chull(df$Dim.1, df$Dim.2), ]))

ggplot(pca, aes(x=Dim.1, y=Dim.2)) + 
  geom_polygon(data=hulls, aes(x=Dim.1, y=Dim.2, group=Cluster), fill=NA, color="black", size=0.5, linetype="dashed") +
  geom_point(aes(colour=factor(Group)), size=10) +
  ggtitle("") + xlab(paste0('PC1 (',v[1],'%)')) + ylab(paste0("PC2 (",v[2],"%)")) +
  geom_text(label=row.names(pca), size=2, color='white', fontface="bold", family = "Times New Roman") + 
  theme_bw() +
  theme(text = element_text(size = 14, 
        family = "Times New Roman"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "right",
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(colour=NA)) +
  scale_color_manual(values = c("CD" = "firebrick3", "HC" = "#238B37", "UC" = "dodgerblue3"),
                     labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК")) 

#Использование V и J сегментов TCR
gu <- geneUsage(data_b, "HomoSapiens.TRAJ", .norm = TRUE, .ambig = c("inc", "exc", "maj"))
gene_names <- gu$Names
gu <- as.matrix(gu[, -which(names(gu) == "Names")])
rownames(gu) <- gene_names
gu[is.na(gu)] <- 0
gu<- gu[, sapply(gu, is.numeric)]
gu <- as.matrix(gu)

pheatmap(
  gu,
  annotation_col = annotations,
  annotation_colors = list(Group = c("CD" = "firebrick3", "HC" = "#238B37", "UC" = "dodgerblue3")),
  cluster_rows = F,
  cluster_cols = T,
  show_colnames = TRUE,
  show_rownames = TRUE,
  border_color = NA,
  border = TRUE,
  main = "",
  fontsize = 14,
  fontsize_row = 12,
  fontsize_col = 12,
  annotation_legend = F,
  fontfamily = "serif",
  color = colorRampPalette(c("white", "#238B37"))(100)  
)

#БК уникальные сегменты
gu <- geneUsage(data_a, "HomoSapiens.TRАJ", .norm = TRUE, .ambig = c("inc", "exc", "maj"))
gu <- gu[rowSums(is.na(gu[, grepl("K", colnames(gu))])) == sum(grepl("K", colnames(gu))), ]
vis(gu, .by = c("Group"), .meta = DATA$meta, .test = F) +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  scale_fill_manual("Группа",
                    values = c("firebrick3", "#238B37", "dodgerblue3"),
                    labels = c("CD" = "БК", "HC" = "К", "UC" = "НЯК")) +  
  theme(text = element_text(size = 18, family = "Times New Roman"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

#Теловая карта с перепредставленными клонотипами
data_list <- list()
df_alpha <- data.frame(aaSeqCDR3 = character())
df_beta <- data.frame(aaSeqCDR3 = character())
for (i in DATA$meta$Sample[DATA$meta$Group != "test" & DATA$meta$Point != "K06" & DATA$meta$Point != "K03" & DATA$meta$Group != "IBD"]) {
  if (grepl("α", i)) {
    data_list[[paste0(i)]] <- read.table(paste0("/home/sikamov/TCR/res/", i, "/", i, ".clones_TRA.tsv"), header = TRUE, sep = "\t")
    data_df <- read.table(paste0("/home/sikamov/TCR/res/", i, "/", i, ".clones_TRA.tsv"), header = TRUE, sep = "\t")
    df <- data_df[,c("aaSeqCDR3", "uniqueMoleculeFraction")]
    df <- top_n(df,5,uniqueMoleculeFraction)
    names(df) <- c("aaSeqCDR3", paste0(i))
    df_alpha <- merge(df, df_alpha, by = "aaSeqCDR3", all = TRUE)
  } else if (grepl("β", i)) {
    data_list[[paste0(i)]] <- read.table(paste0("/home/sikamov/TCR/res/", i, "/", i, ".clones_TRB.tsv"), header = TRUE, sep = "\t")
    data_df <- read.table(paste0("/home/sikamov/TCR/res/", i, "/", i, ".clones_TRB.tsv"), header = TRUE, sep = "\t")
    df <- data_df[, c("aaSeqCDR3", "uniqueMoleculeFraction")]
    df <- top_n(df,5,uniqueMoleculeFraction)
    names(df) <- c("aaSeqCDR3", paste0(i))
    df_beta <- merge(df, df_beta, by = "aaSeqCDR3", all = TRUE)
  }
}

top_a_UMI_ind <- df_alpha$aaSeqCDR3
top_b_UMI_ind <- df_beta$aaSeqCDR3

df <- df_beta
df <- df[!duplicated(df$aaSeqCDR3), ]
rownames(df) <- df$aaSeqCDR3
df$aaSeqCDR3 <- NULL
df <- df[rowSums(!is.na(df)) > 1, ]

df_a <- df[rowSums(is.na(df[, grepl("K", colnames(df))])) == sum(grepl("K", colnames(df))), ]
df_b <- df[rowSums(is.na(df[, grepl("K", colnames(df))])) == sum(grepl("K", colnames(df))), ]
top_a_UMI<-rownames(df_a)
top_b_UMI<-rownames(df_b)

df[is.na(df)] <- 0
df<- df[, sapply(df, is.numeric)]
df <- as.matrix(df)
rownames(df) <- labels_all[rownames(df)]

annotations <- data.frame(Group = meta_b$Group)
rownames(annotations) <- meta_b$Sample
 
pheatmap(
  df,
  annotation_col = annotations,
  annotation_colors = list(Group = c("CD" = "firebrick3", "HC" = "#238B37", "UC" = "dodgerblue3")),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  #border_color = NA,
 # border = TRUE,
  main = "",
  fontsize = 14,
  fontsize_row = 14,
  fontsize_col = 14,
  annotation_legend = TRUE,
  colorRampPalette(c("white", "#238B37"))(100),  
  fontfamily = "serif",
)


#Перепредставленные клонотипы
tc <- trackClonotypes(data_b, top_b_UMI) 
tc <- tc[rowSums(tc[,19:29]==0) == 11, ]
tc[,19:29] <- NULL 

#a
labels <- c(
  "CAESKEKGGSYIPTF" = "CAESKEKGGSYIPTF (ЦМВ) ",
  "CAGIGGTSYGKLTF" = "CAGIGGTSYGKLTF (Вирус гриппа А/ЦМВ) ",
  "CALGENTGNQFYF" = "CALGENTGNQFYF (ЦМВ) ",
  "CAQLIGFGNVLHC" = "CAQLIGFGNVLHC (-) ",
  "CATPMVTSGAGSYQLTF" = "CATPMVTSGAGSYQLTF (-) ",
  "CAVDNTSGTYKYIF" = "CAVDNTSGTYKYIF (ЦМВ) ",
  "CAVDVDSSYKLIF" = "CAVDVDSSYKLIF (ЦМВ) ",
  "CAVEAAGNKLTF" = "CAVEAAGNKLTF (ЦМВ) ",
  "CAVEPPSGTYKYIF" = "CAVEPPSGTYKYIF (ЦМВ)",
  "CAVGAYSGGGADGLTF" = "CAVGAYSGGGADGLTF (ЦМВ) ",
  "CAVYHQAGTALIF" = "CAVYHQAGTALIF (Вирус гриппа А*) ",
  "CAYRTPDGQKLLF" = "CAYRTPDGQKLLF (-) ",
  "CGMIGGSNYKLTF" = "CGMIGGSNYKLTF (ЦМВ) ",
  "CGTEMRNDYKLSF" = "CGTEMRNDYKLSF (-) "
)

#b
labels <- c(
  "CAISPSFYDPDTQYF" = "CAISPSFYDPDTQYF (-) ",
  "CASQDRNRYNEKLFF" = "CASQDRNRYNEKLFF (-) ",
  "CASSFRDSNTQYF" = "CASSFRDSNTQYF (Пшеница (Глютен)) ",
  "CASSFTSAGNNEQFF" = "CASSFTSAGNNEQFF (ЦМВ) ",
  "CASSGPTGGNYEQYF" = "CASSGPTGGNYEQYF (ВИЧ-1/ЦМВ) ",
  "CASSLLVA_GLNEQFF" = "CASSLLVA_GLNEQFF (-) ",
  "CASSLQAGANEQYF" = "CASSLQAGANEQYF (Вирус денге) ",
  "CASSLRGHYFSNTEAFF" = "CASSLRGHYFSNTEAFF (-) ",
  "CASSMFATEAFF" = "CASSMFATEAFF (-) ",
  "CASSQSNQPQHF" = "CASSQSNQPQHF (Человек (SEC24A)) ",
  "CASSRLAGGEDTQYF" = "CASSRLAGGEDTQYF (ЦМВ/ВИЧ-1) ",
  "CATSREIGWEQFF" = "CATSREIGWEQFF (-) ",
  "CSARESWTGGYTF" = "CSARESWTGGYTF (-) ",
  "CSATVDSATNEKLFF" = "CSATVDSATNEKLFF (-) "
)


vis(tc) +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  theme(
    text = element_text(size = 18, family = "Times New Roman"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)) +
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(11, "Spectral"))(nrow(tc)),
    labels = labels,
    name = "Клонотип (Аннотация VDJdb)"
  )

#Клонотипы интереса
interest <- c("CASSFRDSNTQYF","CASSQSNQPQHF")

labels <- c(
  "CASSFRDSNTQYF" = "CASSFRDSNTQYF (Глютен) ",
  "CASSQSNQPQHF" = "CASSQSNQPQHF (SEC24A) "
)

tc <- trackClonotypes(data_b, interest) 
tc <- tc[rowSums(tc[,19:29]==0) == 11, ]
tc <- subset(tc, select = which(colSums(tc != 0) > 0))

vis(tc) +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  theme(
    text = element_text(size = 14, family = "Times New Roman"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top") +
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(11, "Spectral"))(nrow(tc)),
    labels = labels,
    name = "Клонотип (Аннотация VDJdb)"
  )


#Известные клонотипы
top_known <- c("CAVSSNDYKLSF","CAVMNYGGSQGNLIF", "CASSLRNTGELFF","CAASGGGSYIPTF")
top_known <- c("CAVEDSNYQLIW")

#
labels <- c(
  "CAVSSNDYKLSF" = "CAVSSNDYKLSF (Вирус гриппа А)",
  "CAVMNYGGSQGNLIF" = "CAVMNYGGSQGNLIF (ЦМВ/SARS-CoV-2)",
  "CASSLRNTGELFF" = "CASSLRNTGELFF (ЦМВ/ВЭБ)",
  "CAASGGGSYIPTF" = "CAASGGGSYIPTF (ЦМВ)",
  "CAVEDSNYQLIW" = "CAVEDSNYQLIW (ВЭБ/Человек (BST2))"
)

tc <- trackClonotypes(data, top_known) 
tc <- subset(tc, select = which(colSums(tc != 0) > 0))

vis(tc) +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  theme(
    text = element_text(size = 14, family = "Times New Roman"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right") +
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(11, "Spectral"))(nrow(tc)),
    labels = labels,
    name = "Клонотип (Аннотация VDJdb)"
  )


#Индивидуальные профили TCR реперутаров
trackClonotypes(repFilter(DATA, .method = "by.meta", 
                .query = list(Group = include(group), 
                              Point = exclude(point),
                              Patient = include("I"),
                              Chain = include ("α")))$data, list(1,5)) %>%  vis() +
  labs(title = NULL, subtitle = NULL, y = "Доля в репертуаре", x = NULL) +
  theme(
    text = element_text(size = 14, family = "Times New Roman"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)) +
  scale_fill_manual(
    values = colorRampPalette(rev(brewer.pal(11, "Greens")))(5),
    labels = labels_all,
    name = "Клонотип (Аннотация VDJdb)"
  )

#"Blues" "Oranges" "Greens" "Purples" "Reds" "Yellows" "Greys"
#A C E K H G I

# Аннотация клонотипов
#a
labels_all <- c(
  "CAAHARDTGRRALTF" = "CAAHARDTGRRALTF (-))",
  "CAAISGGYNKLIF" = "CAAISGGYNKLIF (Вирус гриппа А)",
  "CAAQGAQKLVF" = "CAAQGAQKLVF (ЦМВ)",
  "CAARGEDAGNMLTF" = "CAARGEDAGNMLTF (-)",
  "CAARKGKLIF" = "CAARKGKLIF (-)",
  "CAASAKNYGQNFVF" = "CAASAKNYGQNFVF (ЦМВ)",
  "CAASARSGGGADGLTF" = "CAASARSGGGADGLTF (ЦМВ*)",
  "CAASEVYNNNDMRF" = "CAASEVYNNNDMRF (-)",
  "CAASGGGSQGNLIF" = "CAASGGGSQGNLIF (ЦМВ)",
  "CAASPTNFGNEKLTF" = "CAASPTNFGNEKLTF (ЦМВ)",
  "CAASRTGFQKLVF" = "CAASRTGFQKLVF (-)",
  "CAAYNNAGNMLTF" = "CAAYNNAGNMLTF (ЦМВ*)",
  "CAESKEKGGSYIPTF" = "CAESKEKGGSYIPTF (ЦМВ)",
  "CAESNKLVF" = "CAESNKLVF (-)",
  "CAETPSGGATNKLIF" = "CAETPSGGATNKLIF (ЦМВ)",
  "CAFIEMNRDDKIIF" = "CAFIEMNRDDKIIF (ЦМВ)",
  "CAFMGYNNNDMRF" = "CAFMGYNNNDMRF (ЦМВ)",
  "CAGDADSGGGADGLTF" = "CAGDADSGGGADGLTF (ЦМВ)",
  "CAGFVRDSNYQLIW" = "CAGFVRDSNYQLIW (ЦМВ)",
  "CAGIGGTSYGKLTF" = "CAGIGGTSYGKLTF (Вирус гриппа А*)",
  "CAGPDGNARLMF" = "CAGPDGNARLMF (ЦМВ)",
  "CAGPTVFGRF" = "CAGPTVFGRF (-)",
  "CAGQPNTGNQFYF" = "CAGQPNTGNQFYF (-)",
  "CAGYGSDKYIF" = "CAGYGSDKYIF (-)",
  "CAIHTGGFKTIF" = "CAIHTGGFKTIF (ЦМВ)",
  "CALEGYGNKLVF" = "CALEGYGNKLVF (ЦМВ)",
  "CALGDSGGSNYKLTF" = "CALGDSGGSNYKLTF (ВЭБ/ЦМВ/SARS-CoV-2)",
  "CALGENTGNQFYF" = "CALGENTGNQFYF (-)",
  "CALGERGGSNYKLTF" = "CALGERGGSNYKLTF (ВЭБ/ЦМВ/SARS-CoV-2)",
  "CALRAPTDSWGKLQF" = "CALRAPTDSWGKLQF (-)",
  "CALRNTGGFKTIF" = "CALRNTGGFKTIF (ЦМВ/ВГС)",
  "CALSDPGGSEKLVF" = "CALSDPGGSEKLVF (ЦМВ)",
  "CALSTNDYKLSF" = "CALSTNDYKLSF (-)",
  "CAMRDPTSGTYKYIF" = "CAMRDPTSGTYKYIF (-)",
  "CAMREGRYGGSQGNLIF" = "CAMREGRYGGSQGNLIF (Вирус гриппа А)",
  "CAMSGLGSARQLTF" = "CAMSGLGSARQLTF (ЦМВ)",
  "CAMSSDNYGQNFVF" = "CAMSSDNYGQNFVF (ЦМВ)",
  "CAPGVRDTGRRALTF" = "CAPGVRDTGRRALTF (-)",
  "CAPWESGSEKLVF" = "CAPWESGSEKLVF (-)",
  "CAQLIGFGNVLHC" = "CAQLIGFGNVLHC (-)",
  "CARERGSYIPTF" = "CARERGSYIPTF (-)",
  "CASLDSWGKLQF" = "CASLDSWGKLQF (ЦМВ*/ВЭБ)",
  "CASNFGNEKLTF" = "CASNFGNEKLTF (ЦМВ)",
  "CASPEAAGNKLTF" = "CASPEAAGNKLTF (ЦМВ)",
  "CATPMVTSGAGSYQLTF" = "CATPMVTSGAGSYQLTF (-)",
  "CATPSGGYNKLIF" = "CATPSGGYNKLIF (Вирус гриппа А/ЦМВ)",
  "CATVKANQAGTALIF" = "CATVKANQAGTALIF (-)",
  "CATYSGGGADGLTF" = "CATYSGGGADGLTF (Вирус гриппа А/Человек (BST2))", #BST2
  "CAVDNTSGTYKYIF" = "CAVDNTSGTYKYIF (ЦМВ)",
  "CAVDVDSSYKLIF" = "CAVDVDSSYKLIF (ЦМВ)",
  "CAVEAAGNKLTF" = "CAVEAAGNKLTF (ЦМВ*)",
  "CAVENFNKFYF" = "CAVENFNKFYF (ЦМВ)",
  "CAVEPPSGTYKYIF" = "CAVEPPSGTYKYIF (ЦМВ)",
  "CAVERASGNTPLVF" = "CAVERASGNTPLVF (-)",
  "CAVGAGGYSTLTF" = "CAVGAGGYSTLTF (-)",
  "CAVGAYSGGGADGLTF" = "CAVGAYSGGGADGLTF (ЦМВ)",
  "CAVNWNSGGYQKVTF" = "CAVNWNSGGYQKVTF (Вирус гриппа А)",
  "CAVQATTGNQFYF" = "CAVQATTGNQFYF (ЦМВ)",
  "CAVQEGSSGSARQLTF" = "CAVQEGSSGSARQLTF (-)",
  "CAVRDSNYQLIW" = "CAVRDSNYQLIW (ВЭБ*/Человек (BST2*)/Mycobacterium tuberculosis*)",
  "CAVREDNTDKLIF" = "CAVREDNTDKLIF (SARS-CoV-2/ЦМВ)",
  "CAVRGFSDGQKLLF" = "CAVRGFSDGQKLLF (-)",
  "CAVRSTGANSKLTF" = "CAVRSTGANSKLTF (ЦМВ)",
  "CAVRTRESSNTGKLIF" = "CAVRTRESSNTGKLIF (-)",
  "CAVRTYLDTGRRALTF" = "CAVRTYLDTGRRALTF (-)",
  "CAVSGYAGNNRKLIW" = "CAVSGYAGNNRKLIW (-)",
  "CAVSPLYNARLMF" = "CAVSPLYNARLMF (-)",
  "CAVTDSNYQLIW" = "CAVTDSNYQLIW (Человек (BST2*)/ВЭБ)",
  "CAVTSKLTF" = "CAVTSKLTF (-)",
  "CAVVTNAGKSTF" = "CAVVTNAGKSTF (ЦМВ)",
  "CAVYHQAGTALIF" = "CAVYHQAGTALIF (Вирус гриппа А*)",
  "CAYIMDSNYQLIW" = "CAYIMDSNYQLIW (Вирус гриппа А/ЦМВ)",
  "CAYLNTDKLIF" = "CAYLNTDKLIF (ЦМВ)",
  "CAYRGGGADGLTF" = "CAYRGGGADGLTF (Человек (BST2)/ЦМВ)",
  "CAYRGRRGSTLGRLYF" = "CAYRGRRGSTLGRLYF (ЦМВ)",
  "CAYRIGAGGTSYGKLTF" = "CAYRIGAGGTSYGKLTF (Вирус гриппа А)",
  "CAYRTPDGQKLLF" = "CAYRTPDGQKLLF (-)",
  "CAYWPDMRF" = "CAYWPDMRF (-)",
  "CGMIGGSNYKLTF" = "CGMIGGSNYKLTF (ЦМВ)",
  "CGTEMRNDYKLSF" = "CGTEMRNDYKLSF (-)",
  "CIVPLRDGIIF" = "CIVPLRDGIIF (-)",
  "CIVRDSYSGAGSYQLTF" = "CIVRDSYSGAGSYQLTF (-)",
  "CLVG*KA_GDKLTF" = "CLVG*KA_GDKLTF (-)",
  "CLVGDVNRDDKIIF" = "CLVGDVNRDDKIIF (-)",
  "CVVIHSNSGYALNF" = "CVVIHSNSGYALNF (-)",
  "CVVIRYGNNRLAF" = "CVVIRYGNNRLAF (-)",
  "CVVNIPSATDKLIF" = "CVVNIPSATDKLIF (-)",
  "CVVRNSGYALNF" = "CVVRNSGYALNF (-)",
  "CVVSDQGFQKLVF" = "CVVSDQGFQKLVF (-)",
  "CVVTLNTGFQKLVF" = "CVVTLNTGFQKLVF (ЦМВ)",
  "CSARILG_GDTEAFF" = "CSARILG_GDTEAFF (-)",
  "CAGQGAQKLVF" = "CAGQGAQKLVF (ЦМВ)",
  "CILRDVNTGFQKLVF" = "CILRDVNTGFQKLVF (-)",
  "CAADFGHDKVIF" = "CAADFGHDKVIF (-)",
  "CAGMDSNYQLIW" = "CAGMDSNYQLIW (Человек (BST2*))",
  "CALSDDTNAGKSTF" = "CALSDDTNAGKSTF (-)",
  "CAVDNAGNNRKLIW" = "CAVDNAGNNRKLIW (-)",
  "CALTPAPSGYSTLTF" = "CALTPAPSGYSTLTF (-)",
  "CILRDISGYNKLIF" = "CILRDISGYNKLIF (-)",
  "CAANAGGSNYKLTF" = "CAANAGGSNYKLTF (ЦМВ)"
)

#b
labels_all <- c(
  "CSARILG_GDTEAFF" = "CSARILG_GDTEAFF (-)",
  "CAIRTGSYNEQFF" = "CAIRTGSYNEQFF (-)",
  "CAISELETTQPQHF" = "CAISELETTQPQHF (-)",
  "CAISEPRGSYEQYF" = "CAISEPRGSYEQYF (-)",
  "CAISESKTGNSPLHF" = "CAISESKTGNSPLHF (-)",
  "CAISPSFYDPDTQYF" = "CAISPSFYDPDTQYF (-)",
  "CALRGALNSDEQFF" = "CALRGALNSDEQFF (-)",
  "CASGSPETSGSFF" = "CASGSPETSGSFF (-)",
  "CASIRKQPQHF" = "CASIRKQPQHF (-)",
  "CASKAGGGGHTDTQYF" = "CASKAGGGGHTDTQYF (-)",
  "CASKEQGENEKLFF" = "CASKEQGENEKLFF (-)",
  "CASNLETLAKNIQYF" = "CASNLETLAKNIQYF (-)",
  "CASQAQGSYEQYF" = "CASQAQGSYEQYF (-)", 
  "CASQDRNRYNEKLFF" = "CASQDRNRYNEKLFF (-)",
  "CASRFSENEQFF" = "CASRFSENEQFF (-)",
  "CASRGLSTDTQYF" = "CASRGLSTDTQYF (Пшеница (Глютен))",  ###??
  "CASRPGTDQYNEQFF" = "CASRPGTDQYNEQFF (ЦМВ/SARS-CoV-2)",
  "CASRSGEDEQFF" = "CASRSGEDEQFF (-)",
  "CASSAHGEGETQYF" = "CASSAHGEGETQYF (ЦМВ)",
  "CASSAPYNEQFF" = "CASSAPYNEQFF (ЦМВ)",
  "CASSDGLGSYEQYF" = "CASSDGLGSYEQYF (Вирус гриппа А)",
  "CASSDTGIGNQPQHF" = "CASSDTGIGNQPQHF (ВИЧ-1/Вирус гриппа А)",
  "CASSEGTASTDTQYF" = "CASSEGTASTDTQYF (ЦМВ)",
  "CASSFDPPPTKLETQYF" = "CASSFDPPPTKLETQYF (-)",
  "CASSFEQGYTYNEQFF" = "CASSFEQGYTYNEQFF (ЦМВ)",
  "CASSFGGGSPLHF" = "CASSFGGGSPLHF (Человек (MLANA))",
  "CASSFRDSNTQYF" = "CASSFRDSNTQYF (Пшеница (Глютен))",
  "CASSFTSAGNNEQFF" = "CASSFTSAGNNEQFF (ЦМВ)",
  "CASSGDSTDTQYF" = "CASSGDSTDTQYF (Пшеница (Глютен))",
  "CASSGGGSGFYNEQFF" = "CASSGGGSGFYNEQFF (ЦМВ)",
  "CASSGPTGGNYEQYF" = "CASSGPTGGNYEQYF (ВИЧ-1/ЦМВ)",
  "CASSLAGGAGYNEQFF" = "CASSLAGGAGYNEQFF (ЦМВ)",
  "CASSLAGSGELFF" = "CASSLAGSGELFF (ЦМВ/Вирус гриппа А)",
  "CASSLAGTGFYEQYF" = "CASSLAGTGFYEQYF (ЦМВ)",
  "CASSLARANEAFF" = "CASSLARANEAFF (-)",
  "CASSLEIFNEQFF" = "CASSLEIFNEQFF (-)",
  "CASSLFAAVQETQYF" = "CASSLFAAVQETQYF (-)",
  "CASSLGQGQSRQPQHF" = "CASSLGQGQSRQPQHF (ЦМВ)",
  "CASSLGTANNYGYTF" = "CASSLGTANNYGYTF (ВЭБ)",
  "CASSLIGNEQFF" = "CASSLIGNEQFF (ВГС)",
  "CASSLLVA_GLNEQFF" = "CASSLLVA_GLNEQFF (-)",
  "CASSLQAGANEQYF" = "CASSLQAGANEQYF (Вирус денге)",
  "CASSLRGHYFSNTEAFF" = "CASSLRGHYFSNTEAFF (-)",
  "CASSLSSGGTSYEQYF" = "CASSLSSGGTSYEQYF (-)",
  "CASSLSVGGTDYEQYF" = "CASSLSVGGTDYEQYF (-)",
  "CASSLVSRPGASEQFF" = "CASSLVSRPGASEQFF (-)",
  "CASSLWTRGSDEQYF" = "CASSLWTRGSDEQYF (ЦМВ/ВИЧ-1)",
  "CASSLYGGVQETQYF" = "CASSLYGGVQETQYF (ВИЧ-1)",
  "CASSMFATEAFF" = "CASSMFATEAFF (-)",
  "CASSPGTSGNNEQFF" = "CASSPGTSGNNEQFF (ВИЧ-1)",
  "CASSPGWRPSGANVLTF" = "CASSPGWRPSGANVLTF (-)",
  "CASSPLGQPNLQETQYF" = "CASSPLGQPNLQETQYF (-)",
  "CASSPPIAGDTQYF" = "CASSPPIAGDTQYF (-)",    #Plasmodium berghei
  "CASSPPSGSYEQYF" = "CASSPPSGSYEQYF (Вирус гриппа А)",
  "CASSPPTGGQETQYF" = "CASSPPTGGQETQYF (ЦМВ)",
  "CASSPRGQTPRNYEQYF" = "CASSPRGQTPRNYEQYF (-)",
  "CASSPRLAGAHSTDTQYF" = "CASSPRLAGAHSTDTQYF (-)",
  "CASSPSG_SG*EQFF" = "CASSPSG_SG*EQFF (-)",
  "CASSPSGLAGRYF" = "CASSPSGLAGRYF (-)",
  "CASSQDSNVQEAFF" = "CASSQDSNVQEAFF (-)",
  "CASSQGTYEQYF" = "CASSQGTYEQYF (ЦМВ)",
  "CASSQSNQPQHF" = "CASSQSNQPQHF (Человек (SEC24A))",
  "CASSQVDLVFSGANVLTF" = "CASSQVDLVFSGANVLTF (-)",
  "CASSRGGGRNTIYF" = "CASSRGGGRNTIYF (ВИЧ-1)",
  "CASSRLAGGEDTQYF" = "CASSRLAGGEDTQYF (ЦМВ/ВИЧ-1)",
  "CASSRLTGFQYF" = "CASSRLTGFQYF (-)",
  "CASSRLYGPSDEAFF" = "CASSRLYGPSDEAFF (-)",
  "CASSSAYTIYF" = "CASSSAYTIYF (-)",
  "CASSSGTTEETQYF" = "CASSSGTTEETQYF (ВИЧ-1/ЦМВ)",
  "CASSSPWDREVMNTEAFF" = "CASSSPWDREVMNTEAFF (-)",
  "CASSSRDYNNNEQFF" = "CASSSRDYNNNEQFF (Trypanosoma cruzi)",
  "CASSSRGTGELFF" = "CASSSRGTGELFF (ЦМВ)",
  "CASSSVGYANTGELFF" = "CASSSVGYANTGELFF (ЦМВ/ВЭБ)",
  "CASSVERGYSYNEQFF" = "CASSVERGYSYNEQFF (ЦМВ)",
  "CATSREIGWEQFF" = "CATSREIGWEQFF (-) ",
  "CSARESWTGGYTF" = "CSARESWTGGYTF (-) ",
  "CSATVDSATNEKLFF" = "CSATVDSATNEKLFF (-) ",
  "CASSWTGNQPQHF" = "CASSWTGNQPQHF (-) ",
  "CASSYSDSGGNQPQHF" = "CASSYSDSGGNQPQHF (-) ",
  "CSVDLTTDTQYF" = "CSVDLTTDTQYF (-) ",
  "CAIRDSNTEAFF" = "CAIRDSNTEAFF (-) ",
  "CASRRGGEGSNTEAFF" = "CASRRGGEGSNTEAFF (-) ",
  "CASSGILAGNEQYF" = "CASSGILAGNEQYF (-) ",
  "CATSDGTSDRTDTQYF" = "CATSDGTSDRTDTQYF (ВИЧ-1) ",
  "CAWGMTSGRYTGELFF" = "CAWGMTSGRYTGELFF (-) ",
  "CASSESDRSSYEQYF" = "CASSESDRSSYEQYF (ЦМВ/ВИЧ-1) ",
  "CASSLAFDEKLFF" = "CASSLAFDEKLFF (Вирус гриппа А) ",
  "CASSLGLSRELFF" = "CASSLGLSRELFF (ЦМВ) ",
  "CSARNGLNEKLFF" = "CSARNGLNEKLFF (-) ",
  "CASQAWVGIHEQFF" = "CASQAWVGIHEQFF (-) ",
  "CASSDRGGFGEKLFF" = "CASSDRGGFGEKLFFs (-) ",
  "CASSVRKGGNTEAFF" = "CASSVRKGGNTEAFF (ЦМВ) ",
  "CATFDGNTGELFF" = "CATFDGNTGELFF (SARS-CoV-2) ",
  "CSASVRGRDGSYEQYF" = "CSASVRGRDGSYEQYF (-) ",
  "CASSLRAGVDYEQYF" = "CASSLRAGVDYEQYF (ЦМВ) ",
  "CASSPLIVNTEAFF" = "CASSPLIVNTEAFF (ЦМВ) ",
  "CASSQSGGGMNTEAFF" = "CASSQSGGGMNTEAFF (-) ",
  "CASTQNSNQPQHF" = "CASTQNSNQPQHF (Вирус гриппа А) "
)


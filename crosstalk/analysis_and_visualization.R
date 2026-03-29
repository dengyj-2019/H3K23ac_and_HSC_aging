library(CellChat)
library(Seurat)
library(openxlsx)
library(future)
library(future.apply)
library(Matrix)
library(openxlsx)

setwd('/data1/02.private/dengyj/analysis/HSC/crosstalk/young')
source('support/calculation.R')
source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

database <- CellChatDB.mouse$interaction[, c('interaction_name', 'ligand', 'receptor')]
database$interaction_name <- paste0(database$ligand, '~', database$receptor)
complex_info <- CellChatDB.mouse$complex
complex_names <- rownames(complex_info)
complex_info <- lapply(complex_names, function(complex){
  component <- as.character(complex_info[rownames(complex_info) %in% complex,
                                         grep('subunit', colnames(complex_info))])
  component <- component[component != ""]
})
names(complex_info) <- complex_names

combined_niche <- readRDS('/data1/02.private/dengyj/analysis/HSC/Integrated/aging_niche_combined/combined_niche.rds')
GSE210584_mature_hema <- readRDS('/data1/02.private/dengyj/analysis/HSC/Integrated/GSE210584/GSE210584_mature_hema.rds')

load('/data1/02.private/dengyj/analysis/HSC/HSC_2M/HSC_2M.RData')


HSC_2M_filtered$age_group <- 'Young'


combined_BM <-merge(HSC_2M_filtered, list(combined_niche[, combined_niche$age_group %in% c('Young')],
                                          GSE210584_mature_hema[, GSE210584_mature_hema$age_group %in% c('Young')] ))
combined_BM$identity <- Idents(combined_BM)

young_crosstalk <- pre_prossessing(expm1(combined_BM$RNA@data), database = database,
                                   identity = combined_BM$identity, condition = combined_BM$age_group, 
                                   complex_info = complex_info)

young_crosstalk <- get_permutation2(young_crosstalk,  multiprocess=F,iteration = 1000,  block_size = 50, slot_used = 'identity',
                                    ligand_pct_threshold = 0.03, receptor_pct_threshold = 0.03, style = 'product')


saveRDS(young_crosstalk, file ='young_crosstalk.rds')




#####plotting
library(Seurat)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(patchwork)
library(aplot)
library(openxlsx)
library(ComplexHeatmap)
library(gridBase)
library(circlize)
library(svglite)
library(Matrix)


lbls_func <- function(x){
    x <- gsub('MLR-qHSC', 'Meg3+-HSC', x)
    x <- gsub('MK-qHSC', 'Pf4+-HSC', x)
    x <- gsub('LR-qHSC', 'Clec4d+-HSC', x)
    x <- gsub('LR-aHSC', 'Rhd+Mt2lo-HSC', x)
    x <- gsub('R-aHSC', 'Rhd+Mt2hi-HSC', x)
    x <- gsub('L-aHSC', 'Pebp1hi-HSC', x)
    x <- gsub('Other', 'Imm_GRAN', x)
}

age_colors <- c('2M'='#42827d','8M'='#dfcc28','28M'='#df886d')

cluster_colors <- c('MLR-qHSC'='#91D0D4', 'MK_1'='#587DB2', 'MK_2'='#00aaff', 'MK-qHSC' = '#587DB2', 
                    'B-primed progenitor1'='#097cc5','B-primed progenitor2'='#f02463',
                    'Procr+MPP'='#fa9041','9'='#097cc5',
'TSP-like'='#00b659','Lymphoid-primed progenitor'='#a86bac','Myeloid-primed progenitor'='#ffd800',
'IFN cell'='#476d87','Histone cell'='#fa92ba', 
'LR-aHSC'='#E9A342', 'R-aHSC'='#A35592', 'L-aHSC'='#42827E',
 #'DC-qHSC'='#B7814C', 
                    'LR-qHSC'='#DBACC2', 
                    'LR/R-aHSC' = '#fa9041',
                    '8' = '#00b659',
                   "Nr4a1+MLR"="#4279C0","Cd34+MLR"="#A8B0E3",#"Histone+"="#71906D",
                    "IFN_MLR"="#e6ba77",
"Pf4_lo"='#cdc97d', #'#cecc90',"#c3a8a0",
                    "Pf4_hi"="#ff927c","Car1+"="#ee587d","Txnip+MLR"="#42BACA", 'Histone+'='#1036bb',#'#9672ac',
                    'Other'='gray80', 'From_8M0'='#53A85F', 
 'From_8M1'='#F3B1A0', 'From_8M2'='#57C3F3', 'From_8M3'='#E95C59', 'From_28M1'='#3c97c4',
  'From_28M0'='#c4c04f', 'From_28M2'='#ff7e00', 'From_28M3'='#d73e4b', 
                    'CITE_seq1'='#9fbe8c', 'CITE_seq2' = '#ceba74', 
                   "MSC"="#E5D2DD","EC-sinusoidal"="#53A85F","EC-arteriolar"='#625D9E',#"#F1BB72",
"EC-arterial"="#D6E7A3","osteoprogenitor"="#57C3F3","OLC-BM"="#476D87",
"CL-PC"="#E95C59","OLC-PSC"='#CCC9E6',#'#F3B1A0',#"#E59CC4",
                    "Chondro"='#712820',#"#AB3282",
                    "FL-PC"="#23452F",
"Pericyte"="#BD956A")


stromal_class<-c('MSC','EC-arterial','EC-arteriolar','EC-sinusoidal','osteoprogenitor','OLC-BM','CL-PC','OLC-PSC','Chondro','FL-PC','Pericyte')
Immunit_class<-c('Bcell','DC','early_Ery','Granulocyte','late_Ery','MK','Tcell')

levels(young_crosstalk@identity)

tbl <- table(young_crosstalk@identity)

cluster_filtered <- names(tbl)[tbl <= 30]


niche_idents <- niche_idents_lvls <- c('MSC',  'EC-sinusoidal','EC-arteriolar', 'EC-arterial', 
 'osteoprogenitor', 'OLC-BM', 'CL-PC', 'OLC-PSC','Chondro', 
  'FL-PC', 'Pericyte', 'Tcell', 'Bcell', 
                                                             'DC', 'Granulocyte', 
                                                             'early_Ery', 'late_Ery', 
                                                             'MK')

HSC_idents <- HSC_idents_lvls <- 
c('MLR-qHSC', 'MK-qHSC', 'LR-qHSC', 'LR-aHSC', 'R-aHSC', 'L-aHSC')                                                             


idents_result <- get_sig_LR_result(young_crosstalk, slot_used = 'identity', ####胞间互作阈值设定为0.01
                                   signal = 'intercellular',inter_p_threshold = 0.01)

HSC_Rec_idents_result <- idents_result[idents_result$ligand_cluster %in% setdiff(niche_idents, cluster_filtered) & 
                       idents_result$receptor_cluster %in% setdiff(HSC_idents, cluster_filtered), ]

HSC_Rec_idents_result_intercellular <- get_plot_df(object = young_crosstalk, clusters_pairs = 
                                             unique(HSC_Rec_idents_result$clu_pairs),
                                do_scale = T,do_log = F,
                              LR_pairs = unique(HSC_Rec_idents_result$interaction_name), 
                             data_to_plot  = 'intercellular')

HSC_Rec_idents_result_pnull <- get_plot_df(object = young_crosstalk, clusters_pairs = 
                                             unique(HSC_Rec_idents_result$clu_pairs),
                                do_scale = F,do_log = T,
                              LR_pairs = unique(HSC_Rec_idents_result$interaction_name), 
                             data_to_plot  = 'pnull')



HSC_Rec_idents_result_plot_df <- merge(HSC_Rec_idents_result_intercellular, HSC_Rec_idents_result_pnull)

HSC_Rec_idents_result_plot_df$clusters_pairs <- 
factor(HSC_Rec_idents_result_plot_df$clusters_pairs, 
      levels = intersect(unlist(lapply(niche_idents_lvls, function(x){paste0(x, '~', HSC_idents_lvls)})), 
                         unique(HSC_Rec_idents_result_plot_df$clusters_pairs)))

HSC_Rec_idents_result_plot_df_filtered <- HSC_Rec_idents_result_plot_df[HSC_Rec_idents_result_plot_df$pnull >= 2, ]


plot_df <- 
data.frame(
    table(HSC_Rec_idents_result_plot_df_filtered$ligand_cluster, HSC_Rec_idents_result_plot_df_filtered$receptor_cluster))
colnames(plot_df) <- c('ligand_cluster', 'receptor_cluster', 'Freq')


options(repr.plot.width=5, repr.plot.height=4)
tmp_plot_df <- plot_df[plot_df$ligand_cluster %in% stromal_class, ]
tmp_plot_df$receptor_cluster <- lbls_func(tmp_plot_df$receptor_cluster)


options(repr.plot.width=6, repr.plot.height=3.5)
ggplot(tmp_plot_df, 
       aes(receptor_cluster, ligand_cluster, fill  = Freq))+
scale_fill_gradientn(colours = c("white", '#5c163e'))+
geom_tile()+
scale_y_discrete(limits = c('EC-arterial','EC-arteriolar','MSC','FL-PC',
                           'CL-PC', 'Chondro','osteoprogenitor',
                               'EC-sinusoidal','Pericyte','OLC-PSC','OLC-BM'))+
scale_x_discrete(limits = rev(c('Meg3+-HSC', 'Pf4+-HSC', 'Clec4d+-HSC', 'Rhd+Mt2lo-HSC', 
              'Rhd+Mt2hi-HSC', 'Pebp1hi-HSC')))+
d_theme_w(size = 10)+
coord_flip()+
labs(x = 'Receptor cluster', y = 'Ligand cluster', fill = 'Counts')+
theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('crosstalk_count.svg', width = 6, height = 3.5)




tmp_df <- HSC_Rec_idents_result_plot_df_filtered[HSC_Rec_idents_result_plot_df_filtered$ligand_cluster %in% 
                                                 stromal_class, ]
tmp_li <- apply((apply(table(tmp_df$receptor_cluster, 
                   tmp_df$interaction_name), 2, function(x){x}) > 0) * 1, 2, 
      function(x){names(x)}[x > 0])
crosstalk_used <- names(tmp_li)[sapply(tmp_li, length) <= 2]
plot_data2 <- HSC_Rec_idents_result_plot_df_filtered[
    HSC_Rec_idents_result_plot_df_filtered$ligand_cluster %in% 
    stromal_class #& HSC_Rec_idents_result_plot_df_filtered$interaction_name %in% crosstalk_used
    , ]
receptor_lvls <- c('MLR-qHSC', 'MK-qHSC', 'LR-qHSC', 
                                         'LR-aHSC', 'R-aHSC', 'L-aHSC','Other')
ligand_lvls <- c('MSC','EC-arterial','EC-arteriolar','EC-sinusoidal','osteoprogenitor',
                 'OLC-BM','CL-PC','OLC-PSC','Chondro','FL-PC','Pericyte')
plot_data2$receptor_cluster_order <- match(plot_data2$receptor_cluster, receptor_lvls)
plot_data2$ligand_cluster_order <- match(plot_data2$ligand_cluster, ligand_lvls)
plot_data2 <- plot_data2[order(plot_data2$receptor_cluster_order, plot_data2$ligand_cluster_order), ]
# plot_data2$clusters_pairs <- factor(plot_data2$clusters_pairs, levels = unique(plot_data2$clusters_pairs))
label_plot_data2 <- data.frame(clusters_pairs = as.character(plot_data2$clusters_pairs))
label_plot_data2$ligand_cluster <- sapply(strsplit(label_plot_data2$clusters_pairs, split = '~'), function(x){x[1]})
label_plot_data2$receptor_cluster <- sapply(strsplit(label_plot_data2$clusters_pairs, split = '~'), function(x){x[2]})

label_plot_data2$ligand_cluster<- factor(label_plot_data2$ligand_cluster, 
                                      levels = unique(label_plot_data2$ligand_cluster) )

label_plot_data2$receptor_cluster<- factor(label_plot_data2$receptor_cluster, 
                                        levels = unique(label_plot_data2$receptor_cluster))
plot_data2$ligand <- sapply(strsplit(plot_data2$interaction_name, split = '~'), 
                         function(x){x[1]})
plot_data2$receptor <- sapply(strsplit(plot_data2$interaction_name, split = '~'), 
                         function(x){x[2]})

lmls_x <- unique(plot_data2$clusters_pairs)
options(repr.plot.width=15, repr.plot.height=0.5)
label_plot_ligand <- ggplot(label_plot_data2)+
geom_tile(color = NA,mapping = aes(clusters_pairs, 1, fill = ligand_cluster))+
scale_fill_manual(values =cluster_colors)+
scale_x_discrete(limits = lmls_x)+
theme(legend.text = element_text(size = 12.5, family='sans'),legend.direction = 'vertical',
           legend.title = element_text(size = 12.5, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
      legend.position = 'bottom',
     text=element_blank(), axis.ticks=element_blank())

label_plot_receptor <- ggplot(label_plot_data2)+
geom_tile(color = NA,mapping = aes(clusters_pairs, 1, fill = receptor_cluster))+
scale_fill_manual(values =cluster_colors, label = lbls_func)+
scale_x_discrete(limits = lmls_x)+
theme(legend.text = element_text(size = 12.5, family='sans'),legend.direction = 'vertical',
           legend.title = element_text(size = 12.5, family='sans'),
     panel.background = element_rect(color = NA, fill = 'white'),
      legend.position = 'bottom',
     text=element_blank(), axis.ticks=element_blank())

p <- ggplot(plot_data2, aes(clusters_pairs,interaction_name, 
                                                   fill = intercellular, size = pnull))+
geom_point(shape=21, stroke=0.2)+
scale_x_discrete(limits = lmls_x)+
scale_y_discrete(limits = rev(unique(plot_data2$interaction_name)))+

scale_fill_gradientn(colours = c('#DAECEB', '#009397'))+
theme(panel.border = element_rect(color ='black', fill =NA),legend.direction = 'horizontal',
      panel.grid.major = element_line(color = 'grey'),
      panel.background = element_rect(color =NA, fill ='white'),
      axis.text.y = element_text(size = 8.5, face = 'bold'), 
      axis.text.x = element_blank())+#element_text(size = 6, angle =45, hjust = 1, face ='bold'))+
labs(fill = 'Interaction', x = NULL, y= NULL)+
scale_size_continuous(range = c(1,3))+
guides(size = guide_legend(title = "-log10 P",nrow = 3, title.position = 'top',
                           override.aes = list(fill = 'black')), 
       fill = guide_colorbar(title = "Intensity",title.position = 'top'))
options(repr.plot.width=12, repr.plot.height=16)
wrap_plots(list(p, label_plot_ligand, label_plot_receptor))+
plot_layout(ncol = 1, guides = 'collect', heights = c(10,0.5,0.5))&
theme(plot.margin = margin(b = 0, l = 0.1,unit = 'cm'))



circlize_data <- plot_data2[plot_data2$interaction_name %in%
                           c('Angpt1~Tek', 'Jag1~Notch1','Jam2~Jam3', 
                               'Esam~Esam','Col2a1~Sdc4',
                               'Col9a1~Sdc4','Ccl2~Ackr1', 
                               'Cd34~Selp', 'Vwf~ITGA2B_ITGB3', 
                               'L1cam~ITGA4_ITGB7','Cxcl12~Cxcr4','Angptl4~ITGA5_ITGB1',
                               'Wnt2~FZD7_LRP6','Wnt11~Fzd7','Wnt5b~Fzd7',
                               'Wnt5a~Fzd7','Icam1~Itgal',  'Ccl8~Ackr1', 'Ccl7~Ackr1', 'Cxcl5~Ackr1', 'Fn1~ITGA5_ITGB1'
                           ), ]


circlize_data <- circlize_data[(circlize_data$interaction_name %in% c('App~Cd74')& 
                               circlize_data$ligand_cluster %in% c('EC-arterial', 'EC-arteriolar')) | 
                               !circlize_data$interaction_name %in% c('App~Cd74'), ]

circlize_data <- circlize_data[,c(#'interaction_name', 
                                  'ligand_cluster', 
                                  'receptor_cluster', 'ligand', 'receptor')]
circlize_data$value <- 1
circlize_data$ligand_combn <- paste0(circlize_data$ligand, '~', circlize_data$ligand_cluster)
circlize_data$receptor_combn <- paste0(circlize_data$receptor, '~', circlize_data$receptor_cluster)

receptor_lvls <- c('MLR-qHSC', 'MK-qHSC', 'LR-qHSC', 
                                         'LR-aHSC', 'R-aHSC', 'L-aHSC','Other')
ligand_lvls <- c('MSC','EC-arterial','EC-arteriolar','EC-sinusoidal','osteoprogenitor',
                 'OLC-BM','CL-PC','OLC-PSC','Chondro','FL-PC','Pericyte')

ligand_order <- circlize_data$ligand_combn[
    order(match(circlize_data$ligand_cluster, ligand_lvls))]


receptor_order <- circlize_data$receptor_combn[
    order(match(circlize_data$receptor_cluster, receptor_lvls))]


order_used <- unique(c(ligand_order, receptor_order))

group_names <- c(circlize_data$ligand_combn, circlize_data$receptor_combn)

group <- sapply(strsplit(group_names, split = '~'), function(x){x[2]})
names(group) <- group_names

col_used <- sapply(strsplit(group_names, split = '~'), function(x){
    cluster_colors[x[2]]
})
names(col_used) <- group_names

clu_used <- unique(c(circlize_data$ligand_cluster, circlize_data$receptor_cluster))

ligand_clu_used <- unique(circlize_data$ligand_cluster)


receptor_clu_used <- unique(circlize_data$receptor_cluster)

lgd1 = Legend(labels =  lbls_func(ligand_clu_used), title = "Ligand Cluster", 
              title_gp = gpar(fontsize = 12.5, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(0.6, "cm"),
              ncol = 1,title_position = 'topleft',
              labels_gp = gpar(fontsize = 12.5, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[ligand_clu_used])))

lgd2 = Legend(labels =  lbls_func(receptor_clu_used), title = "Receptor Cluster", 
              title_gp = gpar(fontsize = 12.5, fontfamily='sans'),
              #grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), 
              row_gap = unit(1, "mm"),
              title_gap = unit(1, "mm"),
              type = "points", pch = 15,size = unit(0.6, "cm"),
              ncol = 1,title_position = 'topleft',
              labels_gp = gpar(fontsize = 12.5, fontfamily = 'sans'), 
              legend_gp = gpar(col = c(cluster_colors[receptor_clu_used])))

pd = packLegend(lgd1, lgd2, direction = "vertical", row_gap = unit(1,'cm'))


options(repr.plot.width=10, repr.plot.height=10)

plot_func <- function(){
circos.clear()
circos.par(gap.degree = 1, canvas.xlim = c(-0.8, 1.5), # 默认大约是 c(-1, 1)
canvas.ylim = c(-1, 1))
chordDiagramFromDataFrame(circlize_data[, c('ligand_combn', 'receptor_combn', 'value')], order = order_used,
             preAllocateTracks = list(track.height = max(strwidth(c(circlize_data$ligandm, 
                                                                    circlize_data$receptor)))*0.8),
                          small.gap = 0.1,
             grid.col = col_used,#group = group,big.gap = 2,
             directional = 1, direction.type = c("arrows"),annotationTrack = "grid" ,
            link.arr.type = "big.arrow")
    circos.track(track.index = 1, 
                 panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ylim[1], gsub('~.*', '', CELL_META$sector.index), 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    }, bg.border = NA)
draw(pd,  x = unit(25, "cm"), y = unit(13, "cm"), just = c("right", "center"))
circos.clear()    
}
plot_func()


svglite("crosstalk_circlize.svg", width = 10, height = 10) # 
plot_func()
dev.off()
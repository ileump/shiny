##------------------
## @knitr 10_setup

knitr::opts_chunk$set(echo = F, warning = F, message = F)


library(magrittr)
library(plyr)
library(ggplot2)
library(grid)
library(extrafont)

source('functions.R')

##--------------------
## @knitr 20_load_data

## Expected variables in the RData file 
## input_params, g_data_params, my_options
load(params$tmp_data)
input <- input_params
g_data <- g_data_params

input$dpi <- as.numeric(input$dpi)

knitr::opts_chunk$set(
    dpi = input$dpi
)

##--------------------------
## @knitr 30-write_summary_1


##---------------------------------
## @knitr 33-metabolites_statistics

p <- g_data$metabolites_statistics

print(p)

##---------------------------------
## @knitr 35-ome_summary

p <- g_data$ome_summary[['ome_summary']]

print(p)


##---------------------------------
## @knitr 36-ome_summary_text

out <- g_data$ome_summary_text

cat(out, '\n\n')

##--------------------
## @knitr 40-pca_scree

pca <- g_data$pca
pca_scree <- g_data$pca_scree

if (!is.null(pca)) {
    grid.newpage()
    grid.draw(pca_scree)
}


##--------------------
## @knitr 45-pca_scree_text

out <- g_data$pca_scree_text

if (!is.null(out)) {
    cat(out, '\n\n')
}


##------------------
## @knitr 50-pca_ind

pca <- g_data$pca
pca_ind <- g_data$pca_ind

if (!is.null(pca)) {
    for (i in 1:length(pca_ind)) {
        gg <- grid_fam(pca_ind[[i]], input$font_family)
        grid.newpage()
        grid.draw(gg)
    }
}



##------------------
## @knitr 60-heatmap

if (!is.null(g_data$heatmap)) {
    # grid.newpage()
    # grid.draw(grid_fam(g_data$pheatmap, input$font_family))
    # print(g_data$heatmap)
    draw(g_data$heatmap, 
         merge_legends = T,
         align_heatmap_legend = 'heatmap_top',
         legend_grouping = 'original'
         )
}

##--------------------------
## @knitr 65-heatmap_subplot
size_params <- input$heatmapSubplot_params

if (is.null(size_params) | length(size_params) == 0) {
    ## To avoid print NULL to report
    cat('')
} else {
    if (length(g_data$heatmap_subplot) > 0) {
        for (i in 1:length(g_data$heatmap_subplot)) {
            w <- size_params[[paste0('width', i)]] / input$dpi
            h <- size_params[[paste0('height', i)]] / input$dpi
            
            sub_chunk <- paste0(
                '\n```{r heatmap_subplot_', i, ', echo=F,', "results='asis',", 'fig.width=', 
                round(w, 3), 
                ', fig.height=',
                round(h, 3),
                '}\n',
                "
# grid.newpage()
# grid.draw(grid_fam(g_data$pheatmap_subplot[[i]], input$font_family))
# g_data$heatmap_subplot[[i]]
draw(g_data$heatmap_subplot[[i]], 
     merge_legends = T,
     align_heatmap_legend = 'heatmap_top',
     legend_grouping = 'original'
     )

```\n
        "
            )
            
            cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
        }
    }
}



##------------------
## @knitr 70-oplsda

oplsda_list <- g_data$oplsda
font_family <- input$font_family
oplsda_label <- input$oplsda_label
parCexMetricN <- 0.8

# encoding <- ifelse(input$file_utf8, 'UTF-8', '')

if(length(oplsda_list) > 0) {
    for(i in 1:length(oplsda_list)) {
        oplsda <- oplsda_list[[i]]
        cat('\n\n### ', paste0(levels(attr(oplsda, 'suppLs')$y),
                            collapse = '-'), '\n\n')
        
        if (inherits(oplsda, 'error')) {
            cat(oplsda$message)
        } else {
            
            if(ropls::getSummaryDF(oplsda)['ort'] > 0) {
                
                sub_chunk <- paste0(
                    '\n```{r oplsda_subchunk_', i, ', echo=F,', "results='asis',", 'fig.width=', 
                    round(input$oplsda_width * 2 / input$dpi, 3), 
                    ', fig.height=',
                    round(input$oplsda_height * 3 / input$dpi, 3),
                    '}\n',
                    "
predicted <- ropls::predict(oplsda)
# set levels to be the same as predicted
reference <- factor(
    attr(oplsda, 'suppLs')$yMCN[, 1],
    levels = levels(predicted)
)

par(mfrow = c(3, 2), family = font_family)
ropls::plot(oplsda, typeVc = 'permutation', parDevNewL = F)
ropls::plot(oplsda, typeVc = 'overview', parDevNewL = F)
ropls::plot(oplsda, typeVc = 'outlier', parDevNewL = F)
if(oplsda_label) {
  ropls::plot(oplsda, typeVc = 'x-score', parDevNewL = F,
              parCexMetricN = parCexMetricN)
} else {
  ropls::plot(oplsda, typeVc = 'x-score', parDevNewL = F, 
              parCexMetricN = parCexMetricN,
              parLabVc = rep('\u25CF', nrow(oplsda@scoreMN)))
}
                
pROC::plot.roc(pROC::roc(attr(oplsda, 'suppLs')$y,
                         attr(oplsda, 'scoreMN')[, 1]),
               print.auc = T, print.auc.cex = 2,
               main = 'ROC')

cat('\n\n')
cat(sprintf('Sensitivity: %.2f; Specificity: %.2f; Negative level: %s; Positive level: %s',
        caret::sensitivity(table(predicted, reference)),
        caret::specificity(table(predicted, reference)),
        levels(predicted)[1],
        levels(predicted)[2]
))
                    
```\n
")
                cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
                
            } else {
                
                sub_chunk <- paste0(
                    '\n```{r oplsda_subchunk_', i, ', echo=F,', "results='asis',", 'fig.width=', 
                    round(input$oplsda_width * 2 / input$dpi, 3), 
                    ', fig.height=',
                    round(input$oplsda_height * 2 / input$dpi, 3),
                    '}\n',
                    "
predicted <- ropls::predict(oplsda)
# set levels to be the same as predicted
reference <- factor(
    attr(oplsda, 'suppLs')$yMCN[, 1],
    levels = levels(predicted)
)

## OPLS-DA with 0 orthogonal comp
par(mfrow = c(2, 2), family = font_family)
ropls::plot(oplsda, typeVc = 'permutation', parDevNewL = F)
ropls::plot(oplsda, typeVc = 'overview', parDevNewL = F)
# ropls::plot(oplsda, typeVc = 'outlier', parDevNewL = F)
# ropls::plot(oplsda, typeVc = 'x-score', parDevNewL = F)

## plot p1 only
d <- data.frame(
    x = jitter(as.numeric(attr(oplsda, 'suppLs')$y),
               amount = 0.3),
    g = attr(oplsda, 'suppLs')$y,
    y = ropls::getScoreMN(oplsda)[, 1],
    labels = rownames(ropls::getScoreMN(oplsda))
)
    
plot(d$x, d$y, type = 'n', xaxt = 'n',
     xlim = c(0, 3),
     xlab = 'Group', ylab = 't1', main = 'Predictive score only')
    
axis(1, at = c(1, 2), 
     labels=levels(attr(oplsda, 'suppLs')$y))
    
if (oplsda_label) {
  text(d$x, d$y, d$labels, 
       col = plyr::mapvalues(
         d$g,
         levels(d$g),
         c('blue', 'red')
       ) %>% as.character)
} else {
  text(d$x, d$y, '\u25CF', 
       col = plyr::mapvalues(
         d$g,
         levels(d$g),
         c('blue', 'red')
       ) %>% as.character)
}

pROC::plot.roc(pROC::roc(attr(oplsda, 'suppLs')$y,
                         attr(oplsda, 'scoreMN')[, 1]),
               print.auc = T, print.auc.cex = 2,
               main = 'ROC')

cat('\n\n')
cat(sprintf('Sensitivity: %.2f; Specificity: %.2f; Negative level: %s; Positive level: %s',
        caret::sensitivity(table(predicted, reference)),
        caret::specificity(table(predicted, reference)),
        levels(predicted)[1],
        levels(predicted)[2]
))
```\n
")

            cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))

            }
        }
    }
}


##------------------
## @knitr 75-oplsda_text

out <- g_data$oplsda_text

cat(out, '\n\n')

##----------------------
## @knitr 80-oplsda_vip



##--------------------------
## @knitr 90-hypothesis_test
## export selected hypothesis tests


    

##-------------------
## @knitr 100-volcano

p_list <- g_data$volcano
if (length(p_list) > 0)
    l_ply(p_list, print)

##-------------------
## @knitr 105-volcano_text

out <- g_data$volcano_text

cat(out, '\n\n')


##-------------------
## @knitr 108-differential_metabolites_text

out <- g_data$differential_metabolites_text

cat(out, '\n\n')

##-------------------
## @knitr 110-boxplot


## add jitter
add_jitter <- input$boxplot_add_jitter

## add point
add_point <- input$boxplot_add_point

## show significance
add_sig <- input$boxplot_add_sig

## show violin
add_violin <- input$boxplot_add_violin



## notch
notch <- input$boxplot_notch

## x label angle
x_angle <- input$boxplot_x_angle
h_just <- input$boxplot_hjust
v_just <- input$boxplot_vjust


## only produce the plots that are shown in shiny
## different from shiny, which produces all the plots at one go
shown_id <- g_data$boxplot_shown_id

## unit
unit <- input$unit

## ggplot theme
ggtheme <- g_data$ggtheme

font_family <- input$font_family


d <- g_data$d1[, shown_id, drop = F]
g <- g_data$group
dh <- g_data$hypothesis_test[shown_id, , drop = F]


posthoc <- NULL

if (add_sig){
    posthoc <- input$boxplot_posthoc
}

if (input$pair_var == 'None') {
    pair <- NA
} else {
    pair <- as.factor(g_data$pair)
}

p_list <- batch_boxplot(d, g, dh, pair, 
                        posthoc = posthoc, unit = unit, 
                        add_point = add_point, add_sig = add_sig, 
                        add_violin = add_violin, add_jitter = add_jitter,
                        notch, x_angle, h_just, v_just,
                        ggtheme = ggtheme, font_family = font_family)


ncols <- input$boxplot_ncol

nrows <- ceiling(length(p_list) / ncols)


if (ncols == 1) {
    p_list <- p_list[shown_id]
    if (length(p_list) > 0)
        l_ply(p_list, print)
} else {
    grid.newpage()
    tmp <- grid.draw(gridExtra::marrangeGrob(p_list, nrow=nrows, ncol=ncols, top = ''))
}

##-------------------
## @knitr 105-boxplot_text

out <- g_data$boxplot_text

cat(out, '\n\n')


##-------------------
## @knitr 120-barplot

## show significance
add_sig <- input$barplot_add_sig

## add jitter
add_jitter <- input$barplot_add_jitter

## add point
add_point <- input$barplot_add_point

## x label angle
x_angle <- input$barplot_x_angle
h_just <- input$barplot_hjust
v_just <- input$barplot_vjust

## confidence level
conf.int <- input$barplot_ci

## only produce the plots that are shown in shiny
## different from shiny, which produces all the plots at one go
shown_id <- g_data$barplot_shown_id

## unit
unit <- input$unit

## ggplot theme
ggtheme <- g_data$ggtheme

font_family <- input$font_family


d <- g_data$d1[, shown_id, drop = F]
g <- g_data$group
dh <- g_data$hypothesis_test[shown_id, , drop = F]


posthoc <- NULL

if (add_sig){
    posthoc <- input$barplot_posthoc
}

if (input$pair_var == 'None') {
    pair <- NA
} else {
    pair <- as.factor(g_data$pair)
}

p_list <- batch_barplot(d, g, dh, pair,
                        posthoc = posthoc, unit = unit, 
                        add_point = add_point, add_sig = add_sig, 
                        add_jitter = add_jitter,
                        x_angle, h_just, v_just, 
                        as.numeric(conf.int),
                        ggtheme = ggtheme, font_family = font_family)


ncols <- input$barplot_ncol

nrows <- ceiling(length(p_list) / ncols)


if (ncols == 1) {
    p_list <- p_list[shown_id]
    if (length(p_list) > 0)
        l_ply(p_list, print)
} else {
    grid.newpage()
    tmp <- grid.draw(gridExtra::marrangeGrob(p_list, nrow=nrows, ncol=ncols, top = NULL))
}

##-------------------
## @knitr 125-barplot_text

out <- g_data$barplot_text

cat(out, '\n\n')

##-----------------------
## @knitr 130-sessionInfo

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
#------------------------------------
# Games-Howell Post Hoc Test in R
# https://stats.stackexchange.com/questions/83941/games-howell-post-hoc-test-in-r
#-------------------------------------

## obselete
posthoc.tgh <- function(y, x, method=c("games-howell", "tukey"), digits=2) {
    ### Based on http://www.psych.yorku.ca/cribbie/6130/games_howell.R
    method <- tolower(method);
    tryCatch(method <- match.arg(method), error=function(err) {
        stop("Argument for 'method' not valid!");
    });
    
    res <- list(input = list(x=x, y=y, method=method, digits=digits));
    
    res$intermediate <- list(x = factor(x[complete.cases(x,y)]),
                             y = y[complete.cases(x,y)]);
    res$intermediate$n <- tapply(y, x, length);
    res$intermediate$groups <- length(res$intermediate$n);
    res$intermediate$df <- sum(res$intermediate$n) - res$intermediate$groups;
    res$intermediate$means <- tapply(y, x, mean);
    res$intermediate$variances <- tapply(y, x, var);
    
    res$intermediate$pairNames <- combn(levels(res$intermediate$x),
                                        2, paste0, collapse=":");
    
    res$intermediate$descriptives <- cbind(res$intermediate$n,
                                           res$intermediate$means,
                                           res$intermediate$variances);
    rownames(res$intermediate$descriptives) <- levels(res$intermediate$x);
    colnames(res$intermediate$descriptives) <- c('n', 'means', 'variances');
    
    ### Start on Tukey
    res$intermediate$errorVariance <-
        sum((res$intermediate$n-1) * res$intermediate$variances) /
        res$intermediate$df;
    res$intermediate$t <- combn(res$intermediate$groups, 2, function(ij) {
        abs(diff(res$intermediate$means[ij]))/
            sqrt(res$intermediate$errorVariance*sum(1/res$intermediate$n[ij]));
    } );
    res$intermediate$p.tukey <- ptukey(res$intermediate$t*sqrt(2),
                                       res$intermediate$groups,
                                       res$intermediate$df,
                                       lower.tail=FALSE);
    res$output <- list();
    res$output$tukey <- cbind(res$intermediate$t,
                              res$intermediate$df,
                              res$intermediate$p.tukey)                                     
    rownames(res$output$tukey) <- res$intermediate$pairNames;
    colnames(res$output$tukey) <- c('t', 'df', 'p');
    
    ### Start on Games-Howell
    res$intermediate$df.corrected <- combn(res$intermediate$groups, 2, function(ij) {               
        sum(res$intermediate$variances[ij] /
                res$intermediate$n[ij])^2 / 
            sum((res$intermediate$variances[ij] /
                     res$intermediate$n[ij])^2 / 
                    (res$intermediate$n[ij]-1));
    } );
    res$intermediate$t.corrected <- combn(res$intermediate$groups, 2, function(ij) {               
        abs(diff(res$intermediate$means[ij]))/
            sqrt(sum(res$intermediate$variances[ij] /
                         res$intermediate$n[ij]));
    } );    
    res$intermediate$p.gameshowell <- ptukey(res$intermediate$t.corrected*sqrt(2),
                                             res$intermediate$groups,
                                             res$intermediate$df.corrected,
                                             lower.tail=FALSE)  
    res$output$games.howell <- cbind(res$intermediate$t.corrected,
                                     res$intermediate$df.corrected,
                                     res$intermediate$p.gameshowell);
    rownames(res$output$games.howell) <- res$intermediate$pairNames;
    colnames(res$output$games.howell) <- c('t', 'df', 'p');
    
    ### Set class and return object
    class(res) <- 'posthocTukeyGamesHowell';
    return(res);
    
}

print.posthocTukeyGamesHowell <- function(x, digits=x$input$digits, ...) {
    print(x$intermediate$descriptives, digits=digits);
    cat('\n');
    if (x$input$method == 'tukey') {
        print(x$output$tukey);
    }
    else if (x$input$method == 'games-howell') {
        print(x$output$games.howell, digits=digits);
    }
}

##--------------------------------------------------
## Hypothesis test
##--------------------------------------------------
hypothesisTest <- function(d1, g, pair = NA, method = union('none', p.adjust.methods)) {
    ## d1: data
    ## g: group
    method <- match.arg(method)
    
    if (!(length(pair) == 1 && is.na(pair))) {
        pair <- as.factor(pair)
        testthat::expect_equal(
            length(g), length(pair)
        )
        testthat::expect(
            all(table(g, pair) == 1),
            "Each pair must be in two groups."
        )
    }
    if (!is.data.frame(d1)) d1 <- as.data.frame(d1)
    if (!is.factor(g)) g <- as.factor(g)
    
    testthat::expect_equal(
        nrow(d1), length(g)
    )
    
    if (nlevels(g) == 1) {
        message('There is only 1 group. Return NULL.')
        return(NULL)
    }
    if (is.na(pair)) {
        ## not paired design
        test_result <- sapply(d1, function(x) {
            # g_mean <- aggregate(x, by = list(Group = g), mean)
            
            ## combination of group pairs
            g_pair <- combn(levels(g), 2)
            
            g_fold <- apply(g_pair, 2, function(pair.i) {
                mean(x[g == pair.i[2]], na.rm = T) /
                    mean(x[g == pair.i[1]], na.rm = T)
            })
            
            names(g_fold) <- apply(g_pair, 2, function(x) {
                paste0('Fold: ', paste0(rev(x), collapse = '/'))
            })
            
            ## welch's ANOVA
            ## var.equal = FALSE is the default, just to make it visible
            aov.out <- tryCatch(
                oneway.test(x ~ g, var.equal = FALSE),
                error = function(e) {
                    ## Number of replicate is 1
                    print(e)
                    list(p.value = NA)
                }
            )
            
            ## non-parametric
            ## wilcox when #groups == 2
            ## kruskal-wallis when #groups > 2
            if (nlevels(g) > 2)
                nonparam.out <- kruskal.test(x ~ g)
            else
                nonparam.out <- wilcox.test(x ~ g)
            
            if (nlevels(g) > 2) {
                ## posthoc.tgh in functions.R
                # posthoc.out1 <- posthoc.tgh(x, g)
                ## do not adjust at post-hoc test level
                posthoc.out <- rosetta::posthocTGH(
                    x, g, 
                    p.adjust = 'none', 
                    formatPvalue = F
                )
                
                ## replace - by :
                posthoc.names <- apply(
                    combn(levels(g),2), 2, 
                    function(x) {paste0(rev(x), collapse=':')})
                
                posthoc.out$output$tukey <- as.matrix(posthoc.out$output$tukey)
                rownames(posthoc.out$output$tukey) <- posthoc.names
                posthoc.out$output$games.howell <- as.matrix(posthoc.out$output$games.howell)
                rownames(posthoc.out$output$games.howell) <- posthoc.names
                
                ## dunn's test
                ## p value adjust
                ## do not adjust at post-hoc test level
                dunn.out <- FSA::dunnTest(x ~ g, method = 'none')
                dunn.out$res[, 'Comparison'] <- stringr::str_replace(
                    dunn.out$res[, 'Comparison'],
                    ' - ', ':'
                ) ## make it same as Tukey
                ## dunnTest does not honor level order
                dunn.out$res <- dunn.out$res[
                    match(
                        combn(levels(g), 2, function(x) { paste0(sort(x), collapse=":") }),
                        dunn.out$res[, 'Comparison']
                    ),
                    ]
                ## reverse the order
                paste0Rev <- function(x) {paste0(rev(x), collapse = ':')}
                dunn.out$res$Comparison <- combn(levels(g), 2, paste0Rev)
                
                res <- c(`parametric pvalue` = aov.out$p.value,
                         posthoc.out$output$tukey[, 'p'],
                         posthoc.out$output$games.howell[, 'p'],
                         `non-parametric pvalue` = nonparam.out$p.value,
                         dunn.out$res[, 'P.adj'],
                         g_fold)
                
                ## Need to set names because column of data.frame
                ## is not named
                # names(res)[2:(1+nrow(posthoc.out$output$tukey))] <-
                #     stringr::rownames(posthoc.out$output$tukey)
                
                ## number of pair-wise comparison
                n_comb <- nrow(posthoc.out$output$tukey)
                names(res)[2:(1 + 2 * n_comb)] <- paste0(
                    rep(c('TukeyHSD: ', 'Games-Howell: '), each = n_comb),
                    names(res)[2:(1 + 2 * n_comb)]
                )
                
                names(res)[(3 + 2 * n_comb):(2 + 3 * n_comb)] <- paste0(
                    rep('Dunn: ', n_comb),
                    dunn.out$res$Comparison
                )
                
            } else {
                res <- c(`parametric pvalue` = aov.out$p.value, 
                         `non-parametric pvalue` = nonparam.out$p.value,
                         g_fold)
            }
            
            res 
        }) %>% t
    } else {
        ## paired design
        ## paired t-test
        test_result <- sapply(d1, function(x) {
            order_pair <- order(pair)
            x_order <- x[order_pair]
            g_order <- g[order_pair]
            g1 <- which(g_order == levels(g)[1])
            g2 <- which(g_order == levels(g)[2])
            
            message('Fold change is the median of all pair-wise fold changes for paired design.')
            g_fold <- median(x_order[g2] / x_order[g1], na.rm = T)
            
            paired.t.out <- tryCatch(
                t.test(x_order[g1], x_order[g2], paired = T),
                error = function(e) {
                    print(e)
                    list(p.value = NA)
                }
            )
            
            paired.wilcox.out <- wilcox.test(x_order[g1], x_order[g2], paired = T)
            c(paired.t.out$p.value, 
              paired.wilcox.out$p.value,
              g_fold)
        }) %>% t
        
        colnames(test_result) <- c('parametric pvalue', 'non-parametric pvalue', paste0('Fold: ', levels(g)[2], '/', levels(g)[1]))
        
    }
    
    if (method != 'none') {
        test_result_adjusted <- as.data.frame(test_result) %>%
            tibble::rownames_to_column('rn') %>%
            dplyr::mutate(
                dplyr::across(!starts_with('Fold') & ! starts_with('rn'),
                              ~ p.adjust(.x, method = method)
                              )) %>%
            tibble::column_to_rownames('rn') %>%
            as.matrix
    } else {
        test_result_adjusted <- NULL
    }
    return(list(unadjusted = test_result, adjusted = test_result_adjusted))
}

##--------------------------------------------------
## Volcano plot
##--------------------------------------------------
plot_volcano <- function(d, sig_lvl = 0.05, fold_cutoff = 1.5, 
                         max_label = 10, title = '', 
                         ggtheme = theme_bw(), font_family = 'Arial',
                         force = 1, expand = 0.05,
                         min.segment.length = 0.5) {
    ## d is data.frame with 3 columns: pvalue, fold, label
    ## sig_lvl is the significance level
    ## fold_cutoff is fold cutoff
    assertthat::assert_that(ncol(d) == 3)
    
    colnames(d) <- c('pvalue', 'fold', 'label')
    
    d <- d %>% dplyr::mutate(
        # pvalue and fold change both must pass threshold
        significant = ifelse(
            pvalue >= sig_lvl, 'Not significant',
            ifelse(
                fold > fold_cutoff, 
                paste0('P<', sig_lvl, '&FC>', fold_cutoff),
                ifelse(
                    fold < 1/fold_cutoff,
                    paste0('P<', sig_lvl, '&FC<1/', fold_cutoff),
                    'Not significant'
                ))),
        label = ifelse(significant != 'Not significant', label, '')
    ) %>% dplyr::filter(
        !is.na(pvalue) & !is.na(fold) & is.finite(fold) & fold > 0
    )
    
    
    if (max_label > 0) {
        ## number of significant variables
        n_sig <- length(which(d$label != ''))
        if (n_sig > max_label) {
            ## ids of top max_label smallest P value
            ## order by P value followed by fold change
            id_label <- intersect(
                # order(rank(d[, 'pvalue']) + rank(-abs(log2(d[, 'fold'])))),
                order(d[, 'pvalue']),
                which(d$significant != 'Not significant')
            )[1:max_label]
            d[-id_label, 'label'] <- ''
        }
    } else if(max_label == 0) {
        d[, 'label'] <- ''
    }
    
    if (fold_cutoff == 1) {
        ## change significant to P > 0.05 and P < 0.05
        ## change color_values accordingly
        d <- d %>% dplyr::mutate(
            significant = ifelse(
                significant == 'Not significant',
                paste0('P\u2265', sig_lvl),
                paste0('P<', sig_lvl)
            )
        )
        color_values <- c('gray40', 'red') %>%
            `names<-`(c(paste0('P\u2265', sig_lvl),
                        paste0('P<', sig_lvl)))
    } else {
        color_values <- c('gray40', 'red', 'green') %>%
            `names<-`(c('Not significant', 
                        paste0('P<', sig_lvl, '&FC>', fold_cutoff),
                        paste0('P<', sig_lvl, '&FC<1/', fold_cutoff)))
    }
    
    p <- ggplot(d, aes(x = log2(fold), y = -log10(pvalue))) +
        geom_point(aes(color = significant)) +
        scale_color_manual(values = color_values) +
        ggrepel::geom_text_repel(aes(label = label), family = font_family,
                                 force = force, max.iter = 2000) +
        geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.5) +
        geom_hline(yintercept = -log10(sig_lvl), linetype = 'dashed', alpha = 0.5) +
        scale_x_continuous(
            limits = c(-max(abs(log2(d$fold))), max(abs(log2(d$fold)))),
            expand = expansion(mult = c(expand, expand), add = c(0, 0))) +
        scale_y_continuous(expand = expansion(mult = c(0, expand), add = c(0, 0))) +
        guides(color = guide_legend(title = 'Significant')) +
        labs(title = title, x = 'log2(fold change)', y = '-log10(P value)') +
        ggtheme
    
    ## add vertical corresponds to fold change cut-off
    if (fold_cutoff > 1) {
        p <- p +
            geom_vline(xintercept = log2(fold_cutoff), linetype = 'dashed', alpha = 0.5) +
            geom_vline(xintercept = -log2(fold_cutoff), linetype = 'dashed', alpha = 0.5)
    }
    
    p
}

volcanoText <- function(d, sig_lvl = 0.05, fold_cutoff = 1.5, max_label = 10) {
    ## d is data.frame with 3 columns: pvalue, fold, label
    ## sig_lvl is the significance level
    ## fold_cutoff is fold cutoff
    assertthat::assert_that(ncol(d) == 3)
    
    colnames(d) <- c('pvalue', 'fold', 'label')
    
    d <- d %>% dplyr::mutate(
        # pvalue and fold change both must pass threshold
        significant = ifelse(
            pvalue >= sig_lvl, 'Not significant',
            ifelse(
                fold > fold_cutoff, 
                paste0('P<', sig_lvl, '&FC>', fold_cutoff),
                ifelse(
                    fold < 1/fold_cutoff,
                    paste0('P<', sig_lvl, '&FC<1/', fold_cutoff),
                    'Not significant'
                ))),
        label = ifelse(significant != 'Not significant', label, '')
    ) %>% dplyr::filter(
        !is.na(pvalue) & !is.na(fold) & is.finite(fold) & fold > 0
    )
    
    
    if (max_label > 0) {
        ## number of significant variables
        n_sig <- length(which(d$label != ''))
        if (n_sig > max_label) {
            ## ids of top max_label smallest P value
            ## order by P value followed by fold change
            id_label <- intersect(
                # order(rank(d[, 'pvalue']) + rank(-abs(log2(d[, 'fold'])))),
                order(d[, 'pvalue']),
                which(d$significant != 'Not significant')
            )[1:max_label]
            d[-id_label, 'label'] <- ''
        }
    } else if(max_label == 0) {
        d[, 'label'] <- ''
    }
    
    if (fold_cutoff == 1) {
        ## change significant to P > 0.05 and P < 0.05
        ## change color_values accordingly
        d <- d %>% dplyr::mutate(
            significant = ifelse(
                significant == 'Not significant',
                paste0('P\u2265', sig_lvl),
                paste0('P<', sig_lvl)
            )
        )
    }
    
    list(
        up = subset(d, label != '' & fold > 1)[, 'label'],
        down = subset(d, label != '' & fold < 1)[, 'label']
    )
}

plotDifferentialMetabolites <- function(d, sig_lvl = 0.05, fold_cutoff = 1.5, 
                                        max_label = 10, title = '',
                                        ggtheme = theme_bw(), font_family = 'Arial',
                                        force = 1, expand = 0.05,
                                        min.segment.length = 0.5) {
    ## d is data.frame with 3 columns: pvalue, fold, label
    ## sig_lvl is the significance level
    ## fold_cutoff is fold cutoff
    assertthat::assert_that(ncol(d) == 3)
    
    colnames(d) <- c('pvalue', 'fold', 'label')
    
    d <- d %>% dplyr::mutate(
        # pvalue and fold change both must pass threshold
        significant = ifelse(
            pvalue >= sig_lvl, 'Not significant',
            ifelse(
                fold > fold_cutoff,
                paste0('P<', sig_lvl, '&FC>', fold_cutoff),
                ifelse(
                    fold < 1/fold_cutoff,
                    paste0('P<', sig_lvl, '&FC<1/', fold_cutoff),
                    'Not significant'
                )))
    ) %>% dplyr::filter(
        !is.na(pvalue) & !is.na(fold) & is.finite(fold) & fold > 0
    )
    
    if (nrow(d) == 0) return(NULL)
    
    ## sort by P value
    d <- d[order(d$pvalue, decreasing = F), ]
    ## select top n
    d <- d[1:min(max_label, nrow(d)), ]
    d <- d %>% mutate(
        yval = nrow(d):1
    )
    
    ylim.min <- 0.5
    ylim.max <- nrow(d)
    expand.y <- 1
    
    
    ## labels in forest plot
    # labels, could be extended to show more information
    table_plot <-
        ggplot(d) +
        ggtheme +
        geom_text(aes(label = label, x = 0, y = yval - 0.2), hjust = 0,
                  family = font_family) +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank()
        ) +
        scale_y_continuous(expand = c(0, expand.y), limits = c(ylim.min, ylim.max)) +
        xlim(0, 3)
    
    
    # Forest plot
    forest1 <-
        ggplot(d, aes(x = yval, y = log2(fold))) +
        ggtheme +
        geom_bar(stat = 'identity', position = 'dodge') +
        geom_hline(yintercept = 0) +
        theme(
            axis.line.x = element_line(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.border = element_blank()
        )  +
        coord_flip() +
        scale_x_continuous(expand = c(0, expand.y), limits = c(ylim.min, ylim.max + 0.5)) +
        ylab("log2(fold change)")
    
    sig1 <- ggplot(d, aes(y = yval - 0.2)) +
        ggtheme +
        geom_text(aes(label = ifelse(
            pvalue > 0.001,
            round(pvalue, digits = 3),
            format(pvalue, digits = 3, scientific = T)), x = 0), hjust = 0) +
        # geom_text(label = 'sig', x = 0.1, y = nrow(d1) + 1.5) +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank()
        ) +
        xlim(0, 2) +
        scale_y_continuous(expand = c(0, expand.y), limits = c(ylim.min, ylim.max))
    
    grid.newpage()
    pushViewport(
        viewport(
            layout=grid.layout(2, 1,
                               widths=unit(c(1),
                                           'null'),
                               heights=unit(c(1.2, 15),
                                            "null"))))
    
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
    g <- gridExtra:::gtable_cbind(
        ggplotGrob(table_plot),
        ggplotGrob(forest1),
        ggplotGrob(sig1),
        size = "max")
    panels <- g$layout$r[grep("panel", g$layout$name)]
    w <- c(2.3, 2, 0.6)
    ws <- sum(w)
    g$widths[panels] <- unit(w, "null")
    grid.draw(g)
    
    popViewport()
    
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row=1))
    grid.lines(x = unit(c(0.01, 1), 'npc'),
               y = unit(-0.1, 'npc'))
    grid.lines(x = unit(c(0.01, 1), 'npc'),
               y = unit(0.8, 'npc'),
               gp = gpar(lwd = 2))
    grid.text(x = unit(0.05, 'npc'), just = 'left',
              y = unit(0.4, 'npc'),
              label = paste0('Top metabolites in ', title))
    grid.text(x = unit(0.65, 'npc'),
              y = unit(0.4, 'npc'),
              label = 'Fold change')
    grid.text(x = unit(0.9, 'npc'),
              y = unit(0.4, 'npc'),
              label = 'P-value')
    
}

differentialMetabolitesText <- function(d, sig_lvl = 0.05, fold_cutoff = 1.5, 
                                        max_label = 10) {
    ## d is data.frame with 3 columns: pvalue, fold, label
    ## sig_lvl is the significance level
    ## fold_cutoff is fold cutoff
    assertthat::assert_that(ncol(d) == 3)
    
    colnames(d) <- c('pvalue', 'fold', 'label')
    
    d <- d %>% dplyr::mutate(
        # pvalue and fold change both must pass threshold
        significant = ifelse(
            pvalue >= sig_lvl, 'Not significant',
            ifelse(
                fold > fold_cutoff,
                paste0('P<', sig_lvl, '&FC>', fold_cutoff),
                ifelse(
                    fold < 1/fold_cutoff,
                    paste0('P<', sig_lvl, '&FC<1/', fold_cutoff),
                    'Not significant'
                )))
    ) %>% dplyr::filter(
        !is.na(pvalue) & !is.na(fold) & is.finite(fold) & fold > 0
    )
    
    if (nrow(d) == 0) return(NULL)
    
    ## sort by P value
    d <- d[order(d$pvalue, decreasing = F), ]
    ## select top n
    d <- d[1:min(max_label, nrow(d)), ]
    d <- d %>% mutate(
        yval = nrow(d):1
    )
    
    list(
        all = subset(d, label != '')[, 'label'],
        up = subset(d, label != '' & fold > 1 & pvalue < sig_lvl)[, 'label'],
        down = subset(d, label != '' & fold < 1 & pvalue < sig_lvl)[, 'label']
    )
}


##--------------------------------------------------
## subchunkify
## set figure width and height dynamically within a chunk
## https://stackoverflow.com/questions/15365829/dynamic-height-and-width-for-knitr-plots
##--------------------------------------------------
# 
# subchunkify <- function(g, fig_height=7, fig_width=5) {
#     g_deparsed <- paste0(deparse(
#         function() {g}
#     ), collapse = '')
#     
#     sub_chunk <- paste0("
#                         `","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=",
#                         fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
#                         "\n(", 
#                         g_deparsed
#                         , ")()",
#                         "\n`","``
#                         ")
#     
#     cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
                        # }

# And use the function like this, defining your own figure sizes:
#     
# ```{r echo=FALSE, results='asis'}
# g <- ggplot(economics, aes(date, unemploy)) + 
#     geom_line()
# subchunkify(g, 10, 3)
# subchunkify(g, 7, 7)
# ```

# Or let the data define the sizes:
#     
#     ```{r echo=FALSE, results='asis'}
# g <- ggplot(economics, aes(date, unemploy)) + 
#     geom_line()
# for (i in seq(2, 5)) {
#     subchunkify(g, i / 2, i)
# }

# g <- ggplot(economics, aes(date, unemploy)) + 
#     geom_line()
# 
# cat('<h2>A Small Square Plot</h2>')
# subchunkify(g, 3, 3)


##--------------------------------------------------
## boxplot 
##--------------------------------------------------
batch_boxplot <- function(d, g, dh, posthoc, unit, pair = NA,
                          add_point = T, add_sig = T, add_jitter = T,
                          add_violin = F, notch = F, 
                          x_angle = 0, h_just = 0.5, v_just = 0, 
                          ggtheme = theme_bw(), font_family = 'Arial') {
    # d: data.frame
    # g: group, length is the same as row number of d
    # dh: hypothesis test resutls
    d.colnames <- colnames(d)
    d <- cbind(d, group = g)
    
    p_list <- lapply(d.colnames, function(var.i) {
        ## add violin
        if (add_violin) {
            p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                geom_violin() +
                geom_boxplot(width = 0.2, notch = notch)
            
            ## add violin and jitter
            if (add_point) {
                p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                    geom_violin() + 
                    geom_boxplot(width = 0.2, outlier.size = 0, notch = notch) +
                    geom_jitter(width = 0.1)
            }
       
        } else if (add_point) {
            
            if (add_jitter) {
                if (!is.na(pair)) {
                    p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                        geom_boxplot(outlier.size = 0, notch = notch) +
                        geom_jitter(aes(group = pair), position = position_dodge(0.2)) +
                        geom_line(aes(group = pair), color = 'black', alpha = 0.6, position = position_dodge(0.2))
                } else {
                    p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                        geom_boxplot(outlier.size = 0, notch = notch) +
                        geom_jitter(width = 0.1)
                }
            } else {
                ##no jitter
                if (!is.na(pair)) {
                    p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                        geom_boxplot(outlier.size = 0, notch = notch) +
                        geom_point() + 
                        geom_line(aes(group = pair), color = 'black', alpha = 0.6)
                } else {
                    p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                        geom_boxplot(outlier.size = 0, notch = notch) +
                        geom_point()
                }
            } 
            
        } else {
            ## no voilin no jitter
            p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                 geom_boxplot(notch = notch) 
        }
            
        if (add_sig) {
            if (nlevels(g) == 1) {
                
            } else if(nlevels(g) > 2) {
                
                p_matrix <- matrix(
                    1, nrow = nlevels(g), ncol = nlevels(g),
                    dimnames = list(levels(g), levels(g))    
                )
                ## multcompLetters work on lower.tri
                p_matrix[lower.tri(p_matrix)] <- dh[var.i, stringr::str_detect(colnames(dh), posthoc)]
                
                tryCatch(
                    {
                        sig_letters <- multcompView::multcompLetters(
                            p_matrix, threshold = 0.05
                        )$Letters
                        
                        d_sig <- data.frame(group = names(sig_letters), 
                                            y = 1.05 * max(d[, var.i], na.rm = T) - 0.05 * min(d[, var.i], na.rm = T),
                                            label = sig_letters)                    
                        
                        suppressWarnings(
                            p <- p + geom_text(
                                data = d_sig,
                                aes(x = group, y = y, label = label), family = font_family
                            ) + scale_y_continuous(
                                expand = c(0.1, 0)
                            )
                        )
                    },
                    error = function(e) {
                        print(paste(var.i, ':', e$message))
                    }
                )
            } else {
                p_col <- ifelse(posthoc == 'Parametric',
                                'parametric pvalue', 'non-parametric pvalue')
                d_sig <- data.frame(
                    start = levels(g)[1],
                    end = levels(g)[2],
                    y = 1.1 * max(d[, var.i], na.rm = T) - 0.1 * min(d[, var.i], na.rm = T),
                    label = formatC(dh[var.i, p_col], digits = 2)
                )
                
                suppressWarnings(
                    p <- p + ggsignif::geom_signif(
                        data = d_sig,
                        aes(xmin = start, xmax = end,
                            annotations = label, y_position = y),
                        manual = T,
                        tip_length = min(
                            0.01 * (max(d[, var.i], na.rm = T) - min(d[, var.i], na.rm = T)),
                            0.01),
                        family = font_family
                    ) + scale_y_continuous(expand = c(0.15, 0))
                )
            }
            
        }
        
        ## y-axis label
        ylab <- ifelse(unit == ' ',
                       'Concentration',
                       parse(text = paste0('Concentration (', unit, ')'))
        )
        
        p <- p + labs(x = '', title = var.i, 
                      y = ylab) +
            ggtheme +
            theme(axis.text.x = element_text(angle = x_angle, vjust = v_just, hjust = h_just))
        
        p
    })
    return(p_list)
}

##--------------------------------------------------
## barplot 
##--------------------------------------------------

batch_barplot <- function(d, g, dh, 
                          posthoc, unit, pair = NA,
                          add_point = T, add_jitter = T, add_sig = T, 
                          x_angle = 0, h_just = 0.5, v_just = 0, 
                          conf.int = 0.95, 
                          ggtheme = theme_bw(), font_family = 'Arial') {
    # d: data.frame
    # g: group, length is the same as row number of d
    # dh: hypothesis test resutls
    # conf.int: confidence interval
    d.colnames <- colnames(d)
    d <- cbind(d, group = g)
    mult <- qnorm((1 + conf.int)/2)
    if (abs(conf.int - 0.68) < 0.01)
        mult = 1
    
    p_list <- lapply(d.colnames, function(var.i) {
        
        p <- ggplot(
            d,
            aes(x = group,
                y = eval(parse(
                    text = paste0('`', var.i, '`')
                )))) +
            stat_summary(
                fun = mean,
                geom = 'bar',
                fill = 'gray80',
                color = 'black'
            ) +
            stat_summary(
                fun.data = mean_se,
                fun.args = list(mult = mult),
                geom = 'errorbar',
                width = 0.2
            )
        
        if (add_point) {
            if (add_jitter) {
                if (!is.na(pair)) {
                    p <- p +
                         geom_jitter(aes(group = pair), position = position_dodge(0.2)) +
                         geom_line(aes(group = pair), color = 'black', alpha = 0.6, position = position_dodge(0.2))
                } else {
                    p <- p +
                         geom_jitter(width = 0.1)
                }
            } else {
                ##no jitter
                if (!is.na(pair)) {
                    p <- p +
                         geom_point() +
                         geom_line(aes(group = pair), color = 'black', alpha = 0.6)
                } else {
                    p <- p + geom_point(aes(group = g))
                }
            }
        } else {
            ## no point
            p <- p
        }
        
        if (add_sig) {
            if (nlevels(g) == 1) {
                
            } else if(nlevels(g) > 2) {
                p_matrix <- matrix(
                    1, nrow = nlevels(g), ncol = nlevels(g),
                    dimnames = list(levels(g), levels(g))    
                )
                ## multcompLetters work on lower.tri
                p_matrix[lower.tri(p_matrix)] <- dh[var.i, stringr::str_detect(colnames(dh), posthoc)]
                
                tryCatch(
                    {
                        sig_letters <- multcompView::multcompLetters(
                            p_matrix, threshold = 0.05
                        )$Letters
                        ## compute y limits
                        ymax <- vaggregate(d[, var.i], g, mean_se) %>%
                            unlist %>% max
                        d_sig <- data.frame(group = names(sig_letters), 
                                            # y = ymax * 1.05,
                                            # y = 1.05 * max(d[, var.i]) - 0.05 * min(d[, var.i]),
                                            # y = 1.1 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                                            y = 1.05 * max(d[, var.i], na.rm = T) - 0.05 * min(d[, var.i], na.rm = T),
                                            label = sig_letters)                    
                        
                        suppressWarnings(
                            p <- p + geom_text(
                                data = d_sig,
                                aes(x = group, y = y, label = label),
                                family = font_family
                            ) + scale_y_continuous(
                                expand = c(0.1, 0)
                            )
                        )
                    },
                    error = function(e) {
                        print(paste(var.i, ':', e$message))
                    }
                )
                
            } else {
                ## 2 groups
                p_col <- ifelse(posthoc == 'Parametric',
                                'parametric pvalue', 'non-parametric pvalue')
                d_sig <- data.frame(
                    start = levels(g)[1],
                    end = levels(g)[2],
                    # y = 1.1 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                    y = 1.05 * max(d[, var.i], na.rm = T) - 0.05 * min(d[, var.i], na.rm = T),
                    label = formatC(dh[var.i, p_col], digits = 2)
                )
                
                suppressWarnings(
                    p <- p + ggsignif::geom_signif(
                        data = d_sig,
                        aes(xmin = start, xmax = end, y_position = y,
                            annotations = label),
                        manual = T,
                        tip_length = min(
                            0.01 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                            0.01),
                        family = font_family
                    ) + scale_y_continuous(expand = c(0.15, 0))
                )
            }
            
        }
        
        ## y-axis label
        ylab <- ifelse(unit == ' ',
                       'Concentration',
                       parse(text = paste0('Concentration (', unit, ')'))
        )
        
        p <- p + labs(x = '', title = var.i, 
                      y = ylab) +
            ggtheme +
            theme(axis.text.x = element_text(angle = x_angle, vjust = v_just, hjust = h_just))
        
        p
    })

    return(p_list)
}

##--------------------------------------------------
## grid_fam
## change font family of grid plot 
## 20190530
##--------------------------------------------------
grid_fam <- function(p, font_family="Arial") 
{
    print('grid_fam')
    if (class(p)[1] == 'pheatmap') {
        ## return gtable if p is of clas pheatmap
        ## colname and rowname
        # print('1')
        g <- p$gtable
        # print('2')
        px <- which(g$layout$name %in% c('col_names', 'row_names'))
        # print('3')
        for (i in px) {
            g$grobs[[i]]$gp$fontfamily <- font_family
            g$grobs[[i]]$gp$fontface <- 1L
            g$grobs[[i]]$gp$font <- 1L
        }
        # print('4')
        ## legend
        px <- which(g$layout$name == 'legend')
        # print('5')
        id <- grep('text', names(g$grobs[[px]]$children))
        # print('6')
        for (i in id) {
            g$grobs[[px]]$children[[i]]$gp$fontfamily <- font_family
            g$grobs[[px]]$children[[i]]$gp$fontface <- 'plain'
            g$grobs[[px]]$children[[i]]$gp$font <- structure(1L, names = 'plain')
        }
        # print('7')
        return(invisible(g))
    }
    
    ## https://stackoverflow.com/questions/17012518/why-does-this-r-ggplot2-code-bring-up-a-blank-display-device
    ## pdf(NULL) prevents ggplotGrob from generating Rplots.pdf
    # pdf(NULL)
    # suppressWarnings(g <- ggplotGrob(p))
    # dev.off()
    g <- ggplotGrob(p)
    # dev.off()
    
    px <- which(g$layout$name=="panel")
    
    id <- grep("text", names(g$grobs[[px]]$children))
    
    for(i in id)  {
        n <- length(g$grobs[[px]]$children[[i]]$gp$fontfamily)
        g$grobs[[px]]$children[[i]]$gp$fontfamily <- rep(font_family, n)
        g$grobs[[px]]$children[[i]]$gp$fontface <- rep('plain', n)
        g$grobs[[px]]$children[[i]]$gp$font <- structure(rep(1L, n), 
                                                         names = rep('plain', n))
    }
    # grid::grid.newpage()
    # grid::grid.draw(g)
    print('end grid_fam')
    invisible(g)
}


##--------------------------------------------------
## OPLS-DA splot
## 20190606
## Calculate p1 and pcorr1 for splot
##--------------------------------------------------
splotCal <- function(opls.obj, x) {
    ## opls.obj: ropls::opls object
    ## x: normalized data passed to opls
    ## normalized data
    # s <- as.matrix(mSetObj$dataSet$norm)
    s <- as.matrix(x)
    ## p1
    # T <- as.matrix(mSetObj$analSet$oplsda$scoreMN)
    T <- as.matrix(opls.obj@scoreMN)
    ## p1 is pearson correlation between p1 and variable
    p1 <- c()
    for (i in 1:ncol(s)) {
        scov <- cov(s[, i], T)
        p1 <- matrix(c(p1, scov), ncol = 1)
    }
    ## pcorr1 is pearson correlation divided by standard deviation of each variable
    pcorr1 <- c()
    for (i in 1:nrow(p1)) {
        den <- apply(T, 2, sd) * sd(s[, i])
        corr1 <- p1[i, ]/den
        pcorr1 <- matrix(c(pcorr1, corr1), ncol = 1)
    }
    data.frame(p1, pcorr1, vip = attr(opls.obj, 'vipVn'))
}

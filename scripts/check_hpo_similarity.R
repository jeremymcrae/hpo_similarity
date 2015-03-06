# script to check the distribution of P values from HPO similarity testing

library(qqman)
library(Cairo)


HPO_PATH = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/results/hpo_similarity.max.txt"
ENRICHMENT_PATH = "/nfs/users/nfs_j/jm33/apps/mupit/results/de_novos.ddd_4k.with_diagnosed.ddd_only.enrichment_results.txt"
DDG2P_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt"

open_hpo_similarity_results <- function(path) {
    hpo = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # adjust any P values at 0 to the next lowest P value, so that these points
    # can be included in a QQ plot.
    p_values = unique(hpo$hpo_similarity_p_value)
    p_values = p_values[p_values!= 0]
    min_p = min(p_values)
    
    hpo$hpo_similarity_p_value[hpo$hpo_similarity_p_value == 0] = min_p
    
    return(hpo)
}

open_enrichment_results <- function(path) {
    # open the P values from enrichment testing, so we can check that the P values
    # from de novo enrichment are independent of the P values from testing for HPO
    # similarity
    enrichment = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    enrichment = enrichment[, c("hgnc", "p_lof", "p_func")]
    
    return(enrichment)
}

open_ddg2p_genes <- function(path) {
    # open the DDG2P, so we can check for enrichment in the more significant genes
    ddg2p = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    ddg2p = ddg2p[!grepl("Possible DD Gene", ddg2p$DDD_category), ]
    
    return(unique(ddg2p$ddg2p_gene_name))
}

check_correlation_with_enrichment <- function(hpo) {
    # test for correlation between the HPO similarity P value and the enrichment P values
    to_lof = cor(log10(hpo$hpo_similarity_p_value), log10(hpo$p_lof), use="pairwise.complete.obs")
    to_func = cor(log10(hpo$hpo_similarity_p_value), log10(hpo$p_func), use="pairwise.complete.obs")
    
    to_lof = signif(to_lof, 3)
    to_func = signif(to_func, 3)
    
    cat(paste("HPO similarity vs lof enrichment correlation", to_lof, "\n"))
    cat(paste("HPO similarity vs func enrichment correlation", to_func, "\n"))
}

plot_qq <- function(hpo, path) {
    
    path = paste(basename(path), ".pdf", sep="")
    
    # plot a QQ plot, and  HPO similarity P values vs enrichment P values
    Cairo(file=path, type="pdf", height=15, width=15, units="cm")
    qq(hpo$hpo_similarity_p_value, main="QQ plot of HPO similarity tests",
        xlim=c(0, 6), ylim=c(0, 6), las=1)
    
    plot(-log10(hpo$hpo_similarity_p_value), -log10(hpo$p_lof),
        main="HPO similarity vs LoF enrichment",
        xlab="-log10 P HPO similarity", ylab="-log10 P LoF enrichment", las=1)
    plot(-log10(hpo$hpo_similarity_p_value), -log10(hpo$p_func),
        main="HPO similarity vs Func enrichment",
        xlab="-log10 P HPO similarity", ylab="-log10 P Func enrichment", las=1)
    dev.off()
}

test_ddg2p_enrichment <- function(hpo, ddg2p_genes) {
    
    hpo$in_ddg2p = hpo$hgnc %in% ddg2p_genes
    
    low_p_genes = hpo[hpo$hpo_similarity_p_value < 0.05, ]
    high_p_genes = hpo[hpo$hpo_similarity_p_value >= 0.05, ]
    
    # count the number of DDG2P genes in the two subsets
    low_p_true = sum(low_p_genes$in_ddg2p)
    low_p_false = sum(!low_p_genes$in_ddg2p)
    high_p_true = sum(high_p_genes$in_ddg2p)
    high_p_false = sum(!high_p_genes$in_ddg2p)
    
    # construct a matrix for the fisher's exact test
    ddg2p_counts =
        matrix(c(low_p_true, high_p_true, low_p_false, high_p_false),
            nrow = 2,
            dimnames = list(c("low_p", "high_p"),
                            c("ddg2p", "non_ddg2p")))
    
    # test for enrichment of DDG2P genes in the more significant genes
    fisher_p = fisher.test(ddg2p_counts)$p.value
    
    # find the DDG2P enrichment in the most significant genes
    ratio = (low_p_true/sum(low_p_true + low_p_false))/(sum(hpo$in_ddg2p)/nrow(hpo))
    
    fisher_p = signif(fisher_p, 3)
    ratio = signif(ratio, 3)
    
    cat(paste("The more significant genes (HPO similarity P<0.05) show a ", ratio,
        "-fold enrichment of DDG2P genes compared to all genes (P=", fisher_p,")\n", sep=""))
}


main <- function() {
    hpo = open_hpo_similarity_results(HPO_PATH)
    enrichment = open_enrichment_results(ENRICHMENT_PATH)
    ddg2p_genes = open_ddg2p_genes(DDG2P_PATH)
    hpo = merge(hpo, enrichment, by="hgnc")
    
    check_correlation_with_enrichment(hpo)
    plot_qq(hpo, HPO_PATH)
    test_ddg2p_enrichment(hpo, ddg2p_genes)
}

main()

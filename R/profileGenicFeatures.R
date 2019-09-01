#' Calculate metadata profiles across genic features
#' 
#' This function calculates metadata profiles across genic features using a procedure similar to that used in \href{https://www.sciencedirect.com/science/article/pii/S1097276519303533?via%3Dihub#fig1}{Viphakone et al., 2019}.
#' @param genicRegions A named list of \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} objects.
#' @param sampleObject A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} or \code{\link[GenomicAlignments:GAlignments-class]{GenomicAlignments}} object.
#' @param tx2gene A \code{"\linkS4class{data.frame}"} object. The first column must be of transcript identifiers, while the second must be of gene identifiers. Additional columns will be discarded.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. Required if \code{genicRegions} or \code{tx2gene} are not provided. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param bins An ordered integer vector (must be greater equal that the length of \code{genicRegions}). The vector order determines the number of bins for each region. If more bins than regions are provided, the additional will be ignored.
#' @param weightCol A single character string. This must be the name of an integer column in the \code{sampleObject} object.
#' @param ignoreStrand When set to 'TRUE', the strand information in \code{sampleObject} is ignored. This does not affect the features in \code{genicRegions}.
#' @param dropEmpty When set to 'TRUE', the transcripts with no signal in any of their sub-regions will be discarded. When set to 'FALSE', all values of these regions will be set to 0.
#' @param collapseBy A character string indicating which summarising method to use. One of "region" (default), "gene", "transcript", "recursive": can be abbreviated.
#'      Profiles are always calculated per transcript (tx_id), and reported as such when this parameter is set to 'transcript'.
#'      When set to 'gene', profiles are collapsed (using sum, mean and sd) by gene_id, while with 'region' they are collapsed in the same way but by region_id.
#'      When set to 'recursive', profiles are first collapse by gene_id and then by region_id, and sum, mean (pooled) and sd (pooled) are reported. 
#'      Please, be aware that each method will result in a different output, as the number and names of the columns by the chosen method.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
#' @return A \code{\link[data.table:data.table-class]{data.table}} of the normalised binned coverage across the genic features. Names and number of columns are determined by \code{collapseBy}.
#'
#' @import IRanges
#' @import BiocGenerics
#' @importFrom S4Vectors mcols
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom data.table ":=" ".N" ".GRP" ".I" data.table as.data.table setkey setkeyv setnames rbindlist
#' @importFrom stats setNames
#' @export
#' @examples
#' library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
#' TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
#' 
#' query = extractGenicFeatures(TxDb)
#' 
#' library("pasillaBamSubset")
#' library("GenomicAlignments")
#' 
#' fl1 <- untreated1_chr4()
#' subject = readGAlignments(fl1)
#'
#' profile = profileGenicFeatures(genicRegions=query, sampleObject=subject, TxDb=TxDb)
#'
#' library("ggplot2")
#'
#' ggplot(profile, aes(x=bin, y=Mean, col=region_id)) + 
#'    geom_ribbon(aes(ymin=Mean-Sd, ymax=Mean+Sd, fill=region_id), alpha=0.5) +
#'    geom_line() +
#'    geom_vline(xintercept=c(20.5, 120.5), linetype="dashed") +
#'    scale_x_continuous("Relative position") +
#'    scale_y_continuous("Average normalised signal")

profileGenicFeatures <- function(genicRegions=NULL, sampleObject=NULL, bins=c(20, 100, 70, 100), TxDb=NULL, tx2gene=NULL, weightCol=NULL, ignoreStrand=FALSE, dropEmpty=TRUE, collapseBy=c("region", "gene", "transcript", "recursive"), verbose=TRUE){

        stopifnot( !is.null(genicRegions) | !is.null(TxDb) )
        stopifnot( !is.null(genicRegions) & length(genicRegions)<=length(bins) ) # Early call
        stopifnot( !is.null(tx2gene) | !is.null(TxDb) )
        stopifnot( !is.null(weightCol)  & weightCol%in%colnames(mcols(sampleObject)) )
        stopifnot( !any(is.na(seqlengths(sampleObject))) | !is.null(TxDb) )

        collapseBy = match.arg(collapseBy)

        if( is.null(tx2gene) ){
            tx2gene = unlist(transcriptsBy(TxDb, by="gene"), use.names=TRUE)
            tx2gene = data.table(tx_id = tx2gene$tx_name, gene_id = names(tx2gene))
        }
        tx2gene = data.table(tx2gene[, 1:2])
        setnames(tx2gene, c("tx_id", "gene_id"))

        if( is.null(genicRegions) ){
            genicRegions = extractGenicFeatures(TxDb, tx2gene=tx2gene, bins=bins)
        }
        stopifnot( !is.null(genicRegions) & length(genicRegions)<=length(bins) ) # Late call

        GenomeInfoDb::seqlevelsStyle(sampleObject) = "UCSC"
        sampleObject = sortSeqlevels(sampleObject)
        if( any(is.na(seqlengths(sampleObject))) | any(!(seqlevels(GenomeInfoDb::keepStandardChromosomes(TxDb)) %in% seqlevels(sampleObject))) ){
            sampleObject = GenomeInfoDb::keepStandardChromosomes(sampleObject, pruning.mode="coarse")
            if(  length(seqlevels(sampleObject)) < length(seqlevels(GenomeInfoDb::keepStandardChromosomes(TxDb))) )
                seqlevels(sampleObject) = seqlevels(GenomeInfoDb::keepStandardChromosomes(TxDb))
            seqinfo(sampleObject) = seqinfo(GenomeInfoDb::keepStandardChromosomes(TxDb))
        }

        profiles = lapply(seq_along(genicRegions),
            function(i){
                region_name = names(genicRegions[i])
                region = genicRegions[[i]]
                nbin = bins[i]

                tmp = subsetByOverlaps(sampleObject, region, ignore.strand=ignoreStrand)

                weightVec = which(colnames(mcols(tmp))==weightCol)
                if( length(weightVec) > 0 )
                    tmp = rep(tmp, mcols(tmp)[, weightVec])

                tmp_list = lapply(setNames(c("+", "-"), c("Plus", "Minus")),
                    function(std){
                        if( any(strand(tmp)==std) & any(any(strand(region)==std)) ){
                            tmp_std = coverage(tmp[strand(tmp)==std,])

                            # Sort regions (crucial for '-' strand)
                            region_std = sort(unlist(region[strand(region)==std,]))
                            region_std = GRangesList(split(region_std, names(region_std)))
                            if( verbose )
                                message(paste(region_name, "regions on", std, "strand:", length(region_std)))
                            region_std = unlist(region_std, use.names=FALSE)
                            region_std$tx_id = names(region_std)

                            region_dt = data.table(tx_id = region_std$tx_id)
                            region_dt[, exon_index := .I]
                            region_dt[, tx_index := .GRP, by="tx_id"]
                            setkey(region_dt, tx_id)
                            setkey(tx2gene, tx_id)
                            region_dt[tx2gene, gene_id := gene_id]
                            region_dt[, gene_index := .GRP, by="gene_id"]
                            setkey(region_dt, exon_index)

                            tmp_std = as.data.table(tmp_std[region_std], key="group")
                            if( std=="+" ){
                                tmp_std[, pos := 1:.N, by="group"]
                            }else{
                                tmp_std[, pos := .N:1, by="group"]
                            }
                            setkey(tmp_std, group)

                            tmp_std = tmp_std[region_dt, nomatch=0]
                            tmp_std = tmp_std[order(tx_index, group, pos)]
                            tmp_std[, pos := 1:.N, by="tx_index"]
                            tmp_std[, bin := findInterval(pos, seq(0.5, max(pos)+0.5, length.out=nbin+1)), by="tx_id"]

                            tmp_std[, value := as.numeric(value)]
                            tmp_std[, value := value/(max(pos)/nbin), by=c("gene_id", "tx_id")] # Normalise by transcript bin width
                            tmp_std = tmp_std[, list(value = sum(value)), by=c("gene_id", "tx_id", "bin")]
                            
                            # tmp_std[, norm := value/sum(value), by=c("gene_id", "tx_id")] # Normalise by transcript total count
                            # tmp_std = tmp_std[!is.na(norm),]
                        }
                        return(tmp_std)
                    })

                tmp_list = rbindlist(tmp_list, idcol="strand")

                tmp_list[, bin := as.numeric(bin) + sum(bins[0:(i-1)])]

                return(tmp_list)
            })

        names(profiles) = names(genicRegions)
        profiles = rbindlist(profiles, idcol="region_id")

        profiles = profiles[, list(region_id, strand, bin, value = value/sum(value)), by=c("gene_id", "tx_id")]

        if( !dropEmpty )
            profiles[!is.finite(value), value := 0]
        profiles = profiles[is.finite(value),]

        profiles[, region_id := factor(region_id, levels=names(genicRegions))]

        switch(collapseBy,
             region = profiles[,
                            list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin")],
            gene = profiles[,
                            list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin", "gene_id")],
            transcript = profiles,
            recursive = (profiles[,
                            list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin", "gene_id")][,
                            list(Sum = sum(Mean),   PMean = mean(Mean), PSd = sqrt(sum((Mean-mean(Mean)^2))/(.N-1))), by=c("region_id", "bin")]))
}

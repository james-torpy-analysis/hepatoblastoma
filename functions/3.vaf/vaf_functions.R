reverseStrand <- function (x)
{

    old2new <- c("+" = "-", "-" = "+", "*" = "*")
    strand(x) <- old2new[as.character(strand(x))]
    return(x)

}

writeSam <- function (bam_in, selected, sam_out)
{

    ## read bam alignments
    ga <- readGAlignments(bam_in, use.names = TRUE, param = param)
    cols <- c(
        "qname",
        "flag",
        "rname",
        "pos",
        "mapq",
        "cigar",
        "mrnm",
        "mpos",
        "isize",
        "seq",
        "qual")

    ## subset to select reads
    index <- which(mcols(ga)$qname %in% selected)
    sam <- mcols(ga)[index, cols]

    ## write sam
    if (file.exists(sam_out)) file.remove(sam_out)

    system(paste0("samtools view -H ", bam_in, " > ", sam_out))
    write.table(sam,
                file = sam_out,
                append = TRUE,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    
    # convert to bam and index:
    bam_out <- gsub("sam", "bam", sam_out)
    system(paste0("samtools view -bh ", sam_out, " > ", bam_out))
    system(paste0("samtools index ", bam_out))
    
    # remove sams:
    system(paste0("rm ", sam_out, "*"))

}

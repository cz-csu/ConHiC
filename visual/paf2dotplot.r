#!/usr/bin/env Rscript
install.packages("cowplot",
                 repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
# 新增拼图包，如果未安装请先 install.packages("patchwork")
if (!require("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required for combining plots. Please install it: install.packages('patchwork')")
}

option_list <- list(
  make_option(c("-o","--output"), type="character", default="combined_dotplot",
              help="output filename prefix [%default]", 
              dest="output_filename"),
  make_option(c("-p","--plot-size"), type="numeric", default=15,
              help="Total plot width in inches [%default]",
              dest="plot_size"),
  make_option(c("-f","--flip"), action="store_true", default=FALSE,
              help="flip query if most alignments are in reverse complement [%default]",
              dest="flip"),
  make_option(c("-b","--break-point"), action="store_true", default=FALSE,
              help="show break points [%default]",
              dest="break_point"),
  make_option(c("-c","--identity-lower-color"), type="numeric", default=0,
              help="percent identity that lower than this value will be assigned the same color [%default]",
              dest="identity_lower"),
  make_option(c("-s", "--sort-by-refid"), action="store_true", default=FALSE,
              help="sort reference IDs in alphabetical order, default by length [%default]",
              dest="sortbyID"),
  make_option(c("-q", "--min-query-length"), type="numeric", default=400000,
              help="filter queries with total alignments less than cutoff X bp [%default]",
              dest="min_query_aln"),
  make_option(c("-m", "--min-alignment-length"), type="numeric", default=10000,
              help="filter alignments less than cutoff X bp [%default]",
              dest="min_align"),
  make_option(c("-r", "--min-ref-len"), type="numeric", default=1000000,
              help="filter references with length less than cutoff X bp [%default]",
              dest="min_ref_len"),
  make_option(c("-i", "--reference-ids"), type="character", default=NULL,
              help="comma-separated list of reference IDs to keep and order [%default]",
              dest="refIDs"),
  make_option(c("-t", "--titles"), type="character", default=NULL,
              help="comma-separated list of titles for each input file [%default]",
              dest="titles")
)

options(error=traceback)
parser <- OptionParser(usage = "%prog [options] file1.paf file2.paf ...", option_list = option_list)
opts = parse_args(parser, positional_arguments = TRUE) # 允许接收多个位置参数
opt = opts$options
input_files = opts$args

if(length(input_files) <= 0){
  cat(sprintf("Error: missing input file(s)!\n\n"))
  print_help(parser)
  quit()
}

# 辅助函数：从文件名提取标题（去掉前缀）
extract_title_from_filename <- function(filename) {
  # 去掉路径和扩展名
  base_name = tools::file_path_sans_ext(basename(filename))
  
  # 去掉_前面的部分，只保留_后面的部分
  parts <- strsplit(base_name, "_")[[1]]
  if (length(parts) > 1) {
    # 只保留_后面的部分
    return(paste(parts[-1], collapse = "_"))
  } else {
    # 如果没有下划线，返回整个文件名
    return(base_name)
  }
}

# --- 定义核心处理函数 ---
process_single_paf <- function(input_file, opt, plot_title = NULL, all_chromosomes = NULL) {
  
  if (file.access(input_file, mode=4) == -1){
    warning(sprintf("Skipping: input file %s does not exist or cannot be read!", input_file))
    return(NULL)
  }

  # read in alignments
  alignments = read.table(input_file, stringsAsFactors = F, row.names=NULL, fill = T, header = F)[, c(1:12)]
  alignments[, c(2:4, 7:12)] = apply(alignments[, c(2:4, 7:12)], 2, as.double)
  colnames(alignments)[1:12] = c("queryID","queryLen","queryStart","queryEnd","strand","refID","refLen","refStart","refEnd","numResidueMatches","lenAln","mapQ")
  
  # 修改：提取染色体号 - 只保留_后面的部分
  alignments$refID_original <- alignments$refID  # 保存原始名称
  alignments$refID <- sapply(alignments$refID, function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) > 1) {
      # 只保留_后面的部分
      return(paste(parts[-1], collapse = "_"))
    } else {
      return(x)
    }
  })

  # calculate similarity
  alignments$percentID = alignments$numResidueMatches / alignments$lenAln
  if (opt$identity_lower) {
      alignments$percentID[which(alignments$percentID < opt$identity_lower)] <- opt$identity_lower
  }

  queryStartTemp = alignments$queryStart
  alignments$queryStart[which(alignments$strand == "-")] = alignments$queryEnd[which(alignments$strand == "-")]
  alignments$queryEnd[which(alignments$strand == "-")] = queryStartTemp[which(alignments$strand == "-")]
  rm(queryStartTemp)

  # Filtering logic
  if(is.null(opt$refIDs)){
    if (opt$sortbyID){
      # 按字母顺序排序染色体
      refIDsToKeepOrdered = unique(sort(alignments$refID))
    }else{
      chromMax = tapply(alignments$refLen, alignments$refID, max)
      refIDsToKeepOrdered = names(sort(chromMax, decreasing = T))
    }
  }else{
    # 如果用户提供了refIDs参数，也需要处理
    refIDsToKeepOrdered = unlist(strsplit(opt$refIDs, ","))
    # 先提取染色体号
    refIDsToKeepOrdered_processed <- sapply(refIDsToKeepOrdered, function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) > 1) {
        return(paste(parts[-1], collapse = "_"))
      } else {
        return(x)
      }
    })
    alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered_processed),]
  }
  
  # 如果指定了统一的染色体顺序，确保所有子图使用相同的顺序
  if(!is.null(all_chromosomes)) {
    # 只保留all_chromosomes中存在的染色体
    alignments = alignments[which(alignments$refID %in% all_chromosomes),]
    if(nrow(alignments) == 0) return(NULL)
    # 使用统一的染色体顺序
    refIDsToKeepOrdered = all_chromosomes[all_chromosomes %in% unique(alignments$refID)]
  }

  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  # Check if filtering removes everything to avoid crash
  valid_queries = names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]
  if(length(valid_queries) == 0) return(NULL)
  
  alignments = alignments[which(alignments$queryID %in% valid_queries),]
  alignments = alignments[which(alignments$lenAln > opt$min_align),]
  alignments = alignments[which(alignments$refLen > opt$min_ref_len),]
  
  if (nrow(alignments) == 0) return(NULL)

  # Ref processing - 使用统一的染色体顺序
  refIDsToKeepOrdered = refIDsToKeepOrdered[which(refIDsToKeepOrdered %in% alignments$refID)]
  alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered)
  alignments = alignments[with(alignments,order(refID,refStart)),]
  
  # 重新计算chromMax，确保即使某个染色体在当前文件中没有数据也保留位置
  chromMax = numeric(length(refIDsToKeepOrdered))
  names(chromMax) = refIDsToKeepOrdered
  for(chr in refIDsToKeepOrdered) {
    if(chr %in% alignments$refID) {
      chromMax[chr] = max(alignments$refLen[alignments$refID == chr])
    } else {
      # 如果该染色体在当前文件中没有数据，使用最小长度
      chromMax[chr] = opt$min_ref_len
    }
  }

  alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
  alignments$refEnd2 = alignments$refEnd + sapply(as.character(alignments$refID), function(x) ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )

  # Query sorting
  alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID)))
  queryMaxAlnIndex = tapply(alignments$lenAln, alignments$queryID, which.max, simplify = F)
  alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
    function(x, i) alignments$refStart2[which(i == alignments$queryID)][x],
    queryMaxAlnIndex, names(queryMaxAlnIndex)
  ))])

  queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
  if(length(levels(alignments$refID)) > 1){
    queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
  } else {
    queryID_Ref = sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])
  }
  alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])
  queryMax = tapply(alignments$queryLen, alignments$queryID, max)

  if(opt$flip){
    queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
    queryRevComp = names(queryRevComp)[which(queryRevComp)]
    if(length(queryRevComp) > 0) {
        alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
        alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1
    }
  }

  alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
  alignments$queryEnd2 = alignments$queryEnd + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

  if (opt$break_point) {
    alignments$break_col = rep(0, length(alignments$percentID));
  }
  
  # --- Plotting ---
  gp = ggplot(alignments) + 
    theme_bw() + 
    theme(
      text = element_text(size = 8),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(angle = 15, size=6),
      axis.text.x = element_text(hjust = 1, angle = 45, size=6),
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
      legend.position = "none"  # 所有子图都不显示图例
    ) +
    scale_color_distiller(palette = "Spectral")

  # 为子图添加标题
  if (!is.null(plot_title)) {
    # 如果用户提供了标题，使用用户提供的标题
    gp = gp + ggtitle(plot_title)
  } else {
    # 如果没有提供标题，使用提取的文件名（去掉前缀）
    clean_title <- extract_title_from_filename(input_file)
    gp = gp + ggtitle(clean_title)
  }

  # X axis - 使用统一的染色体顺序
  if (length(unique(alignments$refID)) == 1){
    reflen = unique(alignments$refLen)
    gp = gp + scale_x_continuous(expand = c(0, 0), limits = c(0, reflen), breaks = pretty(c(0, reflen), n=3), labels = function(x) paste0(round(x/1e6), "Mb")) +
      xlab(unique(alignments$refID))
  }else{
    # 使用统一的染色体顺序和标签
    chromosome_labels <- levels(alignments$refID)
    
    gp = gp + 
      theme(panel.grid.major.x=element_blank()) +
      geom_vline(xintercept = cumsum(as.numeric(chromMax)), col="#ebebeb") + 
      scale_x_continuous(
        expand = c(0, 0), 
        limits = c(0, sum(as.numeric(chromMax))), 
        breaks = cumsum(as.numeric(chromMax)) - chromMax/2,
        labels = chromosome_labels
      ) + 
      xlab("Chromosome")
  }

  # Y axis
  if (length(unique(alignments$queryID)) == 1){
    queryLen = unique(alignments$queryLen)
    gp = gp + scale_y_continuous(expand = c(0, 0), limits = c(0, queryLen), breaks = pretty(c(0, queryLen), n=3), labels = function(x) paste0(round(x/1e6), "Mb")) +
      ylab(unique(alignments$queryID))
  }else{
      gp = gp + 
      theme(panel.grid.major.y=element_blank()) +
      geom_hline(yintercept = cumsum(as.numeric(queryMax)), col="#ebebeb") + 
      scale_y_continuous(
        expand = c(0, 0), 
        limits = c(0, sum(as.numeric(queryMax))), 
        breaks = cumsum(as.numeric(queryMax)) - queryMax/2,
        labels = substr(levels(alignments$queryID), start = 1, stop = 15)
      ) + 
      ylab("Query")
  }

  gp = gp + geom_segment(aes(x = refStart2, xend = refEnd2, y = queryStart2, yend = queryEnd2, color = percentID))
  
  if (opt$break_point) {
    gp = gp + geom_point(aes(x = refStart2, y = queryStart2, color = break_col), size = 0.1) +
      geom_point(aes(x = refEnd2, y = queryEnd2, color = break_col), size = 0.1)
  }
  
  return(gp)
}

# --- Main Loop ---

plot_list <- list()
cat(sprintf("Processing %d input files...\n", length(input_files)))

# 处理标题参数
titles <- NULL
if (!is.null(opt$titles)) {
  titles <- unlist(strsplit(opt$titles, ","))
  if (length(titles) != length(input_files)) {
    warning(sprintf("Number of titles (%d) does not match number of input files (%d). Using default titles.", length(titles), length(input_files)))
    titles <- NULL
  }
}

# 第一步：先扫描所有文件，获取统一的染色体集合
cat("Scanning all files to get unified chromosome set...\n")
all_chromosomes <- c()
for (f in input_files) {
  if (file.access(f, mode=4) == -1) next
  
  alignments = tryCatch({
    temp = read.table(f, stringsAsFactors = F, row.names=NULL, fill = T, header = F)[, c(1:12)]
    colnames(temp)[1:12] = c("queryID","queryLen","queryStart","queryEnd","strand","refID","refLen","refStart","refEnd","numResidueMatches","lenAln","mapQ")
    
    # 提取染色体号
    temp$refID <- sapply(temp$refID, function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) > 1) {
        return(paste(parts[-1], collapse = "_"))
      } else {
        return(x)
      }
    })
    
    unique(temp$refID)
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", f, e$message))
    return(NULL)
  })
  
  if(!is.null(alignments)) {
    all_chromosomes <- unique(c(all_chromosomes, alignments))
  }
}

# 确保染色体按数字顺序排序（如果是数字的话）
if(length(all_chromosomes) > 0) {
  # 尝试按数字排序
  chrom_nums <- suppressWarnings(as.numeric(all_chromosomes))
  if(!any(is.na(chrom_nums))) {
    # 都是数字，按数字排序
    all_chromosomes <- all_chromosomes[order(chrom_nums)]
  } else {
    # 不是纯数字，按字母排序
    all_chromosomes <- sort(all_chromosomes)
  }
  cat(sprintf("Found chromosomes: %s\n", paste(all_chromosomes, collapse = ", ")))
}

# 第二步：使用统一的染色体顺序处理每个文件
for (i in seq_along(input_files)) {
    f <- input_files[i]
    cat(sprintf("  -> Processing %s\n", f))
    
    # 获取对应标题
    plot_title <- if (!is.null(titles)) titles[i] else NULL
    
    p <- process_single_paf(f, opt, plot_title, all_chromosomes)
    if (!is.null(p)) {
        plot_list[[f]] <- p
    } else {
        cat(sprintf("     Warning: No valid data in %s\n", f))
    }
}

if (length(plot_list) == 0) {
    stop("No plots generated.")
}

cat("Combining plots...\n")

# 方法：创建一个特殊的图作为图例
# 首先，从第一个图中提取颜色范围信息
if (length(plot_list) > 0) {
  # 创建一个空的ggplot对象用于生成图例
  legend_plot <- ggplot(data.frame(x = 0:1, y = 0:1, col = 0:1), 
                       aes(x = x, y = y, color = col)) +
    geom_point() +
    scale_color_distiller(palette = "Spectral", 
                         name = "Identity (%)",
                         guide = guide_colorbar(
                           title.position = "top",
                           title.hjust = 0.5,
                           barwidth = unit(4, "cm"),
                           barheight = unit(0.3, "cm")
                         )) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          legend.box.margin = margin(t = 5, r = 0, b = 5, l = 0))
  
  # 提取图例
  legend <- cowplot::get_legend(legend_plot)
  
  # 移除所有子图的图例
  for (i in 1:length(plot_list)) {
    plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  }
  
  # 使用Patchwork组合图片
  combined_plots <- wrap_plots(plot_list, ncol = 6)
  
  # 创建最终布局：主图在上，图例在下
  final_plot <- plot_spacer() / combined_plots / legend +
    plot_layout(heights = c(0.05, 1, 0.1))  # 调整各部分高度比例
  
} else {
  final_plot <- wrap_plots(plot_list, ncol = 6)
}

# 计算合适的画布大小
out_w <- opt$plot_size
n_rows <- ceiling(length(plot_list) / 6)
# 调整高度计算
out_h <- (out_w / 6) * n_rows * 0.9 + 1.0  # 增加高度以适应底部图例

ggsave(filename = paste0(opt$output_filename, ".pdf"), plot = final_plot, 
       width = out_w, height = out_h, units = "in", dpi = 300, limitsize = FALSE)

cat(sprintf("Done! Saved to %s.pdf\n", opt$output_filename))
options(warn=0)
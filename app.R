# Sys.setlocale("LC_ALL", "English_United States.1252")
# Sys.setenv(LANG = "en_US.UTF-8")
# options(encoding = "UTF-8")
# options(timeout = 600)
# options(rsconnect.http.timeout = 600)
# rsconnect::deployApp()

# app.R --------------------------------------------------------------------
# Volcano Explorer (2 pages):
#   1) Load & Process
#   2) Volcano explorer
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyWidgets)
  library(shinythemes)
  library(DT)
  library(vroom)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(plotly)
  library(viridisLite)   
  library(tibble)
  library(shinyBS)
  library(shinycssloaders)
  library(limma)
  library(RColorBrewer)
  library(zip)
})

options(shiny.maxRequestSize = 1024 * 1024^2)

# ----------------------------- Helpers -----------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

read_msdial_robust <- function(path) {
  max_cols <- NA
  try({
    n_fields <- utils::count.fields(path, sep = ",")
    max_cols <- max(n_fields, na.rm = TRUE)
  }, silent = TRUE)
  
  if (!is.na(max_cols)) {
    col_names <- paste0("V", seq_len(max_cols))
    df <- utils::read.csv(path, header = FALSE, col.names = col_names, 
                          stringsAsFactors = FALSE, colClasses = "character", na.strings = "")
  } else {
    df <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE, 
                          colClasses = "character", fill = TRUE)
  }
  
  hdr_i <- NA
  for (i in 1:min(50, nrow(df))) {
    row_txt <- tolower(as.character(unlist(df[i, ])))
    row_txt <- gsub("[^a-z0-9]", "", row_txt)
    if ("averagemz" %in% row_txt || "alignmentid" %in% row_txt) {
      hdr_i <- i
      break
    }
  }
  
  if (!is.na(hdr_i)) {
    new_names <- trimws(as.character(unlist(df[hdr_i, , drop = TRUE])))
    mask_bad <- is.na(new_names) | new_names == "" | new_names == "NA"
    if (any(mask_bad)) {
      new_names[mask_bad] <- paste0("Unknown_", seq_len(sum(mask_bad)))
    }
    new_names <- make.unique(new_names, sep = "_")
    names(df) <- new_names
    
    if (hdr_i < nrow(df)) {
      df <- df[(hdr_i + 1):nrow(df), , drop = FALSE]
    } else {
      df <- df[0, , drop = FALSE]
    }
  }
  df
}

clean_mzmine_export <- function(df) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) > 0) {
    last <- df[[ncol(df)]]
    if (all(is.na(last)) || all(trimws(as.character(last)) == "")) {
      if (grepl("^Unnamed", names(df)[ncol(df)]) || names(df)[ncol(df)] == "") {
        df <- df[, -ncol(df), drop = FALSE]
      }
    }
  }
  df
}

multi_sample_idx <- function(cols, kws) {
  kws <- as.character(kws)
  kws <- kws[nzchar(kws)]
  if (!length(kws)) return(integer(0))
  hits <- Reduce(`|`, lapply(kws, function(k) grepl(k, cols, fixed = TRUE)))
  which(hits)
}

clean_sample_names <- function(x) {
  x <- gsub(" Peak area$", "", x, ignore.case = TRUE)
  x <- gsub("\\.(mzML|mzXML|raw|cdf)$", "", x, ignore.case = TRUE)
  x <- trimws(x)
  x
}

make_label_table <- function(sample_names, labels) {
  tibble::tibble(
    Sample = as.character(sample_names),
    Label  = trimws(as.character(labels))
  )
}

labels_from_sample_names_or_raw <- function(sample_names,
                                            token_sep = "_",
                                            token_index = 2,
                                            clean_names = TRUE) {
  sn <- if (isTRUE(clean_names)) clean_sample_names(sample_names) else sample_names
  sep <- token_sep %||% "_"
  idx <- as.integer(token_index %||% 2)

  parts <- strsplit(sn, sep, fixed = TRUE)
  ok <- vapply(parts, function(z) length(z) >= idx, logical(1))

  labs <- vapply(seq_along(parts), function(i) {
    if (ok[i] && nzchar(parts[[i]][[idx]])) {
      parts[[i]][[idx]]
    } else {
      sn[i]
    }
  }, character(1))

  labs
}

stop_if_one_group <- function(labs) {
  labs <- trimws(as.character(labs))
  labs <- labs[nzchar(labs)]
  u <- unique(labs)

  if (length(u) < 2) {
    shiny::showNotification(
      paste0("Need at least 2 different groups in 'Label' to run statistics."),
      type = "error",
      duration = 8
    )
    return(TRUE)  
  }
  FALSE
}

read_onecol_csv <- function(path) {
  v <- suppressWarnings(vroom::vroom(path, col_names = FALSE, delim = ","))
  as.character(v[[1]])
}

parse_suffix_list <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

clean_sample_names_for_match <- function(x, enabled = FALSE, remove_suffixes = NULL) {
  x0 <- as.character(x)

  # normalize even without suffix cleaning
  out <- x0
  out <- gsub("\u00A0", " ", out, fixed = TRUE)
  out <- gsub("[[:space:]]+", " ", out)
  out <- trimws(out)
  out <- gsub('^"|"$', "", out)

  if (!isTRUE(enabled)) {
    return(out)
  }

  default_suffixes <- c(
    " Peak area", " Peak Area", "Peak area", "Peak Area",
    " Peak height", " Peak Height", "Peak height", "Peak Height",
    "_Area", "_Height",
    " Area", " Height",
    ".mzML", ".mzXML", ".raw", ".RAW",
    ".cdf", ".CDF",
    ".mzData", ".mzdata",
    ".wiff", ".WIFF",
    ".d", ".D"
  )

  suffixes <- unique(c(parse_suffix_list(remove_suffixes), default_suffixes))
  suffixes <- suffixes[nzchar(suffixes)]

  strip_one_suffix <- function(v, sfx) {
    n <- nchar(sfx)
    if (!is.finite(n) || n < 1) return(v)

    hit <- nchar(v) >= n &
      tolower(substr(v, nchar(v) - n + 1, nchar(v))) == tolower(sfx)

    v[hit] <- substr(v[hit], 1, nchar(v[hit]) - n)
    v <- gsub("[[:space:]]+", " ", v)
    trimws(v)
  }

  # repeat because names can end as: ".mzML Peak area"
  # first pass removes "Peak area", second pass removes ".mzML"
  for (pass in seq_len(20)) {
    old <- out

    for (sfx in suffixes) {
      out <- strip_one_suffix(out, sfx)
    }

    if (identical(old, out)) break
  }

  out[!nzchar(out)] <- x0[!nzchar(out)]
  out
}

guess_metadata_sample_col <- function(cols) {
  candidates <- c(
    "Sample", "sample",
    "SampleName", "sample_name",
    "Filename", "FileName", "filename",
    "File", "Name", "Injection", "Run"
  )

  hit <- candidates[candidates %in% cols]
  if (length(hit)) hit[1] else cols[1]
}

guess_metadata_label_col <- function(cols, sample_col = NULL) {
  cols2 <- setdiff(cols, sample_col)

  candidates <- c(
    "Condition", "condition",
    "Label", "label",
    "Group", "group",
    "Treatment", "treatment",
    "Class", "class"
  )

  hit <- candidates[candidates %in% cols2]
  if (length(hit)) hit[1] else cols2[1]
}

read_metadata_csv <- function(upload, context = "metadata labels") {
  req(upload)

  ext <- tolower(tools::file_ext(upload$name))
  validate(
    need(ext == "csv", paste0("Upload a .csv file for ", context, "."))
  )

  as.data.frame(
    vroom::vroom(
      upload$datapath,
      delim = ",",
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

metadata_labels_by_sample <- function(upload,
                                      sample_names,
                                      sample_col,
                                      label_col,
                                      clean_enabled = FALSE,
                                      remove_suffixes = NULL,
                                      context = "metadata labels") {
  meta <- read_metadata_csv(upload, context = context)

  validate(
    need(sample_col %in% names(meta), "Selected metadata sample-name column was not found."),
    need(label_col %in% names(meta), "Selected metadata label column was not found.")
  )

  app_key <- clean_sample_names_for_match(
    sample_names,
    enabled = clean_enabled,
    remove_suffixes = remove_suffixes
  )

  meta_key <- clean_sample_names_for_match(
    meta[[sample_col]],
    enabled = clean_enabled,
    remove_suffixes = remove_suffixes
  )

  validate(
    need(!anyDuplicated(app_key),
         "Detected sample names are duplicated after optional cleaning."),
    need(!anyDuplicated(meta_key),
         "Metadata sample names are duplicated after optional cleaning.")
  )

  idx <- match(app_key, meta_key)

  if (any(is.na(idx))) {
    missing_samples <- app_key[is.na(idx)]

    validate(
      need(
        FALSE,
        paste0(
          "Metadata file is missing these sample names: ",
          paste(head(missing_samples, 10), collapse = ", "),
          if (length(missing_samples) > 10) " ..." else ""
        )
      )
    )
  }

  labs <- trimws(as.character(meta[[label_col]][idx]))

  validate(
    need(!any(is.na(labs) | labs == ""),
         "Selected metadata label column contains empty values.")
  )

  labs
}

finite_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  c(min(x), max(x))
}

guess_col <- function(cols, candidates) {
  if (!length(cols)) return(NULL)
  cols_l <- tolower(cols)
  cols_norm <- gsub("[^a-z0-9]", "", cols_l)
  cand_norm <- gsub("[^a-z0-9]", "", tolower(candidates))
  
  for (cand in cand_norm) {
    j <- which(cols_norm == cand)
    if (length(j)) return(cols[j[1]])
  }
  for (cand in cand_norm) {
    j <- which(grepl(cand, cols_norm, fixed = TRUE))
    if (length(j)) return(cols[j[1]])
  }
  NULL 
}

safe_ttest_p <- function(x, g, paired = FALSE, var.equal = FALSE) {
  x <- as.numeric(x)
  g <- as.factor(g)
  if (length(unique(g)) != 2) return(1)
  a <- x[g == levels(g)[1]]
  b <- x[g == levels(g)[2]]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(1)
  if (stats::sd(a) == 0 && stats::sd(b) == 0) return(1)
  out <- try(stats::t.test(x ~ g, paired = paired, var.equal = var.equal)$p.value, silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) 1 else as.numeric(out)
}

safe_wilcox_p <- function(x, g, paired = FALSE) {
  x <- as.numeric(x)
  g <- as.factor(g)
  if (length(unique(g)) != 2) return(1)
  a <- x[g == levels(g)[1]]
  b <- x[g == levels(g)[2]]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 1 || length(b) < 1) return(1)
  out <- try(stats::wilcox.test(x ~ g, paired = paired)$p.value, silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) 1 else as.numeric(out)
}

parse_feature_table_to_matrix <- function(raw_df,
                                         row_id_col,
                                         mz_col,
                                         rt_col,
                                         sample_keywords) {
  raw_df <- clean_mzmine_export(raw_df)
  cols <- names(raw_df)

  validate(
    need(row_id_col %in% cols, "Row ID column not found."),
    need(mz_col %in% cols,     "m/z column not found."),
    need(rt_col %in% cols,     "rt column not found.")
  )

  sidx <- multi_sample_idx(cols, sample_keywords)
  validate(need(length(sidx) > 0,
                sprintf("No sample columns matched keywords: %s",
                        paste(sample_keywords, collapse = ", "))))

  sample_cols <- cols[sidx]

  # samples x features
  # Using data.table::transpose exactly as original
  mat <- as.data.frame(data.table::transpose(raw_df[, sample_cols, drop = FALSE]),
                       check.names = FALSE, stringsAsFactors = FALSE)
  rownames(mat) <- sample_cols

  # feature definition
  id <- as.character(raw_df[[row_id_col]])
  mz <- suppressWarnings(as.numeric(raw_df[[mz_col]]))
  rt <- suppressWarnings(as.numeric(raw_df[[rt_col]]))
  
  # safety for NA
  mz[is.na(mz)] <- 0
  rt[is.na(rt)] <- 0

  feat_raw <- paste0(round(mz, 4), "@", round(rt, 2))
  Feature <- make.unique(feat_raw, sep = "_")
  colnames(mat) <- Feature

  fmap <- tibble(
    id = id,
    mz = mz,
    RT = rt,
    Feature = Feature
  )

  mat[] <- lapply(mat, function(z) suppressWarnings(as.numeric(z)))
  mat[is.na(mat)] <- 0

  list(mat = mat, fmap = fmap, raw = raw_df)
}

impute_lod_random <- function(X,
                              noise_mode = c("quantile", "manual"),
                              noise_quantile = 0.25,
                              noise_manual = 50,
                              sd_val = 30,
                              seed = 1234) {
  noise_mode <- match.arg(noise_mode)
  X <- as.matrix(X)
  X[X == 0] <- NA

  nz <- 1:min(X, na.rm = T) 
  if (!length(nz)) return(X * 0)

  if (noise_mode == "quantile") {
    noise <- as.numeric(stats::quantile(nz, probs = noise_quantile, na.rm = TRUE))
  } else {
    noise <- as.numeric(noise_manual)
  }
  noise <- ifelse(is.finite(noise) && noise > 0, noise, 1)

  sd_val <- as.numeric(sd_val)
  sd_val <- ifelse(is.finite(sd_val) && sd_val >= 0, sd_val, 0)

  set.seed(seed)
  X[is.na(X)] <- 0
  idx <- X == 0
  imp <- abs(stats::rnorm(sum(idx), mean = noise, sd = sd_val))
  X[idx] <- imp
  X
}

compute_stats_long <- function(df_used,
                               test = c("Student", "Wilcoxon", "limma"),
                               adj  = c(
                                 "BH", "holm", "hochberg", "hommel",
                                 "bonferroni", "BY", "fdr", "none"),
                               paired = FALSE,
                               eqvar = FALSE,
                               pseudocount = 1.1,
                               log2_test = FALSE,
                               scale_data = FALSE,
                               ref_group = NULL,
                               comparisons = NULL) {
  test <- match.arg(test)
  adj  <- match.arg(adj)

  validate(need("Label" %in% names(df_used), "Internal error: Label column missing."))

  feats <- setdiff(colnames(df_used), "Label")
  validate(need(length(feats) > 0, "No feature columns detected."))

  gr <- as.factor(df_used$Label)
  lev <- levels(gr)
  validate(need(length(lev) >= 2, "Need at least 2 Label groups to run statistics."))

  if (!is.null(comparisons)) {

  comparisons <- as.data.frame(
    comparisons,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(
      all(c("Group_num", "Group_den") %in% names(comparisons)),
      paste0(
        "Manual comparisons must contain ",
        "Group_num and Group_den columns."
      )
    )
  )

  comparisons <- comparisons %>%
    dplyr::transmute(
      Group_num = trimws(as.character(Group_num)),
      Group_den = trimws(as.character(Group_den))
    ) %>%
    dplyr::filter(
      Group_num %in% lev,
      Group_den %in% lev,
      Group_num != Group_den
    ) %>%
    dplyr::distinct()

  validate(
    need(
      nrow(comparisons) > 0,
      "Select at least one valid manual comparison."
    )
  )

  comb <- as.matrix(
    comparisons[, c("Group_num", "Group_den"), drop = FALSE]
  )

} else if (!is.null(ref_group) && ref_group %in% lev) {

  others <- setdiff(lev, ref_group)

  # First group / second group
  comb <- cbind(
    rep(ref_group, length(others)),
    others
  )

} else {

  # Fallback when neither method is supplied
  comb <- t(utils::combn(lev, 2))
}

  out_list <- vector("list", nrow(comb))

  for (i in seq_len(nrow(comb))) {
    gnum <- comb[i, 1] # Reference (Denominator)
    gden <- comb[i, 2] # Comparison (Numerator)
    comp <- paste0(gnum, " / ", gden)

    sub <- df_used[df_used$Label %in% c(gden, gnum), c("Label", feats), drop = FALSE]
    sub$Label <- factor(as.character(sub$Label), levels = c(gden, gnum))

    # choose scale for hypothesis tests
    if (isTRUE(log2_test)) {
      sub_test <- sub
      sub_test[feats] <- lapply(sub_test[feats], function(z) log2(as.numeric(z) + pseudocount))
    } else {
      sub_test <- sub
    }

    if (isTRUE(scale_data)) {
      sub_test[feats] <- lapply(sub_test[feats], function(z) {
        vec <- as.numeric(z)
        s <- stats::sd(vec, na.rm = TRUE)
        # Prevent division by zero if variance is 0
        if (is.na(s) || s == 0) return(rep(0, length(vec))) 
        as.numeric(scale(vec, center = TRUE, scale = TRUE))
      })
    }
    
    # p-values (raw or log2 scale, depending on toggle)
    if (test == "Student") {
      p <- vapply(
        feats,
        function(f) safe_ttest_p(sub_test[[f]], sub_test$Label, paired = paired, var.equal = eqvar),
        numeric(1)
      )
    } else if (test == "Wilcoxon") {
      p <- vapply(
        feats,
        function(f) safe_wilcox_p(sub_test[[f]], sub_test$Label, paired = paired),
        numeric(1)
      )
    } else if (test == "limma") {
      # limma requires transposed matrix: features in rows, samples in columns
      emat <- t(as.matrix(sub_test[feats]))
      
      # Create design matrix for the two groups
      design <- stats::model.matrix(~ sub_test$Label)
      
      # Run Moderated t-test pipeline
      fit <- limma::lmFit(emat, design)
      fit <- limma::eBayes(fit)
      
      # Extract p-values for the comparison (second coefficient)
      p <- fit$p.value[, 2]
      names(p) <- feats # Ensure order matches
    }

    padj <- if (adj == "none") p else stats::p.adjust(p, method = adj)

    # group means on RAW scale (keep as-is)
    Xnum <- sub[sub$Label == gnum, feats, drop = FALSE]
    Xden <- sub[sub$Label == gden, feats, drop = FALSE]
    mean_num_raw <- colMeans(Xnum, na.rm = TRUE)
    mean_den_raw <- colMeans(Xden, na.rm = TRUE)

    # FC (keep as-is)
    mean_num_log2 <- log2(colMeans(as.matrix(Xnum), na.rm = TRUE) + pseudocount)
    mean_den_log2 <- log2(colMeans(as.matrix(Xden), na.rm = TRUE) + pseudocount)
    FC <- mean_num_log2 - mean_den_log2

    Mean <- 0.5 * (mean_num_raw + mean_den_raw)

    dd <- tibble(
      Groups = comp,
      Group_num = gnum,
      Group_den = gden,
      Feature = feats,
      `Adj.p-value` = as.numeric(padj),
      Mean = as.numeric(Mean),
      mean_num = as.numeric(mean_num_raw),
      mean_den = as.numeric(mean_den_raw),
      FC = as.numeric(FC),
      TestScale = if (isTRUE(log2_test)) "log2" else "raw"
    )

    dd$`Adj.p-value.log` <- -log10(pmax(dd$`Adj.p-value`, .Machine$double.xmin))
    dd$Significant <- (dd$`Adj.p-value` <= 0.05) & (abs(dd$FC) >= 1)

    out_list[[i]] <- dd
  }

  dplyr::bind_rows(out_list)
}

volcano_to_wide_if_needed <- function(volc) {

  volc <- as.data.frame(
    volc,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"Groups" %in% names(volc) || nrow(volc) == 0) {
    return(volc)
  }

  static_cols <- intersect(
    c(
      "Feature",
      "id",
      "mz",
      "RT",
      "NPC#class",
      "ClassyFire#class"
    ),
    names(volc)
  )

  n_comp <- dplyr::n_distinct(volc$Groups)

  if (n_comp == 1) {

  group_num_name <- unique(as.character(volc$Group_num))
  group_den_name <- unique(as.character(volc$Group_den))

  if (
    length(group_num_name) == 1 &&
    length(group_den_name) == 1
  ) {
    names(volc)[names(volc) == "mean_num"] <-
      paste0("Mean_", group_num_name)

    names(volc)[names(volc) == "mean_den"] <-
      paste0("Mean_", group_den_name)
  }

  # These columns are not needed in the exported volcano table
  volc <- volc %>%
    dplyr::select(
      -dplyr::any_of(c("Mean", "TestScale",  "Adj.p-value.log"))
    )

  return(volc)
}

  comparison_metrics <- intersect(
  c(
    "FC",
    "Adj.p-value",
    "Significant"
  ),
  names(volc)
)

  comparison_wide <- volc %>%
    dplyr::select(
      dplyr::all_of(static_cols),
      Groups,
      dplyr::all_of(comparison_metrics)
    ) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(static_cols),
      names_from = Groups,
      values_from = dplyr::all_of(comparison_metrics),
      names_glue = "{.value}__{Groups}",
      names_repair = "minimal"
    )

  group_means <- dplyr::bind_rows(

    volc %>%
      dplyr::transmute(
        dplyr::across(dplyr::all_of(static_cols)),
        group_name = as.character(Group_num),
        group_mean = as.numeric(mean_num)
      ),

    volc %>%
      dplyr::transmute(
        dplyr::across(dplyr::all_of(static_cols)),
        group_name = as.character(Group_den),
        group_mean = as.numeric(mean_den)
      )

  ) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(static_cols),
      names_from = group_name,
      values_from = group_mean,
      names_glue = "Mean_{group_name}",
      names_repair = "minimal"
    )

  dplyr::left_join(
    comparison_wide,
    group_means,
    by = static_cols
  )
}

make_dark2_color_map <- function(conditions) {
  conditions <- unique(as.character(conditions))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]

  if (!length(conditions)) {
    return(tibble::tibble(
      Condition = character(),
      Colour = character()
    ))
  }

  base_cols <- RColorBrewer::brewer.pal(8, "Dark2")

  cols <- if (length(conditions) <= 8) {
    base_cols[seq_along(conditions)]
  } else {
    grDevices::colorRampPalette(base_cols)(length(conditions))
  }

  tibble::tibble(
    Condition = conditions,
    Colour = cols
  )
}

palette_choices <- c(
  "Dark2", "Set1", "Set2", "Paired", "Accent", "Pastel1", "Pastel2",
  "viridis", "plasma", "magma", "inferno", "cividis"
)

make_palette <- function(pal = "Dark2", n = 8) {
  pal <- pal %||% "Dark2"
  n <- max(1, as.integer(n))

  if (pal %in% c("viridis", "plasma", "magma", "inferno", "cividis")) {
    return(get(pal, asNamespace("viridisLite"))(n))
  }

  if (pal %in% rownames(RColorBrewer::brewer.pal.info)) {
    maxn <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
    base <- RColorBrewer::brewer.pal(max(3, min(maxn, n)), pal)

    if (n > length(base)) {
      base <- grDevices::colorRampPalette(base)(n)
    }

    return(base[seq_len(n)])
  }

  viridisLite::viridis(n)
}

make_autoplotter_data <- function(df_used, sample_names) {
  df_used <- as.data.frame(df_used, check.names = FALSE, stringsAsFactors = FALSE)

  feature_cols <- setdiff(names(df_used), "Label")

  out <- df_used[, feature_cols, drop = FALSE]

  out <- cbind(
    Sample = as.character(sample_names),
    out
  )

  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

make_autoplotter_metadata <- function(df_used, sample_names) {
  df_used <- as.data.frame(df_used, check.names = FALSE, stringsAsFactors = FALSE)

  labs <- as.character(df_used$Label)

  cmap <- make_dark2_color_map(labs)

  meas_rep <- ave(
    seq_along(labs),
    labs,
    FUN = seq_along
  )

  first_condition_row <- !duplicated(labs)

  plot_order_full <- match(labs, cmap$Condition)
  colour_full <- cmap$Colour[match(labs, cmap$Condition)]

  tibble::tibble(
    Filename = as.character(sample_names),
    Condition = labs,
    MeasRep = as.integer(meas_rep),
    ExpRep = NA_character_,

    # Filled only once per unique Condition
    PlottingOrder = ifelse(first_condition_row, plot_order_full, NA_integer_),
    Colour = ifelse(first_condition_row, colour_full, NA_character_),

    Relative_Correction = NA_real_,
    Absolute_Correction = NA_real_
  )
}

make_autoplotter_name_map <- function(fmap, volcano = NULL) {
  fmap <- as.data.frame(fmap, check.names = FALSE, stringsAsFactors = FALSE)

  out <- tibble::tibble(
    Name = as.character(fmap$Feature)
  )

  if (!is.null(volcano) && all(c("Feature", "NPC#class", "ClassyFire#class") %in% names(volcano))) {
    ann <- volcano %>%
      dplyr::select(
        Feature,
        `NPC#class`,
        `ClassyFire#class`
      ) %>%
      dplyr::distinct(Feature, .keep_all = TRUE) %>%
      dplyr::mutate(
        `NPC#class` = dplyr::if_else(
          is.na(`NPC#class`) | `NPC#class` == "Not provided",
          "",
          as.character(`NPC#class`)
        ),
        `ClassyFire#class` = dplyr::if_else(
          is.na(`ClassyFire#class`) | `ClassyFire#class` == "Not provided",
          "",
          as.character(`ClassyFire#class`)
        )
      )

    out <- out %>%
      dplyr::left_join(ann, by = c("Name" = "Feature"))

  } else {
    out$`NPC#class` <- ""
    out$`ClassyFire#class` <- ""
  }

  out$`NPC#class`[is.na(out$`NPC#class`)] <- ""
  out$`ClassyFire#class`[is.na(out$`ClassyFire#class`)] <- ""

  out
}

# ----------------------------- UI -----------------------------------------

ui <- fluidPage(
  useShinyjs(),

tags$head(tags$style(HTML("
  .app-footer { position: fixed; left:0; right:0; bottom:0; 
                text-align:center; font-size:12px; opacity:0.75;
                padding:8px; background: rgba(255,255,255,0.8);
                border-top: 1px solid #ddd; z-index: 9999; }
  body { padding-bottom: 45px; }
"))),

tags$head(tags$style(HTML("
  /* Editable Labels table: prevent white-on-white editing issue */

  #labels_table.html-widget.datatables {
    background-color: transparent !important;
  }

  #labels_table .dataTables_wrapper,
  #labels_table table.dataTable,
  #labels_table .dataTables_scroll,
  #labels_table .dataTables_scrollHead,
  #labels_table .dataTables_scrollBody {
    background-color: #ffffff !important;
    color: #2c3e50 !important;
  }

  #labels_table table.dataTable th,
  #labels_table table.dataTable td {
    background-color: #ffffff !important;
    color: #2c3e50 !important;
  }

  /* Cell when focused / double-clicked / edited */
  #labels_table table.dataTable tbody td.focus,
  #labels_table table.dataTable tbody td:focus,
  #labels_table table.dataTable tbody tr.selected td,
  #labels_table table.dataTable tbody td.selected {
    background-color: #ffffff !important;
    color: #000000 !important;
    box-shadow: inset 0 0 0 2px #66CDAA !important;
  }

  /* Input box created during editing */
  #labels_table input,
  #labels_table textarea,
  #labels_table .dataTables_wrapper input,
  #labels_table .dataTables_wrapper textarea {
    background-color: #ffffff !important;
    color: #000000 !important;
    -webkit-text-fill-color: #000000 !important;
    caret-color: #000000 !important;
    border: 1px solid #66CDAA !important;
  }

  /* Keep search / info / pagination readable if shown */
  #labels_table .dataTables_length,
  #labels_table .dataTables_filter,
  #labels_table .dataTables_info,
  #labels_table .dataTables_paginate {
    color: #000000 !important;
    font-weight: bold;
    padding: 5px;
  }
  
  /* Colored pickerInput buttons for Volcano filters */

.bootstrap-select > .dropdown-toggle[data-id='sel_feat'] {
  background-color: #66CDAA !important;
  border-color: #45b894 !important;
  color: white !important;
  font-weight: bold !important;
}

.bootstrap-select > .dropdown-toggle[data-id='npc_filter_values'] {
  background-color: #18bc9c !important;
  border-color: #13a085 !important;
  color: white !important;
  font-weight: bold !important;
}

.bootstrap-select > .dropdown-toggle[data-id='classyfire_filter_values'] {
  background-color: #18bc9c !important;
  border-color: #13a085 !important;
  color: white !important;
  font-weight: bold !important;
}

/* Highlight selected options inside dropdown */
.bootstrap-select .dropdown-menu li.selected a {
  background-color: #dff7ef !important;
  color: #000000 !important;
  font-weight: bold !important;
}

.bootstrap-select .dropdown-menu li.selected a span.check-mark {
  color: #18bc9c !important;
}
  
"))),

tags$head(
  tags$title("Metabocano"),
  tags$link(rel = "icon", type = "image/png",
            href = "https://raw.githubusercontent.com/plyush1993/Metabocano/main/inst/www/sticker.png")
),

tags$head(
    tags$style(HTML("
      /* Increase max-width to prevent wrapping and align text left */
      .tooltip-inner {
        max-width: none !important;
        white-space: nowrap;
        text-align: left !important;
        font-size: 18px;
      }
    "))
  ),

tags$head(tags$style(HTML("
  /* make disabled download links truly inactive */
  a.shiny-download-link.disabled, 
  .shiny-download-link.disabled {
    pointer-events: none !important;
    opacity: 0.5 !important;
    cursor: not-allowed !important;
  }
"))),

div(
  class = "app-footer",
  HTML('Created by: Ivan Plyushchenko &nbsp;|&nbsp;
       <a href="https://github.com/plyush1993/Metabocano" target="_blank">GitHub repository</a>')
),
  
  div(
  style = "
    width: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    margin-bottom: 20px;
  ",
  
  tags$img(
    src = 'https://raw.githubusercontent.com/plyush1993/Metabocano/main/inst/www/sticker.png',
    height = '150px',
    style = 'margin-right: 20px;'
  ),
  
  div(
    style = '
      font-size: 32px;
      font-weight: 900;
      color: #66CDAA; 
      text-align: center;
    ',
    "Enhanced Interactive Volcano Plot for Metabolomics Studies"
  )
), 
  
 theme = shinytheme("flatly"), 
  setBackgroundColor(color = c("#FFFFFF", "#FFFFFF", "#67CFAC61"), gradient = "linear", direction = "bottom"),

  tags$head(
    tags$style(HTML("
      .nav-tabs > li > a {
        font-size: 20px !important;
        font-weight: bold !important;
        padding: 12px 18px !important;
      }
      .nav-tabs > li.active > a {
        font-size: 22px !important;
      }
    "))
  ),

  tags$head(tags$style(HTML("
    .shiny-output-error-validation { color:#b00020 !important; font-size: 18px !important; font-weight:800 !important; padding:10px; }
    .highlight { background:#fff; border:2px solid #000; color:#000; padding:8px; font-size: 18px; border-radius:8px; font-weight:bold; }
    .nav-tabs>li>a { font-size: 18px; padding: 10px 14px; }
    .small-note { font-size: 13px; opacity: 0.85; }
  "))),

tags$head(tags$style(HTML("
    /* Existing Footer and layout styles */
    .app-footer { position: fixed; left:0; right:0; bottom:0; 
                  text-align:center; font-size:12px; opacity:0.75;
                  padding:8px; background: rgba(255,255,255,0.8);
                  border-top: 1px solid #ddd; z-index: 9999; }
    body { padding-bottom: 45px; }

    /* --- NEW: Thicker Upload Progress Bar --- */
    .progress.shiny-file-input-progress {
      height: 20px !important;
      margin-top: 10px !important;
      border-radius: 5px !important;
    }
    
    .progress.shiny-file-input-progress .progress-bar {
      line-height: 20px !important;
      font-size: 14px !important;
      font-weight: bold !important;
      background-color: #66CDAA !important; /* Matches your app's theme color */
    }
    /* ---------------------------------------- */

    .tooltip-inner {
      max-width: none !important;
      white-space: nowrap;
      text-align: left !important;
      font-size: 18px;
    }
  "))),

  tabsetPanel(id = "tabs",
    tabPanel("1) Load & Process", value = "load",
      sidebarLayout(
        sidebarPanel(
          h3(class = "highlight", "Upload"),
          selectInput(
            "software_tool",
            "Software tool:",
            choices = c("mzMine" = "mzmine", "xcms" = "xcms", "MS-DIAL" = "msdial", "Default" = "default"),
            selected = "mzmine"
          ),
          fileInput("file_data", "Upload feature table (.csv)", accept = ".csv"),
          helpText(HTML("<i class='fa fa-info-circle'></i> Need data to test? Download example datasets from our <a href='https://github.com/plyush1993/Metabocano' target='_blank'>GitHub</a>.")),
          uiOutput("col_pickers"),

          selectizeInput(
            "sample_keywords",
            "Sample column keywords (pick/add multiple):",
            choices  = c(".mzML", ".mzXML", ".raw", " Peak area", "Area"),
            selected = c(".mzML", ".mzXML"),
            multiple = TRUE,
            options  = list(create = TRUE, createOnBlur = TRUE,
                            placeholder = "Type to add keyword and press Enter")
          ),

          tags$hr(),
          h3(class = "highlight", "Labels"),
          helpText(HTML("<i class='fa fa-info-circle'></i> Avoid double underscores (`__`) in group labels")),

radioButtons(
  "label_source",
  "Label source:",
  choices = c(
    "From sample names (by token)" = "token",
    "From metadata CSV (match by sample name)" = "metadata",
    "From uploaded labels CSV (1 column, no header)" = "csv",
    "Manual editable table" = "manual"
  ),
  selected = "token"
),

conditionalPanel(
  condition = "input.label_source == 'token' || input.label_source == 'manual'",
  textInput("token_sep", "Token separator", value = "_"),
  numericInput("token_index", "Token index (1-based)", value = 2, min = 1, step = 1)
),

conditionalPanel(
  condition = "input.label_source == 'metadata'",

  fileInput(
    "file_metadata_labels",
    "Upload metadata CSV with column names",
    accept = ".csv"
  ),

  uiOutput("metadata_sample_col_ui"),
  uiOutput("metadata_label_col_ui"),

  checkboxInput(
    "metadata_clean_sample_names",
    "Clean sample names only for metadata matching",
    value = TRUE
  ),

  conditionalPanel(
    condition = "input.metadata_clean_sample_names == true",

    selectizeInput(
      "metadata_remove_suffixes",
      "Remove suffixes/extensions:",
      choices = c(
        ".mzML", ".mzXML", ".raw", ".RAW",
        ".cdf", ".CDF", ".mzData", ".mzdata",
        ".wiff", ".WIFF", ".d", ".D",
        " Peak area", " Peak Area",
        " Peak height", " Peak Height",
        "_Area", "_Height",
        " Area", " Height"
      ),
      selected = c(
        " Peak area", " Peak height",
        "_Area", "_Height",
        " Area", " Height"
      ),
      multiple = TRUE,
      options = list(
        create = TRUE,
        createOnBlur = TRUE,
        placeholder = "Type custom suffix and press Enter"
      )
    )
  ),

  div(
    class = "small-note",
    "Metadata rows are matched by sample name, not by row order. Cleaning is used for matching."
  )
),

conditionalPanel(
  condition = "input.label_source == 'csv'",
  fileInput("file_labels", "Upload labels CSV", accept = ".csv")
),

conditionalPanel(
  condition = "input.label_source == 'manual'",
  div(
    style = "display: inline-flex; align-items: center; gap: 6px; margin-bottom: 10px;",

    actionButton(
      "fill_manual_labels",
      label = tags$span(
        HTML("Fill editable table from<br>current token labels"),
        style = "line-height: 1.1;"
      ),
      class = "btn-success",
      style = "
        font-size: 12px;
        padding: 4px 8px;
        line-height: 1.1;
        width: 150px;
        white-space: normal;
      "
    )
  ),

  div(
    class = "small-note",
    "Double-click cells in the Label column to edit group names. The Sample column is locked."
  )
),

checkboxInput("show_labels_table", "Show labels table", TRUE),
tags$hr(),

          h3(class = "highlight", "Join with Annotation"),
      div(
        style = "display: flex; align-items: center; margin-bottom: 15px;",
        
        materialSwitch(
          inputId = "use_sirius", 
          label = "Join SIRIUS summary", 
          value = FALSE, 
          status = "success", 
          width = "auto"
        ),
        
        actionButton(
          inputId = "btn5", 
          label = "?", 
          class = "btn-xs", 
          style = "font-weight: bold; margin-left: 10px; margin-top: -20px;"
        )
      ),
      bsTooltip(
        id = "btn5", 
        title = "<b>Join SIRIUS csv output with processed table (by <em>Row ID</em> column)</b>", 
        placement = "right", 
        trigger = "click", 
        options = list(container = "body")
      ),

      conditionalPanel(
        condition = "input.use_sirius",
        fileInput("file_sirius", "Upload SIRIUS .csv output", accept = ".csv"),
        uiOutput("sirius_pickers")
      ),

          tags$hr(),
          h3(class = "highlight", "Imputation by Noise"),
          radioButtons("do_mvi", "Imputation:", c("No"="no", "Yes"="yes"), selected = "yes", inline = TRUE),
          conditionalPanel(
            condition = "input.do_mvi == 'yes'",
            conditionalPanel(
              condition = "input.do_mvi == 'yes'",
              radioButtons("noise_mode", "Noise:",
                           c("Quantile in the range 1:min"="quantile", "Manual value"="manual"),
                           selected = "quantile"),
              conditionalPanel(condition = "input.noise_mode == 'quantile'",
                               numericInput("noise_quantile", "Quantile (0-1)", value = 0.25, min = 0, max = 1, step = 0.01)),
              conditionalPanel(condition = "input.noise_mode == 'manual'",
                               numericInput("noise_manual", "Noise value", value = 50, min = 0, step = 1)),
              numericInput("noise_sd", "SD for random values", value = 30, min = 0, step = 1)
            )
          ),

          tags$hr(),
          h3(class = "highlight", "Statistics"),
          radioButtons("comparison_mode", "Comparison selection:",
            choices = c(
              "Reference group vs all others" = "reference",
              "Choose comparisons manually" = "manual"), selected = "reference"),
          uiOutput("comparison_picker"),
          selectInput(
            "test_type", "Test:",
            c("Student", "Wilcoxon", "limma (Moderated t-test)" = "limma"), selected = "Student"),
          selectInput("p_adjust", "p-adjust:", c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"), selected = "BH"),
          conditionalPanel(
            condition = "input.test_type == 'Student' || input.test_type == 'Wilcoxon'",
            checkboxInput("paired", "Paired test", FALSE)
          ),
          conditionalPanel(
            condition = "input.test_type == 'Student'",
            checkboxInput("eqvar", "Equal variances (Student t-test)", FALSE)
          ),
          materialSwitch(
          "log2_test",
          "Log Transformation",
          value = FALSE,
          status = "success"
          ),
          materialSwitch(
            "standard_scaling",
            "Auto Scaling",
            value = FALSE,
            status = "success"
          ),
          tags$hr(),
          actionButton("run_proc", "Run preprocessing", class = "btn btn-success"),
          tags$br(), tags$br(),
          downloadButton("dl_volcano", "Volcano table csv", class = "btn-info"),
          actionButton("btn1", "?"),
          bsTooltip("btn1", 
          title = "<b>Download table with all calculated statistical values. Can be merged <em>GNPS-derived .cys file</em> in Cytoscape by <em>id</em> column.</b>", "right", trigger = "click", options = list(container = "body")),

          tags$br(),tags$br(),
          downloadButton("dl_matrix", "MetaboAnalyst-ready csv", class = "btn-info"),
          actionButton("btn2", "?"),
          bsTooltip("btn2", 
          title = "<b>Download peak table after MVI and with <em>Label</em> column.</b><br>Suitable as input in MetaboAnalyst (www.metaboanalyst.ca/).", "right", trigger = "click", options = list(container = "body")),
        
        tags$br(),tags$br(),
        downloadButton(
          "dl_autoplotter_zip",
          "AutoPlotter-ready ZIP",
          class = "btn-info"
        ),
        actionButton("btn_auto", "?"),
        bsTooltip(
          "btn_auto",
          title = "<b>Download AutoPlotter-ready ZIP archive.</b><br>Suitable as input as <em>Compounds in Columns</em> in Metabolite AutoPlotter (https://mpietzke.shinyapps.io/AutoPlotter/).",
          placement = "right",
          trigger = "click",
          options = list(container = "body"))
        ),

        mainPanel(
          uiOutput("raw_header"),
          DTOutput("raw_preview"),
          tags$hr(),
          uiOutput("labels_header"),
          uiOutput("label_upload_warning"),
          conditionalPanel(
  condition = "input.show_labels_table || input.label_source == 'manual'",
  DTOutput("labels_table")
),
          tags$hr(),
          uiOutput("proc_summary")
        )
      )
    ),

    tabPanel("2) Volcano explorer", value = "volcano",
      sidebarLayout(
        sidebarPanel(uiOutput("volcano_sidebar")),
        mainPanel(uiOutput("volcano_main"))
      )
    ),

tabPanel("App) SIRIUS & GNPS stats", value = "sirius_gnps",
  sidebarLayout(
    sidebarPanel(uiOutput("sirius_gnps_sidebar")),
    mainPanel(uiOutput("sirius_gnps_main"))
  )
)

  )
)

# ----------------------------- Server -------------------------------------

server <- function(input, output, session) {

  output$label_upload_warning <- renderUI({
  src <- input$label_source %||% "token"

  need_file <- FALSE

  if (identical(src, "csv") && is.null(input$file_labels)) {
    need_file <- TRUE
  }

  if (identical(src, "metadata") && is.null(input$file_metadata_labels)) {
    need_file <- TRUE
  }

  if (!need_file) return(NULL)

  div(
    class = "alert alert-warning",
    style = "
      margin-top: 10px;
      margin-bottom: 12px;
      font-size: 16px;
      font-weight: 700;
      border: 2px solid #f0ad4e;
      border-radius: 8px;
    ",
    "Please upload a CSV file first!"
  )
})
  
  session$onFlushed(function() {
    shinyjs::disable("run_proc")
    shinyjs::disable("dl_volcano")
    shinyjs::disable("dl_matrix")
    shinyjs::disable("dl_autoplotter_zip")
  }, once = TRUE)

  observe({
    shinyjs::toggleState("run_proc", condition = !is.null(input$file_data))
  })

 observe({
  ids <- c(
    "dl_volcano",
    "dl_matrix",
    "dl_autoplotter_zip"
  )

  if (procReady()) {
    lapply(ids, shinyjs::enable)
  } else {
    lapply(ids, shinyjs::disable)
  }
})

  observeEvent(input$file_data, {
    rv$raw <- NULL
    rv$mat <- NULL
    rv$fmap <- NULL
    rv$labels <- NULL
    rv$df_used <- NULL
    rv$volcano <- NULL
  }, ignoreInit = TRUE)

  observeEvent(
  list(
    input$label_source,
    input$file_labels,
    input$file_metadata_labels,
    input$metadata_sample_col,
    input$metadata_label_col,
    input$metadata_clean_sample_names,
    input$metadata_remove_suffixes,
    input$token_sep,
    input$token_index
  ),
  {
    rv$labels <- NULL
    rv$df_used <- NULL
    rv$volcano <- NULL
  },
  ignoreInit = TRUE
)
  
  rv <- reactiveValues(
    raw = NULL,
    mat = NULL,
    fmap = NULL,
    labels = NULL,
    df_used = NULL,
    volcano = NULL
  )

  procReady <- reactive({
    !is.null(rv$volcano) && nrow(rv$volcano) > 0
  })
  
  dataset_name <- reactive({
  nm <- input$file_data$name %||% "dataset.csv"
  tools::file_path_sans_ext(basename(nm))
})

  # ---- Load raw data (Robust Switch) ----
  raw_df <- reactive({
    req(input$file_data)
    ext <- tools::file_ext(input$file_data$name)
    validate(need(tolower(ext) == "csv", "Please upload a .csv file"))
    
    tool <- input$software_tool %||% "mzmine"
    
    if (tool == "msdial") {
      # Use custom robust reader for MS-DIAL
      df <- read_msdial_robust(input$file_data$datapath)
      clean_mzmine_export(df) # cleanup empty trailing cols
    } else {
      # Use fast vroom reader for standard CSVs (xcms/mzmine/default)
      clean_mzmine_export(vroom::vroom(input$file_data$datapath, delim = ","))
    }
  }) %>% bindCache(input$file_data$name, input$software_tool)

  # ---- Column Pickers (Smart Defaults) ----
  output$col_pickers <- renderUI({
    req(raw_df())
    cols <- names(raw_df())
    tool <- input$software_tool %||% "mzmine"
    
    rid_cand <- switch(tool,
                       xcms   = c("...1", "X...1", "feature_id", "feature", "id", "row id"),
                       msdial = c("alignment id", "alignmentid", "spot id"),
                       default = c("Feature", "row id", "id", "feature_id"),
                       c("row id", "id", "feature_id")
    )
    mz_cand <- switch(tool,
                      xcms   = c("mzmed", "mz", "m/z", "mzmin", "mzmax"),
                      msdial = c("average mz", "averagemz", "mz"),
                      default = c("mz", "m/z", "mass", "average mz"),
                      c("row m/z", "row mz", "mz")
    )
    rt_cand <- switch(tool,
                      xcms   = c("rtmed", "rt", "rtmin"),
                      msdial = c("average rt(min)", "average rt", "averagertmin", "rt"),
                      default = c("rt", "retention time", "time"),
                      c("row retention time", "row rt", "rt")
    )
    
    def_rid <- guess_col(cols, rid_cand)
    def_mz  <- guess_col(cols, mz_cand)
    def_rt  <- guess_col(cols, rt_cand)
    
    choices_rid <- c("<Auto-generate>", cols)
    sel_rid <- if (is.null(def_rid)) "<Auto-generate>" else def_rid
    
    tagList(
      selectInput("row_id_col", "Row ID column:", choices = choices_rid, selected = sel_rid),
      selectInput("mz_col",     "m/z column:",   choices = cols, selected = def_mz  %||% cols[1]),
      selectInput("rt_col",     "rt column:",    choices = cols, selected = def_rt  %||% cols[1])
    )
  })

  output$raw_header <- renderUI({
    req(raw_df())
    h3(sprintf("Raw dataset: %d rows × %d columns", nrow(raw_df()), ncol(raw_df())))
  })

  output$raw_preview <- renderDT({
    req(raw_df())
    datatable(head(raw_df(), 20), options = list(scrollX = TRUE, pageLength = 8))
  })

  # ---- Build matrix + fmap ----
  built <- reactive({
    req(raw_df(), input$row_id_col, input$mz_col, input$rt_col)
    
    df <- raw_df()
    rid_col <- input$row_id_col
    
    if (rid_col == "<Auto-generate>") {
      df$FeatureID_Auto <- as.numeric(seq_len(nrow(df)))
      rid_col <- "FeatureID_Auto"
    }
    
    kws <- input$sample_keywords
    if (is.null(kws) || length(kws) == 0) kws <- c(".") 

    parse_feature_table_to_matrix(
      df, 
      row_id_col = rid_col,
      mz_col     = input$mz_col,
      rt_col     = input$rt_col,
      sample_keywords = kws
    )
  })

  sample_names <- reactive({
    req(built())
    rownames(built()$mat)
  })

  manual_labels <- reactiveVal(NULL)

auto_label_table <- reactive({
  req(sample_names())

  make_label_table(
    sample_names(),
    labels_from_sample_names_or_raw(
      sample_names(),
      token_sep = input$token_sep %||% "_",
      token_index = input$token_index %||% 2,
      clean_names = FALSE
    )
  )
})

observeEvent(sample_names(), {
  req(auto_label_table())
  manual_labels(auto_label_table())
}, ignoreInit = FALSE)

observeEvent(input$fill_manual_labels, {
  req(auto_label_table())

  manual_labels(auto_label_table())

  rv$labels <- NULL
  rv$df_used <- NULL
  rv$volcano <- NULL

  showNotification(
    "Editable label table was filled from current token labels.",
    type = "message",
    duration = 3
  )
}, ignoreInit = TRUE)

observeEvent(input$labels_table_cell_edit, {
  req(input$label_source == "manual")

  info <- input$labels_table_cell_edit

  tbl <- manual_labels()
  req(tbl)

  row_i <- as.integer(info$row)

  if (!is.finite(row_i) || row_i < 1 || row_i > nrow(tbl)) {
    showNotification("Edited row is outside label table.", type = "error", duration = 3)
    return(NULL)
  }

  tbl$Label[row_i] <- trimws(as.character(info$value))

  manual_labels(tbl)

  rv$labels <- NULL
  rv$df_used <- NULL
  rv$volcano <- NULL

  showNotification(
    paste0("Label updated: ", tbl$Sample[row_i], " -> ", tbl$Label[row_i]),
    type = "message",
    duration = 2
  )
}, ignoreInit = TRUE)
  
metadata_raw <- reactive({
  req(input$file_metadata_labels)
  read_metadata_csv(input$file_metadata_labels, context = "metadata labels")
})


output$metadata_sample_col_ui <- renderUI({
  req(metadata_raw())

  cols <- names(metadata_raw())

  selectInput(
    "metadata_sample_col",
    "Metadata sample-name column:",
    choices = cols,
    selected = guess_metadata_sample_col(cols)
  )
})


output$metadata_label_col_ui <- renderUI({
  req(metadata_raw())

  cols <- names(metadata_raw())
  sample_col <- input$metadata_sample_col %||% guess_metadata_sample_col(cols)
  choices <- setdiff(cols, sample_col)

  validate(
    need(length(choices) > 0, "Metadata file has no column available for labels.")
  )

  selectizeInput(
    "metadata_label_col",
    "Metadata column to use as Label:",
    choices = choices,
    selected = guess_metadata_label_col(cols, sample_col),
    multiple = FALSE
  )
})


metadata_labels <- reactive({
  req(
    input$file_metadata_labels,
    input$metadata_sample_col,
    input$metadata_label_col,
    sample_names()
  )

  metadata_labels_by_sample(
    upload = input$file_metadata_labels,
    sample_names = sample_names(),
    sample_col = input$metadata_sample_col,
    label_col = input$metadata_label_col,
    clean_enabled = isTRUE(input$metadata_clean_sample_names),
    remove_suffixes = input$metadata_remove_suffixes %||% character(0),
    context = "metadata labels"
  )
})

  # ---- Labels ----
  labels_vec <- reactive({
  req(sample_names())

  src <- input$label_source %||% "token"

  if (identical(src, "csv")) {

  req(input$file_labels)

  v <- read_onecol_csv(input$file_labels$datapath)
  v <- trimws(v)

  validate(
    need(
      length(v) == length(sample_names()),
      sprintf(
        "Labels count (%d) must match #samples (%d).",
        length(v), length(sample_names())
      )
    )
  )

  v

} else if (identical(src, "metadata")) {

  metadata_labels()

} else if (identical(src, "manual")) {

    tbl <- manual_labels()
    req(tbl)

    validate(
      need(
        nrow(tbl) == length(sample_names()),
        "Manual label table must match the number of samples."
      ),
      need(
        identical(as.character(tbl$Sample), as.character(sample_names())),
        "Manual label table does not match current sample names. Click 'Fill editable table from current token labels'."
      ),
      need(
        !any(is.na(tbl$Label) | trimws(tbl$Label) == ""),
        "All samples must have labels."
      )
    )

    trimws(as.character(tbl$Label))

  } else {

    sn <- sample_names()
    sn2 <- sn
    sep <- input$token_sep %||% "_"
    idx <- as.integer(input$token_index %||% 2)

    parts <- strsplit(sn2, sep, fixed = TRUE)
    ok <- vapply(parts, function(z) length(z) >= idx, logical(1))

    validate(
      need(
        all(ok),
        sprintf("Token %d missing in some sample names. Adjust separator/index.", idx)
      )
    )

    labs <- vapply(parts, function(z) z[[idx]], character(1))

    validate(
      need(all(nzchar(labs)), "Parsed empty labels — adjust separator/index.")
    )

    labs
  }
})

  output$labels_header <- renderUI({
    req(sample_names(), labels_vec())
    h3(sprintf("Labels ready: %d samples", length(labels_vec())))
  })

  output$labels_table <- renderDT({
  req(sample_names())

  src <- input$label_source %||% "token"

  tbl <- if (identical(src, "manual")) {
    manual_labels()
  } else {
    req(labels_vec())
    make_label_table(sample_names(), labels_vec())
  }

  req(tbl)

  datatable(
    tbl,
    editable = if (identical(src, "manual")) {
      list(
        target = "cell",
        disable = list(columns = c(0)) # lock Sample column
      )
    } else {
      FALSE
    },
    options = list(
      pageLength = 8,
      scrollX = TRUE,
      ordering = FALSE,
      searching = FALSE
    ),
    rownames = FALSE
  )
}, server = FALSE)
  
  comparison_pairs <- reactive({
  req(labels_vec())

  levs <- sort(
    unique(
      trimws(as.character(labels_vec()))
    )
  )

  levs <- levs[nzchar(levs)]

  validate(
    need(
      length(levs) >= 2,
      "Need at least 2 groups to define comparisons."
    )
  )

  # All directional comparisons:
  # A / B and B / A are separate options
  tidyr::expand_grid(
    Group_num = levs,
    Group_den = levs
  ) %>%
    dplyr::filter(Group_num != Group_den) %>%
    dplyr::mutate(
      Comparison_ID = sprintf(
        "comparison_%04d",
        dplyr::row_number()
      ),
      Comparison = paste0(
        Group_num,
        " / ",
        Group_den
      )
    ) %>%
    dplyr::select(
      Comparison_ID,
      Comparison,
      Group_num,
      Group_den
    )
})


output$comparison_picker <- renderUI({
  req(labels_vec())

  levs <- sort(
    unique(
      trimws(as.character(labels_vec()))
    )
  )

  pairs <- comparison_pairs()

  # Preserve the currently selected reference group
  old_ref <- isolate(input$ref_group)

  if (is.null(old_ref) || !old_ref %in% levs) {
    old_ref <- levs[1]
  }

  # Preserve valid manual choices when UI is rebuilt
  old_manual <- isolate(
    input$manual_comparisons
  ) %||% character(0)

  old_manual <- intersect(
    old_manual,
    pairs$Comparison_ID
  )

  tagList(

    conditionalPanel(
      condition = "input.comparison_mode == 'reference'",

      selectInput(
        "ref_group",
        "Reference Group:",
        choices = levs,
        selected = old_ref
      )
    ),

    conditionalPanel(
      condition = "input.comparison_mode == 'manual'",

      pickerInput(
        inputId = "manual_comparisons",
        label = "Select comparisons:",
        choices = stats::setNames(
          pairs$Comparison_ID,
          pairs$Comparison
        ),
        selected = old_manual,
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE,
          `none-selected-text` =
            "Select one or more comparisons",
          `selected-text-format` = "count > 2",
          `count-selected-text` =
            "{0} comparison(s) selected",
          `style` = "btn-success"
        )
      )
    )
  )
})

selected_manual_comparisons <- reactive({
  req(
    identical(
      input$comparison_mode %||% "reference",
      "manual"
    )
  )

  selected_ids <- input$manual_comparisons %||%
    character(0)

  validate(
    need(
      length(selected_ids) > 0,
      "Select at least one manual comparison."
    )
  )

  pairs <- comparison_pairs()

  selected_rows <- match(
    selected_ids,
    pairs$Comparison_ID
  )

  validate(
    need(
      !any(is.na(selected_rows)),
      paste0(
        "One or more selected comparisons ",
        "are no longer available."
      )
    )
  )

  pairs[
    selected_rows,
    c("Group_num", "Group_den"),
    drop = FALSE
  ]
})

  # ---- SIRIUS pickers ----
  sirius_df <- reactive({
    req(input$use_sirius)
    req(input$file_sirius)
    ext <- tools::file_ext(input$file_sirius$name)
    validate(need(tolower(ext) == "csv", "SIRIUS file must be .csv"))
    vroom::vroom(input$file_sirius$datapath, delim = ",") %>% as.data.frame(check.names = FALSE)
  })

  output$sirius_pickers <- renderUI({
    req(sirius_df())
    cols <- names(sirius_df())
    tagList(
      selectInput("sirius_idcol", "SIRIUS Feature ID column:",
                  choices = cols, selected = if ("mappingFeatureId" %in% cols) "mappingFeatureId" else cols[1]),
      selectInput("sirius_npcol", "NPC column:",
                  choices = cols, selected = if ("NPC#class" %in% cols) "NPC#class" else cols[1]),
      selectInput("sirius_cfcol", "ClassyFire column:",
                  choices = cols, selected = if ("ClassyFire#class" %in% cols) "ClassyFire#class" else cols[1])
    )
  })
  
  # ---- SIRIUS & GNPS stats tab ----

guess_col_ci <- function(cols, candidates, default = NULL) {
  if (is.null(cols) || !length(cols)) return(default)

  cols_norm <- gsub("[^a-z0-9]", "", tolower(cols))
  cand_norm <- gsub("[^a-z0-9]", "", tolower(candidates))

  for (cand in cand_norm) {
    hit <- which(cols_norm == cand)
    if (length(hit)) return(cols[hit[1]])
  }

  for (cand in cand_norm) {
    hit <- which(grepl(cand, cols_norm, fixed = TRUE))
    if (length(hit)) return(cols[hit[1]])
  }

  if (!is.null(default)) default else cols[1]
}

clean_stats_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | x %in% c("NA", "NaN", "null", "Not provided")] <- NA_character_
  x
}

gnps_pairs_df <- reactive({
  req(input$file_gnps_pairs)

  ext <- tolower(tools::file_ext(input$file_gnps_pairs$name))

  validate(
    need(
      ext %in% c("tsv", "txt", "csv"),
      "Upload GNPS pairs as .tsv, .txt, or .csv."
    )
  )

  delim <- if (ext == "csv") "," else "\t"

  as.data.frame(
    vroom::vroom(
      input$file_gnps_pairs$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})

output$sirius_gnps_sidebar <- renderUI({
  if (!isTRUE(input$use_sirius) || is.null(input$file_sirius)) {
    return(
      div(
        class = "highlight",
        "Switch on 'Join SIRIUS summary' in the Load & Process tab and upload the SIRIUS .csv file."
      )
    )
  }

  s <- sirius_df()
  s_cols <- names(s)

  default_id <- if (!is.null(input$sirius_idcol) && input$sirius_idcol %in% s_cols) {
    input$sirius_idcol
  } else {
    guess_col_ci(
      s_cols,
      c("mappingFeatureId", "id", "featureId", "row ID"),
      s_cols[1]
    )
  }

  default_ann <- guess_col_ci(
    s_cols,
    c("NPC#class", "ClassyFire#class", "NPC class", "ClassyFire class", "class", "superclass"),
    s_cols[1]
  )

  tagList(
    h3(class = "highlight", "SIRIUS annotation frequency"),

    selectInput(
      "stats_sirius_id_col",
      "SIRIUS ID column:",
      choices = s_cols,
      selected = default_id
    ),

    selectInput(
      "stats_sirius_col",
      "SIRIUS column for frequency statistics:",
      choices = s_cols,
      selected = default_ann
    ),

    uiOutput("stats_class_picker"),

    tags$hr(),

    h3(class = "highlight", "Optional GNPS network pairs"),

    fileInput(
      "file_gnps_pairs",
      "Upload GNPS Pairs List file (.tsv/.txt/.csv)",
      accept = c(".tsv", ".txt", ".csv")
    ),

    uiOutput("gnps_pickers_ui"),

    div(
      class = "small-note",
      "GNPS pairs are converted from ClusterID1/ClusterID2 into one ClusterID column. By default, no peak-table pruning is applied. If a peak-table ID column is selected, SIRIUS IDs and GNPS ClusterIDs are restricted to IDs present in the uploaded peak table."
    )
  )
})

output$gnps_pickers_ui <- renderUI({
  req(gnps_pairs_df())

  g_cols <- names(gnps_pairs_df())

  raw_cols <- tryCatch(names(raw_df()), error = function(e) character(0))

  default_peak <- if (!is.null(input$row_id_col) && input$row_id_col %in% raw_cols) {
    input$row_id_col
  } else {
    guess_col_ci(
      raw_cols,
      c("row ID", "row id", "id", "feature_id"),
      default = "__none__"
    )
  }

  tagList(
    selectInput(
      "gnps_cluster1_col",
      "GNPS ClusterID1 column:",
      choices = g_cols,
      selected = guess_col_ci(g_cols, c("ClusterID1", "CLUSTERID1"), g_cols[1])
    ),

    selectInput(
      "gnps_cluster2_col",
      "GNPS ClusterID2 column:",
      choices = g_cols,
      selected = guess_col_ci(g_cols, c("ClusterID2", "CLUSTERID2"), g_cols[min(2, length(g_cols))])
    ),

    selectInput(
      "gnps_component_col",
      "GNPS ComponentIndex column:",
      choices = g_cols,
      selected = guess_col_ci(g_cols, c("ComponentIndex", "component"), g_cols[1])
    ),

    selectInput(
      "gnps_peak_id_col",
      "Optional peak-table pruning by Peak ID:",
      choices = c("Do not use peak-table ID matching" = "__none__", raw_cols),
      selected = "__none__"
    )
  )
})

sirius_stats_data <- reactive({
  req(sirius_df(), input$stats_sirius_id_col, input$stats_sirius_col)

  s <- as.data.frame(
    sirius_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(input$stats_sirius_id_col %in% names(s), "Selected SIRIUS ID column was not found."),
    need(input$stats_sirius_col %in% names(s), "Selected SIRIUS statistics column was not found.")
  )

  out <- tibble::tibble(
    SIRIUS_ID = trimws(as.character(s[[input$stats_sirius_id_col]])),
    Annotation = clean_stats_value(s[[input$stats_sirius_col]])
  ) %>%
    dplyr::filter(nzchar(SIRIUS_ID), !is.na(Annotation))

  keep_ids <- selected_peak_ids_for_gnps()

  if (!is.null(keep_ids) && length(keep_ids)) {
    out <- out %>%
      dplyr::filter(SIRIUS_ID %in% keep_ids)
  }

  out
})

output$stats_class_picker <- renderUI({
  req(sirius_stats_data())

  vals <- sirius_stats_data() %>%
    dplyr::count(Annotation, sort = TRUE, name = "Frequency") %>%
    dplyr::pull(Annotation)

  if (!length(vals)) {
    return(
      div(
        class = "small-note",
        "No non-empty values detected for the selected SIRIUS column."
      )
    )
  }

  pickerInput(
    "stats_selected_class",
    "Specific class/value for GNPS ComponentIndex statistics:",
    choices = vals,
    selected = vals[1],
    multiple = FALSE,
    options = list(
      `live-search` = TRUE,
      `style` = "btn-success"
    )
  )
})

sirius_frequency_table <- reactive({
  req(sirius_stats_data())

  n_total <- nrow(sirius_stats_data())

  sirius_stats_data() %>%
    dplyr::group_by(Annotation) %>%
    dplyr::summarise(
      Frequency = dplyr::n(),
      Unique_IDs = dplyr::n_distinct(SIRIUS_ID),
      Percent = round(100 * Frequency / n_total, 2),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Frequency), Annotation)
})

output$sirius_frequency_table <- DT::renderDT({
  datatable(
    sirius_frequency_table(),
    rownames = FALSE,
    class = "compact stripe hover nowrap",
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      order = list(list(1, "desc"))
    )
  )
})

selected_peak_ids_for_gnps <- reactive({
  if (is.null(input$gnps_peak_id_col) || identical(input$gnps_peak_id_col, "__none__")) {
    return(NULL)
  }

  req(raw_df())

  raw <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(input$gnps_peak_id_col %in% names(raw), "Selected peak-table ID column was not found.")
  )

  ids <- trimws(as.character(raw[[input$gnps_peak_id_col]]))
  ids[nzchar(ids)]
})

gnps_component_map <- reactive({
  req(
    gnps_pairs_df(),
    input$gnps_cluster1_col,
    input$gnps_cluster2_col,
    input$gnps_component_col
  )

  g <- as.data.frame(
    gnps_pairs_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(input$gnps_cluster1_col %in% names(g), "GNPS ClusterID1 column was not found."),
    need(input$gnps_cluster2_col %in% names(g), "GNPS ClusterID2 column was not found."),
    need(input$gnps_component_col %in% names(g), "GNPS ComponentIndex column was not found.")
  )

  out <- g %>%
    dplyr::transmute(
      ComponentIndex = as.character(.data[[input$gnps_component_col]]),
      ClusterID1 = as.character(.data[[input$gnps_cluster1_col]]),
      ClusterID2 = as.character(.data[[input$gnps_cluster2_col]])
    ) %>%
    tidyr::pivot_longer(
      cols = c("ClusterID1", "ClusterID2"),
      names_to = "Cluster_side",
      values_to = "ClusterID"
    ) %>%
    dplyr::mutate(
      ClusterID = trimws(as.character(ClusterID)),
      ComponentIndex = trimws(as.character(ComponentIndex))
    ) %>%
    dplyr::filter(nzchar(ClusterID), nzchar(ComponentIndex)) %>%
    dplyr::distinct(ComponentIndex, ClusterID)

  keep_ids <- selected_peak_ids_for_gnps()

  if (!is.null(keep_ids) && length(keep_ids)) {
    out <- out %>%
      dplyr::filter(ClusterID %in% keep_ids)
  }

  out
})

selected_class_component_stats <- reactive({
  req(sirius_stats_data(), input$stats_selected_class)

  ids <- sirius_stats_data() %>%
    dplyr::filter(Annotation == input$stats_selected_class) %>%
    dplyr::distinct(SIRIUS_ID)

  if (!nrow(ids)) {
    return(
      tibble::tibble(
        ComponentIndex = character(),
        Points = integer(),
        ClusterIDs = character()
      )
    )
  }

  if (is.null(input$file_gnps_pairs)) {
    return(
      tibble::tibble(
        ComponentIndex = "GNPS pairs not uploaded",
        Points = dplyr::n_distinct(ids$SIRIUS_ID),
        ClusterIDs = paste(sort(unique(ids$SIRIUS_ID)), collapse = ", ")
      )
    )
  }

  comp <- gnps_component_map()

  ids %>%
    dplyr::left_join(comp, by = c("SIRIUS_ID" = "ClusterID")) %>%
    dplyr::mutate(
      ComponentIndex = dplyr::if_else(
        is.na(ComponentIndex) | !nzchar(ComponentIndex),
        "No ComponentIndex match",
        ComponentIndex
      )
    ) %>%
    dplyr::group_by(ComponentIndex) %>%
    dplyr::summarise(
      Points = dplyr::n_distinct(SIRIUS_ID),
      ClusterIDs = paste(sort(unique(SIRIUS_ID)), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Points), ComponentIndex)
})

output$selected_class_summary <- renderUI({
  req(sirius_stats_data(), input$stats_selected_class)

  ids_all <- sirius_stats_data() %>%
    dplyr::filter(Annotation == input$stats_selected_class) %>%
    dplyr::distinct(SIRIUS_ID)

  comp_stats <- selected_class_component_stats()

  mapped <- if (!is.null(input$file_gnps_pairs)) {
    sum(comp_stats$Points[comp_stats$ComponentIndex != "No ComponentIndex match"], na.rm = TRUE)
  } else {
    NA_integer_
  }

  div(
    style = "
      background:#ffffffcc;
      padding:10px;
      border-radius:10px;
      border:1px solid #ddd;
      margin-bottom:10px;
    ",
    h4("Selected class summary"),
    tags$ul(
      tags$li(strong("Column: "), input$stats_sirius_col),
      tags$li(strong("Selected value: "), input$stats_selected_class),
      tags$li(strong("Unique SIRIUS IDs / points: "), dplyr::n_distinct(ids_all$SIRIUS_ID)),
      tags$li(
        strong("SIRIUS entries: "),
        nrow(
          sirius_stats_data() %>%
            dplyr::filter(Annotation == input$stats_selected_class)
        )
      ),
      if (!is.na(mapped)) {
        tags$li(strong("Points with ComponentIndex match: "), mapped)
      },
      if (!is.null(input$file_gnps_pairs)) {
        tags$li(strong("Number of ComponentIndex groups shown: "), nrow(comp_stats))
      }
    )
  )
})

output$selected_class_component_table <- DT::renderDT({
  datatable(
    selected_class_component_stats(),
    rownames = FALSE,
    class = "compact stripe hover nowrap",
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      order = list(list(1, "desc"))
    )
  )
})

output$sirius_gnps_main <- renderUI({
  if (!isTRUE(input$use_sirius) || is.null(input$file_sirius)) {
    return(
      div(
        class = "highlight",
        "No SIRIUS annotation table uploaded yet. Go to Load & Process -> Join with Annotation -> upload SIRIUS summary."
      )
    )
  }

  tagList(
    h3("Frequency by selected SIRIUS column"),
    withSpinner(DTOutput("sirius_frequency_table"), type = 8, color = "#66CDAA"),
    tags$hr(),
    uiOutput("selected_class_summary"),
    h3("Selected class distribution by GNPS ComponentIndex"),
    withSpinner(DTOutput("selected_class_component_table"), type = 8, color = "#66CDAA")
  )
})

  # ---- Process button ----
  observeEvent(input$run_proc, {
  req(built(), labels_vec())
  labs_pre <- labels_vec()
  if (stop_if_one_group(labs_pre)) {
    return(NULL)
  }
  comparison_mode <- input$comparison_mode %||% "reference"
  manual_pairs <- NULL
  if (identical(comparison_mode, "manual")) {
    if (
      is.null(input$manual_comparisons) ||
      length(input$manual_comparisons) == 0
    ) {
      showNotification(
        "Select at least one manual comparison.",
        type = "error",
        duration = 6
      )
      return(NULL)
    }
    manual_pairs <- selected_manual_comparisons()
  } else {
    req(input$ref_group)
  }
  withProgress(message = "Processing...", value = 0, {
      incProgress(0.15, detail = "Building matrix")
      X <- built()$mat
      fmap <- built()$fmap

      incProgress(0.20, detail = "Adding labels")
      labs <- labels_vec()
      validate(need(length(labs) == nrow(X),
                    sprintf("Labels length (%d) must match #samples (%d).",
                            length(labs), nrow(X))))

      df_raw <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
      df_used <- cbind(Label = labs, df_raw)
      df_used$Label <- as.factor(df_used$Label)

      incProgress(0.25, detail = "Imputation (if enabled)")
      if (identical(input$do_mvi, "yes")) {
          X0 <- df_used[, -1, drop = FALSE]
          Xm <- impute_lod_random(
            X0,
            noise_mode     = input$noise_mode %||% "quantile",
            noise_quantile = input$noise_quantile %||% 0.25,
            noise_manual   = input$noise_manual %||% 50,
            sd_val         = input$noise_sd %||% 30,
            seed           = 1234
          )
        df_used <- as.data.frame(cbind(Label = df_used$Label, as.data.frame(Xm, check.names = FALSE)),
                                 check.names = FALSE, stringsAsFactors = FALSE)
        df_used$Label <- as.factor(df_used$Label)
      }

      incProgress(0.30, detail = "Running...")
      volc <- compute_stats_long(
        df_used,
        test = input$test_type %||% "Student",
        adj = input$p_adjust %||% "BH",
        paired = isTRUE(input$paired),
        eqvar = isTRUE(input$eqvar),
        pseudocount = 1.1,
        log2_test = isTRUE(input$log2_test),
        scale_data = isTRUE(input$standard_scaling),
      
        ref_group = if (
          identical(comparison_mode, "reference")
        ) {
          input$ref_group
        } else {
          NULL
        },
      
        comparisons = manual_pairs
      )

      incProgress(0.05, detail = "Joining mz/rt/id")
      volc <- volc %>% left_join(fmap, by = "Feature")

      # default annotation columns
      volc$`NPC#class` <- "Not provided"
      volc$`ClassyFire#class` <- "Not provided"

      incProgress(0.05, detail = "Joining SIRIUS (optional)")
      if (isTRUE(input$use_sirius) && !is.null(input$file_sirius)) {
        s <- sirius_df()
        req(input$sirius_idcol, input$sirius_npcol, input$sirius_cfcol)

        ss <- s %>%
          transmute(
            id = as.character(.data[[input$sirius_idcol]]),
            `NPC#class` = as.character(.data[[input$sirius_npcol]]),
            `ClassyFire#class` = as.character(.data[[input$sirius_cfcol]])
          )

        volc <- volc %>%
          mutate(id = as.character(id)) %>%
          left_join(ss, by = "id", suffix = c("", ".sirius")) %>%
          mutate(
            `NPC#class` = dplyr::coalesce(.data[["NPC#class.sirius"]], .data[["NPC#class"]]),
            `ClassyFire#class` = dplyr::coalesce(.data[["ClassyFire#class.sirius"]], .data[["ClassyFire#class"]])
          ) %>%
          select(-any_of(c("NPC#class.sirius", "ClassyFire#class.sirius")))
      }

      volc <- volc %>%
        mutate(
          mz = round(mz, 6),
          RT = round(RT, 4),
          Mean = signif(Mean, 4),
          FC = round(FC, 4),
          `Adj.p-value` = as.numeric(`Adj.p-value`),
          `Adj.p-value.log` = round(`Adj.p-value.log`, 5)
        )

      rv$raw <- built()$raw
      rv$mat <- X
      rv$fmap <- fmap
      rv$labels <- labs
      rv$df_used <- df_used
      rv$volcano <- volc
    })

    showNotification("Processing finished. Switching to Volcano explorer...", type = "message", duration = 4)
    updateTabsetPanel(session, "tabs", selected = "volcano")
  }, ignoreInit = TRUE)

  output$proc_summary <- renderUI({
    if (!procReady()) {
      div(class = "highlight", "Not processed yet. Upload -> Set settings -> click 'Run preprocessing'.")
    } else {
      div(
        style="background:#ffffffcc; padding:10px; border-radius:10px; font-size:20px; border:1px solid #ddd;",
        h4("Processing summary"),
        tags$ul(
          tags$li(sprintf("Samples: %d", nrow(rv$df_used))),
          tags$li(sprintf("Features: %d", ncol(rv$df_used) - 1)),
          tags$li(sprintf("Comparisons: %d", length(unique(rv$volcano$Groups)))),
          tags$li(sprintf("Rows in volcano table: %d", nrow(rv$volcano)))
        ),
        div(class="small-note", "")
      )
    }
  })

  # ---------------- Volcano tab UI ----------------

      observeEvent(input$reset_filters, {
      req(procReady())
    
      # switches
      updateMaterialSwitch(session, "sig_only", value = FALSE)
      updateMaterialSwitch(session, "use_fdr_filter", value = FALSE)
      updateMaterialSwitch(session, "use_fc_filter",  value = FALSE)
      updateMaterialSwitch(session, "use_npc_filter", value = FALSE)
      updateMaterialSwitch(session, "use_classyfire_filter", value = FALSE)
      
      # pickers / radios
      updatePickerInput(session, "sel_feat", selected = character(0))
      updatePickerInput(session, "npc_filter_values", selected = character(0))
      updatePickerInput(session, "classyfire_filter_values", selected = character(0))
      updateRadioButtons(session, "color_by", selected = "Groups")
      updateRadioButtons(session, "present_as", selected = "Boxplot")
      updateRadioButtons(session, "fc_dir", selected = "both")
    
      # sliders
      dd <- rv$volcano
      mzr <- finite_range(dd$mz)
      rtr <- finite_range(dd$RT)
      mr  <- finite_range(log10(dd$Mean + 1.1))
    
      if (!is.null(mzr)) {
        updateSliderInput(
          session, "mz_range",
          value = c(
            floor(mzr[1] * 10000) / 10000,
            ceiling(mzr[2] * 10000) / 10000
          )
        )
      }
      
      if (!is.null(rtr)) {
        updateSliderInput(
          session, "rt_range",
          value = c(
            floor(rtr[1] * 1000) / 1000,
            ceiling(rtr[2] * 1000) / 1000
          )
        )
      }
      
      if (!is.null(mr)) {
        updateSliderInput(
          session, "intensity_range",
          value = c(
            floor(mr[1] * 1000) / 1000,
            ceiling(mr[2] * 1000) / 1000
          )
        )
      }
    
      updateNumericInput(session, "sig_p_cutoff", value = 0.05)
      updateSliderInput(session, "fc_thr", value = 1)
      updateSelectInput(session, "volcano_palette", selected = "Set1")
      updateSelectInput(session, "box_palette", selected = "Dark2")
    })
  
  annotation_filter_choices <- reactive({
  req(procReady())

  dd <- rv$volcano

  clean_choices <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x) & x != "Not provided"]
    sort(unique(x))
  }

  list(
    npc = clean_choices(dd$`NPC#class`),
    classyfire = clean_choices(dd$`ClassyFire#class`)
  )
})
  
  output$volcano_sidebar <- renderUI({
    if (!procReady()) {
      return(div(class="highlight",
                 "No processed dataset yet. Go to '1) Load & Process' and click 'Run preprocessing'."))
    }

    tagList(
      actionButton("reset_filters", "Reset filters", class = "btn btn-warning"),
      br(),
      br(),
      materialSwitch("sig_only", "Significant only", value = FALSE, status = "success"),

tags$hr(),
h4(class = "highlight", "Annotation filters"),

materialSwitch(
  "use_npc_filter",
  "Filter by NPC class",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_npc_filter == true",
  if (length(annotation_filter_choices()$npc) > 0) {
    pickerInput(
      inputId = "npc_filter_values",
      label = "NPC class(es):",
      choices = annotation_filter_choices()$npc,
      selected = character(0),
      multiple = TRUE,
      options = list(
  `actions-box` = TRUE,
  `live-search` = TRUE,
  `none-selected-text` = "Select NPC class(es)",
  `style` = "btn-success",
  `selected-text-format` = "count > 1",
  `count-selected-text` = "{0} NPC class(es) selected"
)
    )
  } else {
    div(
      class = "small-note",
      "No NPC classes detected. Upload SIRIUS annotation and rerun preprocessing."
    )
  }
),

materialSwitch(
  "use_classyfire_filter",
  "Filter by ClassyFire class",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_classyfire_filter == true",
  if (length(annotation_filter_choices()$classyfire) > 0) {
    pickerInput(
      inputId = "classyfire_filter_values",
      label = "ClassyFire class(es):",
      choices = annotation_filter_choices()$classyfire,
      selected = character(0),
      multiple = TRUE,
      options = list(
  `actions-box` = TRUE,
  `live-search` = TRUE,
  `none-selected-text` = "Select ClassyFire class(es)",
  `style` = "btn-success",
  `selected-text-format` = "count > 1",
  `count-selected-text` = "{0} ClassyFire class(es) selected"
)
    )
  } else {
    div(
      class = "small-note",
      "No ClassyFire classes detected. Upload SIRIUS annotation and rerun preprocessing."
    )
  }
),

tags$hr(),

pickerInput(
  inputId = "sel_feat",
  label   = "Select/deselect features (optional):",
  choices = sort(unique(rv$volcano$Feature)),
  options = list(
  `actions-box` = TRUE,
  `live-search` = TRUE,
  `style` = "btn-success",
  `selected-text-format` = "count > 2",
  `count-selected-text` = "{0} feature(s) selected"
),
  multiple = TRUE
),
      br(),
      radioButtons(
        "color_by",
        "Color points by:",
        choices = c("Groups" = "Groups", "Mean" = "Mean"),
        selected = "Groups",
        inline = TRUE
      ),

      radioButtons(
        "present_as",
        "Click plot shows:",
        choices = c("Boxplot" = "Boxplot", "Scatterplot" = "Scatterplot"),
        selected = "Boxplot",
        inline = TRUE
      ),

selectInput(
  "volcano_palette",
  "Volcano palette:",
  choices = palette_choices,
  selected = "Set1"
),

selectInput(
  "box_palette",
  "Boxplot/scatter palette:",
  choices = palette_choices,
  selected = "Dark2"
),

      uiOutput("volcano_sliders")
    )
  })

  output$volcano_main <- renderUI({
    if (!procReady()) return(NULL)
    tagList(
      withSpinner(plotlyOutput("volcano_plot", height = "520px"), type = 8, color = "#66CDAA"),
      div(style = "height:8px;"),
      plotlyOutput("feature_plot", height = "260px")
    )
  })

  output$volcano_sliders <- renderUI({
  req(procReady())
  dd <- rv$volcano

  mzr <- finite_range(dd$mz)
  rtr <- finite_range(dd$RT)
  mr  <- finite_range(log10(dd$Mean + 1.1))
  fcmax <- max(abs(dd$FC), na.rm = TRUE)

  validate(need(!is.null(mzr) && !is.null(rtr) && !is.null(mr),
                "No finite mz/RT/Mean values available for sliders."))

  mz_min <- floor(mzr[1] * 10000) / 10000
  mz_max <- ceiling(mzr[2] * 10000) / 10000

  rt_min <- floor(rtr[1] * 1000) / 1000
  rt_max <- ceiling(rtr[2] * 1000) / 1000

  int_min <- floor(mr[1] * 1000) / 1000
  int_max <- ceiling(mr[2] * 1000) / 1000

  tagList(
    sliderInput("mz_range", "m/z:",
                min = mz_min,
                max = mz_max,
                value = c(mz_min, mz_max),
                step = 0.0001),

    sliderInput("rt_range", "RT:",
                min = rt_min,
                max = rt_max,
                value = c(rt_min, rt_max),
                step = 0.001),

    sliderInput("intensity_range", "Mean log10(Intensity):",
                min = int_min,
                max = int_max,
                value = c(int_min, int_max),
                step = 0.001),
      
          tags$hr(),

h4(class = "highlight", "Significance thresholds"),

numericInput(
  "sig_p_cutoff",
  "Adj.p-value cutoff:",
  value = 0.05,
  min = 0,
  max = 1,
  step = 0.001
),

sliderInput(
  "fc_thr",
  "FC threshold (|log2FC| ≥):",
  min = 0,
  max = max(1, round(fcmax, 1)),
  value = 1,
  step = 0.1
),

materialSwitch(
  "use_fdr_filter",
  "Filter by Adj.p-value threshold",
  value = FALSE,
  status = "success"
),

materialSwitch(
  "use_fc_filter",
  "Filter by Fold-Change threshold",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_fc_filter == true",
  radioButtons(
    "fc_dir",
    "Direction:",
    choices = c("Both sides" = "both", "Up only" = "up", "Down only" = "down"),
    selected = "both",
    inline = TRUE
  )
)
    )
  })

  # ---- Filtered volcano data
  filtered_volcano <- reactive({
    req(procReady(), input$mz_range, input$rt_range, input$intensity_range)
  
    dd <- rv$volcano
    pcut <- suppressWarnings(as.numeric(input$sig_p_cutoff %||% 0.05))
    if (!is.finite(pcut) || pcut <= 0 || pcut > 1) pcut <- 0.05
    
    fcut <- suppressWarnings(as.numeric(input$fc_thr %||% 1))
    if (!is.finite(fcut) || fcut < 0) fcut <- 1
  
    if (!is.null(input$sel_feat) && length(input$sel_feat) > 0) {
  dd <- dd %>% dplyr::filter(Feature %in% input$sel_feat)
}
    if (isTRUE(input$sig_only)) {
  dd <- dd %>%
    dplyr::filter(
      is.finite(`Adj.p-value`),
      `Adj.p-value` <= pcut,
      is.finite(FC),
      abs(FC) >= fcut
    )
}
  
    if (isTRUE(input$use_npc_filter)) {
  validate(
    need(
      !is.null(input$npc_filter_values) && length(input$npc_filter_values) > 0,
      "NPC filter is enabled. Select at least one NPC class."
    )
  )

  dd <- dd %>%
    dplyr::filter(`NPC#class` %in% input$npc_filter_values)
}

if (isTRUE(input$use_classyfire_filter)) {
  validate(
    need(
      !is.null(input$classyfire_filter_values) && length(input$classyfire_filter_values) > 0,
      "ClassyFire filter is enabled. Select at least one ClassyFire class."
    )
  )

  dd <- dd %>%
    dplyr::filter(`ClassyFire#class` %in% input$classyfire_filter_values)
}
    
    dd <- dd %>%
      dplyr::filter(
        is.finite(mz), is.finite(RT), is.finite(Mean),
        mz >= input$mz_range[1], mz <= input$mz_range[2],
        RT >= input$rt_range[1], RT <= input$rt_range[2],
        log10(Mean + 1.1) >= input$intensity_range[1],
        log10(Mean + 1.1) <= input$intensity_range[2]
      )
  
    if (isTRUE(input$use_fdr_filter)) {
  dd <- dd %>%
    dplyr::filter(
      is.finite(`Adj.p-value`),
      `Adj.p-value` <= pcut
    )
}
    
    if (isTRUE(input$use_fc_filter)) {
      thr <- fcut
      if (input$fc_dir == "both") {
        dd <- dd %>% dplyr::filter(abs(FC) >= thr)
      } else if (input$fc_dir == "up") {
        dd <- dd %>% dplyr::filter(FC >= thr)
      } else {
        dd <- dd %>% dplyr::filter(FC <= -thr)
      }
    }
  
    dd
  })

  # ---- Volcano plot
  output$volcano_plot <- renderPlotly({
    req(filtered_volcano())
    dd <- filtered_volcano()
    validate(need(nrow(dd) > 0, "No points left after filtering."))

    dd$key <- paste(dd$Groups, dd$Feature, sep = "__")

    hover_txt <- paste0(
      "Groups: ", dd$Groups,
      "<br>FC: ", dd$FC,
      "<br>FDR: ", format(dd$`Adj.p-value`, digits = 3, scientific = TRUE),
      "<br>Test scale: ", dd$TestScale,
      "<br>Feature: ", dd$Feature,
      "<br>ID: ", dd$id,
      "<br>m/z: ", dd$mz,
      "<br>RT: ", dd$RT,
      "<br>Mean Intensity: ", format(dd$Mean, big.mark = ",", scientific = FALSE),
      "<br>NPC: ", dd$`NPC#class`,
      "<br>ClassyFire: ", dd$`ClassyFire#class`
    )

    fc_line <- suppressWarnings(as.numeric(input$fc_thr %||% 1))
    if (!is.finite(fc_line) || fc_line < 0) fc_line <- 1
    
    p_thr <- suppressWarnings(as.numeric(input$sig_p_cutoff %||% 0.05))
    if (!is.finite(p_thr) || p_thr <= 0 || p_thr > 1) p_thr <- 0.05
    
    ythr <- -log10(pmax(p_thr, .Machine$double.xmin))
    
    shapes <- list(
  list(
    type = "line",
    x0 = -fc_line,
    x1 = -fc_line,
    xref = "x",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    line = list(dash = "dot")
  ),
  list(
    type = "line",
    x0 = fc_line,
    x1 = fc_line,
    xref = "x",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    line = list(dash = "dot")
  ),
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = ythr,
    y1 = ythr,
    yref = "y",
    line = list(dash = "dot")
  )
)

    if (input$color_by == "Groups") {
      plot_ly(
        data = dd,
        x = ~FC, y = ~`Adj.p-value.log`,
        colors = make_palette(input$volcano_palette %||% "Set1", length(unique(dd$Groups))),
        color = ~Groups,
        type = "scatter", mode = "markers",
        text = hover_txt, hoverinfo = "text",
        key = ~key,
        marker = list(size = 12, opacity = 0.85, line = list(color = "black", width = 1)),
        source = "volcano"
      ) %>%
        layout(
          shapes = shapes,
          xaxis = list(title = "log2(FC)"),
          yaxis = list(title = "-log10(FDR)"),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    } else {
      plot_ly(
        data = dd,
        x = ~FC, y = ~`Adj.p-value.log`,
        color = ~log10(Mean + 1.1),
        symbol = ~Groups,
        colors = make_palette(input$volcano_palette %||% "viridis", 12),
        type = "scatter", mode = "markers",
        text = hover_txt, hoverinfo = "text",
        key = ~key,
        marker = list(size = 12, opacity = 0.9, line = list(color = "black", width = 1)),
        source = "volcano"
      ) %>%
        layout(
          shapes = shapes,
          xaxis = list(title = "log2(FC)"),
          yaxis = list(title = "-log10(FDR)"),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    }
  })

  # ---- Feature plot on click
  output$feature_plot <- renderPlotly({
    req(procReady(), rv$df_used, rv$mat)

    click <- event_data("plotly_click", source = "volcano")
    if (is.null(click) || is.null(click$key)) return(NULL)

    key <- click$key[[1]]
    parts <- str_split_fixed(key, "__", 2)
    comp  <- parts[1, 1]
    feat  <- parts[1, 2]

    df_used <- rv$df_used
    validate(need(feat %in% colnames(df_used), "Clicked feature not found in matrix."))

    row <- rv$volcano %>% filter(Groups == comp, Feature == feat) %>% slice(1)
    s_names <- rownames(rv$mat)
    
    yy <- as.numeric(df_used[[feat]])
    xx <- as.character(df_used$Label)

    box_cols <- make_palette(input$box_palette %||% "Dark2", length(unique(xx)))
    
    if (input$present_as == "Boxplot") {
      plot_ly(
        x = xx, y = yy,
        colors = box_cols,
        type = "box",
        color = xx,
        boxpoints = FALSE,
        marker = list(opacity = 0.8)
      ) %>% hide_legend() %>% 
        layout(
          title = list(text = paste0(feat, "<br><span style='font-size:12px;'>")),
          xaxis = list(title = "Group"),
          yaxis = list(title = "Intensity")
        )
    } else {
      plot_ly(
        x = xx, y = yy,
        text = s_names,
        colors = box_cols,
        type = "scatter", mode = "markers",
        color = xx,
        hovertemplate = paste(
          "<b>Sample:</b> %{text}<br>",
          "<b>Group:</b> %{x}<br>",
          "<b>Intensity:</b> %{y}",
          "<extra></extra>" 
        ),
        marker = list(size = 20, opacity = 0.85, symbol = "diamond", line = list(color = "black", width = 2))
      ) %>% hide_legend() %>% 
        layout(
          title = list(text = paste0(feat, "<br><span style='font-size:12px;'>")),
          xaxis = list(title = "Group"),
          yaxis = list(title = "Intensity")
        )
    }
  })

  # ---- Downloads
 output$dl_volcano <- downloadHandler(
  filename = function() {
    req(rv$volcano)

    n_comp <- length(unique(rv$volcano$Groups))

    if (n_comp > 1) {
      paste0(dataset_name(), "_volcano_table_wide.csv")
    } else {
      paste0(dataset_name(), "_volcano_table.csv")
    }
  },

  content = function(file) {
    req(rv$volcano)

    out <- volcano_to_wide_if_needed(rv$volcano)

    data.table::fwrite(out, file, na = "")
  }
)

output$dl_matrix <- downloadHandler(
  filename = function() {
    paste0(dataset_name(), "_MetaboAnalyst_table.csv")
  },
  content = function(file) {
    req(rv$df_used, rv$mat)

    out <- as.data.frame(rv$df_used, check.names = FALSE, stringsAsFactors = FALSE)
    out <- cbind(
      Sample = rownames(rv$mat),
      out
    )

    data.table::fwrite(out, file, na = "")
  }
)

output$dl_autoplotter_zip <- downloadHandler(
  filename = function() {
    nm <- input$file_data$name %||% "dataset.csv"
    paste0(tools::file_path_sans_ext(basename(nm)), "_AutoPlotter.zip")
  },

  content = function(file) {
    req(rv$df_used, rv$mat, rv$fmap)

    zip_dir <- tempfile("autoplotter_")
    dir.create(zip_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(zip_dir, recursive = TRUE, force = TRUE), add = TRUE)

    data_file <- file.path(zip_dir, "data_table.csv")
    meta_file <- file.path(zip_dir, "metadata.csv")
    name_file <- file.path(zip_dir, "name_map.csv")

    autoplotter_data <- make_autoplotter_data(
      df_used = rv$df_used,
      sample_names = rownames(rv$mat)
    )

    autoplotter_metadata <- make_autoplotter_metadata(
      df_used = rv$df_used,
      sample_names = rownames(rv$mat)
    )

    autoplotter_name_map <- make_autoplotter_name_map(
      fmap = rv$fmap,
      volcano = rv$volcano
    )

    data.table::fwrite(autoplotter_data, data_file, na = "")
    data.table::fwrite(autoplotter_metadata, meta_file, na = "")
    data.table::fwrite(autoplotter_name_map, name_file, na = "")

    tmp_zip <- tempfile(fileext = ".zip")
    on.exit(unlink(tmp_zip, force = TRUE), add = TRUE)

    zip::zipr(
      zipfile = tmp_zip,
      files = c(data_file, meta_file, name_file),
      root = zip_dir
    )

    ok <- file.copy(tmp_zip, file, overwrite = TRUE)
    if (!ok) {
      stop("Failed to copy AutoPlotter ZIP archive to download file.")
    }
  },

  contentType = "application/zip"
)

}

shinyApp(ui, server)
#.....................................................
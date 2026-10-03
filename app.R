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
  library(igraph)
  library(ComplexHeatmap)
  library(InteractiveComplexHeatmap)
  library(ggrepel)
})

options(shiny.maxRequestSize = 1024 * 1024^2)

# ----------------------------- Helpers -----------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

clean_missing_text <- function(x) {

  x <- trimws(
    as.character(x)
  )

  bad <- (
    is.na(x) |
      !nzchar(x) |
      tolower(x) %in% c(
        "na",
        "nan",
        "null",
        "not provided"
      )
  )

  x[bad] <- NA_character_

  x
}

make_prefixed_colmap <- function(cols, prefix) {

  cols <- as.character(
    cols %||% character(0)
  )

  cols <- cols[
    !is.na(cols) &
      nzchar(cols)
  ]

  if (!length(cols)) {
    return(
      stats::setNames(
        character(0),
        character(0)
      )
    )
  }

  output_names <- make.unique(
    paste0(prefix, cols),
    sep = "_"
  )

  stats::setNames(
    output_names,
    cols
  )
}

format_extra_value <- function(x) {

  if (is.null(x) || !length(x)) {
    return("NA")
  }

  value <- x[[1]]

  if (length(value) == 0 || is.na(value)) {
    return("NA")
  }

  if (is.numeric(value)) {

    if (!is.finite(value)) {
      return("NA")
    }

    return(
      format(
        value,
        digits = 12,
        scientific = FALSE,
        trim = TRUE
      )
    )
  }

  value <- trimws(
    as.character(value)
  )

  if (!nzchar(value)) {
    "NA"
  } else {
    value
  }
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

if (!all(ok)) {
  msg <- sprintf(
    "Warning: Token %d missing in some sample names. Falling back to the full sample name.",
    idx
  )

  warning(msg)

  if (!is.null(shiny::getDefaultReactiveDomain())) {
    shiny::showNotification(
      msg,
      type = "warning",
      duration = 8
    )
  }
}

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

  suffixes <- unique(parse_suffix_list(remove_suffixes))
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

safe_ttest_p <- function(
    x,
    g,
    paired = FALSE,
    var.equal = FALSE
) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  g <- droplevels(
    as.factor(g)
  )

  # Remove observations without a valid group
  valid_group <- !is.na(g)

  x <- x[valid_group]
  g <- droplevels(
    g[valid_group]
  )

  if (nlevels(g) != 2) {
    return(NA_real_)
  }

  group_levels <- levels(g)

  a <- x[
    g == group_levels[1]
  ]

  b <- x[
    g == group_levels[2]
  ]

  if (isTRUE(paired)) {

    # Paired samples must have equal lengths
    if (length(a) != length(b)) {
      return(NA_real_)
    }

    # Remove missing values pairwise
    complete_pairs <- is.finite(a) &
      is.finite(b)

    a <- a[complete_pairs]
    b <- b[complete_pairs]

    if (length(a) < 2) {
      return(NA_real_)
    }

    p_value <- tryCatch(

      stats::t.test(
        x = a,
        y = b,
        paired = TRUE
      )$p.value,

      error = function(e) {
        NA_real_
      }
    )

  } else {

    # For unpaired tests, remove missing values independently
    a <- a[
      is.finite(a)
    ]

    b <- b[
      is.finite(b)
    ]

    if (
      length(a) < 2 ||
      length(b) < 2
    ) {
      return(NA_real_)
    }

    p_value <- tryCatch(

      stats::t.test(
        x = a,
        y = b,
        paired = FALSE,
        var.equal = isTRUE(var.equal)
      )$p.value,

      error = function(e) {
        NA_real_
      }
    )
  }

  if (
    length(p_value) != 1 ||
    !is.finite(p_value)
  ) {
    return(NA_real_)
  }

  as.numeric(p_value)
}

safe_wilcox_p <- function(
    x,
    g,
    paired = FALSE
) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  g <- droplevels(
    as.factor(g)
  )

  # Remove observations without a valid group
  valid_group <- !is.na(g)

  x <- x[valid_group]
  g <- droplevels(
    g[valid_group]
  )

  if (nlevels(g) != 2) {
    return(NA_real_)
  }

  group_levels <- levels(g)

  a <- x[
    g == group_levels[1]
  ]

  b <- x[
    g == group_levels[2]
  ]

  if (isTRUE(paired)) {

    # Paired samples must have equal lengths
    if (length(a) != length(b)) {
      return(NA_real_)
    }

    # Remove missing values pairwise
    complete_pairs <- is.finite(a) &
      is.finite(b)

    a <- a[complete_pairs]
    b <- b[complete_pairs]

    if (length(a) < 1) {
      return(NA_real_)
    }

    p_value <- tryCatch(

      suppressWarnings(
        stats::wilcox.test(
          x = a,
          y = b,
          paired = TRUE,

          # Explicit approximation gives more consistent
          # behavior across different R versions
          exact = FALSE,
          correct = TRUE
        )$p.value
      ),

      error = function(e) {
        NA_real_
      }
    )

  } else {

    # For unpaired tests, remove missing values independently
    a <- a[
      is.finite(a)
    ]

    b <- b[
      is.finite(b)
    ]

    if (
      length(a) < 1 ||
      length(b) < 1
    ) {
      return(NA_real_)
    }

    p_value <- tryCatch(

      suppressWarnings(
        stats::wilcox.test(
          x = a,
          y = b,
          paired = FALSE,

          # Explicit approximation gives more consistent
          # behavior across different R versions
          exact = FALSE,
          correct = TRUE
        )$p.value
      ),

      error = function(e) {
        NA_real_
      }
    )
  }

  if (
    length(p_value) != 1 ||
    !is.finite(p_value)
  ) {
    return(NA_real_)
  }

  as.numeric(p_value)
}

parse_feature_table_to_matrix <- function(
  raw_df,
  feature_id_source,
  annotation_id_col,
  mz_col,
  rt_col,
  sample_keywords = NULL,
  sample_cols = NULL,
  mz_rt_sep = "@"
) {

  raw_df <- clean_mzmine_export(raw_df)

  cols <- names(raw_df)

  validate(

    need(
        annotation_id_col %in% cols,
        "Annotation matching ID column not found."
      ),
    
    need(
      mz_col %in% cols,
      "m/z column not found."
    ),

    need(
      rt_col %in% cols,
      "RT column not found."
    ),

    need(
      feature_id_source %in% c("combine_mz_rt", "auto") ||
        feature_id_source %in% cols,
      "Selected Feature ID source was not found."
    )
  )


  # -----------------------------
  # Sample columns
  # -----------------------------
  if (
    !is.null(sample_cols) &&
    length(sample_cols) > 0
  ) {

    sample_cols <- intersect(
      sample_cols,
      cols
    )

    validate(
      need(
        length(sample_cols) > 0,
        "No selected sample columns were found in the peak table."
      )
    )

  } else {

    sidx <- multi_sample_idx(
      cols,
      sample_keywords
    )

    validate(
      need(
        length(sidx) > 0,
        sprintf(
          "No sample columns matched keywords: %s",
          paste(
            sample_keywords,
            collapse = ", "
          )
        )
      )
    )

    sample_cols <- cols[sidx]
  }


  # -----------------------------
  # Samples x features matrix
  # -----------------------------
  mat <- as.data.frame(
    data.table::transpose(
      raw_df[
        ,
        sample_cols,
        drop = FALSE
      ]
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  rownames(mat) <- sample_cols


  # -----------------------------
  # m/z and RT
  # -----------------------------
  mz <- suppressWarnings(
    as.numeric(
      raw_df[[mz_col]]
    )
  )

  rt <- suppressWarnings(
    as.numeric(
      raw_df[[rt_col]]
    )
  )


  # -----------------------------
  # Feature ID
  # -----------------------------
  if (
    identical(
      feature_id_source,
      "combine_mz_rt"
    )
  ) {

    sep <- as.character(
      mz_rt_sep %||% "@"
    )

    if (!nzchar(sep)) {
      sep <- "@"
    }

    mz_text <- ifelse(
      is.na(mz),
      "NA",
      as.character(
        round(mz, 4)
      )
    )

    rt_text <- ifelse(
      is.na(rt),
      "NA",
      as.character(
        round(rt, 2)
      )
    )

    Feature <- paste0(
      mz_text,
      sep,
      rt_text
    )

  } else if (
    identical(
      feature_id_source,
      "auto"
    )
  ) {

    Feature <- paste0(
      "feat_",
      seq_len(
        nrow(raw_df)
      )
    )

  } else {

    Feature <- trimws(
  as.character(
    raw_df[[feature_id_source]]
  )
)

    # Protect against missing/empty IDs
    bad_id <- is.na(Feature) |
      !nzchar(Feature)

    if (any(bad_id)) {

      Feature[bad_id] <- paste0(
        "feat_",
        which(bad_id)
      )
    }
  }


  # Ensure unique feature names
  Feature <- make.unique(
    Feature,
    sep = "_"
  )

  colnames(mat) <- Feature

  id <- trimws(
  as.character(
    raw_df[[annotation_id_col]]
  )
)
  
  # Keep "id" because the rest of the current app
  # uses this column for annotation/network matching.
  # It now corresponds to the selected Feature ID.
  fmap <- tibble::tibble(
  id = id,
  mz = mz,
  RT = rt,
  Feature = Feature
)


  # Numeric intensities
  mat[] <- lapply(
    mat,
    function(z) {
      suppressWarnings(
        as.numeric(z)
      )
    }
  )

  mat[is.na(mat)] <- 0


  list(
    mat = mat,
    fmap = fmap,
    raw = raw_df
  )
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
    gnum <- comb[i, 1] # Numerator
    gden <- comb[i, 2] # Denominator
    comp <- paste0(gnum, " / ", gden)

    if (
  isTRUE(paired) &&
  test %in% c("Student", "Wilcoxon")
) {

  n_num <- sum(as.character(gr) == gnum, na.rm = TRUE)
  n_den <- sum(as.character(gr) == gden, na.rm = TRUE)

  minimum_pairs <- if (test == "Student") 2L else 1L

  validate(
    need(
      n_num == n_den,
      paste0(
        "Paired test: ", comp,
        " has unequal group sizes (",
        n_num, " and ", n_den, ")."
      )
    ),
    need(
      min(n_num, n_den) >= minimum_pairs,
      paste0(
        "Paired test: ", comp,
        " needs at least ", minimum_pairs,
        " sample pair(s)."
      )
    )
  )
}
    
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
    dd$Significant_default <- (dd$`Adj.p-value` <= 0.05) & (abs(dd$FC) >= 1)

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

  # Columns that vary between statistical comparisons
comparison_specific_cols <- c(
  "Groups",
  "Group_num",
  "Group_den",
  "Adj.p-value",
  "Mean",
  "mean_num",
  "mean_den",
  "FC",
  "TestScale",
  "Adj.p-value.log",
  "Significant_default"
)

# Everything else is feature-level metadata and should
# remain in the wide report, including Peak_* and SIRIUS_*.
static_cols <- setdiff(
  names(volc),
  comparison_specific_cols
)

# Keep the most important columns first
preferred_static_order <- c(
  "Feature",
  "id",
  "mz",
  "RT",
  "NPC#class",
  "ClassyFire#class",
  "GNPS_annotation",
  "Other_annotation"
)

static_cols <- c(
  intersect(
    preferred_static_order,
    static_cols
  ),

  setdiff(
    static_cols,
    preferred_static_order
  )
)

  n_comp <- dplyr::n_distinct(volc$Groups)

  if (n_comp == 1) {

  group_num_name <- unique(as.character(volc$Group_num))
  group_den_name <- unique(as.character(volc$Group_den))
  comparison_name <- unique(as.character(volc$Groups))

  if (
    length(group_num_name) == 1 &&
    length(group_den_name) == 1
  ) {
    names(volc)[names(volc) == "mean_num"] <-
      paste0("Mean_", group_num_name)

    names(volc)[names(volc) == "mean_den"] <-
      paste0("Mean_", group_den_name)
  }

  if (length(comparison_name) == 1) {
    names(volc)[names(volc) == "Mean"] <-
      paste0("Mean__", comparison_name)
  }

  volc <- volc %>%
    dplyr::select(
      -dplyr::any_of(
        c(
          "Group_num",
          "Group_den",
          "TestScale",
          "Adj.p-value.log"
        )
      )
    )

  return(volc)
}

  comparison_metrics <- intersect(
  c(
    "FC",
    "Adj.p-value",
    "Mean",
    "Significant_default"
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

  fmap <- as.data.frame(
    fmap,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Keep names and ordering consistent with the AutoPlotter data.
  out <- tibble::tibble(
    Name = as.character(fmap$Feature)
  )

  if (
    is.null(volcano) ||
    !"Feature" %in% names(volcano) ||
    nrow(volcano) == 0
  ) {
    return(out)
  }

  volcano <- as.data.frame(
    volcano,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Exclude comparison-specific statistics.
  # Keep all feature metadata and annotation columns.
  statistical_cols <- c(
    "Groups",
    "Group_num",
    "Group_den",
    "Adj.p-value",
    "Mean",
    "mean_num",
    "mean_den",
    "FC",
    "TestScale",
    "Adj.p-value.log",
    "Significant_default"
  )

  annotation_cols <- setdiff(
    names(volcano),
    c("Feature", statistical_cols)
  )

  ann <- volcano %>%
    dplyr::select(
      Feature,
      dplyr::all_of(annotation_cols)
    ) %>%
    dplyr::mutate(
      Feature = as.character(Feature),
      dplyr::across(
        dplyr::where(
          function(x) is.character(x) || is.factor(x)
        ),
        clean_missing_text
      )
    ) %>%
    # Annotations repeat across comparisons:
    # retain one row per feature.
    dplyr::distinct(
      Feature,
      .keep_all = TRUE
    )

  out %>%
    dplyr::left_join(
      ann,
      by = c("Name" = "Feature")
    )
}

volcano_main_ui <- function() {

  tagList(

        div(
      style = paste(
        "background:rgba(255,255,255,0.95);",
        "border:1px solid #ddd;",
        "border-radius:8px;",
        "padding:12px 16px;",
        "margin-bottom:15px;"
      ),

      div(
        style = "font-size:17px;font-weight:600;",
        textOutput("volcano_feature_count", inline = TRUE)
      ),

      tags$details(
  style = "margin-top:8px;",

  tags$summary(
    style = "cursor:pointer;color:#228B22;",
    "Applied filters"
  ),

  uiOutput("volcano_applied_filters"),

  tags$details(
    style = "margin-top:12px;",

    tags$summary(
      style = "cursor:pointer;color:#228B22;",
      "Download filtered peak table"
    ),

    div(
      style = "margin-top:10px;",

      p(
        class = "small-note",
        paste(
          "Export the original peak table with only",
          "features retained by the current filters."
        )
      ),

      downloadButton(
        outputId = "dl_filtered_feature_table",
        label = "Download filtered CSV",
        class = "btn-success"
      )
    )
  )
)
    ),
    
    # ========================================================
    # VOLCANO VIEW
    # ========================================================

    conditionalPanel(
      condition = "!input.show_interactive_heatmap",

      withSpinner(
        plotlyOutput(
          "volcano_plot",
          height = "520px"
        ),
        type = 8,
        color = "#66CDAA"
      ),

      div(
        style = "height:8px;"
      ),

      uiOutput(
        "selected_feature_panel"
      )
    ),


    # ========================================================
    # INTERACTIVE HEATMAP VIEW
    # ========================================================

    conditionalPanel(
      condition = "input.show_interactive_heatmap",

      div(

        style = "
          background: rgba(255,255,255,0.90);
          border: 1px solid #ddd;
          border-radius: 8px;
          padding: 12px;
          margin-bottom: 15px;
        ",

        h4(
          class = "highlight",
          "Interactive heatmap"
        ),

        uiOutput(
          "heatmap_filter_summary"
        ),


        # ----------------------------------------------------
        # Heatmap-specific settings
        # ----------------------------------------------------

        fluidRow(

          column(
            width = 4,

            selectInput(
              "hm_scale",
              "Scaling:",
              choices = c(
                "Unit variance, no centering" = "uv",
                "Z-score" = "zscore",
                "None" = "none"
              ),
              selected = "uv"
            )
          ),

          column(
            width = 4,

            selectInput(
              "hm_distance",
              "Clustering distance:",
              choices = c(
                "Euclidean" = "euclidean",
                "Manhattan" = "manhattan",
                "Correlation" = "correlation"
              ),
              selected = "euclidean"
            )
          ),

          column(
            width = 4,

            selectInput(
              "hm_method",
              "Clustering method:",
              choices = c(
                "Ward.D2" = "ward.D2",
                "Complete" = "complete",
                "Average" = "average"
              ),
              selected = "ward.D2"
            )
          )
        ),


        fluidRow(

          column(
            width = 3,

            checkboxInput(
              "hm_cluster_samples",
              "Cluster samples",
              value = TRUE
            )
          ),

          column(
            width = 3,

            checkboxInput(
              "hm_cluster_features",
              "Cluster features",
              value = TRUE
            )
          ),

          column(
            width = 3,

            checkboxInput(
              "hm_show_samples",
              "Show sample names",
              value = FALSE
            )
          ),

          column(
            width = 3,

            checkboxInput(
              "hm_show_features",
              "Show feature names",
              value = FALSE
            ),
            
            checkboxInput(
                "hm_show_borders",
                "Show cell borders",
                value = TRUE
              )
          )
        ),


        fluidRow(

          column(
            width = 4,

            selectInput(
              "hm_palette",
              "Heatmap palette:",
              choices = c(
                "Viridis" = "viridis",
                "Magma" = "magma",
                "Blue - White - Red" = "bwr"
              ),
              selected = "bwr"
            )
          ),

          column(
            width = 4,

            selectInput(
              "hm_group_palette",
              "Group annotation palette:",
              choices = palette_choices,
              selected = "Dark2"
            )
          )
        )
      ),


  tagList(
  conditionalPanel(
    condition = "output.heatmap_has_data === 'yes'",

    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
      heatmap_id = "metabocano_heatmap",
      layout = "1-(2|3)",
      width1 = 700,
      height1 = 550,
      width2 = 350,
      height2 = 300,
      action = "click",
      cursor = TRUE,
      output_ui = shiny::uiOutput("heatmap_feature_info")
    )
  ),

  conditionalPanel(
    condition = "output.heatmap_has_data !== 'yes'",

    div(
      class = "small-note",
      paste(
        "No heatmap is available for the current selection.",
        "Check the filters and preprocessing status."
      )
    )
  )
)
    )
  )
}

add_volcano_top_labels <- function(
    p, dd, n = 0L, width_px = 800, height_px = 520,
    label_col = "Feature"
) {

  n <- suppressWarnings(as.integer(n))
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    return(p)
  }

  # Keep valid plotted coordinates.
  d <- as.data.frame(dd, stringsAsFactors = FALSE)

  d <- d[
    is.finite(d$FC) & is.finite(d$plot_y),
    ,
    drop = FALSE
  ]

  if (!nrow(d)) return(p)

  # Rank only valid FDR values.
  d$.score <- NA_real_

  valid <- is.finite(d$`Adj.p-value`) &
    d$`Adj.p-value` >= 0 &
    d$`Adj.p-value` <= 1

  d$.score[valid] <-
    -log10(
      pmax(d$`Adj.p-value`[valid], .Machine$double.xmin)
    ) * abs(d$FC[valid])

  ranked <- which(is.finite(d$.score))

  ranked <- ranked[
    order(
      -d$.score[ranked],
      as.character(d$Feature[ranked]),
      ranked
    )
  ]

  # One label per feature.
  # With multiple comparisons, use its highest-scoring point.
  ranked <- ranked[
    !duplicated(as.character(d$Feature[ranked]))
  ]

  selected <- head(ranked, n)

  if (!length(selected)) return(p)

  # Put labeled points first so ggrepel's label numbering
  # corresponds to these rows.
  d <- d[
    c(selected, setdiff(seq_len(nrow(d)), selected)),
    ,
    drop = FALSE
  ]

  n_labels <- length(selected)

  # Read the selected label column.
if (
  length(label_col) != 1L ||
  is.na(label_col) ||
  !label_col %in% names(d)
) {
  label_col <- "Feature"
}

label_values <- clean_missing_text(
  as.character(d[[label_col]])
)

# Missing annotations fall back to the Feature name.
missing_label <- is.na(label_values) |
  !nzchar(trimws(label_values))

label_values[missing_label] <-
  as.character(d$Feature[missing_label])

d$.label <- ""
d$.label[seq_len(n_labels)] <-
  label_values[seq_len(n_labels)]

  padded_range <- function(x) {
    r <- range(x, finite = TRUE)
    span <- diff(r)

    if (span == 0) {
      span <- max(abs(r), 1)
    }

    r + c(-1, 1) * span * 0.08
  }

  xr <- padded_range(d$FC)
  yr <- padded_range(d$plot_y)

  # Off-screen plot used only to calculate label positions.
  label_plot <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = FC,
      y = plot_y,
      label = .label
    )
  ) +
    ggrepel::geom_text_repel(
      size = 3.2,
      family = "sans",
      seed = 123,
      max.overlaps = Inf,
      max.time = 1,
      max.iter = 10000,
      box.padding = 0.5,
      point.padding = 0.3,
      point.size = 3,
      min.segment.length = 0
    ) +
    ggplot2::coord_cartesian(
      xlim = xr,
      ylim = yr,
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = grid::unit(rep(0, 4), "pt")
    )

  grDevices::pdf(
    file = NULL,
    width = max(300, width_px - 160) / 96,
    height = max(250, height_px - 120) / 96
  )

  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::grid.draw(ggplot2::ggplotGrob(label_plot))
  grid::grid.force()

  # Enter the panel viewport to interpret ggrepel coordinates.
  tree <- grid::grid.ls(
    viewports = TRUE,
    print = FALSE
  )

  panel_vp <- tree$name[
    tree$type == "vpListing" &
      grepl("^panel([.-]|$)", tree$name)
  ]

  if (!length(panel_vp)) {
    stop("Could not locate the ggrepel panel viewport.")
  }

  grid::seekViewport(panel_vp[[1]])

  annotations <- lapply(seq_len(n_labels), function(i) {

    label_grob <- grid::grid.get(
      paste0("textrepelgrob", i),
      grep = FALSE,
      global = TRUE
    )

    if (is.null(label_grob)) {
      stop("Could not retrieve a ggrepel label position.")
    }

    nx <- grid::convertX(
      label_grob$x, "npc", valueOnly = TRUE
    )

    ny <- grid::convertY(
      label_grob$y, "npc", valueOnly = TRUE
    )

    list(
      # Connector points to the original feature.
      x = d$FC[i],
      y = d$plot_y[i],
      xref = "x",
      yref = "y",

      # Repelled label position.
      ax = xr[1] + nx * diff(xr),
      ay = yr[1] + ny * diff(yr),
      axref = "x",
      ayref = "y",

      text = as.character(
        htmltools::htmlEscape(d$.label[i])
      ),
      showarrow = TRUE,
      arrowhead = 0,
      arrowwidth = 0.8,
      arrowcolor = "grey50",
      xanchor = "center",
      yanchor = "middle",
      font = list(
        size = 12,
        color = "black",
        family = "Arial"
      ),
      captureevents = FALSE
    )
  })

  plotly::layout(
    p,
    annotations = annotations,
    xaxis = list(range = xr),
    yaxis = list(range = yr)
  )
}

annotation_panel <- function(
    switch_id, label, tooltip_id, tooltip_text, ...
) {
  div(
    style = paste(
      "background:rgba(255,255,255,0.95);",
      "border:1px solid #d9e2dc;",
      "border-left:4px solid #5cb85c;",
      "border-radius:8px;",
      "padding:14px;",
      "margin-bottom:12px;"
    ),

    div(
      style = paste(
        "display:flex;",
        "align-items:flex-start;",
        "justify-content:space-between;",
        "gap:8px;"
      ),

      shinyWidgets::materialSwitch(
        inputId = switch_id,
        label = label,
        value = FALSE,
        status = "success",
        width = "auto"
      ),

      actionButton(
        inputId = tooltip_id,
        label = "?",
        class = "btn-xs",
        style = "font-weight:bold;flex-shrink:0;"
      )
    ),

    shinyBS::bsTooltip(
      id = tooltip_id,
      title = tooltip_text,
      placement = "right",
      trigger = "click",
      options = list(container = "body")
    ),

    conditionalPanel(
      condition = paste0("input.", switch_id, " == true"),

      tags$details(
        open = NA,

        tags$summary(
          style = paste(
            "cursor:pointer;",
            "color:#27823b;",
            "font-weight:600;",
            "padding:6px 0;"
          ),
          "Upload and column settings"
        ),

        div(
          style = paste(
            "border-top:1px solid #e5e5e5;",
            "padding-top:12px;",
            "margin-top:6px;"
          ),
          ...
        )
      )
    )
  )
}

# ----------------------------- UI -----------------------------------------

ui <- fluidPage(
  useShinyjs(),

  tags$head(
  tags$script(HTML("
    (function () {

      function fixFileInputs() {
        document.querySelectorAll('input[type=file]').forEach(
          function (input) {

            const button = input.closest('.btn-file');

            if (!button) return;

            // Keep the hidden input inside its Browse button.
            button.style.setProperty(
              'position', 'relative', 'important'
            );

            input.style.setProperty(
              'position', 'absolute', 'important'
            );
            input.style.setProperty(
              'top', '0', 'important'
            );
            input.style.setProperty(
              'left', '0', 'important'
            );
            input.style.setProperty(
              'width', '1px', 'important'
            );
            input.style.setProperty(
              'height', '1px', 'important'
            );
            input.style.setProperty(
              'opacity', '0', 'important'
            );
            input.style.setProperty(
              'pointer-events', 'none', 'important'
            );
          }
        );
      }

      // Inputs present when the page first loads.
      $(fixFileInputs);

      // Inputs subsequently created by renderUI().
      $(document).on('shiny:bound', fixFileInputs);

    })();
  "))
),
  
tags$head(tags$style(HTML("
  .app-footer { position: fixed; left:0; right:0; bottom:0; 
                text-align:center; font-size:12px; opacity:0.75;
                padding:8px; background: rgba(255,255,255,0.8);
                border-top: 1px solid #ddd; z-index: 9999; }
  body { padding-bottom: 45px; }
"))),

tags$head(
  tags$script(
    HTML("
      $(document).on(
        'click',
        '#selected_feature_info_btn',
        function() {

          var txt = $(this).attr('data-copy');

          if (!txt || !navigator.clipboard) {
            return;
          }

          navigator.clipboard.writeText(txt).then(function() {
            Shiny.setInputValue(
              'selected_feature_info_copied',
              Date.now(),
              {priority: 'event'}
            );
          });
        }
      );
    ")
  )
),

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

div(class = "app-footer", HTML('
    <span class="footer-text">by Plyushchenko I.V.</span>
    <span class="footer-sep">&nbsp;|&nbsp;</span>
    <span class="footer-text">GPLv3</span>
    <span class="footer-sep">&nbsp;|&nbsp;</span>
     <a id="latest-release-link"
     class="footer-link"
     href="https://github.com/plyush1993/metabocano/releases/latest"
     target="_blank">v. </a>
    
    <script>
    fetch("https://api.github.com/repos/plyush1993/metabocano/releases/latest")
      .then(function(response) {
        if (!response.ok) throw new Error("GitHub release request failed");
        return response.json();
      })
      .then(function(data) {
        var link = document.getElementById("latest-release-link");
        if (link && data.tag_name) {
          link.textContent = "v. " + data.tag_name;
          if (data.html_url) {
            link.href = data.html_url;
          }
        }
      })
      .catch(function(error) {
        var link = document.getElementById("latest-release-link");
        if (link) {
          link.textContent = "Latest release";
          link.href = "https://github.com/plyush1993/metabocano/releases/latest";
        }
      });
  </script>
    
  ')),
  
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
  choices = c(
    "mzMine" = "mzmine",
    "xcms" = "xcms",
    "MS-DIAL" = "msdial",
    "Generic" = "default"
  ),
  selected = "mzmine"
),

fileInput(
  "file_data",
  "Upload feature table (.csv)",
  accept = ".csv"
),

helpText(
  HTML(
    "<i class='fa fa-info-circle'></i> Need data to test? Download example datasets from our <a href='https://github.com/plyush1993/Metabocano' target='_blank'>GitHub</a>."
  )
),

uiOutput("upload_tab_error"),

uiOutput("col_pickers"),

          radioButtons(
  "sample_mode",
  "How to define sample columns?",
  choices = c(
    "Auto-detect numeric sample columns" = "auto",
    "By keyword match" = "kws",
    "Pick columns manually" = "manual"
  ),
  selected = "kws"
),

conditionalPanel(
  condition = "input.sample_mode == 'kws'",

  selectizeInput(
    "sample_keywords",
    "Sample column keywords (pick/add multiple):",
    choices = c(
      ".mzML", ".mzXML", ".raw", ".d", ".wiff", ".lcd",
      "Peak area", "Area", "_Area"
    ),
    selected = c(".mzML", ".mzXML"),
    multiple = TRUE,
    options = list(
      create = TRUE,
      createOnBlur = TRUE,
      placeholder = "Type to add keyword and press Enter"
    )
  )
),

conditionalPanel(
  condition = "input.sample_mode == 'manual'",
  uiOutput("manual_sample_cols_ui")
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
        ".mzML", ".mzXML", ".raw", ".RAW", ".lcd",
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

          h3(class = "highlight", "Join with Annotation"),

annotation_panel(
  switch_id = "use_peak_extra_cols",
  label = "Additional peak-table columns",
  tooltip_id = "btnAD",

  tooltip_text = paste0(
    "Selected columns will be added to the downloaded volcano table ",
    "and displayed after clicking a volcano point."
  ),

  uiOutput("peak_extra_cols_ui")
),

annotation_panel(
  switch_id = "use_sirius",
  label = "SIRIUS / CANOPUS annotation",
  tooltip_id = "btn5",

  tooltip_text = paste0(
    "<b>Join SIRIUS annotations with the processed table.</b><br>",
    "The selected peak-table Feature ID column is matched ",
    "to the selected SIRIUS mapping ID column.<br>",
    "Default mapping ID: <em>mappingFeatureId</em><br>",
    "Default NPC column: <em>NPC#class</em><br>",
    "Default ClassyFire column: <em>ClassyFire#class</em>"
  ),

  fileInput(
    "file_sirius",
    "Upload SIRIUS output (.csv/.tsv/.txt)",
    accept = c(".csv", ".tsv", ".txt")
  ),

  uiOutput("sirius_pickers")
),

annotation_panel(
  switch_id = "use_gnps_annotation",
  label = "GNPS library annotation",
  tooltip_id = "btn_gnps_annotation",

  tooltip_text = paste0(
    "<b>Join GNPS library annotations with the processed table.</b><br>",
    "The selected peak-table Feature ID column is matched ",
    "to the selected GNPS ID column.<br>",
    "Default GNPS ID: <em>#Scan#</em><br>",
    "Default annotation: <em>Compound_Name</em>"
  ),

  fileInput(
    "file_gnps_annotation",
    "Upload GNPS library results (.tsv/.txt/.csv)",
    accept = c(".tsv", ".txt", ".csv")
  ),

  uiOutput("gnps_annotation_pickers")
),

annotation_panel(
  switch_id = "use_main_gnps_pairs",
  label = "GNPS network / ComponentIndex",
  tooltip_id = "btn_main_gnps_pairs",

  tooltip_text = paste0(
    "<b>Enable filtering by GNPS network component.</b><br>",
    "Select the peak-table matching ID column, both pairs-file ",
    "node ID columns, and the component column.<br>",
    "Default pairs columns: <em>CLUSTERID1</em>, ",
    "<em>CLUSTERID2</em>, and <em>ComponentIndex</em>.<br>",
    "Both endpoints are assigned to their component."
  ),

  fileInput(
    "file_main_gnps_pairs",
    "Upload GNPS network pairs (.tsv/.txt/.csv)",
    accept = c(".tsv", ".txt", ".csv")
  ),

  uiOutput("main_gnps_pairs_pickers")
),

annotation_panel(
  switch_id = "use_other_annotation",
  label = "Other annotation source",
  tooltip_id = "btn_other_annotation",

  tooltip_text = paste0(
    "<b>Join annotations from an external table.</b><br>",
    "Choose a peak-table ID column and the corresponding ",
    "ID column in the annotation file.<br>",
    "Choose one primary annotation column and optionally ",
    "add additional columns."
  ),

  fileInput(
    "file_other_annotation",
    "Upload annotation table (.csv/.tsv/.txt)",
    accept = c(".csv", ".tsv", ".txt")
  ),

  uiOutput("other_annotation_pickers")
),

          tags$hr(),
          h3(class = "highlight", "Imputation by Noise"),
          radioButtons("do_mvi", "Imputation:", c("No"="no", "Yes"="yes"), selected = "no", inline = TRUE),
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
            checkboxInput(
  "paired",
  "Paired test — samples paired by order",
  FALSE
),

conditionalPanel(
  condition = "input.paired == true",

  helpText(
    paste(
      "Samples are paired by their order within each group.",
      "Check every pair below before running preprocessing.",
      "Sample names are not used to identify matching subjects."
    )
  ),

  tags$details(
    tags$summary("Show sample pairs"),
    div(
      style = "overflow-x: auto;",
      tableOutput("paired_sample_preview")
    )
  )
)
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
          downloadButton("dl_annotation", "Feature table csv", class = "btn-info"),
          actionButton("btn_annotation", "?"),
          bsTooltip("btn_annotation",
            title = paste0(
              "<b>Download table with one row per feature with annotation information.</b><br>",
              "Always includes <em>Feature ID</em>, <em>Annotation matching ID</em>, ",
              "<em>m/z</em>, and <em>RT</em>.<br>",
              "After preprocessing, all joined Peak table, SIRIUS, GNPS, ",
              "and Other Annotation columns are also included."
            ),
            placement = "right",
            trigger = "click",
            options = list(container = "body")
          ),
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
        mainPanel(
  conditionalPanel(
    condition = "output.volcano_ready === 'yes'",
    volcano_main_ui()
  )
)
      )
    ),

navbarMenu(
  title = "3) Other Utils",

  tabPanel(
    title = "SIRIUS & GNPS merging",
    value = "sirius_gnps",

    sidebarLayout(
      sidebarPanel(
        uiOutput("sirius_gnps_sidebar")
      ),
      mainPanel(
        uiOutput("sirius_gnps_main")
      )
    )
  )
)

  )
)

# ----------------------------- Server -------------------------------------

server <- function(input, output, session) {

  upload_error <- reactiveVal(NULL)


output$upload_tab_error <- renderUI({

  msg <- upload_error()

  if (is.null(msg)) {
    return(NULL)
  }

  div(
    style = "
      color: #a94442;
      background-color: #f2dede;
      border: 1px solid #ebccd1;
      padding: 12px;
      margin-bottom: 12px;
      border-radius: 5px;
      font-size: 15px;
      font-weight: bold;
      text-align: center;
    ",
    icon("exclamation-triangle"),
    " ",
    msg
  )
})

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
  
  observeEvent(
  input$selected_feature_info_copied,
  {
    showNotification(
      "m/z and RT copied as: mz,rt",
      type = "message",
      duration = 2
    )
  },
  ignoreInit = TRUE
)
  
  session$onFlushed(function() {
    shinyjs::disable("run_proc")
    shinyjs::disable("dl_annotation")
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

 observe({

  annotation_ready <-
    !is.null(input$file_data) &&

    !is.null(input$feature_id_source) &&

    !is.null(input$annotation_id_col) &&

    !is.null(input$mz_col) &&
    !identical(input$mz_col, "None") &&

    !is.null(input$rt_col) &&
    !identical(input$rt_col, "None")


  shinyjs::toggleState(
    "dl_annotation",
    condition = annotation_ready
  )
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
  
  processing_input_ids <- c(
  "software_tool",
  "feature_id_source",
  "annotation_id_col",
  "mz_col",
  "rt_col",
  "mz_rt_sep",

  "sample_mode",
  "sample_keywords",
  "sample_cols_manual",

  "comparison_mode",
  "ref_group",
  "manual_comparisons",
  "test_type",
  "p_adjust",
  "paired",
  "eqvar",
  "log2_test",
  "standard_scaling",

  "do_mvi",
  "noise_mode",
  "noise_quantile",
  "noise_manual",
  "noise_sd",

  "use_peak_extra_cols",
  "peak_extra_cols",

  "use_sirius",
  "file_sirius",
  "sirius_idcol",
  "sirius_npcol",
  "sirius_cfcol",
  "use_sirius_extra_cols",
  "sirius_extra_cols",

  "use_gnps_annotation",
  "file_gnps_annotation",
  "gnps_annotation_idcol",
  "gnps_annotation_col",
  "use_gnps_extra_cols",
  "gnps_extra_cols",

  "use_other_annotation",
  "file_other_annotation",
  "other_peak_id_col",
  "other_annotation_idcol",
  "other_annotation_col",
  "use_other_extra_cols",
  "other_extra_cols"
)

observeEvent(
  lapply(processing_input_ids, function(id) input[[id]]),

  {
    if (is.null(rv$volcano)) {
      return(invisible(NULL))
    }

    rv$raw <- NULL
    rv$mat <- NULL
    rv$fmap <- NULL
    rv$labels <- NULL
    rv$df_used <- NULL
    rv$volcano <- NULL

    showNotification(
      "Processing settings changed. Run preprocessing again.",
      type = "warning",
      duration = 6
    )
  },

  ignoreInit = TRUE,
  priority = 100
)
  
  dataset_name <- reactive({
  nm <- input$file_data$name %||% "dataset.csv"
  tools::file_path_sans_ext(basename(nm))
})

  # ---- Load raw data (Robust Switch) ----
raw_df <- reactive({

  req(input$file_data)

  tryCatch({

    ext <- tolower(
      tools::file_ext(input$file_data$name)
    )

    if (!identical(ext, "csv")) {
      stop("Please upload a .csv file.")
    }

    tool <- input$software_tool %||% "mzmine"

    if (tool == "msdial") {

      df <- read_msdial_robust(
        input$file_data$datapath
      )

    } else {

      df <- vroom::vroom(
        input$file_data$datapath,
        delim = ",",
        show_col_types = FALSE
      )
    }

    df <- clean_mzmine_export(df)

    if (nrow(df) == 0 || ncol(df) == 0) {
      stop("The uploaded table is empty.")
    }

    upload_error(NULL)

    df

  }, error = function(e) {

    msg <- paste0(
      "Parsing error: ",
      conditionMessage(e)
    )

    upload_error(msg)

    showNotification(
      msg,
      type = "error",
      duration = 6
    )

    validate(
      need(FALSE, msg)
    )
  })

}) 

  output$volcano_label_column_ui <- renderUI({

  req(procReady(), rv$volcano)

  excluded <- c(
    "Groups",
    "Group_num",
    "Group_den",
    "Adj.p-value",
    "Mean",
    "mean_num",
    "mean_den",
    "FC",
    "TestScale",
    "Adj.p-value.log",
    "Significant_default",
    "key",
    "plot_y",
    "FC_status"
  )

  available <- setdiff(
    names(rv$volcano),
    excluded
  )

  preferred <- c(
    "Feature",
    "id",
    "GNPS_annotation",
    "Other_annotation",
    "NPC#class",
    "ClassyFire#class"
  )

  available <- c(
    intersect(preferred, available),
    setdiff(available, preferred)
  )

  req(length(available) > 0)

  selected <- isolate(input$volcano_label_column)

  if (
    is.null(selected) ||
    !selected %in% available
  ) {
    selected <- available[[1]]
  }

  selectInput(
    "volcano_label_column",
    "Label text:",
    choices = available,
    selected = selected
  )
})
  
  # ---- Column Pickers (Smart Defaults) ----
output$col_pickers <- renderUI({

  req(raw_df())

  cols <- names(raw_df())
  tool <- input$software_tool %||% "mzmine"

  # Smart defaults for m/z
  mz_cand <- switch(
    tool,

    xcms = c(
      "mzmed",
      "mz",
      "m/z",
      "mzmin",
      "mzmax"
    ),

    msdial = c(
      "average mz",
      "averagemz",
      "mz"
    ),

    default = c(
      "mz",
      "m/z",
      "mass",
      "average mz"
    ),

    c(
      "row m/z",
      "row mz",
      "mz"
    )
  )

  # Smart defaults for RT
  rt_cand <- switch(
    tool,

    xcms = c(
      "rtmed",
      "rt",
      "rtmin"
    ),

    msdial = c(
      "average rt(min)",
      "average rt",
      "averagertmin",
      "rt"
    ),

    default = c(
      "rt",
      "retention time",
      "time"
    ),

    c(
      "row retention time",
      "row rt",
      "rt"
    )
  )

  def_mz <- guess_col(
    cols,
    mz_cand
  )

  def_rt <- guess_col(
    cols,
    rt_cand
  )

  ann_id_cand <- switch(
  tool,

  xcms = c(
    "...1",
    "X...1",
    "feature_id",
    "feature",
    "id"
  ),

  msdial = c(
    "alignment id",
    "alignmentid",
    "spot id"
  ),

  default = c(
    "row id",
    "feature_id",
    "feature id",
    "id"
  ),

  c(
    "row id",
    "id",
    "feature_id"
  )
)

def_ann_id <- guess_col(
  cols,
  ann_id_cand
)

if (is.null(def_ann_id)) {
  def_ann_id <- cols[1]
}
  
  choices_mzrt <- c(
    "None",
    cols
  )

  tagList(

    selectInput(
      "feature_id_source",
      "Feature ID:",
      choices = c(
        "Combine m/z and RT" = "combine_mz_rt",
        "Auto-generate (feat_1)" = "auto",
        stats::setNames(cols, cols)
      ),
      selected = "combine_mz_rt"
    ),
    
    selectInput(
  "annotation_id_col",
  "Annotation matching ID:",
  choices = cols,
  selected = def_ann_id
),

    fluidRow(

      column(
        6,
        selectInput(
          "mz_col",
          "m/z column:",
          choices = choices_mzrt,
          selected = def_mz %||% "None"
        )
      ),

      column(
        6,
        selectInput(
          "rt_col",
          "RT column:",
          choices = choices_mzrt,
          selected = def_rt %||% "None"
        )
      )
    ),

    conditionalPanel(
      condition = "input.feature_id_source == 'combine_mz_rt'",

      textInput(
        "mz_rt_sep",
        "m/z–RT separator:",
        value = "@"
      )
    )
  )
})

  output$manual_sample_cols_ui <- renderUI({

  req(raw_df())

  selectizeInput(
    "sample_cols_manual",
    "Pick sample columns:",
    choices = names(raw_df()),
    selected = NULL,
    multiple = TRUE,
    options = list(
      placeholder = "Select first and last sample columns"
    )
  )
})
  
  observeEvent(input$sample_cols_manual, {

  req(raw_df())

  sel <- input$sample_cols_manual

  if (is.null(sel) || length(sel) < 2)
    return()

  cols <- names(raw_df())

  pos <- match(sel, cols)
  pos <- pos[!is.na(pos)]

  if (length(pos) < 2)
    return()

  # Select everything between leftmost and rightmost selection
  range_cols <- cols[min(pos):max(pos)]

  if (!setequal(sel, range_cols)) {

    updateSelectizeInput(
      session,
      "sample_cols_manual",
      selected = range_cols
    )
  }

}, ignoreInit = TRUE)
  
  output$raw_header <- renderUI({
  req(raw_df())

  h3(
    sprintf(
      "Raw dataset: %d Features × %d Samples",
      nrow(raw_df()),
      length(sample_cols_selected())
    )
  )
})

  output$raw_preview <- renderDT({
    req(raw_df())
    datatable(head(raw_df(), 20), options = list(scrollX = TRUE, pageLength = 8))
  })

  sample_cols_selected <- reactive({

  req(
  raw_df(),
  input$feature_id_source,
  input$mz_col,
  input$rt_col
)

  # Remove an old mapping error when inputs are reevaluated
  upload_error(NULL)

  df <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cols <- names(df)

  mode <- input$sample_mode %||% "kws"

  feature_id_col <- if (
  !is.null(input$feature_id_source) &&
  !input$feature_id_source %in%
    c(
      "combine_mz_rt",
      "auto"
    )
) {
  input$feature_id_source
} else {
  character(0)
}

meta <- unique(
  c(
    feature_id_col,
    input$annotation_id_col,
    input$mz_col,
    input$rt_col,
    "None"
  )
)

  meta <- meta[
    !is.na(meta) &
      nzchar(meta)
  ]


  sample_error <- function(msg) {

    upload_error(msg)

    showNotification(
      msg,
      type = "error",
      duration = 8
    )

    validate(
      need(FALSE, msg)
    )
  }


  # ---------------------
  # MANUAL
  # ---------------------
  if (mode == "manual") {

    selected <- input$sample_cols_manual %||%
      character(0)

    if (!length(selected)) {

      sample_error(
        "Pick at least one sample column."
      )
    }

    sc <- intersect(
      selected,
      cols
    )

    sc <- setdiff(
      sc,
      meta
    )

    if (!length(sc)) {

      sample_error(
        "Selected sample columns were not found or contain only Feature ID, m/z, or RT columns."
      )
    }

    return(sc)
  }


  # ---------------------
  # KEYWORDS
  # ---------------------
  if (mode == "kws") {

    kws <- input$sample_keywords %||%
      character(0)

    kws <- as.character(kws)

    kws <- kws[
      !is.na(kws) &
        nzchar(kws)
    ]

    if (!length(kws)) {

      sample_error(
        "Add at least one sample-column keyword."
      )
    }

    idx <- multi_sample_idx(
      cols,
      kws
    )

    if (!length(idx)) {

      sample_error(
        paste0(
          "No sample columns matched the keywords: ",
          paste(kws, collapse = ", ")
        )
      )
    }

    sc <- cols[idx]

    sc <- setdiff(
      sc,
      meta
    )

    if (!length(sc)) {

      sample_error(
        "Keyword matches contain only Feature ID, m/z, or RT columns. Use Manual or Auto."
      )
    }

    return(sc)
  }


  # ---------------------
  # AUTO
  # ---------------------
  cand <- setdiff(
    cols,
    meta
  )

  cand <- cand[
    !grepl(
      "^row\\b",
      cand,
      ignore.case = TRUE
    )
  ]

  if (!length(cand)) {

    sample_error(
      "No candidate sample columns were found."
    )
  }

  prop_num <- vapply(
    df[cand],
    function(x) {

      x2 <- suppressWarnings(
        as.numeric(
          as.character(x)
        )
      )

      mean(
        is.finite(x2),
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  sc <- cand[
    prop_num >= 0.7
  ]

  if (!length(sc)) {

    sample_error(
      "Auto-detect found no numeric sample columns. Switch to Manual or Keywords."
    )
  }

  sc
})
  
  output$peak_extra_cols_ui <- renderUI({

  req(raw_df())

  df <- raw_df()
  cols <- names(df)

  # Identify sample-intensity columns so they are not offered
  # as additional feature metadata columns.
  sample_cols <- sample_cols_selected()

  # Exclude columns already used as essential feature metadata
  feature_id_col <- if (
  !is.null(input$feature_id_source) &&
  !input$feature_id_source %in%
    c(
      "combine_mz_rt",
      "auto"
    )
) {
  input$feature_id_source
} else {
  character(0)
}

excluded_cols <- unique(
  c(
    feature_id_col,
    input$annotation_id_col %||% character(0),
    input$mz_col %||% character(0),
    input$rt_col %||% character(0),
    sample_cols
  )
)

  choices <- setdiff(
    cols,
    excluded_cols
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional non-sample peak-table columns were detected."
      )
    )
  }

  current_selection <- isolate(
    input$peak_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "peak_extra_cols",
    label = "Peak-table columns to include:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more peak-table columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} peak-table column(s) selected",
      `style` = "btn-success"
    )
  )
})
  
  # ---- Build matrix + fmap ----
built <- reactive({

  req(
  raw_df(),
  input$feature_id_source,
  input$annotation_id_col,
  input$mz_col,
  input$rt_col
)

  validate(

    need(
      !identical(
        input$mz_col,
        "None"
      ) &&
        input$mz_col %in%
        names(raw_df()),
      "Select a valid m/z column."
    ),

    need(
      !identical(
        input$rt_col,
        "None"
      ) &&
        input$rt_col %in%
        names(raw_df()),
      "Select a valid RT column."
    )
  )

  df <- raw_df()

  sc <- sample_cols_selected()

  parse_feature_table_to_matrix(
  raw_df = df,
  feature_id_source =
    input$feature_id_source %||%
    "combine_mz_rt",
  annotation_id_col = input$annotation_id_col,
  mz_col = input$mz_col,
  rt_col = input$rt_col,
  sample_cols = sc,
  mz_rt_sep =
    input$mz_rt_sep %||%
    "@"
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

  labels_from_sample_names_or_raw(
    sample_names(),
    token_sep = input$token_sep %||% "_",
    token_index = input$token_index %||% 2,
    clean_names = FALSE
  )
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

  output$paired_sample_preview <- renderTable({

  req(
    isTRUE(input$paired),
    input$test_type %in% c("Student", "Wilcoxon")
  )

  labs <- as.character(labels_vec())
  samples <- sample_names()

  validate(
    need(
      length(labs) == length(samples),
      "Sample names and labels do not match."
    )
  )

  if (identical(input$comparison_mode, "manual")) {

    comparisons <- selected_manual_comparisons()

  } else {

    req(input$ref_group)

    others <- setdiff(
      sort(unique(labs)),
      input$ref_group
    )

    comparisons <- data.frame(
      Group_num = rep(input$ref_group, length(others)),
      Group_den = others,
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(seq_len(nrow(comparisons)), function(i) {

    group_a <- comparisons$Group_num[i]
    group_b <- comparisons$Group_den[i]

    a <- samples[which(labs == group_a)]
    b <- samples[which(labs == group_b)]

    # Padding makes an unmatched sample visible.
    n <- max(length(a), length(b))

    data.frame(
      Comparison = rep(paste(group_a, "/", group_b), n),
      Pair = seq_len(n),
      Numerator_sample = a[seq_len(n)],
      Denominator_sample = b[seq_len(n)],
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)

}, striped = TRUE, bordered = TRUE, na = "UNMATCHED")

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

  ext <- tolower(
    tools::file_ext(
      input$file_sirius$name
    )
  )

  validate(
    need(
      ext %in% c("csv", "tsv", "txt"),
      "SIRIUS file must be .csv, .tsv, or .txt."
    )
  )

  delim <- if (identical(ext, "csv")) {
    ","
  } else {
    "\t"
  }

  as.data.frame(
    vroom::vroom(
      input$file_sirius$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})

  output$sirius_pickers <- renderUI({

  req(sirius_df())

  cols <- names(
    sirius_df()
  )

  tagList(

    selectInput(
      "sirius_idcol",
      "SIRIUS Feature ID column:",
      choices = cols,
      selected = if (
        "mappingFeatureId" %in% cols
      ) {
        "mappingFeatureId"
      } else {
        cols[1]
      }
    ),

    selectInput(
      "sirius_npcol",
      "NPC column:",
      choices = cols,
      selected = if (
        "NPC#class" %in% cols
      ) {
        "NPC#class"
      } else {
        cols[1]
      }
    ),

    selectInput(
      "sirius_cfcol",
      "ClassyFire column:",
      choices = cols,
      selected = if (
        "ClassyFire#class" %in% cols
      ) {
        "ClassyFire#class"
      } else {
        cols[1]
      }
    ),

    materialSwitch(
      inputId = "use_sirius_extra_cols",
      label = "Add additional SIRIUS columns",
      value = FALSE,
      status = "success",
      width = "auto"
    ),

    conditionalPanel(
      condition = "input.use_sirius_extra_cols == true",

      uiOutput(
        "sirius_extra_cols_ui"
      )
    )
  )
})
  
  output$sirius_extra_cols_ui <- renderUI({

  req(
    sirius_df(),
    input$sirius_idcol,
    input$sirius_npcol,
    input$sirius_cfcol
  )

  cols <- names(
    sirius_df()
  )

  choices <- setdiff(
    cols,
    c(
      input$sirius_idcol,
      input$sirius_npcol,
      input$sirius_cfcol
    )
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional SIRIUS columns are available."
      )
    )
  }

  current_selection <- isolate(
    input$sirius_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "sirius_extra_cols",
    label = "Additional SIRIUS columns:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more SIRIUS columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} SIRIUS column(s) selected",
      `style` = "btn-success"
    )
  )
})
  
gnps_annotation_df <- reactive({

  req(input$use_gnps_annotation)
  req(input$file_gnps_annotation)

  ext <- tolower(
    tools::file_ext(
      input$file_gnps_annotation$name
    )
  )

  validate(
    need(
      ext %in% c("tsv", "txt", "csv"),
      "GNPS annotation file must be .tsv, .txt, or .csv."
    )
  )

  delim <- if (identical(ext, "csv")) {
    ","
  } else {
    "\t"
  }

  as.data.frame(
    vroom::vroom(
      input$file_gnps_annotation$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})

output$gnps_annotation_pickers <- renderUI({

  req(gnps_annotation_df())

  cols <- names(
    gnps_annotation_df()
  )

  validate(
    need(
      length(cols) > 0,
      "No columns were detected in the GNPS file."
    )
  )

  # Default GNPS ID column
  default_id <- guess_col(
    cols,
    c(
      "#Scan#",
      "Scan",
      "scan",
      "ClusterIndex",
      "Cluster ID",
      "row ID",
      "id"
    )
  ) %||% cols[1]

  # Default primary annotation column
  default_annotation <- guess_col(
    cols,
    c(
      "Compound_Name",
      "Compound_name",
      "Compound name",
      "CompoundName",
      "Library compound name",
      "Annotation",
      "Name"
    )
  ) %||% cols[1]

  tagList(

    selectInput(
      "gnps_annotation_idcol",
      "GNPS ID column:",
      choices = cols,
      selected = default_id
    ),

    selectInput(
      "gnps_annotation_col",
      "Primary GNPS annotation column:",
      choices = cols,
      selected = default_annotation
    ),

    materialSwitch(
      inputId = "use_gnps_extra_cols",
      label = "Add additional GNPS columns",
      value = FALSE,
      status = "success",
      width = "auto"
    ),

    conditionalPanel(
      condition = "input.use_gnps_extra_cols == true",

      uiOutput(
        "gnps_extra_cols_ui"
      )
    )
  )
})
  
output$gnps_extra_cols_ui <- renderUI({

  req(
    gnps_annotation_df(),
    input$gnps_annotation_idcol,
    input$gnps_annotation_col
  )

  cols <- names(
    gnps_annotation_df()
  )

  # Do not offer the join ID or primary annotation again
  choices <- setdiff(
    cols,
    c(
      input$gnps_annotation_idcol,
      input$gnps_annotation_col
    )
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional GNPS columns are available."
      )
    )
  }

  current_selection <- isolate(
    input$gnps_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "gnps_extra_cols",
    label = "Additional GNPS columns:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more GNPS columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} GNPS column(s) selected",
      `style` = "btn-success"
    )
  )
})

# ============================================================
# Other annotation source
# ============================================================

other_annotation_df <- reactive({

  req(input$use_other_annotation)
  req(input$file_other_annotation)

  ext <- tolower(
    tools::file_ext(
      input$file_other_annotation$name
    )
  )

  validate(
    need(
      ext %in% c("csv", "tsv", "txt"),
      "Other annotation file must be .csv, .tsv, or .txt."
    )
  )

  delim <- if (identical(ext, "csv")) {
    ","
  } else {
    "\t"
  }

  as.data.frame(
    vroom::vroom(
      input$file_other_annotation$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})


output$other_annotation_pickers <- renderUI({

  req(
    raw_df(),
    other_annotation_df()
  )

  peak_cols <- names(raw_df())
  ann_cols  <- names(other_annotation_df())

  validate(
    need(
      length(peak_cols) > 0,
      "No peak-table columns were detected."
    ),
    need(
      length(ann_cols) > 0,
      "No columns were detected in the annotation file."
    )
  )

  # ---------------------------------
  # Default peak-table matching ID
  # ---------------------------------
  default_peak_id <- guess_col(
    peak_cols,
    c(
      "row ID",
      "row id",
      "Row ID",
      "alignment id",
      "Alignment ID",
      "feature_id",
      "feature id",
      "id"
    )
  ) %||% peak_cols[1]


  # ---------------------------------
  # Default annotation-file ID
  # ---------------------------------
  default_ann_id <- guess_col(
    ann_cols,
    c(
      "row ID",
      "row id",
      "Row ID",
      "feature_id",
      "feature id",
      "mappingFeatureId",
      "#Scan#",
      "Scan",
      "id"
    )
  ) %||% ann_cols[1]


  # ---------------------------------
  # Default annotation column
  # ---------------------------------
  default_annotation <- guess_col(
    ann_cols,
    c(
      "Annotation",
      "annotation",
      "Compound_Name",
      "Compound_name",
      "Compound name",
      "CompoundName",
      "Name",
      "Identification",
      "Compound"
    )
  ) %||% ann_cols[
    min(
      2,
      length(ann_cols)
    )
  ]


  tagList(

    selectInput(
      "other_peak_id_col",
      "Peak-table ID column:",
      choices = peak_cols,
      selected = default_peak_id
    ),

    selectInput(
      "other_annotation_idcol",
      "Annotation-file ID column:",
      choices = ann_cols,
      selected = default_ann_id
    ),

    selectInput(
      "other_annotation_col",
      "Primary annotation column:",
      choices = ann_cols,
      selected = default_annotation
    ),

    materialSwitch(
      inputId = "use_other_extra_cols",
      label = "Add additional annotation columns",
      value = FALSE,
      status = "success",
      width = "auto"
    ),

    conditionalPanel(
      condition = "input.use_other_extra_cols == true",

      uiOutput(
        "other_extra_cols_ui"
      )
    )
  )
})

# ============================================================
# Annotation table for direct download
# Does NOT require Run preprocessing
# ============================================================

annotation_export_table <- reactive({

  req(
    raw_df(),
    input$feature_id_source,
    input$annotation_id_col,
    input$mz_col,
    input$rt_col
  )

  raw_peak <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(
      input$annotation_id_col %in% names(raw_peak),
      "Annotation matching ID column was not found."
    ),
    need(
      !identical(input$mz_col, "None") &&
        input$mz_col %in% names(raw_peak),
      "Select a valid m/z column."
    ),
    need(
      !identical(input$rt_col, "None") &&
        input$rt_col %in% names(raw_peak),
      "Select a valid RT column."
    )
  )


  # ==========================================================
  # Basic feature information
  # ==========================================================

  mz <- suppressWarnings(
    as.numeric(
      raw_peak[[input$mz_col]]
    )
  )

  rt <- suppressWarnings(
    as.numeric(
      raw_peak[[input$rt_col]]
    )
  )


  # Same Feature ID logic used by the main app
  if (
    identical(
      input$feature_id_source,
      "combine_mz_rt"
    )
  ) {

    sep <- input$mz_rt_sep %||% "@"

    if (!nzchar(sep)) {
      sep <- "@"
    }

    mz_text <- ifelse(
      is.na(mz),
      "NA",
      as.character(round(mz, 4))
    )

    rt_text <- ifelse(
      is.na(rt),
      "NA",
      as.character(round(rt, 2))
    )

    Feature <- paste0(
      mz_text,
      sep,
      rt_text
    )

  } else if (
    identical(
      input$feature_id_source,
      "auto"
    )
  ) {

    Feature <- paste0(
      "feat_",
      seq_len(nrow(raw_peak))
    )

  } else {

    validate(
      need(
        input$feature_id_source %in% names(raw_peak),
        "Selected Feature ID column was not found."
      )
    )

    Feature <- trimws(
      as.character(
        raw_peak[[input$feature_id_source]]
      )
    )

    bad_id <- is.na(Feature) |
      !nzchar(Feature)

    if (any(bad_id)) {
      Feature[bad_id] <- paste0(
        "feat_",
        which(bad_id)
      )
    }
  }

  Feature <- make.unique(
    Feature,
    sep = "_"
  )


  annotation_id <- trimws(
    as.character(
      raw_peak[[input$annotation_id_col]]
    )
  )


  out <- data.frame(
    Feature_ID = Feature,
    Annotation_Feature_ID = annotation_id,
    mz = round(mz, 6),
    RT = round(rt, 4),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


  # Internal joining columns
  out$.Feature_internal <- Feature
  out$.annotation_id_internal <- annotation_id


  # ==========================================================
  # Selected additional peak-table columns
  # ==========================================================

  if (isTRUE(input$use_peak_extra_cols)) {

    selected_peak_cols <- intersect(
      input$peak_extra_cols %||% character(0),
      names(raw_peak)
    )

    if (length(selected_peak_cols)) {

      peak_colmap <- make_prefixed_colmap(
        selected_peak_cols,
        prefix = "Peak_"
      )

      for (original_name in names(peak_colmap)) {

        output_name <- peak_colmap[[original_name]]

        out[[output_name]] <-
          raw_peak[[original_name]]
      }
    }
  }


  # ==========================================================
  # SIRIUS
  # SIRIUS is matched against Annotation Feature ID
  # ==========================================================

  if (
    isTRUE(input$use_sirius) &&
    !is.null(input$file_sirius) &&
    !is.null(input$sirius_idcol) &&
    !is.null(input$sirius_npcol) &&
    !is.null(input$sirius_cfcol)
  ) {

    s <- as.data.frame(
      sirius_df(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    validate(
      need(
        input$sirius_idcol %in% names(s),
        "Selected SIRIUS Feature ID column was not found."
      ),
      need(
        input$sirius_npcol %in% names(s),
        "Selected SIRIUS NPC column was not found."
      ),
      need(
        input$sirius_cfcol %in% names(s),
        "Selected SIRIUS ClassyFire column was not found."
      )
    )


    selected_sirius_extra <- character(0)

    if (isTRUE(input$use_sirius_extra_cols)) {

      selected_sirius_extra <- intersect(
        input$sirius_extra_cols %||% character(0),
        names(s)
      )

      selected_sirius_extra <- setdiff(
        selected_sirius_extra,
        c(
          input$sirius_idcol,
          input$sirius_npcol,
          input$sirius_cfcol
        )
      )
    }

    selected_sirius_extra <- unique(c(
  selected_sirius_extra,
  grep(
    "probability",
    names(s),
    ignore.case = TRUE,
    value = TRUE
  )
))
    
    sirius_colmap <- make_prefixed_colmap(
      selected_sirius_extra,
      prefix = "SIRIUS_"
    )


    ss <- data.frame(
      .annotation_id_internal = trimws(
        as.character(
          s[[input$sirius_idcol]]
        )
      ),

      `NPC#class` = clean_missing_text(
        s[[input$sirius_npcol]]
      ),

      `ClassyFire#class` = clean_missing_text(
        s[[input$sirius_cfcol]]
      ),

      check.names = FALSE,
      stringsAsFactors = FALSE
    )


    if (length(sirius_colmap)) {

      for (original_name in names(sirius_colmap)) {

        output_name <-
          sirius_colmap[[original_name]]

        value <- s[[original_name]]

        if (
          is.character(value) ||
          is.factor(value)
        ) {
          value <- clean_missing_text(value)
        }

        ss[[output_name]] <- value
      }
    }


    ss <- ss %>%
      dplyr::filter(
        !is.na(.annotation_id_internal),
        nzchar(.annotation_id_internal)
      ) %>%
      dplyr::distinct(
        .annotation_id_internal,
        .keep_all = TRUE
      )


    out <- out %>%
      dplyr::left_join(
        ss,
        by = ".annotation_id_internal"
      )
  }


  # ==========================================================
  # GNPS
  # GNPS is also matched against Annotation Feature ID
  # ==========================================================

  if (
    isTRUE(input$use_gnps_annotation) &&
    !is.null(input$file_gnps_annotation) &&
    !is.null(input$gnps_annotation_idcol) &&
    !is.null(input$gnps_annotation_col)
  ) {

    g <- as.data.frame(
      gnps_annotation_df(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )


    validate(
      need(
        input$gnps_annotation_idcol %in% names(g),
        "Selected GNPS ID column was not found."
      ),
      need(
        input$gnps_annotation_col %in% names(g),
        "Selected GNPS annotation column was not found."
      )
    )


    gnps_primary <- data.frame(
      .annotation_id_internal = trimws(
        as.character(
          g[[input$gnps_annotation_idcol]]
        )
      ),

      GNPS_annotation = clean_missing_text(
        g[[input$gnps_annotation_col]]
      ),

      check.names = FALSE,
      stringsAsFactors = FALSE
    ) %>%

      dplyr::filter(
        !is.na(.annotation_id_internal),
        nzchar(.annotation_id_internal)
      ) %>%

      dplyr::group_by(
        .annotation_id_internal
      ) %>%

      dplyr::summarise(

        GNPS_annotation = {

          values <- unique(
            GNPS_annotation[
              !is.na(GNPS_annotation)
            ]
          )

          if (length(values)) {
            paste(
              values,
              collapse = " | "
            )
          } else {
            NA_character_
          }
        },

        .groups = "drop"
      )


    # Additional GNPS columns
    selected_gnps_extra <- character(0)

    if (isTRUE(input$use_gnps_extra_cols)) {

      selected_gnps_extra <- intersect(
        input$gnps_extra_cols %||% character(0),
        names(g)
      )

      selected_gnps_extra <- setdiff(
        selected_gnps_extra,
        c(
          input$gnps_annotation_idcol,
          input$gnps_annotation_col
        )
      )
    }


    gnps_colmap <- make_prefixed_colmap(
      selected_gnps_extra,
      prefix = "GNPS_"
    )


    gnps_join <- gnps_primary


    if (length(gnps_colmap)) {

      gnps_extra <- data.frame(
        .annotation_id_internal = trimws(
          as.character(
            g[[input$gnps_annotation_idcol]]
          )
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )


      for (original_name in names(gnps_colmap)) {

        output_name <-
          gnps_colmap[[original_name]]

        value <- g[[original_name]]

        if (
          is.character(value) ||
          is.factor(value)
        ) {
          value <- clean_missing_text(value)
        }

        gnps_extra[[output_name]] <- value
      }


      gnps_extra <- gnps_extra %>%
        dplyr::filter(
          !is.na(.annotation_id_internal),
          nzchar(.annotation_id_internal)
        ) %>%
        dplyr::distinct(
          .annotation_id_internal,
          .keep_all = TRUE
        )


      gnps_join <- gnps_join %>%
        dplyr::left_join(
          gnps_extra,
          by = ".annotation_id_internal"
        )
    }


    out <- out %>%
      dplyr::left_join(
        gnps_join,
        by = ".annotation_id_internal"
      )
  }


  # ==========================================================
  # OTHER ANNOTATION SOURCE
  # Can use an arbitrary peak-table ID column
  # ==========================================================

  if (
    isTRUE(input$use_other_annotation) &&
    !is.null(input$file_other_annotation) &&
    !is.null(input$other_peak_id_col) &&
    !is.null(input$other_annotation_idcol) &&
    !is.null(input$other_annotation_col)
  ) {

    ann <- as.data.frame(
      other_annotation_df(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )


    validate(
      need(
        input$other_peak_id_col %in%
          names(raw_peak),
        "Selected peak-table ID column for Other Annotation Source was not found."
      ),
      need(
        input$other_annotation_idcol %in%
          names(ann),
        "Selected annotation-file ID column was not found."
      ),
      need(
        input$other_annotation_col %in%
          names(ann),
        "Selected primary annotation column was not found."
      )
    )


    peak_other_map <- data.frame(
      .Feature_internal = Feature,

      .other_join_id = trimws(
        as.character(
          raw_peak[[input$other_peak_id_col]]
        )
      ),

      check.names = FALSE,
      stringsAsFactors = FALSE
    )


    other_primary <- data.frame(
      .other_join_id = trimws(
        as.character(
          ann[[input$other_annotation_idcol]]
        )
      ),

      Other_annotation = clean_missing_text(
        ann[[input$other_annotation_col]]
      ),

      check.names = FALSE,
      stringsAsFactors = FALSE
    ) %>%

      dplyr::filter(
        !is.na(.other_join_id),
        nzchar(.other_join_id)
      ) %>%

      dplyr::group_by(
        .other_join_id
      ) %>%

      dplyr::summarise(

        Other_annotation = {

          values <- unique(
            Other_annotation[
              !is.na(Other_annotation)
            ]
          )

          if (length(values)) {
            paste(
              values,
              collapse = " | "
            )
          } else {
            NA_character_
          }
        },

        .groups = "drop"
      )


    selected_other_extra <- character(0)

    if (isTRUE(input$use_other_extra_cols)) {

      selected_other_extra <- intersect(
        input$other_extra_cols %||%
          character(0),
        names(ann)
      )

      selected_other_extra <- setdiff(
        selected_other_extra,
        c(
          input$other_annotation_idcol,
          input$other_annotation_col
        )
      )
    }


    other_colmap <- make_prefixed_colmap(
      selected_other_extra,
      prefix = "Other_"
    )


    other_join <- other_primary


    if (length(other_colmap)) {

      other_extra <- data.frame(

        .other_join_id = trimws(
          as.character(
            ann[[input$other_annotation_idcol]]
          )
        ),

        check.names = FALSE,
        stringsAsFactors = FALSE
      )


      for (original_name in names(other_colmap)) {

        output_name <-
          other_colmap[[original_name]]

        value <- ann[[original_name]]

        if (
          is.character(value) ||
          is.factor(value)
        ) {
          value <- clean_missing_text(value)
        }

        other_extra[[output_name]] <- value
      }


      other_extra <- other_extra %>%
        dplyr::filter(
          !is.na(.other_join_id),
          nzchar(.other_join_id)
        ) %>%
        dplyr::distinct(
          .other_join_id,
          .keep_all = TRUE
        )


      other_join <- other_join %>%
        dplyr::left_join(
          other_extra,
          by = ".other_join_id"
        )
    }


    out <- out %>%

      dplyr::left_join(
        peak_other_map,
        by = ".Feature_internal"
      ) %>%

      dplyr::left_join(
        other_join,
        by = ".other_join_id"
      ) %>%

      dplyr::select(
        -dplyr::any_of(
          ".other_join_id"
        )
      )
  }


  # Remove internal helper columns
  out <- out %>%
    dplyr::select(
      -dplyr::any_of(
        c(
          ".Feature_internal",
          ".annotation_id_internal"
        )
      )
    )


  out
})

output$other_extra_cols_ui <- renderUI({

  req(
    other_annotation_df(),
    input$other_annotation_idcol,
    input$other_annotation_col
  )

  cols <- names(
    other_annotation_df()
  )

  choices <- setdiff(
    cols,
    c(
      input$other_annotation_idcol,
      input$other_annotation_col
    )
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional annotation columns are available."
      )
    )
  }

  current_selection <- isolate(
    input$other_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "other_extra_cols",
    label = "Additional annotation columns:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more annotation columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} annotation column(s) selected",
      `style` = "btn-success"
    )
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

materialSwitch(
  inputId = "prune_sirius_to_peak",
  label = "Restrict SIRIUS statistics to uploaded peak table",
  value = TRUE,
  status = "success",
  width = "auto"
),

uiOutput(
  "sirius_pruning_notice"
),

uiOutput(
  "sirius_peak_id_picker"
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
    uiOutput("network_component_picker"),
    
    div(
      class = "small-note",
      "GNPS pairs are converted from ClusterID1/ClusterID2 into one ClusterID column. By default, no peak-table pruning is applied. If a peak-table ID column is selected, SIRIUS IDs and GNPS ClusterIDs are restricted to IDs present in the uploaded peak table."
    )
  )
})

output$sirius_pruning_notice <- renderUI({

  if (!isTRUE(input$prune_sirius_to_peak)) {
    return(NULL)
  }

  if (!is.null(input$file_data)) {
    return(NULL)
  }

})

output$sirius_peak_id_picker <- renderUI({

    if (!isTRUE(input$prune_sirius_to_peak)) {
    return(NULL)
  }

  if (is.null(input$file_data)) {
    return(NULL)
  }
  
  req(raw_df())

  raw_cols <- names(
    raw_df()
  )

  validate(
    need(
      length(raw_cols) > 0,
      "No peak-table columns were detected."
    )
  )

selected_annotation_id <-
  input$annotation_id_col %||%
  ""

default_peak_id <- if (
  selected_annotation_id %in% raw_cols
) {

  selected_annotation_id

} else {

  guess_col_ci(
    raw_cols,
    c(
      "row ID",
      "row id",
      "id",
      "feature_id",
      "feature id"
    ),
    default = raw_cols[1]
  )
}

pickerInput(
  inputId = "sirius_peak_id_col",
  label = "Peak-table column used for SIRIUS pruning:",

  choices = c(
    "Auto-generated feature ID (row number)" = "__auto__",
    raw_cols
  ),

  selected = default_peak_id,
  multiple = FALSE,

  options = list(
    `live-search` = TRUE,
    `size` = 10,
    `style` = "btn-success",
    `none-selected-text` = "Choose a peak-table ID column"
  )
)
})

output$gnps_pickers_ui <- renderUI({

  req(
    gnps_pairs_df()
  )

  g_cols <- names(
    gnps_pairs_df()
  )

  tagList(

    selectInput(
      "gnps_cluster1_col",
      "GNPS ClusterID1 column:",
      choices = g_cols,
      selected = guess_col_ci(
        g_cols,
        c(
          "ClusterID1",
          "CLUSTERID1"
        ),
        g_cols[1]
      )
    ),

    selectInput(
      "gnps_cluster2_col",
      "GNPS ClusterID2 column:",
      choices = g_cols,
      selected = guess_col_ci(
        g_cols,
        c(
          "ClusterID2",
          "CLUSTERID2"
        ),
        g_cols[
          min(
            2,
            length(g_cols)
          )
        ]
      )
    ),

    selectInput(
      "gnps_component_col",
      "GNPS ComponentIndex column:",
      choices = g_cols,
      selected = guess_col_ci(
        g_cols,
        c(
          "ComponentIndex",
          "component"
        ),
        g_cols[1]
      )
    )
  )
})

sirius_stats_data <- reactive({
  req(sirius_df(), input$stats_sirius_id_col, input$stats_sirius_col)

  validate(
    need(
      !isTRUE(input$prune_sirius_to_peak) ||
        !is.null(input$file_data),
      paste0(
        "SIRIUS pruning is enabled, but no peak table is uploaded. ",
        "Upload a peak table or switch off SIRIUS pruning."
      )
    )
  )
  
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

  keep_ids <- selected_peak_ids_for_sirius()

  if (!is.null(keep_ids)) {
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

selected_peak_ids_for_sirius <- reactive({

  # Switch off means: use all SIRIUS rows.
  if (!isTRUE(input$prune_sirius_to_peak)) {
    return(NULL)
  }

  req(
    raw_df(),
    input$sirius_peak_id_col
  )

  raw <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  selected_col <- as.character(
    input$sirius_peak_id_col
  )

  ids <- if (
    identical(
      selected_col,
      "__auto__"
    )
  ) {

    as.character(
      seq_len(
        nrow(raw)
      )
    )

  } else {

    validate(
      need(
        selected_col %in% names(raw),
        "Selected peak-table feature ID column was not found."
      )
    )

    trimws(
      as.character(
        raw[[selected_col]]
      )
    )
  }

  ids <- ids[
    !is.na(ids) &
      nzchar(ids)
  ]

  unique(
    ids
  )
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

  out
})

gnps_component_edges <- reactive({

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
    need(
      input$gnps_cluster1_col %in% names(g),
      "GNPS ClusterID1 column was not found."
    ),
    need(
      input$gnps_cluster2_col %in% names(g),
      "GNPS ClusterID2 column was not found."
    ),
    need(
      input$gnps_component_col %in% names(g),
      "GNPS ComponentIndex column was not found."
    )
  )

  edges <- g %>%
    dplyr::transmute(

      ComponentIndex = trimws(
        as.character(
          .data[[input$gnps_component_col]]
        )
      ),

      ClusterID1 = trimws(
        as.character(
          .data[[input$gnps_cluster1_col]]
        )
      ),

      ClusterID2 = trimws(
        as.character(
          .data[[input$gnps_cluster2_col]]
        )
      )
    ) %>%

    dplyr::filter(
      !is.na(ComponentIndex),
      !is.na(ClusterID1),
      !is.na(ClusterID2),
      nzchar(ComponentIndex),
      nzchar(ClusterID1),
      nzchar(ClusterID2),
      ClusterID1 != ClusterID2
    ) %>%

    dplyr::distinct(
      ComponentIndex,
      ClusterID1,
      ClusterID2
    )

  # Optional pruning using the selected peak-table ID column

  edges
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

output$network_component_picker <- renderUI({

  req(
    input$file_gnps_pairs,
    input$stats_selected_class,
    selected_class_component_stats(),
    gnps_component_edges()
  )

  class_stats <- selected_class_component_stats()

  # Remove rows that cannot be displayed as a network
  class_stats <- class_stats %>%
    dplyr::filter(
      !ComponentIndex %in% c(
        "No ComponentIndex match",
        "GNPS pairs not uploaded"
      )
    )

  if (!nrow(class_stats)) {

    return(
      tagList(
        tags$hr(),

        div(
          class = "small-note",
          style = "
            padding: 8px;
            border: 1px solid #dddddd;
            border-radius: 6px;
          ",
          paste0(
            "The selected SIRIUS value does not occur in a ",
            "displayable GNPS component."
          )
        )
      )
    )
  }

  edge_data <- gnps_component_edges()

  # Count total nodes in each available component
  component_totals <- edge_data %>%

    tidyr::pivot_longer(
      cols = c(
        "ClusterID1",
        "ClusterID2"
      ),
      values_to = "ClusterID"
    ) %>%

    dplyr::distinct(
      ComponentIndex,
      ClusterID
    ) %>%

    dplyr::count(
      ComponentIndex,
      name = "Total_nodes"
    )

  component_options <- class_stats %>%

    dplyr::left_join(
      component_totals,
      by = "ComponentIndex"
    ) %>%

    # Exclude components with no surviving edge after pruning
    dplyr::filter(
      !is.na(Total_nodes),
      Total_nodes >= 2
    ) %>%

    dplyr::mutate(
      Display_name = paste0(
        "Component ",
        ComponentIndex,
        " — ",
        Points,
        " selected / ",
        Total_nodes,
        " total"
      )
    ) %>%

    dplyr::arrange(
      dplyr::desc(Points),
      dplyr::desc(Total_nodes),
      ComponentIndex
    )

  if (!nrow(component_options)) {

    return(
      tagList(
        tags$hr(),

        div(
          class = "small-note",
          paste0(
            "No GNPS component containing the selected class ",
"has displayable edges."
          )
        )
      )
    )
  }

  old_component <- isolate(
    input$network_component
  )

  selected_component <- if (
    !is.null(old_component) &&
    old_component %in% component_options$ComponentIndex
  ) {
    old_component
  } else {
    component_options$ComponentIndex[1]
  }

  tagList(

    tags$hr(),

    h3(
      class = "highlight",
      "Interactive component network"
    ),

    pickerInput(
      inputId = "network_component",
      label = "GNPS component to display:",

      choices = stats::setNames(
        component_options$ComponentIndex,
        component_options$Display_name
      ),

      selected = selected_component,
      multiple = FALSE,

      options = list(
        `live-search` = TRUE,
        `style` = "btn-success"
      )
    ),

    radioButtons(
      inputId = "network_label_mode",
      label = "Node labels:",
      choices = c(
        "Selected class only" = "selected",
        "All nodes" = "all",
        "No labels" = "none"
      ),
      selected = "selected"
    ),

    uiOutput(
      "network_volcano_options"
    )
  )
})

output$network_volcano_options <- renderUI({

  if (
    is.null(rv$volcano) ||
    !is.data.frame(rv$volcano) ||
    !nrow(rv$volcano) ||
    !"id" %in% names(rv$volcano)
  ) {

    return(
      div(
        class = "small-note",
        style = "margin-top: 8px;",
        paste0(
          "Run preprocessing to add statistical results from all ",
          "processed comparisons to network-node hover and click details."
        )
      )
    )
  }

  tagList(

    materialSwitch(
      inputId = "network_use_volcano",
      label = "Add all processed statistical comparisons to network nodes",
      value = TRUE,
      status = "success",
      width = "auto"
    ),

    div(
      class = "small-note",
      style = "margin-top: 4px;",
      paste0(
        "For each network node, FC, adjusted p-value, -log10(FDR), ",
        "overall mean, group means, test scale, and significance ",
        "are shown for every processed comparison."
      )
    )
  )
})

selected_component_network_data <- reactive({

  req(
    input$file_gnps_pairs,
    input$network_component,
    input$stats_selected_class,
    gnps_component_edges()
  )

  selected_component <- as.character(
    input$network_component
  )

  component_edges <- gnps_component_edges() %>%

    dplyr::filter(
      ComponentIndex == selected_component
    ) %>%

    dplyr::select(
      ClusterID1,
      ClusterID2
    ) %>%

    dplyr::distinct()

  validate(
    need(
      nrow(component_edges) > 0,
      "The selected component contains no displayable GNPS edges."
    )
  )

  # Build and simplify graph
  graph_object <- igraph::graph_from_data_frame(
    component_edges,
    directed = FALSE
  )

  graph_object <- igraph::simplify(
    graph_object,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )

  node_count <- igraph::vcount(
    graph_object
  )

  validate(
    need(
      node_count <= 500,
      paste0(
        "This component contains ",
        node_count,
        " nodes. Interactive display is limited to 500 nodes."
      )
    )
  )

  # Reproducible force-directed layout
  set.seed(1234)

  coordinates <- igraph::layout_with_fr(
    graph_object,
    niter = 500
  )

  node_names <- igraph::V(
    graph_object
  )$name

  node_data <- tibble::tibble(
    ClusterID = as.character(
      node_names
    ),

    x = as.numeric(
      coordinates[, 1]
    ),

    y = as.numeric(
      coordinates[, 2]
    ),

    Degree = as.integer(
      igraph::degree(
        graph_object
      )
    )
  )

  # Aggregate all values from the selected SIRIUS column
  sirius_node_annotations <- sirius_stats_data() %>%

    dplyr::group_by(
      SIRIUS_ID
    ) %>%

    dplyr::summarise(

      SIRIUS_values = paste(
        sort(
          unique(
            Annotation[
              !is.na(Annotation)
            ]
          )
        ),
        collapse = " | "
      ),

      Selected_class = any(
        Annotation ==
          input$stats_selected_class
      ),

      .groups = "drop"
    ) %>%

    dplyr::rename(
      ClusterID = SIRIUS_ID
    )

  node_data <- node_data %>%

    dplyr::left_join(
      sirius_node_annotations,
      by = "ClusterID"
    ) %>%

    dplyr::mutate(

      Selected_class = dplyr::coalesce(
        Selected_class,
        FALSE
      ),

      Has_SIRIUS = (
        !is.na(SIRIUS_values) &
          nzchar(SIRIUS_values)
      ),

      Node_type = dplyr::case_when(
        Selected_class ~ "Selected class",
        Has_SIRIUS ~ "Other SIRIUS annotation",
        TRUE ~ "No SIRIUS annotation"
      )
    )

  volcano_added <- FALSE

format_network_number <- function(
    x,
    digits = 4
) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  out <- rep(
    "NA",
    length(x)
  )

  ok <- is.finite(x)

  out[ok] <- format(
    signif(
      x[ok],
      digits
    ),
    scientific = FALSE,
    trim = TRUE
  )

  out
}


# Empty structure in case volcano statistics are unavailable
volcano_stats_long <- tibble::tibble(

  ClusterID = character(),
  Comparison = character(),

  Group_num = character(),
  Group_den = character(),

  FC = numeric(),
  Adj_p = numeric(),
  Adj_p_log = numeric(),

  Mean = numeric(),
  mean_num = numeric(),
  mean_den = numeric(),

  TestScale = character(),
  Significant_default = character(),

  GNPS_annotation = character(),

  Comparison_hover = character()
)


if (
  isTRUE(input$network_use_volcano) &&
  !is.null(rv$volcano) &&
  is.data.frame(rv$volcano) &&
  nrow(rv$volcano) &&
  "id" %in% names(rv$volcano)
) {

  volc <- as.data.frame(
    rv$volcano,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


  n_volc <- nrow(volc)


  get_chr <- function(column_name) {

    if (!column_name %in% names(volc)) {

      return(
        rep(
          NA_character_,
          n_volc
        )
      )
    }

    as.character(
      volc[[column_name]]
    )
  }


  get_num <- function(column_name) {

    if (!column_name %in% names(volc)) {

      return(
        rep(
          NA_real_,
          n_volc
        )
      )
    }

    suppressWarnings(
      as.numeric(
        volc[[column_name]]
      )
    )
  }


  comparison_values <- get_chr(
    "Groups"
  )

  comparison_values[
    is.na(comparison_values) |
      !nzchar(
        trimws(
          comparison_values
        )
      )
  ] <- "Comparison"


  gnps_values <- if (
    "GNPS_annotation" %in% names(volc)
  ) {

    clean_missing_text(
      volc$GNPS_annotation
    )

  } else {

    rep(
      NA_character_,
      n_volc
    )
  }


  volcano_stats_long <- tibble::tibble(

    ClusterID = trimws(
      as.character(
        volc$id
      )
    ),

    Comparison = comparison_values,

    Group_num = get_chr(
      "Group_num"
    ),

    Group_den = get_chr(
      "Group_den"
    ),

    FC = get_num(
      "FC"
    ),

    Adj_p = get_num(
      "Adj.p-value"
    ),

    Adj_p_log = get_num(
      "Adj.p-value.log"
    ),

    Mean = get_num(
      "Mean"
    ),

    mean_num = get_num(
      "mean_num"
    ),

    mean_den = get_num(
      "mean_den"
    ),

    TestScale = get_chr(
      "TestScale"
    ),

    Significant_default = get_chr(
      "Significant_default"
    ),

    GNPS_annotation = gnps_values
  ) %>%

    dplyr::filter(
      !is.na(ClusterID),
      nzchar(ClusterID)
    ) %>%

    dplyr::distinct(
      ClusterID,
      Comparison,
      .keep_all = TRUE
    )


  if (nrow(volcano_stats_long)) {

    group_num_label <- volcano_stats_long$Group_num
    group_den_label <- volcano_stats_long$Group_den


    group_num_label[
      is.na(group_num_label) |
        !nzchar(group_num_label)
    ] <- "Group 1"


    group_den_label[
      is.na(group_den_label) |
        !nzchar(group_den_label)
    ] <- "Group 2"


    test_scale_label <- volcano_stats_long$TestScale

    test_scale_label[
      is.na(test_scale_label) |
        !nzchar(test_scale_label)
    ] <- "NA"


    significant_label <-
      volcano_stats_long$Significant_default

    significant_label[
      is.na(significant_label) |
        !nzchar(significant_label)
    ] <- "NA"


    volcano_stats_long$Comparison_hover <- paste0(

      "<br><br><b>Comparison: ",
      htmltools::htmlEscape(
        volcano_stats_long$Comparison
      ),
      "</b>",

      "<br>FC: ",
      format_network_number(
        volcano_stats_long$FC
      ),

      "<br>Adjusted p-value: ",
      format_network_number(
        volcano_stats_long$Adj_p
      ),

      "<br>-log10(FDR): ",
      format_network_number(
        volcano_stats_long$Adj_p_log
      ),

      "<br>Mean intensity: ",
      format_network_number(
        volcano_stats_long$Mean
      ),

      "<br>Mean [",
      htmltools::htmlEscape(
        group_num_label
      ),
      "]: ",
      format_network_number(
        volcano_stats_long$mean_num
      ),

      "<br>Mean [",
      htmltools::htmlEscape(
        group_den_label
      ),
      "]: ",
      format_network_number(
        volcano_stats_long$mean_den
      ),

      "<br>Test scale: ",
      htmltools::htmlEscape(
        test_scale_label
      ),

      "<br>Significant: ",
      htmltools::htmlEscape(
        significant_label
      )
    )


    statistics_by_node <- volcano_stats_long %>%

      dplyr::group_by(
        ClusterID
      ) %>%

      dplyr::summarise(

        Statistical_results = paste0(
          Comparison_hover,
          collapse = ""
        ),

        .groups = "drop"
      )


    node_data <- node_data %>%
  dplyr::left_join(
    statistics_by_node,
    by = "ClusterID"
  )

    volcano_added <- TRUE
  }
}


if (!"Statistical_results" %in% names(node_data)) {
  node_data$Statistical_results <- NA_character_
}


# Load feature annotations independently of the statistics switch.
annotation_cols <- c(
  "GNPS_annotation",
  "Other_annotation"
)

if (
  !is.null(rv$volcano) &&
  is.data.frame(rv$volcano) &&
  nrow(rv$volcano) > 0 &&
  "id" %in% names(rv$volcano)
) {
  annotations <- as.data.frame(
    rv$volcano,
    check.names = FALSE
  )

  for (column in annotation_cols) {
    if (!column %in% names(annotations)) {
      annotations[[column]] <- NA_character_
    }
  }

  annotations_by_node <- annotations %>%
    dplyr::mutate(
      ClusterID = trimws(as.character(id))
    ) %>%
    dplyr::filter(
      !is.na(ClusterID),
      nzchar(ClusterID)
    ) %>%
    dplyr::group_by(ClusterID) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(annotation_cols),
        ~ {
          values <- clean_missing_text(.x)
          values <- unique(values[!is.na(values)])

          if (length(values)) {
            paste(values, collapse = " | ")
          } else {
            NA_character_
          }
        }
      ),
      .groups = "drop"
    )

  node_data <- dplyr::left_join(
    node_data,
    annotations_by_node,
    by = "ClusterID"
  )
}

for (column in annotation_cols) {
  if (!column %in% names(node_data)) {
    node_data[[column]] <- NA_character_
  }
}


sirius_text <- node_data$SIRIUS_values

sirius_text[
  is.na(sirius_text) |
    !nzchar(sirius_text)
] <- "NA"


gnps_text <- node_data$GNPS_annotation

gnps_text[
  is.na(gnps_text) |
    !nzchar(gnps_text)
] <- "NA"

other_text <- clean_missing_text(
  node_data$Other_annotation
)

other_text[is.na(other_text)] <- "NA"

node_data$Hover <- paste0(

  "<b>Cluster ID:</b> ",
  node_data$ClusterID,

  "<br><b>ComponentIndex:</b> ",
  selected_component,

  "<br><b>Node degree:</b> ",
  node_data$Degree,

  "<br><b>",
  input$stats_sirius_col,
  ":</b> ",
  sirius_text,

  "<br><b>GNPS annotation:</b> ",
htmltools::htmlEscape(gnps_text),

"<br><b>Other annotation:</b> ",
htmltools::htmlEscape(other_text),

  "<br><b>Selected class:</b> ",
  ifelse(
    node_data$Selected_class,
    "Yes",
    "No"
  )
)


if (volcano_added) {

  stats_hover <-
    node_data$Statistical_results

  stats_hover[
    is.na(stats_hover) |
      !nzchar(stats_hover)
  ] <- paste0(
    "<br><br>",
    "No processed statistical results matched this node."
  )


  node_data$Hover <- paste0(

    node_data$Hover,

    "<br><br><b>Processed statistical comparisons</b>",

    stats_hover
  )
}

  label_mode <- input$network_label_mode %||%
    "selected"

  node_data$Node_label <- dplyr::case_when(
    label_mode == "all" ~ node_data$ClusterID,
    label_mode == "selected" &
      node_data$Selected_class ~ node_data$ClusterID,
    TRUE ~ ""
  )

  # Create edge coordinates using the simplified graph
  simplified_edges <- igraph::as_data_frame(
    graph_object,
    what = "edges"
  ) %>%

    dplyr::transmute(
      ClusterID1 = as.character(from),
      ClusterID2 = as.character(to)
    )

  edge_coordinates <- simplified_edges %>%

    dplyr::left_join(
      node_data %>%
        dplyr::select(
          ClusterID,
          x1 = x,
          y1 = y
        ),
      by = c(
        "ClusterID1" = "ClusterID"
      )
    ) %>%

    dplyr::left_join(
      node_data %>%
        dplyr::select(
          ClusterID,
          x2 = x,
          y2 = y
        ),
      by = c(
        "ClusterID2" = "ClusterID"
      )
    )

  list(
  component = selected_component,
  graph = graph_object,
  nodes = node_data,
  edges = edge_coordinates,
  volcano_stats = volcano_stats_long,
  volcano_added = volcano_added,
  node_count = igraph::vcount(graph_object),
  edge_count = igraph::ecount(graph_object)
)
})

output$gnps_component_network <- plotly::renderPlotly({

  network <- selected_component_network_data()

  node_data <- network$nodes
  edge_data <- network$edges

  plot_object <- plotly::plot_ly(
    source = "gnps_component_network"
  )

  # GNPS edges
  plot_object <- plot_object %>%

    plotly::add_segments(
      data = edge_data,

      x = ~x1,
      y = ~y1,
      xend = ~x2,
      yend = ~y2,

      inherit = FALSE,

      line = list(
        color = "rgba(130,130,130,0.55)",
        width = 1
      ),

      hoverinfo = "skip",
      showlegend = FALSE
    )

  node_types <- c(
    "Selected class",
    "Other SIRIUS annotation",
    "No SIRIUS annotation"
  )

  node_colors <- c(
    "Selected class" = "#E74C3C",
    "Other SIRIUS annotation" = "#66CDAA",
    "No SIRIUS annotation" = "#BDBDBD"
  )

  node_sizes <- c(
    "Selected class" = 16,
    "Other SIRIUS annotation" = 11,
    "No SIRIUS annotation" = 9
  )

  for (current_type in node_types) {

    current_nodes <- node_data %>%
      dplyr::filter(
        Node_type == current_type
      )

    if (!nrow(current_nodes)) {
      next
    }

    plot_object <- plot_object %>%

      plotly::add_trace(
        data = current_nodes,

        x = ~x,
        y = ~y,

        type = "scatter",
        mode = "markers+text",

        text = ~Node_label,
        hovertext = ~Hover,
        hoverinfo = "text",

        key = ~ClusterID,

        name = current_type,
        inherit = FALSE,

        textposition = "top center",

        textfont = list(
          size = 10,
          color = "#222222"
        ),

        marker = list(
          size = unname(
            node_sizes[
              current_type
            ]
          ),

          color = unname(
            node_colors[
              current_type
            ]
          ),

          line = list(
            color = "#333333",
            width = 1
          )
        )
      )
  }

  plot_object %>%

    plotly::layout(

      title = list(
        text = paste0(
          "GNPS Component ",
          network$component,
          " — ",
          network$node_count,
          " nodes, ",
          network$edge_count,
          " edges"
        )
      ),

      xaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
      ),

      yaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        scaleanchor = "x",
        scaleratio = 1
      ),

      hovermode = "closest",
      dragmode = "pan",

      legend = list(
        orientation = "h",
        x = 0,
        y = -0.05
      ),

      margin = list(
        l = 20,
        r = 20,
        b = 70,
        t = 70
      )
    ) %>%

    plotly::config(
      displaylogo = FALSE,

      modeBarButtonsToRemove = c(
        "lasso2d",
        "select2d"
      )
    )
})

output$gnps_network_node_details <- renderUI({

  click <- plotly::event_data(
    "plotly_click",
    source = "gnps_component_network"
  )


  if (
    is.null(click) ||
    is.null(click$key) ||
    !length(click$key)
  ) {

    return(
      div(
        class = "small-note",
        style = "margin-top: 8px;",
        "Click a network node to show its information."
      )
    )
  }


  network <- selected_component_network_data()


  clicked_id <- as.character(
    click$key[1]
  )


  row <- network$nodes %>%

    dplyr::filter(
      ClusterID == clicked_id
    ) %>%

    dplyr::slice(1)


  if (!nrow(row)) {
    return(NULL)
  }


  sirius_value <- row$SIRIUS_values[1]


  if (
    is.na(sirius_value) ||
    !nzchar(sirius_value)
  ) {

    sirius_value <- "NA"
  }


  gnps_value <- row$GNPS_annotation[1]


  if (
    is.na(gnps_value) ||
    !nzchar(gnps_value)
  ) {

    gnps_value <- "NA"
  }

  other_value <- clean_missing_text(
  row$Other_annotation[1]
)

if (is.na(other_value)) {
  other_value <- "NA"
}

  format_one <- function(x) {

    x <- suppressWarnings(
      as.numeric(x)
    )

    if (
      !length(x) ||
      !is.finite(x[1])
    ) {
      return("NA")
    }

    format(
      signif(
        x[1],
        5
      ),
      scientific = FALSE,
      trim = TRUE
    )
  }


  stats_rows <- network$volcano_stats %>%

    dplyr::filter(
      ClusterID == clicked_id
    )


  stats_table_ui <- if (!nrow(stats_rows)) {

    div(
      class = "small-note",
      style = "margin-top: 10px;",
      "No processed statistical comparison matched this network node."
    )

  } else {


    table_rows <- lapply(

      seq_len(
        nrow(stats_rows)
      ),

      function(i) {


        current <- stats_rows[
          i,
          ,
          drop = FALSE
        ]


        group_num <- current$Group_num[1]
        group_den <- current$Group_den[1]


        if (
          is.na(group_num) ||
          !nzchar(group_num)
        ) {
          group_num <- "Group 1"
        }


        if (
          is.na(group_den) ||
          !nzchar(group_den)
        ) {
          group_den <- "Group 2"
        }


        test_scale <- current$TestScale[1]

        if (
          is.na(test_scale) ||
          !nzchar(test_scale)
        ) {
          test_scale <- "NA"
        }


        significant_value <-
          current$Significant_default[1]

        if (
          is.na(significant_value) ||
          !nzchar(significant_value)
        ) {
          significant_value <- "NA"
        }


        tags$tr(

          tags$td(
            current$Comparison[1]
          ),

          tags$td(
            group_num
          ),

          tags$td(
            format_one(
              current$mean_num
            )
          ),

          tags$td(
            group_den
          ),

          tags$td(
            format_one(
              current$mean_den
            )
          ),

          tags$td(
            format_one(
              current$FC
            )
          ),

          tags$td(
            format_one(
              current$Adj_p
            )
          ),

          tags$td(
            format_one(
              current$Adj_p_log
            )
          ),

          tags$td(
            format_one(
              current$Mean
            )
          ),

          tags$td(
            test_scale
          ),

          tags$td(
            significant_value
          )
        )
      }
    )


    tagList(

      h4(
        style = "margin-top: 14px;",
        "All processed statistical comparisons"
      ),


      div(
        style = "overflow-x:auto;",


        tags$table(

          class =
            "table table-striped table-bordered table-condensed",


          tags$thead(

            tags$tr(

              tags$th("Comparison"),

              tags$th("Group 1"),

              tags$th("Mean 1"),

              tags$th("Group 2"),

              tags$th("Mean 2"),

              tags$th("FC"),

              tags$th("Adjusted p-value"),

              tags$th("-log10(FDR)"),

              tags$th("Mean intensity"),

              tags$th("Test scale"),

              tags$th("Significant")
            )
          ),


          tags$tbody(
            table_rows
          )
        )
      )
    )
  }


  div(

    style = "
      background:#ffffffcc;
      padding:10px;
      border-radius:10px;
      border:1px solid #dddddd;
      margin-top:10px;
    ",


    h4(
      paste0(
        "Selected network node: ",
        clicked_id
      )
    ),


    tags$ul(

      tags$li(
        strong("ComponentIndex: "),
        network$component
      ),

      tags$li(
        strong("Node degree: "),
        row$Degree[1]
      ),

      tags$li(
        strong(
          paste0(
            input$stats_sirius_col,
            ": "
          )
        ),
        sirius_value
      ),

      tags$li(
        strong("Matches selected class: "),
        ifelse(
          isTRUE(
            row$Selected_class[1]
          ),
          "Yes",
          "No"
        )
      ),

      tags$li(
  strong("GNPS annotation: "),
  gnps_value
),

tags$li(
  strong("Other annotation: "),
  other_value
)
    ),


    stats_table_ui
  )
})

output$gnps_network_section <- renderUI({

  if (is.null(input$file_gnps_pairs)) {
    return(NULL)
  }

  if (
    is.null(input$network_component) ||
    !nzchar(
      as.character(
        input$network_component
      )
    )
  ) {

    return(
      tagList(
        tags$hr(),

        div(
          class = "small-note",
          paste0(
            "Choose a SIRIUS class that occurs in a GNPS ",
            "component to display its network."
          )
        )
      )
    )
  }

  tagList(

    tags$hr(),

    h3(
      "Interactive selected GNPS component"
    ),

    div(
      class = "small-note",
      style = "margin-bottom: 8px;",

      HTML(
        paste0(
          "<b>Red:</b> selected SIRIUS class &nbsp;|&nbsp; ",
          "<b>Green:</b> another SIRIUS value &nbsp;|&nbsp; ",
          "<b>Grey:</b> no value in the selected SIRIUS column"
        )
      )
    ),

    withSpinner(
      plotlyOutput(
        "gnps_component_network",
        height = "650px"
      ),
      type = 8,
      color = "#66CDAA"
    ),

    uiOutput(
      "gnps_network_node_details"
    )
  )
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
    withSpinner(DTOutput("selected_class_component_table"), type = 8, color = "#66CDAA"),
    uiOutput("gnps_network_section")
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

      incProgress(
  0.02,
  detail = "Joining additional peak-table columns"
)

if (isTRUE(input$use_peak_extra_cols)) {

  selected_peak_cols <- intersect(
    input$peak_extra_cols %||% character(0),
    names(built()$raw)
  )

  if (length(selected_peak_cols)) {

    peak_raw <- as.data.frame(
      built()$raw,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    validate(
      need(
        nrow(peak_raw) == nrow(fmap),
        paste0(
          "Peak-table metadata could not be joined because ",
          "the peak table and feature map have different row counts."
        )
      )
    )

    peak_colmap <- make_prefixed_colmap(
      selected_peak_cols,
      prefix = "Peak_"
    )

    peak_extra <- peak_raw[
      ,
      names(peak_colmap),
      drop = FALSE
    ]

    names(peak_extra) <- unname(
      peak_colmap
    )

    peak_extra$Feature <- as.character(
      fmap$Feature
    )

    peak_extra <- peak_extra[
      ,
      c(
        "Feature",
        unname(peak_colmap)
      ),
      drop = FALSE
    ]

    volc <- volc %>%
      dplyr::left_join(
        peak_extra,
        by = "Feature"
      )
  }
}
      
      # default annotation columns
      volc$`NPC#class` <- NA_character_
      volc$`ClassyFire#class` <- NA_character_
      volc$GNPS_annotation <- NA_character_
      volc$Other_annotation <- NA_character_

      incProgress(
      0.05,
      detail = "Joining SIRIUS (optional)"
    )

if (
  isTRUE(input$use_sirius) &&
  !is.null(input$file_sirius)
) {

  s <- as.data.frame(
    sirius_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  req(
    input$sirius_idcol,
    input$sirius_npcol,
    input$sirius_cfcol
  )

  validate(
    need(
      input$sirius_idcol %in% names(s),
      "Selected SIRIUS Feature ID column was not found."
    ),

    need(
      input$sirius_npcol %in% names(s),
      "Selected SIRIUS NPC column was not found."
    ),

    need(
      input$sirius_cfcol %in% names(s),
      "Selected SIRIUS ClassyFire column was not found."
    )
  )

  selected_sirius_extra <- character(0)

  if (isTRUE(input$use_sirius_extra_cols)) {

    selected_sirius_extra <- intersect(
      input$sirius_extra_cols %||% character(0),
      names(s)
    )

    selected_sirius_extra <- setdiff(
      selected_sirius_extra,
      c(
        input$sirius_idcol,
        input$sirius_npcol,
        input$sirius_cfcol
      )
    )
  }

  selected_sirius_extra <- unique(c(
  selected_sirius_extra,
  grep(
    "probability",
    names(s),
    ignore.case = TRUE,
    value = TRUE
  )
))
  
  sirius_colmap <- make_prefixed_colmap(
    selected_sirius_extra,
    prefix = "SIRIUS_"
  )

  # Build the core SIRIUS join table
  ss <- data.frame(
    id = trimws(
      as.character(
        s[[input$sirius_idcol]]
      )
    ),

    `NPC#class` = clean_missing_text(
      s[[input$sirius_npcol]]
    ),

    `ClassyFire#class` = clean_missing_text(
      s[[input$sirius_cfcol]]
    ),

    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Add selected additional SIRIUS columns.
  # Numeric columns remain numeric.
  if (length(sirius_colmap)) {

    for (original_name in names(sirius_colmap)) {

      output_name <- sirius_colmap[[
        original_name
      ]]

      value <- s[[
        original_name
      ]]

      if (
        is.character(value) ||
        is.factor(value)
      ) {
        value <- clean_missing_text(
          value
        )
      }

      ss[[
        output_name
      ]] <- value
    }
  }

  ss <- ss %>%
    dplyr::filter(
      !is.na(id),
      nzchar(id)
    )

  # Prevent duplicate SIRIUS IDs from multiplying volcano rows.
  # The first row is retained if duplicate IDs are present.
  duplicate_count <- sum(
    duplicated(ss$id)
  )

  if (duplicate_count > 0) {

    showNotification(
      paste0(
        "SIRIUS contains ",
        duplicate_count,
        " duplicated mapping ID row(s). ",
        "The first row for each ID was retained."
      ),
      type = "warning",
      duration = 7
    )

    ss <- ss %>%
      dplyr::distinct(
        id,
        .keep_all = TRUE
      )
  }

peak_table_ids <- unique(
  trimws(
    as.character(
      volc$id
    )
  )
)

sirius_ids <- unique(
  trimws(
    as.character(
      ss$id
    )
  )
)

peak_table_ids <- peak_table_ids[
  !is.na(peak_table_ids) &
    nzchar(peak_table_ids)
]

sirius_ids <- sirius_ids[
  !is.na(sirius_ids) &
    nzchar(sirius_ids)
]

matched_sirius_ids <- sum(
  peak_table_ids %in% sirius_ids
)

volc <- volc %>%

  dplyr::mutate(
    id = trimws(
      as.character(id)
    )
  ) %>%

  dplyr::left_join(
    ss,
    by = "id",
    suffix = c("", ".sirius")
  ) %>%

  dplyr::mutate(

    `NPC#class` = dplyr::coalesce(
      .data[["NPC#class.sirius"]],
      .data[["NPC#class"]]
    ),

    `ClassyFire#class` = dplyr::coalesce(
      .data[["ClassyFire#class.sirius"]],
      .data[["ClassyFire#class"]]
    )
  ) %>%

  dplyr::select(
    -dplyr::any_of(
      c(
        "NPC#class.sirius",
        "ClassyFire#class.sirius"
      )
    )
  )

if (matched_sirius_ids == 0) {

  showNotification(
    paste0(
      "SIRIUS summary was uploaded, but no IDs matched. ",
      "Check that the peak-table Row ID column and selected ",
      "SIRIUS Feature ID column contain the same values."
    ),
    type = "warning",
    duration = 8
  )

} else {

  showNotification(
    paste0(
      "SIRIUS annotation joined successfully: ",
      format(
        matched_sirius_ids,
        big.mark = ",",
        scientific = FALSE
      ),
      " of ",
      format(
        length(peak_table_ids),
        big.mark = ",",
        scientific = FALSE
      ),
      " unique peak-table IDs matched."
    ),
    type = "message",
    duration = 6
  )
}
}

incProgress(
  0.05,
  detail = "Joining GNPS annotation (optional)"
)

if (
  isTRUE(input$use_gnps_annotation) &&
  !is.null(input$file_gnps_annotation)
) {

  g <- as.data.frame(
    gnps_annotation_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  req(
    input$gnps_annotation_idcol,
    input$gnps_annotation_col
  )

  validate(
    need(
      input$gnps_annotation_idcol %in% names(g),
      "Selected GNPS ID column was not found."
    ),

    need(
      input$gnps_annotation_col %in% names(g),
      "Selected primary GNPS annotation column was not found."
    )
  )

  selected_gnps_extra <- character(0)

  if (isTRUE(input$use_gnps_extra_cols)) {

    selected_gnps_extra <- intersect(
      input$gnps_extra_cols %||% character(0),
      names(g)
    )

    selected_gnps_extra <- setdiff(
      selected_gnps_extra,
      c(
        input$gnps_annotation_idcol,
        input$gnps_annotation_col
      )
    )
  }

  gnps_colmap <- make_prefixed_colmap(
    selected_gnps_extra,
    prefix = "GNPS_"
  )

  gnps_primary <- data.frame(
    id = trimws(
      as.character(
        g[[input$gnps_annotation_idcol]]
      )
    ),

    GNPS_annotation = clean_missing_text(
      g[[input$gnps_annotation_col]]
    ),

    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%

    dplyr::filter(
      !is.na(id),
      nzchar(id)
    ) %>%

    dplyr::group_by(id) %>%

    dplyr::summarise(

      GNPS_annotation = {

        values <- unique(
          GNPS_annotation[
            !is.na(GNPS_annotation)
          ]
        )

        if (length(values)) {
          paste(
            values,
            collapse = " | "
          )
        } else {
          NA_character_
        }
      },

      .groups = "drop"
    )

  gnps_extra <- NULL

  if (length(gnps_colmap)) {

    gnps_extra <- data.frame(
      id = trimws(
        as.character(
          g[[input$gnps_annotation_idcol]]
        )
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    for (original_name in names(gnps_colmap)) {

      output_name <- gnps_colmap[[
        original_name
      ]]

      value <- g[[
        original_name
      ]]

      # Convert empty character values to real NA
      if (
        is.character(value) ||
        is.factor(value)
      ) {
        value <- clean_missing_text(
          value
        )
      }

      # Numeric GNPS columns remain numeric
      gnps_extra[[
        output_name
      ]] <- value
    }

    gnps_extra <- gnps_extra %>%
      dplyr::filter(
        !is.na(id),
        nzchar(id)
      )

    duplicate_extra_ids <- sum(
      duplicated(gnps_extra$id)
    )

    if (duplicate_extra_ids > 0) {

      showNotification(
        paste0(
          "GNPS contains ",
          duplicate_extra_ids,
          " duplicated ID row(s). ",
          "The primary annotations were combined, while the first row ",
          "was retained for each additional GNPS column."
        ),
        type = "warning",
        duration = 8
      )
    }

    gnps_extra <- gnps_extra %>%
      dplyr::distinct(
        id,
        .keep_all = TRUE
      )
  }

  # Combine primary and additional GNPS columns
  gnps_join <- gnps_primary

  if (!is.null(gnps_extra)) {

    gnps_join <- gnps_join %>%
      dplyr::left_join(
        gnps_extra,
        by = "id"
      )
  }

  peak_table_ids <- unique(
  trimws(
    as.character(
      volc$id
    )
  )
)

peak_table_ids <- peak_table_ids[
  !is.na(peak_table_ids) &
    nzchar(peak_table_ids)
]

gnps_ids <- unique(
  trimws(
    as.character(
      gnps_join$id
    )
  )
)

gnps_ids <- gnps_ids[
  !is.na(gnps_ids) &
    nzchar(gnps_ids)
]

matched_gnps_ids <- sum(
  peak_table_ids %in% gnps_ids
)

  volc <- volc %>%

    dplyr::mutate(
      id = trimws(
        as.character(id)
      )
    ) %>%

    dplyr::left_join(
      gnps_join,
      by = "id",
      suffix = c("", ".gnps")
    ) %>%

    dplyr::mutate(
      GNPS_annotation = dplyr::coalesce(
        .data[["GNPS_annotation.gnps"]],
        .data[["GNPS_annotation"]]
      )
    ) %>%

    dplyr::select(
      -dplyr::any_of(
        "GNPS_annotation.gnps"
      )
    )

  if (matched_gnps_ids == 0) {

    showNotification(
      paste0(
        "GNPS annotation was uploaded, but no IDs matched. ",
        "Check that the peak-table Row ID column and selected ",
        "GNPS ID column contain the same values."
      ),
      type = "warning",
      duration = 8
    )

  } else {

    showNotification(
  paste0(
    "GNPS annotation joined successfully: ",
    format(
      matched_gnps_ids,
      big.mark = ",",
      scientific = FALSE
    ),
    " of ",
    format(
      length(peak_table_ids),
      big.mark = ",",
      scientific = FALSE
    ),
    " unique peak-table IDs matched."
  ),
  type = "message",
  duration = 6
)
  }
}

incProgress(
  0.05,
  detail = "Joining other annotation source (optional)"
)

if (
  isTRUE(input$use_other_annotation) &&
  !is.null(input$file_other_annotation)
) {

  ann <- as.data.frame(
    other_annotation_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  raw_peak <- as.data.frame(
    built()$raw,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  req(
    input$other_peak_id_col,
    input$other_annotation_idcol,
    input$other_annotation_col
  )

  validate(

    need(
      input$other_peak_id_col %in% names(raw_peak),
      "Selected peak-table ID column for Other Annotation Source was not found."
    ),

    need(
      input$other_annotation_idcol %in% names(ann),
      "Selected annotation-file ID column was not found."
    ),

    need(
      input$other_annotation_col %in% names(ann),
      "Selected primary annotation column was not found."
    ),

    need(
      nrow(raw_peak) == nrow(fmap),
      "Peak table and internal feature map have different row counts."
    )
  )


  # --------------------------------------------------
  # Map the selected peak-table ID to internal Feature
  # --------------------------------------------------
  peak_other_map <- tibble::tibble(

    Feature = as.character(
      fmap$Feature
    ),

    .other_join_id = trimws(
      as.character(
        raw_peak[[input$other_peak_id_col]]
      )
    )
  )

  # --------------------------------------------------
  # Primary annotation
  # --------------------------------------------------
  other_primary <- tibble::tibble(

    .other_join_id = trimws(
      as.character(
        ann[[input$other_annotation_idcol]]
      )
    ),

    Other_annotation = clean_missing_text(
      ann[[input$other_annotation_col]]
    )
  ) %>%

    dplyr::filter(
      !is.na(.other_join_id),
      nzchar(.other_join_id)
    ) %>%

    dplyr::group_by(
      .other_join_id
    ) %>%

    dplyr::summarise(

      Other_annotation = {

        values <- unique(
          Other_annotation[
            !is.na(Other_annotation)
          ]
        )

        if (length(values)) {
          paste(
            values,
            collapse = " | "
          )
        } else {
          NA_character_
        }
      },

      .groups = "drop"
    )


  # --------------------------------------------------
  # Optional additional columns
  # --------------------------------------------------
  selected_other_extra <- character(0)

  if (isTRUE(input$use_other_extra_cols)) {

    selected_other_extra <- intersect(
      input$other_extra_cols %||%
        character(0),
      names(ann)
    )

    selected_other_extra <- setdiff(
      selected_other_extra,
      c(
        input$other_annotation_idcol,
        input$other_annotation_col
      )
    )
  }


  other_colmap <- make_prefixed_colmap(
    selected_other_extra,
    prefix = "Other_"
  )

  other_extra <- NULL


  if (length(other_colmap)) {

    other_extra <- data.frame(

      .other_join_id = trimws(
        as.character(
          ann[[input$other_annotation_idcol]]
        )
      ),

      check.names = FALSE,
      stringsAsFactors = FALSE
    )


    for (original_name in names(other_colmap)) {

      output_name <-
        other_colmap[[original_name]]

      value <-
        ann[[original_name]]

      if (
        is.character(value) ||
        is.factor(value)
      ) {

        value <- clean_missing_text(
          value
        )
      }

      other_extra[[output_name]] <- value
    }


    duplicate_other_ids <- sum(
      duplicated(
        other_extra$.other_join_id
      )
    )

    if (duplicate_other_ids > 0) {

      showNotification(
        paste0(
          "Other annotation source contains ",
          duplicate_other_ids,
          " duplicated ID row(s). ",
          "Primary annotations were combined and the first row ",
          "was retained for additional columns."
        ),
        type = "warning",
        duration = 8
      )
    }


    other_extra <- other_extra %>%

      dplyr::filter(
        !is.na(.other_join_id),
        nzchar(.other_join_id)
      ) %>%

      dplyr::distinct(
        .other_join_id,
        .keep_all = TRUE
      )
  }


  other_join <- other_primary


  if (!is.null(other_extra)) {

    other_join <- other_join %>%
      dplyr::left_join(
        other_extra,
        by = ".other_join_id"
      )
  }


  # --------------------------------------------------
  # Matching statistics
  # --------------------------------------------------
  peak_other_ids <- unique(
    peak_other_map$.other_join_id
  )

  peak_other_ids <- peak_other_ids[
    !is.na(peak_other_ids) &
      nzchar(peak_other_ids)
  ]

  annotation_other_ids <- unique(
    other_join$.other_join_id
  )

  matched_other_ids <- sum(
    peak_other_ids %in%
      annotation_other_ids
  )


  # --------------------------------------------------
  # Join into volcano table
  # --------------------------------------------------
  volc <- volc %>%

    dplyr::left_join(
      peak_other_map,
      by = "Feature"
    ) %>%

    dplyr::left_join(
      other_join,
      by = ".other_join_id",
      suffix = c("", ".other")
    ) %>%

    dplyr::mutate(

      Other_annotation =
        dplyr::coalesce(
          .data[["Other_annotation.other"]],
          .data[["Other_annotation"]]
        )
    ) %>%

    dplyr::select(
      -dplyr::any_of(
        c(
          ".other_join_id",
          "Other_annotation.other"
        )
      )
    )

  if (matched_other_ids == 0) {

    showNotification(
      paste0(
        "Other annotation source was uploaded, but no IDs matched. ",
        "Check the selected peak-table ID column and annotation-file ID column."
      ),
      type = "warning",
      duration = 8
    )

  } else {

    showNotification(
      paste0(
        "Other annotation source joined successfully: ",
        format(
          matched_other_ids,
          big.mark = ",",
          scientific = FALSE
        ),
        " of ",
        format(
          length(peak_other_ids),
          big.mark = ",",
          scientific = FALSE
        ),
        " unique peak-table IDs matched."
      ),
      type = "message",
      duration = 6
    )
  }
}

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
      updateMaterialSwitch(session, "use_gnps_annotated_filter",value = FALSE)
      updateMaterialSwitch(session, "use_other_annotated_filter", value = FALSE)
      updateMaterialSwitch(session, "use_npc_probability", value = FALSE)
      updateMaterialSwitch(session, "use_classyfire_probability", value = FALSE)
      updateNumericInput(session, "npc_probability_min", value = 0.9)
      updateNumericInput(session, "classyfire_probability_min", value = 0.9)
      updateMaterialSwitch(session, "use_component_filter", value = FALSE)
      updatePickerInput(session, "volcano_components", selected = character(0))
      
      # pickers / radios
      updatePickerInput(session, "sel_feat", selected = character(0))
      updatePickerInput(session, "npc_filter_values", selected = character(0))
      updatePickerInput(session, "classyfire_filter_values", selected = character(0))
      updateRadioButtons(session, "color_by", selected = "Groups")
      updateRadioButtons(session, "volcano_y_axis", selected = "fdr")
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
  
  # Nested probability controls for one SIRIUS annotation system.
sirius_probability_controls <- function(
    prefix, family, annotation_column
) {

  columns <- names(rv$volcano)

  score_cols <- columns[
    grepl(
      paste0("^SIRIUS_", family, "#"),
      columns,
      ignore.case = TRUE
    ) &
      grepl(
        "probability",
        columns,
        ignore.case = TRUE
      )
  ]

  # Prefer the score corresponding to the selected annotation level.
  expected <- paste0(
    "SIRIUS_",
    annotation_column,
    " Probability"
  )

  preferred <- score_cols[
    tolower(score_cols) %in% tolower(expected)
  ]

  selected <- isolate(
    input[[paste0(prefix, "_probability_col")]]
  )

  if (
    is.null(selected) ||
    !selected %in% score_cols
  ) {
    selected <- if (length(preferred)) {
      preferred[[1]]
    } else if (length(score_cols)) {
      score_cols[[1]]
    } else {
      character(0)
    }
  }

  tagList(

    materialSwitch(
      inputId = paste0("use_", prefix, "_probability"),
      label = paste("Also filter", family, "by probability"),
      value = FALSE,
      status = "success"
    ),

    conditionalPanel(
      condition = paste0(
        "input.use_", prefix, "_probability == true"
      ),

      if (length(score_cols)) {

        tagList(
          selectInput(
            inputId = paste0(prefix, "_probability_col"),
            label = paste(family, "probability column:"),

            choices = stats::setNames(
              score_cols,
              sub("^SIRIUS_", "", score_cols)
            ),

            selected = selected
          ),

          numericInput(
            inputId = paste0(prefix, "_probability_min"),
            label = "Keep probability greater than:",
            value = 0.9,
            min = 0,
            max = 1,
            step = 0.01
          )
        )

      } else {

        div(
          class = "small-note",
          paste(
            "No SIRIUS", family, "probability columns found.",
            "Rerun preprocessing after importing SIRIUS."
          )
        )
      }
    )
  )
}


# Called only inside the corresponding class-filter block.
apply_sirius_probability <- function(dd, prefix, family) {

  if (!isTRUE(input[[paste0("use_", prefix, "_probability")]])) {
    return(dd)
  }

  column <- input[[paste0(prefix, "_probability_col")]]

  allowed <- names(dd)[
    grepl(
      paste0("^SIRIUS_", family, "#"),
      names(dd),
      ignore.case = TRUE
    ) &
      grepl("probability", names(dd), ignore.case = TRUE)
  ]

  validate(
    need(
      length(column) == 1L && column %in% allowed,
      paste("Select a SIRIUS", family, "probability column.")
    )
  )

  cutoff <- suppressWarnings(
    as.numeric(input[[paste0(prefix, "_probability_min")]])
  )

  validate(
    need(
      length(cutoff) == 1L &&
        is.finite(cutoff) &&
        cutoff >= 0 &&
        cutoff <= 1,
      "Probability threshold must be between 0 and 1."
    )
  )

  score <- suppressWarnings(
    as.numeric(as.character(dd[[column]]))
  )

  keep <- is.finite(score) &
    score >= 0 &
    score <= 1 &
    score > cutoff

  dd[keep, , drop = FALSE]
}
  
  annotation_filter_choices <- reactive({
  req(procReady())

  dd <- rv$volcano

  clean_choices <- function(x) {

  x <- clean_missing_text(
    x
  )

  x <- x[
    !is.na(x) &
      nzchar(x)
  ]

  sort(
    unique(x)
  )
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
        choices = c(
  "Groups" = "Groups",
  "Mean" = "Mean",
  "Fold change" = "FC"
),
        selected = "Groups",
        inline = TRUE
      ),
      radioButtons(
        "volcano_y_axis",
        "Y-axis:",
        choices = c(
          "-log10(FDR)" = "fdr",
          "Mean log10(Intensity)" = "mean"
        ),
        selected = "fdr",
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

numericInput(
  "volcano_top_n",
  "Number of top features to label:",
  value = 0,
  min = 0,
  step = 1
),

uiOutput("volcano_label_column_ui"),

uiOutput("volcano_sliders"),
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

    tagList(

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
      ),

      sirius_probability_controls(
        prefix = "npc",
        family = "NPC",
        annotation_column =
          input$sirius_npcol %||% "NPC#class"
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

    tagList(

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
          `count-selected-text` =
            "{0} ClassyFire class(es) selected"
        )
      ),

      sirius_probability_controls(
        prefix = "classyfire",
        family = "ClassyFire",
        annotation_column =
          input$sirius_cfcol %||% "ClassyFire#class"
      )
    )

  } else {

    div(
      class = "small-note",
      "No ClassyFire classes detected. Upload SIRIUS annotation and rerun preprocessing."
    )
  }
),

materialSwitch(
  inputId = "use_gnps_annotated_filter",
  label = "Filter by GNPS annotation",
  value = FALSE,
  status = "success"
),

materialSwitch(
  inputId = "use_other_annotated_filter",
  label = "Filter by Other annotation",
  value = FALSE,
  status = "success"
),

materialSwitch(
  "use_component_filter",
  "Filter by GNPS ComponentIndex",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_component_filter == true",

  conditionalPanel(
    condition = paste(
      "input.use_main_gnps_pairs == true &&",
      "input.file_main_gnps_pairs != null"
    ),

    pickerInput(
      inputId = "volcano_components",
      label = "ComponentIndex(es):",
      choices = character(0),
      selected = character(0),
      multiple = TRUE,

      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `none-selected-text` = "Select ComponentIndex(es)",
        `style` = "btn-success",
        `selected-text-format` = "count > 1",
        `count-selected-text` =
          "{0} ComponentIndex(es) selected"
      )
    )
  ),

  conditionalPanel(
    condition = paste(
      "input.use_main_gnps_pairs != true ||",
      "input.file_main_gnps_pairs == null"
    ),

    div(
      class = "small-note",
      paste(
        "Enable the GNPS ComponentIndex source and upload",
        "a GNPS pairs file in Load & Process."
      )
    )
  )
),

tags$hr(),
h4(class = "highlight", "Interactive heatmap"),
materialSwitch(
  inputId = "show_interactive_heatmap",
  label = "Activate",
  value = FALSE,
  status = "success"
)
    )
  })

  output$volcano_ready <- renderText({
  if (procReady()) "yes" else "no"
})

outputOptions(
  output,
  "volcano_ready",
  suspendWhenHidden = FALSE
)
  
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

  # Optional GNPS network pairs.
# Read the optional pairs file once.
main_gnps_pairs_df <- reactive({
  req(
    isTRUE(input$use_main_gnps_pairs),
    input$file_main_gnps_pairs
  )

  file <- input$file_main_gnps_pairs
  ext <- tolower(tools::file_ext(file$name))

  validate(
    need(
      ext %in% c("tsv", "txt", "csv"),
      "Upload GNPS pairs as TSV, TXT, or CSV."
    )
  )

  g <- data.table::fread(
    file$datapath,
    sep = if (ext == "csv") "," else "\t",
    colClasses = "character",
    data.table = FALSE,
    check.names = FALSE
  )

  validate(
    need(ncol(g) >= 3L, "The pairs file needs at least three columns."),
    need(
      !anyDuplicated(names(g)),
      "The pairs file contains duplicate column names. Rename them first."
    )
  )

  g
})


output$main_gnps_pairs_pickers <- renderUI({
  g <- main_gnps_pairs_df()

  # Select a standard name when present.
  # Otherwise leave the selector visibly unselected.
  default_column <- function(columns, candidates) {
    positions <- match(
      tolower(trimws(candidates)),
      tolower(trimws(columns))
    )

    positions <- positions[!is.na(positions)]

    if (length(positions)) columns[positions[1]] else ""
  }

  pairs_cols <- names(g)

  pair_choices <- c(
    "Select a column..." = "",
    stats::setNames(pairs_cols, pairs_cols)
  )

  peak_ui <- if (is.null(input$file)) {
    # raw_df() handles the actual peak-table input below.
    NULL
  } else {
    NULL
  }

  # Allow pairs-column selection even before a peak table
  # has been uploaded.
  peak <- tryCatch(
    raw_df(),
    shiny.silent.error = function(e) NULL
  )

  if (is.null(peak)) {
    peak_ui <- div(
      class = "small-note",
      "Upload a peak table to select its matching feature ID column."
    )
  } else {
    peak_cols <- names(peak)

    default_peak <- input$annotation_id_col

    if (
      length(default_peak) != 1L ||
      !default_peak %in% peak_cols
    ) {
      default_peak <- default_column(
        peak_cols,
        c(
          "row ID",
          "alignment id",
          "feature_id",
          "feature id",
          "id"
        )
      )
    }

    peak_ui <- selectInput(
      inputId = "component_match_col",
      label = "Peak-table feature ID column:",
      choices = c(
        "Select a column..." = "",
        stats::setNames(peak_cols, peak_cols)
      ),
      selected = default_peak
    )
  }

  tagList(
    peak_ui,

    selectInput(
      inputId = "main_pairs_node1_col",
      label = "Pairs: first node ID column:",
      choices = pair_choices,
      selected = default_column(pairs_cols, "CLUSTERID1")
    ),

    selectInput(
      inputId = "main_pairs_node2_col",
      label = "Pairs: second node ID column:",
      choices = pair_choices,
      selected = default_column(pairs_cols, "CLUSTERID2")
    ),

    selectInput(
      inputId = "main_pairs_component_col",
      label = "Pairs: ComponentIndex column:",
      choices = pair_choices,
      selected = default_column(pairs_cols, "ComponentIndex")
    )
  )
})

main_component_map <- reactive({
  g <- main_gnps_pairs_df()

  columns <- c(
    input$main_pairs_node1_col,
    input$main_pairs_node2_col,
    input$main_pairs_component_col
  )

  validate(
    need(
      length(columns) == 3L &&
        all(nzchar(columns)) &&
        all(columns %in% names(g)),
      "Select both node ID columns and the ComponentIndex column."
    ),
    need(
      length(unique(columns)) == 3L,
      "Choose three different columns for the two node IDs and ComponentIndex."
    )
  )

  component <- trimws(as.character(g[[columns[3]]]))

  dplyr::bind_rows(
    tibble::tibble(
      ClusterID = trimws(as.character(g[[columns[1]]])),
      ComponentIndex = component
    ),
    tibble::tibble(
      ClusterID = trimws(as.character(g[[columns[2]]])),
      ComponentIndex = component
    )
  ) %>%
    dplyr::filter(
      !is.na(ClusterID),
      nzchar(ClusterID),
      !is.na(ComponentIndex),
      nzchar(ComponentIndex)
    ) %>%
    dplyr::distinct(ClusterID, ComponentIndex)
})


# Map the original peak-table column to internal feature names.
main_component_feature_ids <- reactive({
  req(procReady(), rv$fmap)

  peak <- raw_df()
  column <- input$component_match_col

  validate(
    need(
      length(column) == 1L &&
        nzchar(column) &&
        column %in% names(peak),
      "Select the peak-table feature ID column in Load & Process."
    ),
    need(
      nrow(peak) == nrow(rv$fmap),
      "The peak table has changed. Run Process again."
    )
  )

  tibble::tibble(
    Feature = as.character(rv$fmap$Feature),
    ClusterID = trimws(as.character(peak[[column]]))
  ) %>%
    dplyr::filter(
      !is.na(ClusterID),
      nzchar(ClusterID)
    ) %>%
    dplyr::distinct()
})

# Turning off the optional source also disables its filter.
observeEvent(input$use_main_gnps_pairs, {

  if (!isTRUE(input$use_main_gnps_pairs)) {

    shinyWidgets::updateMaterialSwitch(
  session,
  inputId = "use_component_filter",
  value = FALSE
)
  }

}, ignoreInit = TRUE)
  
# Network annotations, independent of whether filtering is enabled.
gnps_network_annotations <- reactive({

  empty <- tibble::tibble(
    Feature = character(),
    GNPS_ClusterID = character(),
    GNPS_ComponentIndex = character()
  )

  if (
    !isTRUE(input$use_main_gnps_pairs) ||
    is.null(input$file_main_gnps_pairs) ||
    !isTRUE(procReady())
  ) {
    return(empty)
  }

  tryCatch({

    ids <- main_component_feature_ids()
    mapping <- main_component_map()

    joined <- merge(
      ids,
      mapping,
      by = "ClusterID",
      all = FALSE,
      sort = FALSE
    )

    collapse_ids <- function(x) {
      x <- unique(as.character(x))
      x <- x[!is.na(x) & nzchar(x)]

      paste(
        stringr::str_sort(x, numeric = TRUE),
        collapse = ", "
      )
    }

    joined %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarise(
        GNPS_ClusterID = collapse_ids(ClusterID),
        GNPS_ComponentIndex = collapse_ids(ComponentIndex),
        .groups = "drop"
      )

  }, shiny.silent.error = function(e) empty)
})


# Add columns without changing row order or duplicating features.
add_gnps_network_annotations <- function(dd) {

  annotations <- gnps_network_annotations()

  index <- match(
    as.character(dd$Feature),
    annotations$Feature
  )

  dd$GNPS_ClusterID <- annotations$GNPS_ClusterID[index]
  dd$GNPS_ComponentIndex <- annotations$GNPS_ComponentIndex[index]

  dd
}

  # ---- Filtered volcano data
  filter_volcano_data <- function(exclude = character(), for_choices = FALSE) {
    req(procReady(), input$mz_range, input$rt_range, input$intensity_range)
  
    dd <- rv$volcano
    
    # Keep only features with a GNPS annotation.
if (isTRUE(input$use_gnps_annotated_filter)) {

  dd <- dd %>%
    dplyr::filter(
      !is.na(clean_missing_text(GNPS_annotation))
    )
}

# Keep only features with an Other annotation.
if (isTRUE(input$use_other_annotated_filter)) {

  dd <- dd %>%
    dplyr::filter(
      !is.na(clean_missing_text(Other_annotation))
    )
}
    
if (
  !"component" %in% exclude &&
  isTRUE(input$use_main_gnps_pairs) &&
  isTRUE(input$use_component_filter) &&
  (!for_choices || length(input$volcano_components) > 0)
) {

  mapping <- main_component_map()

validate(
  need(
    length(input$volcano_components) > 0,
    "Select at least one ComponentIndex."
  )
)

keep_ids <- unique(
  mapping$ClusterID[
    mapping$ComponentIndex %in%
      as.character(input$volcano_components)
  ]
)

feature_ids <- main_component_feature_ids()

keep_features <- feature_ids$Feature[
  feature_ids$ClusterID %in% keep_ids
]

dd <- dd[
  as.character(dd$Feature) %in% keep_features,
  ,
  drop = FALSE
]

  if (!for_choices) {
  validate(
    need(
      nrow(dd) > 0,
      paste(
        "No features match these components.",
        "Check the selected ID column."
      )
    )
  )
}
}
    pcut <- suppressWarnings(as.numeric(input$sig_p_cutoff %||% 0.05))
    if (!is.finite(pcut) || pcut < 0 || pcut > 1) pcut <- 0.05
    
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
  
if (
  !"npc" %in% exclude &&
  isTRUE(input$use_npc_filter) &&
  (!for_choices || length(input$npc_filter_values) > 0)
) {

  validate(
    need(
      length(input$npc_filter_values) > 0,
      "NPC filter is enabled. Select at least one NPC class."
    )
  )

  dd <- dd %>%
    dplyr::filter(
      `NPC#class` %in% input$npc_filter_values
    )

  # Applies only while the NPC class filter is active.
  dd <- apply_sirius_probability(
    dd,
    prefix = "npc",
    family = "NPC"
  )
}


if (
  !"classyfire" %in% exclude &&
  isTRUE(input$use_classyfire_filter) &&
  (!for_choices || length(input$classyfire_filter_values) > 0)
) {

  validate(
    need(
      length(input$classyfire_filter_values) > 0,
      "ClassyFire filter is enabled. Select at least one ClassyFire class."
    )
  )

  dd <- dd %>%
    dplyr::filter(
      `ClassyFire#class` %in% input$classyfire_filter_values
    )

  # Applies only while the ClassyFire class filter is active.
  dd <- apply_sirius_probability(
    dd,
    prefix = "classyfire",
    family = "ClassyFire"
  )
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
  }

  filtered_volcano <- reactive({
    filter_volcano_data()
  })

  # ------------------------------------------------------------
# Dynamic annotation choices and unique-feature counts
# ------------------------------------------------------------

# Remember the last update to avoid repeatedly sending
# identical choices back to the browser.
volcano_picker_cache <- new.env(parent = emptyenv())


update_counted_picker <- function(input_id, pairs) {

  # pairs must contain Feature and Value.
  pairs <- pairs %>%
    dplyr::transmute(
      Feature = as.character(Feature),
      Value = as.character(Value)
    ) %>%
    dplyr::filter(
      !is.na(Feature),
      nzchar(Feature),
      !is.na(Value),
      nzchar(trimws(Value))
    ) %>%
    dplyr::distinct(Feature, Value)

  counts <- pairs %>%
    dplyr::count(Value, name = "n")

  selected <- as.character(input[[input_id]])
  selected <- unique(
    selected[!is.na(selected) & nzchar(selected)]
  )

  # Keep selected values even when other filters remove
  # all their matching features.
  values <- stringr::str_sort(
    unique(c(counts$Value, selected)),
    numeric = TRUE
  )

  numbers <- counts$n[match(values, counts$Value)]
  numbers[is.na(numbers)] <- 0L

 # An empty annotation list is valid.
# Explicitly clear the picker instead of creating a label.
if (length(values) == 0L) {

  choices <- character(0)

} else {

  labels <- paste0(
    values,
    " (",
    numbers,
    ifelse(numbers == 1L, " feature", " features"),
    ifelse(
      numbers == 0L & values %in% selected,
      " — selected",
      ""
    ),
    ")"
  )

  choices <- stats::setNames(values, labels)
}

  state <- list(
    choices = choices,
    selected = selected
  )

  previous <- volcano_picker_cache[[input_id]]

  if (!identical(previous, state)) {
    shinyWidgets::updatePickerInput(
      session = session,
      inputId = input_id,
      choices = choices,
      selected = selected
    )

    volcano_picker_cache[[input_id]] <- state
  }
}


# Clear cached updates when processed data changes.
observeEvent(rv$volcano, {
  keys <- ls(
    envir = volcano_picker_cache,
    all.names = TRUE
  )

  if (length(keys)) {
    rm(
      list = keys,
      envir = volcano_picker_cache
    )
  }
}, priority = 100)


# NPC: apply all filters except the NPC filter itself.
observe({
  req(procReady())

  dd <- filter_volcano_data(
    exclude = "npc",
    for_choices = TRUE
  )

  req("NPC#class" %in% names(dd))

  update_counted_picker(
    input_id = "npc_filter_values",
    pairs = tibble::tibble(
      Feature = dd$Feature,
      Value = dd[["NPC#class"]]
    )
  )
})


# ClassyFire: apply all filters except ClassyFire itself.
observe({
  req(procReady())

  dd <- filter_volcano_data(
    exclude = "classyfire",
    for_choices = TRUE
  )

  req("ClassyFire#class" %in% names(dd))

  update_counted_picker(
    input_id = "classyfire_filter_values",
    pairs = tibble::tibble(
      Feature = dd$Feature,
      Value = dd[["ClassyFire#class"]]
    )
  )
})


# ComponentIndex: apply all filters except ComponentIndex.
observe({
  req(
    procReady(),
    isTRUE(input$use_main_gnps_pairs),
    input$file_main_gnps_pairs
  )

  dd <- filter_volcano_data(
    exclude = "component",
    for_choices = TRUE
  )

  feature_ids <- main_component_feature_ids() %>%
  dplyr::filter(
    Feature %in% as.character(dd$Feature)
  )

  mapping <- main_component_map() %>%
    dplyr::transmute(
      ClusterID = trimws(as.character(ClusterID)),
      Value = as.character(ComponentIndex)
    ) %>%
    dplyr::distinct()

  # Merge supports features belonging to multiple components.
  pairs <- merge(
    feature_ids,
    mapping,
    by = "ClusterID",
    all = FALSE,
    sort = FALSE
  )

  update_counted_picker(
    input_id = "volcano_components",
    pairs = pairs
  )
})
  
  output$volcano_feature_count <- renderText({
  req(procReady())

  dd <- filtered_volcano()

  count_features <- function(x) {
    x <- as.character(x)
    length(unique(x[!is.na(x) & nzchar(x)]))
  }

  current <- count_features(dd$Feature)
  total <- count_features(rv$volcano$Feature)

  paste0(
    "Features: ",
    format(current, big.mark = ",", trim = TRUE),
    " / ",
    format(total, big.mark = ",", trim = TRUE),
    " retained"
  )
})


output$volcano_applied_filters <- renderUI({
  req(procReady())

  steps <- character()

  add_step <- function(text) {
    steps <<- c(steps, text)
  }

  fmt <- function(x) {
    format(signif(as.numeric(x), 5), trim = TRUE)
  }

  show_values <- function(x) {
    paste(as.character(x), collapse = ", ")
  }

  pcut <- suppressWarnings(
    as.numeric(input$sig_p_cutoff %||% 0.05)
  )
  if (!is.finite(pcut) || pcut < 0 || pcut > 1) {
    pcut <- 0.05
  }

  fcut <- suppressWarnings(
    as.numeric(input$fc_thr %||% 1)
  )
  if (!is.finite(fcut) || fcut < 0) {
    fcut <- 1
  }

  if (isTRUE(input$use_gnps_annotated_filter)) {
    add_step("GNPS: retain features with a non-missing annotation.")
  }

  if (isTRUE(input$use_other_annotated_filter)) {
    add_step("Other annotation: retain features with a non-missing annotation.")
  }

  if (
    isTRUE(input$use_main_gnps_pairs) &&
    isTRUE(input$use_component_filter)
  ) {
    add_step(paste0(
      "GNPS ComponentIndex: ",
      if (length(input$volcano_components)) {
        show_values(input$volcano_components)
      } else {
        "selection required"
      },
      "; matching column: ",
      input$component_match_col %||% "not selected",
      "."
    ))
  }

  if (length(input$sel_feat)) {
    add_step(paste0(
      "Manual feature selection: ",
      show_values(input$sel_feat),
      "."
    ))
  }

  if (isTRUE(input$sig_only)) {
    add_step(paste0(
      "Significant features only: adjusted p-value ≤ ",
      fmt(pcut),
      " and |log2(FC)| ≥ ",
      fmt(fcut),
      "."
    ))
  }

  for (prefix in c("npc", "classyfire")) {
    if (!isTRUE(input[[paste0("use_", prefix, "_filter")]])) {
      next
    }

    label <- if (prefix == "npc") "NPC" else "ClassyFire"
    values <- input[[paste0(prefix, "_filter_values")]]

    add_step(paste0(
      label, " classes: ",
      if (length(values)) show_values(values) else "selection required",
      "."
    ))

    if (isTRUE(input[[paste0("use_", prefix, "_probability")]])) {
      score_col <- input[[paste0(prefix, "_probability_col")]]
      cutoff <- input[[paste0(prefix, "_probability_min")]]

      add_step(paste0(
        label, " SIRIUS probability: ",
        sub("^SIRIUS_", "", score_col %||% "column not selected"),
        " > ",
        fmt(cutoff %||% 0.9),
        "."
      ))
    }
  }

  # These ranges are always applied by filtered_volcano().
  ranges <- list(
    "m/z" = input$mz_range,
    "Retention time" = input$rt_range,
    "Intensity [log10(Mean + 1.1)]" = input$intensity_range
  )

  for (label in names(ranges)) {
    limits <- ranges[[label]]

    if (length(limits) == 2L) {
      add_step(paste0(
        label, ": ",
        fmt(limits[1]), " to ", fmt(limits[2]),
        " (inclusive)."
      ))
    }
  }

  add_step("Features require finite m/z, retention time and mean intensity.")

  if (isTRUE(input$use_fdr_filter)) {
    add_step(paste0(
      "Adjusted p-value ≤ ", fmt(pcut), "."
    ))
  }

  if (isTRUE(input$use_fc_filter)) {
    direction <- input$fc_dir %||% "both"

    rule <- switch(
      direction,
      both = paste0("|log2(FC)| ≥ ", fmt(fcut)),
      up = paste0("log2(FC) ≥ ", fmt(fcut)),
      paste0("log2(FC) ≤ ", fmt(-fcut))
    )

    add_step(paste0("Fold-change filter: ", rule, "."))
  }

  tags$ul(
    style = paste(
      "margin-top:10px;",
      "padding-left:22px;",
      "max-height:250px;",
      "overflow-y:auto;"
    ),
    lapply(steps, function(step) tags$li(step))
  )
})
  
  # ============================================================
# Interactive heatmap: current filtered feature count
# ============================================================

output$heatmap_has_data <- renderText({

  if (!isTRUE(procReady())) {
    return("no")
  }

  dd <- tryCatch(
    filtered_volcano(),
    shiny.silent.error = function(e) NULL
  )

  if (is.null(dd) || nrow(dd) == 0L) {
    return("no")
  }

  "yes"
})

outputOptions(
  output,
  "heatmap_has_data",
  suspendWhenHidden = FALSE
)

output$heatmap_filter_summary <- renderUI({

  req(
    procReady(),
    filtered_volcano()
  )

  dd <- filtered_volcano()

  features <- unique(
    as.character(
      dd$Feature
    )
  )

  div(
    class = "small-note",

    HTML(
      paste0(
        "<b>",
        format(
          length(features),
          big.mark = ","
        ),
        "</b> feature(s) retained by the current Volcano filters."
      )
    )
  )
})
  
# ============================================================
# Heatmap data
#
# ROWS    = features/metabolites
# COLUMNS = samples
# LABEL   = sample annotation only
# ============================================================

heatmap_data <- reactive({

  req(
    procReady(),
    rv$df_used,
    rv$mat,
    filtered_volcano()
  )


  # ----------------------------------------------------------
  # Features surviving the current Volcano filters
  # ----------------------------------------------------------

  features <- unique(
    as.character(
      filtered_volcano()$Feature
    )
  )

  features <- intersect(
    features,
    colnames(rv$df_used)
  )


  validate(
    need(
      length(features) > 0,
      "No features remain after filtering."
    )
  )


  # ----------------------------------------------------------
  # Original matrix in Metabocano:
  #
  # samples x features
  # ----------------------------------------------------------

  mat <- as.matrix(
    rv$df_used[
      ,
      features,
      drop = FALSE
    ]
  )

  storage.mode(mat) <- "double"


  # Sample names
  sample_names <- rownames(rv$mat)

  if (
    length(sample_names) ==
      nrow(mat)
  ) {
    rownames(mat) <- sample_names
  }


  # ----------------------------------------------------------
  # TRANSPOSE
  #
  # Now:
  # rows    = features
  # columns = samples
  # ----------------------------------------------------------

  mat <- t(mat)


  # ----------------------------------------------------------
  # Scaling
  #
  # Scale each FEATURE across samples.
  #
  # Because features are rows now, use:
  #
  # t(scale(t(mat), ...))
  # ----------------------------------------------------------

   scale_mode <- input$hm_scale %||% "uv"


  if (identical(scale_mode, "uv")) {

  row_sd <- apply(mat, 1L, stats::sd, na.rm = TRUE)

  scalable <- is.finite(row_sd) & row_sd > 0

  if (any(scalable)) {
    mat[scalable, ] <- sweep(
      mat[scalable, , drop = FALSE],
      MARGIN = 1L,
      STATS = row_sd[scalable],
      FUN = "/"
    )
  }

  # Constant rows cannot be scaled to unit variance.
  mat[!scalable, ] <- 0

} else if (identical(scale_mode, "zscore")) {

    mat <- t(
      scale(
        t(mat),
        center = TRUE,
        scale = TRUE
      )
    )
  }


  # Protect against zero-variance features
  mat[!is.finite(mat)] <- 0


  # ----------------------------------------------------------
  # Sample annotation
  #
  # Label is NOT heatmap data.
  # It only colors the sample columns.
  # ----------------------------------------------------------

  annotation_col <- data.frame(

    Group = factor(
      as.character(
        rv$df_used$Label
      )
    ),

    stringsAsFactors = FALSE
  )


  rownames(annotation_col) <- colnames(mat)


  list(
    mat = mat,
    annotation_col = annotation_col,
    features = features
  )
})

  # ============================================================
# Generate interactive heatmap
# ============================================================

# ============================================================
# Generate / update InteractiveComplexHeatmap
# ============================================================

observeEvent({

  req(
    procReady(),
    isTRUE(input$show_interactive_heatmap)
  )

  list(
    input$show_interactive_heatmap,
    filtered_volcano(),
    gnps_network_annotations(),
    rv$df_used,
    rv$mat,
    input$hm_scale,
    input$hm_distance,
    input$hm_method,
    input$hm_palette,
    input$hm_group_palette,
    input$hm_cluster_features,
    input$hm_show_borders,
    input$hm_cluster_samples,
    input$hm_show_features,
    input$hm_show_samples
  )

}, {

  req(
    isTRUE(input$show_interactive_heatmap),

    input$hm_scale,
    input$hm_distance,
    input$hm_method,
    input$hm_palette,
    input$hm_group_palette
  )


  # ----------------------------------------------------------
  # Prepared data
  # ----------------------------------------------------------

  hm <- heatmap_data()

  mat <- hm$mat
  annotation_col <- hm$annotation_col


  req(
    nrow(mat) > 0,
    ncol(mat) > 0
  )


  # ----------------------------------------------------------
  # Annotation colors
  # ----------------------------------------------------------

  groups <- unique(
    as.character(
      annotation_col$Group
    )
  )

  groups <- groups[
    !is.na(groups) &
      nzchar(groups)
  ]


  group_cols <- make_palette(
    input$hm_group_palette %||% "Dark2",
    length(groups)
  )

  names(group_cols) <- groups


  annotation_colors <- list(
    Group = group_cols
  )


  # ----------------------------------------------------------
  # Heatmap colors
  # ----------------------------------------------------------

  heatmap_colors <- switch(

    input$hm_palette %||% "viridis",

    "magma" =
      viridisLite::magma(100),

    "bwr" =
      grDevices::colorRampPalette(
        c(
          "#2166AC",
          "white",
          "#B2182B"
        )
      )(100),

    viridisLite::viridis(100)
  )


  # ----------------------------------------------------------
  # Clustering
  # ----------------------------------------------------------

  cluster_features <-
    isTRUE(input$hm_cluster_features) &&
    nrow(mat) > 1


  cluster_samples <-
    isTRUE(input$hm_cluster_samples) &&
    ncol(mat) > 1


  distance <-
    input$hm_distance %||%
    "euclidean"


  method <-
    input$hm_method %||%
    "ward.D2"


  if (
    identical(method, "ward.D2") &&
    !identical(distance, "euclidean")
  ) {
    distance <- "euclidean"
  }


  # ----------------------------------------------------------
  # Legend
  # ----------------------------------------------------------

  heatmap_title <- switch(

    input$hm_scale %||% "uv",

    "zscore" = "Z-score",

    "none" = "Intensity",

    "Scaled intensity"
  )


  # ----------------------------------------------------------
  # Construct heatmap
  #
  # ROWS    = features
  # COLUMNS = samples
  # ----------------------------------------------------------

  ht <- ComplexHeatmap::pheatmap(

    mat,

    color = heatmap_colors,


    # sample groups
    annotation_col = annotation_col,

    annotation_colors =
      annotation_colors,

    annotation_names_col =
      TRUE,


    # rows = features
    cluster_rows =
      cluster_features,

    # columns = samples
    cluster_cols =
      cluster_samples,


    show_rownames =
      isTRUE(input$hm_show_features),

    show_colnames =
      isTRUE(input$hm_show_samples),


    border_color = if (isTRUE(input$hm_show_borders)) "grey60" else NA,


    clustering_distance_rows =
      distance,

    clustering_distance_cols =
      distance,

    clustering_method =
      method,


    name =
      heatmap_title
  )


  # ----------------------------------------------------------
  # Connect to the PERMANENT heatmap UI
  # ----------------------------------------------------------

  ht <- local({

  grDevices::pdf(
    file = NULL,
    width = 10,
    height = 8
  )

  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  ComplexHeatmap::draw(ht)
})
  
  # Capture data belonging to this heatmap.
# Click indices refer to the original matrix, before clustering.
heatmap_feature_data <- as.data.frame(
  add_gnps_network_annotations(rv$volcano),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

heatmap_intensities <- rv$df_used

output$heatmap_feature_info <- shiny::renderUI({
  div(
    class = "small-note",
    "Click a heatmap cell to see feature and sample information."
  )
})


# Format values without interpreting annotation text as HTML.
format_heatmap_info <- function(x) {
  if (!length(x)) return("NA")

  if (is.numeric(x)) {
    return(paste(
      ifelse(
        is.na(x),
        "NA",
        format(signif(x, 6), trim = TRUE)
      ),
      collapse = ", "
    ))
  }

  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- "NA"

  paste(x, collapse = ", ")
}


# Two-column table used for annotations and statistics.
heatmap_info_table <- function(row, columns) {
  columns <- intersect(columns, names(row))

  tags$table(
    class = "table table-condensed table-striped",
    style = "width:100%;",

    tags$tbody(
      lapply(columns, function(column) {
        label <- sub(
          "^(Peak_|SIRIUS_|GNPS_|Other_)",
          "",
          column
        )

        if (column == "GNPS_ClusterID") {
          label <- "GNPS ClusterID"
        }
        
        if (column == "GNPS_ComponentIndex") {
          label <- "GNPS ComponentIndex"
        }
        
        if (column == "GNPS_annotation") {
          label <- "GNPS annotation"
        }

        if (column == "Other_annotation") {
          label <- "Other annotation"
        }

        tags$tr(
          tags$th(
            style = "vertical-align:top;white-space:normal;",
            label
          ),
          tags$td(
            style = "overflow-wrap:anywhere;white-space:normal;",
            format_heatmap_info(row[[column]])
          )
        )
      })
    )
  )
}


heatmap_click_action <- function(df, output) {

  if (is.null(df) || nrow(df) == 0L) {
    output$heatmap_feature_info <- shiny::renderUI({
      div(
        class = "small-note",
        "Click inside the heatmap body."
      )
    })
    return(invisible(NULL))
  }

  row_index <- as.integer(
    unlist(df$row_index, use.names = FALSE)
  )

  column_index <- as.integer(
    unlist(df$column_index, use.names = FALSE)
  )

  if (!length(row_index) || !length(column_index)) {
    return(invisible(NULL))
  }

  i <- row_index[1]
  j <- column_index[1]

  if (
    is.na(i) || is.na(j) ||
    i < 1L || i > nrow(mat) ||
    j < 1L || j > ncol(mat)
  ) {
    return(invisible(NULL))
  }

  feature <- rownames(mat)[i]
  sample <- colnames(mat)[j]

  feature_rows <- heatmap_feature_data[
    as.character(heatmap_feature_data$Feature) == feature,
    ,
    drop = FALSE
  ]

  # Comparison-dependent columns are displayed separately.
  stats_columns <- c(
    "Groups",
    "Group_num",
    "Group_den",
    "FC",
    "Adj.p-value",
    "Adj.p-value.log",
    "Mean",
    "mean_num",
    "mean_den",
    "TestScale",
    "Significant_default"
  )

  annotation_columns <- setdiff(
    names(feature_rows),
    stats_columns
  )

  processed_intensity <- if (
    feature %in% names(heatmap_intensities) &&
    j <= nrow(heatmap_intensities)
  ) {
    heatmap_intensities[[feature]][j]
  } else {
    NA_real_
  }

  output$heatmap_feature_info <- shiny::renderUI({

    div(
      style = paste(
        "max-height:500px;",
        "overflow-y:auto;",
        "overflow-wrap:anywhere;",
        "padding:8px;"
      ),

      tags$h4("Selected heatmap cell"),

      tags$p(
        tags$strong("Feature: "),
        feature
      ),

      tags$p(
        tags$strong("Sample: "),
        sample
      ),

      tags$p(
        tags$strong("Sample group: "),
        as.character(annotation_col$Group[j])
      ),

      tags$p(
        tags$strong("Processed intensity before heatmap scaling: "),
        format_heatmap_info(processed_intensity)
      ),

      tags$p(
        tags$strong("Displayed heatmap value: "),
        format_heatmap_info(mat[i, j])
      ),

      if (!nrow(feature_rows)) {
        div(
          class = "small-note",
          "No feature annotations matched this heatmap row."
        )
      } else {
        tagList(

          tags$details(
            tags$summary(
              style = "cursor:pointer;color:#228B22;font-weight:600;",
              "All feature information"
            ),

            heatmap_info_table(
              feature_rows[1, , drop = FALSE],
              annotation_columns
            )
          ),

          tags$details(
            style = "margin-top:12px;",

            tags$summary(
              style = "cursor:pointer;color:#228B22;font-weight:600;",
              "All processed comparisons"
            ),

            lapply(seq_len(nrow(feature_rows)), function(k) {
              current <- feature_rows[k, , drop = FALSE]

              tagList(
                tags$h5(
                  tags$strong(
                    as.character(current$Groups[1])
                  )
                ),

                heatmap_info_table(
                  current,
                  setdiff(stats_columns, "Groups")
                )
              )
            })
          )
        )
      }
    )
  })
}


# Keep brushing useful with the custom information panel.
heatmap_brush_action <- function(df, output) {

  output$heatmap_feature_info <- shiny::renderUI({

    if (is.null(df) || nrow(df) == 0L) {
      return(
        div(class = "small-note", "No heatmap cells selected.")
      )
    }

    rows <- unique(
      as.integer(unlist(df$row_index, use.names = FALSE))
    )

    columns <- unique(
      as.integer(unlist(df$column_index, use.names = FALSE))
    )

    rows <- rows[
      !is.na(rows) & rows >= 1L & rows <= nrow(mat)
    ]

    columns <- columns[
      !is.na(columns) & columns >= 1L & columns <= ncol(mat)
    ]

    tagList(
      tags$p(
        tags$strong("Selected region: "),
        paste0(
          length(rows), " features × ",
          length(columns), " samples."
        )
      ),

      tags$p(
        class = "small-note",
        "Click a cell in either heatmap to inspect its feature information."
      )
    )
  })
}


InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
  input = input,
  output = output,
  session = session,

  ht_list = ht,
  heatmap_id = "metabocano_heatmap",

  click_action = heatmap_click_action,
  brush_action = heatmap_brush_action
)

}, ignoreInit = FALSE)
  
  # ---- Volcano plot
  output$volcano_plot <- renderPlotly({
    req(filtered_volcano())
    dd <- add_gnps_network_annotations(filtered_volcano())
    validate(need(nrow(dd) > 0, "No points left after filtering."))

    dd$key <- paste(dd$Groups, dd$Feature, sep = "__")

    use_mean_y <- identical(
        input$volcano_y_axis %||% "fdr",
        "mean"
      )
      
      dd$plot_y <- if (use_mean_y) {
        log10(dd$Mean + 1.1)
      } else {
        dd$`Adj.p-value.log`
      }
      
      y_title <- if (use_mean_y) {
        "Mean log10(Intensity)"
      } else {
        "-log10(FDR)"
      }
    
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
      "<br>ClassyFire: ", dd$`ClassyFire#class`,
      "<br>GNPS annotation: ", dd$GNPS_annotation,
      "<br>GNPS ClusterID: ", htmltools::htmlEscape(ifelse(is.na(dd$GNPS_ClusterID), "NA", dd$GNPS_ClusterID)),
      "<br>GNPS ComponentIndex: ", htmltools::htmlEscape(ifelse(is.na(dd$GNPS_ComponentIndex),"NA", dd$GNPS_ComponentIndex)),
      "<br>Other annotation: ", dd$Other_annotation
    )

    fc_line <- suppressWarnings(as.numeric(input$fc_thr %||% 1))
    if (!is.finite(fc_line) || fc_line < 0) fc_line <- 1
    
    p_thr <- suppressWarnings(as.numeric(input$sig_p_cutoff %||% 0.05))
    if (!is.finite(p_thr) || p_thr < 0 || p_thr > 1) p_thr <- 0.05
    
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
  )
)

# Add the FDR cutoff line only when FDR is used as y-axis
if (!use_mean_y) {
  shapes <- c(
    shapes,
    list(
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
  )
}

  p <- if (identical(input$color_by, "FC")) {

  # FC contains log2 fold changes.
  # A colored point must pass BOTH thresholds.
  dd$FC_status <- "Other"

  passes_p <- is.finite(dd$`Adj.p-value`) &
    dd$`Adj.p-value` <= p_thr

  up <- passes_p &
    is.finite(dd$FC) &
    dd$FC >= fc_line &
    dd$FC > 0

  down <- passes_p &
    is.finite(dd$FC) &
    dd$FC <= -fc_line &
    dd$FC < 0

  dd$FC_status[which(up)] <- "Upregulated"
  dd$FC_status[which(down)] <- "Downregulated"

  dd$FC_status <- factor(
    dd$FC_status,
    levels = c("Other", "Downregulated", "Upregulated")
  )

  plot_ly(
    data = dd,
    x = ~FC,
    y = ~plot_y,
    color = ~FC_status,
    colors = c(
      "Other" = "grey70",
      "Downregulated" = "blue",
      "Upregulated" = "red"
    ),
    type = "scatter",
    mode = "markers",
    text = hover_txt,
    hoverinfo = "text",
    key = ~key,
    marker = list(
      size = 12,
      opacity = 0.85,
      line = list(color = "black", width = 1)
    ),
    source = "volcano"
  ) %>%
    layout(
      shapes = shapes,
      xaxis = list(title = "log2(FC)"),
      yaxis = list(title = y_title),
      legend = list(
        title = list(text = "Regulation")
      )
    ) %>%
    event_register("plotly_click")

} else if (input$color_by == "Groups") {
      plot_ly(
        data = dd,
        x = ~FC, y = ~plot_y,
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
          yaxis = list(title = y_title),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    } else {
      plot_ly(
        data = dd,
        x = ~FC, y = ~plot_y,
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
          yaxis = list(title = y_title),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    }
  
  add_volcano_top_labels(
  p = p,
  dd = dd,
  n = input$volcano_top_n %||% 0,
  width_px =
    session$clientData$output_volcano_plot_width %||% 800,
  height_px = 520,
  label_col = input$volcano_label_column %||% "Feature"
)
  
  })

  # ---- Feature plot on click
  output$selected_feature_panel <- renderUI({
  req(procReady(), rv$volcano)

  click <- plotly::event_data(
    "plotly_click",
    source = "volcano"
  )

  if (
    is.null(click) ||
    is.null(click$key) ||
    !length(click$key)
  ) {
    return(NULL)
  }

  key <- as.character(click$key[[1]])

  parts <- stringr::str_split_fixed(
    key,
    "__",
    2
  )

  comp <- parts[1, 1]
  feat <- parts[1, 2]

  row <- rv$volcano %>%
    dplyr::filter(
      Groups == comp,
      Feature == feat
    ) %>%
    dplyr::slice(1)

  if (nrow(row) == 0) {
    return(NULL)
  }

  # Prefer the original peak-table ID.
  # Fall back to the internal Feature value.
  feature_text <- feat

  if (
    "id" %in% names(row) &&
    !is.na(row$id[[1]]) &&
    nzchar(trimws(as.character(row$id[[1]])))
  ) {
    feature_text <- as.character(row$id[[1]])
  }

  mz_value <- suppressWarnings(
    as.numeric(row$mz[[1]])
  )

  rt_value <- suppressWarnings(
    as.numeric(row$RT[[1]])
  )

  mz_text <- if (is.finite(mz_value)) {
    format(
      mz_value,
      digits = 12,
      scientific = FALSE,
      trim = TRUE
    )
  } else {
    "NA"
  }

  rt_text <- if (is.finite(rt_value)) {
    format(
      rt_value,
      digits = 8,
      scientific = FALSE,
      trim = TRUE
    )
  } else {
    "NA"
  }

  # Include network IDs in the existing annotation details.
row <- add_gnps_network_annotations(row)

detail_cols <- names(row)

# Readable labels for standard columns.
friendly_labels <- c(
  "Feature" = "Feature name",
  "id" = "Original feature ID",
  "mz" = "m/z",
  "RT" = "Retention time",
  "Groups" = "Comparison",
  "Group_num" = "Numerator group",
  "Group_den" = "Denominator group",
  "FC" = "log2 fold change",
  "Adj.p-value" = "Adjusted p-value (FDR)",
  "Adj.p-value.log" = "-log10(FDR)",
  "Mean" = "Mean intensity",
  "mean_num" = "Mean intensity: numerator",
  "mean_den" = "Mean intensity: denominator",
  "TestScale" = "Statistical test scale",
  "Significant_default" = "Significant at default thresholds",
  "NPC#class" = "NPC class",
  "ClassyFire#class" = "ClassyFire class",
  "GNPS_annotation" = "GNPS annotation",
  "GNPS_ClusterID" = "GNPS ClusterID",
  "GNPS_ComponentIndex" = "GNPS ComponentIndex",
  "Other_annotation" = "Other annotation"
)

additional_details_ui <- tags$details(

  style = "
    margin-top: 14px;
    border-top: 1px solid #dddddd;
    padding-top: 10px;
  ",

  tags$summary(
    style = "
      cursor: pointer;
      font-weight: bold;
      color: #2c3e50;
      padding: 6px 0;
    ",
    "Click to expand"
  ),

  div(
    style = "
      margin-top: 8px;
      max-height: 450px;
      overflow-y: auto;
    ",

    lapply(detail_cols, function(column_name) {

      display_name <- if (
        column_name %in% names(friendly_labels)
      ) {
        unname(friendly_labels[[column_name]])
      } else {
        # Keep additional annotation column names recognizable.
        sub(
          "^(Peak|SIRIUS|GNPS|Other)_",
          "\\1: ",
          column_name
        )
      }

      div(
        style = "
          display: grid;
          grid-template-columns: minmax(160px, 35%) minmax(0, 1fr);
          gap: 10px;
          padding: 6px 0;
          border-bottom: 1px solid #eeeeee;
          overflow-wrap: anywhere;
        ",

        tags$strong(
          paste0(display_name, ":")
        ),

        tags$span(
          style = "white-space: pre-wrap; user-select: text;",
          format_extra_value(row[[column_name]])
        )
      )
    })
  )
)
  
  tagList(
    div(
      style = "
        background-color: white;
        border: 2px solid #66CDAA;
        border-radius: 8px;
        padding: 12px 15px;
        margin-bottom: 12px;
      ",

      tags$h4(
        style = "
          margin-top: 0;
          margin-bottom: 12px;
          font-weight: bold;
          color: #2c3e50;
        ",
        "Feature ID: ",
        tags$span(
          style = "color: #18bc9c;",
          feature_text
        )
      ),

      fluidRow(
        column(
          width = 6,

          tags$label(
            `for` = "selected_mz_text",
            style = "font-weight: bold;",
            "m/z:"
          ),

          tags$input(
            id = "selected_mz_text",
            type = "text",
            class = "form-control",
            value = mz_text,
            readonly = "readonly",

            # Clicking or focusing selects the complete value
            onclick = "this.select();",
            onfocus = "this.select();"
          )
        ),

        column(
          width = 6,

          tags$label(
            `for` = "selected_rt_text",
            style = "font-weight: bold;",
            "RT:"
          ),

          tags$input(
            id = "selected_rt_text",
            type = "text",
            class = "form-control",
            value = rt_text,
            readonly = "readonly",

            onclick = "this.select();",
            onfocus = "this.select();"
          )
        )
      ),

      additional_details_ui
    ),

    plotlyOutput(
      "feature_plot",
      height = "260px"
    )
  )
})
  
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
  output$dl_annotation <- downloadHandler(

  filename = function() {
    paste0(
      dataset_name(),
      "_feature_table.csv"
    )
  },

  content = function(file) {

    out <- annotation_export_table()

    if (
  isTRUE(input$use_main_gnps_pairs) &&
  !is.null(input$file_main_gnps_pairs)
) {

  peak <- raw_df()
  id_col <- input$component_match_col

  validate(
    need(
      length(id_col) == 1L &&
        nzchar(id_col) &&
        id_col %in% names(peak),
      "Select the peak-table feature ID column for GNPS pairs."
    ),
    need(
      nrow(out) == nrow(peak),
      "Annotation export rows do not match the peak table."
    )
  )

  ids <- tibble::tibble(
    export_row = seq_len(nrow(peak)),
    ClusterID = trimws(as.character(peak[[id_col]]))
  ) %>%
    dplyr::filter(
      !is.na(ClusterID),
      nzchar(ClusterID)
    )

  joined <- merge(
    ids,
    main_component_map(),
    by = "ClusterID",
    all = FALSE,
    sort = FALSE
  )

  collapse_ids <- function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]

    paste(
      stringr::str_sort(x, numeric = TRUE),
      collapse = ", "
    )
  }

  network <- joined %>%
    dplyr::group_by(export_row) %>%
    dplyr::summarise(
      GNPS_ClusterID = collapse_ids(ClusterID),
      GNPS_ComponentIndex = collapse_ids(ComponentIndex),
      .groups = "drop"
    )

  index <- match(seq_len(nrow(out)), network$export_row)

  out$GNPS_ClusterID <- network$GNPS_ClusterID[index]
  out$GNPS_ComponentIndex <-
    network$GNPS_ComponentIndex[index]
}
    
    data.table::fwrite(
      out,
      file,
      na = ""
    )
  }
)
  
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

    out <- volcano_to_wide_if_needed(
  add_gnps_network_annotations(rv$volcano)
)

    data.table::fwrite(out, file, na = "")
  }
)

 output$dl_filtered_feature_table <- downloadHandler(

  filename = function() {
    original_name <- input$file_data$name %||% "dataset.csv"

    paste0(
      tools::file_path_sans_ext(basename(original_name)),
      "_metabocano.csv"
    )
  },

  content = function(file) {
    req(procReady(), rv$fmap)

    original <- as.data.frame(
      raw_df(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    validate(
      need(
        nrow(original) == nrow(rv$fmap),
        "The feature table has changed. Run Process again before downloading."
      )
    )

    retained_features <- unique(
      as.character(filtered_volcano()$Feature)
    )

    # The feature map follows the original peak-table row order.
    keep <- as.character(rv$fmap$Feature) %in%
      retained_features

    if (identical(input$software_tool, "msdial")) {

  # Read the original file without treating any row as a header.
  source_file <- input$file_data$datapath

  field_counts <- utils::count.fields(
    source_file,
    sep = ",",
    quote = "\"",
    comment.char = ""
  )

  max_cols <- max(field_counts, na.rm = TRUE)

  validate(
    need(
      is.finite(max_cols) && max_cols > 0,
      "Could not determine the MS-DIAL table structure."
    )
  )

  full_table <- utils::read.csv(
    source_file,
    header = FALSE,
    col.names = paste0("V", seq_len(max_cols)),
    colClasses = "character",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character(0),
    comment.char = "",
    fill = TRUE
  )

  # Locate the same header used by read_msdial_robust().
  header_row <- NA_integer_

  for (i in seq_len(min(50L, nrow(full_table)))) {
    row_text <- tolower(
      as.character(unlist(full_table[i, ], use.names = FALSE))
    )

    row_text <- gsub("[^a-z0-9]", "", row_text)

    if (
      "averagemz" %in% row_text ||
      "alignmentid" %in% row_text
    ) {
      header_row <- i
      break
    }
  }

  validate(
    need(
      !is.na(header_row),
      "Could not locate the MS-DIAL feature header."
    )
  )

  validate(
    need(
      nrow(full_table) - header_row == length(keep),
      paste(
        "MS-DIAL rows do not match the processed feature map.",
        "Run Process again before downloading."
      )
    )
  )

  # Preserve ALL introductory rows and the original header.
  # Filter only the feature rows below that header.
  rows_to_export <- c(
    seq_len(header_row),
    header_row + which(keep)
  )

  out <- full_table[
    rows_to_export,
    ,
    drop = FALSE
  ]

  data.table::fwrite(
    out,
    file,
    col.names = FALSE,
    row.names = FALSE,
    na = "",
    quote = "auto"
  )

} else {

  out <- original[
    keep,
    ,
    drop = FALSE
  ]

  data.table::fwrite(
    out,
    file,
    na = ""
  )
}
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
  volcano = add_gnps_network_annotations(rv$volcano)
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

#.....................................................
shinyApp(ui, server)
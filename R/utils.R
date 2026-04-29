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
  rownames(mat) <- clean_sample_names(sample_cols)

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
                               adj  = c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"),
                               paired = FALSE,
                               eqvar  = FALSE,
                               pseudocount = 1.1,
                               log2_test = FALSE,
                               scale_data = FALSE,
                               ref_group = NULL) {
  test <- match.arg(test)
  adj  <- match.arg(adj)

  validate(need("Label" %in% names(df_used), "Internal error: Label column missing."))

  feats <- setdiff(colnames(df_used), "Label")
  validate(need(length(feats) > 0, "No feature columns detected."))

  gr <- as.factor(df_used$Label)
  lev <- levels(gr)
  validate(need(length(lev) >= 2, "Need at least 2 Label groups to run statistics."))

  if (!is.null(ref_group) && ref_group %in% lev) {
    others <- setdiff(lev, ref_group)
    # Each row: [Control, Treatment]
    comb <- cbind(rep(ref_group, length(others)), others)
  } else {
    # Fallback to all pairs if no ref selected
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

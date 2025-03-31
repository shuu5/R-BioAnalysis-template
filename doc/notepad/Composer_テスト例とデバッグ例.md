# R言語によるテストとデバッグの実践例

#テスト #デバッグ #ロギング #testthat #futile.logger #cli #covr #R #実践例 #生物情報学

> [!NOTE]
> このファイルは[[Composer_ProjectRulesTDD]]から分割された、テストとデバッグの具体例を示すものです。
> 関連ガイド: [[Composer_ProjectRulesTDD]] | [[Composer_ドキュメント例]] | [[Composer_ProjectRules_Cursor実装例]]

## 1. テスト例

### 1.1 基本的なテスト

テストは`testthat`パッケージを使用して記述します。以下は基本的な文字列操作関数のテスト例です。

```r
# 良いテストの例 (tests/testthat/test-string_utils.R)
library(testthat)

context("文字列操作関数")

test_that("capitalize_first関数は文字列の最初の文字を大文字にする", {
  expect_equal(capitalize_first("hello"), "Hello")
  expect_equal(capitalize_first("world"), "World")
  expect_equal(capitalize_first(""), "")
})

test_that("capitalize_first関数はNAを入力するとNAを返す", {
  expect_equal(capitalize_first(NA), NA)
})
```

### 1.2 テスト環境の設定

テスト環境の設定は以下のように行います。これにより、テスト実行時にパッケージ環境が正しく読み込まれます。

```r
# tests/testthat.R の例
# プロジェクト環境を確実に読み込む
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# テストライブラリを読み込む
library(testthat)
library(yourpackage)  # パッケージ名を実際のものに置き換えてください

# テスト実行
test_check("yourpackage")  # パッケージ名を実際のものに置き換えてください
```

### 1.3 モックを使ったテスト

外部依存関係をモック化することで、テストの独立性と再現性を確保します。

```r
# mockeryを使用した外部依存関係のモック例
test_that("get_gene_info関数は外部APIにリクエストを送信し結果を正しく処理する", {
  # モック関数の定義
  mock_response <- list(
    status_code = 200,
    content = '{"symbol": "BRCA1", "name": "breast cancer 1", "location": "17q21.31"}'
  )
  
  # httrパッケージのGET関数をモック化
  mockery::stub(get_gene_info, "httr::GET", function(...) mock_response)
  
  # テスト実行
  result <- get_gene_info("BRCA1")
  
  # 結果の検証
  expect_equal(result$symbol, "BRCA1")
  expect_equal(result$name, "breast cancer 1")
  expect_equal(result$location, "17q21.31")
})
```

### 1.4 例外テスト

エラー条件をテストすることも重要です。以下は例外処理のテスト例です。

```r
test_that("validate_input関数は不正な入力に対してエラーを投げる", {
  # NULLの場合
  expect_error(validate_input(NULL), "Input cannot be NULL")
  
  # 数値以外の場合
  expect_error(validate_input("string"), "Input must be numeric")
  
  # 負の値の場合
  expect_error(validate_input(-10), "Input must be positive")
})
```

### 1.5 カバレッジ測定

テストカバレッジは`covr`パッケージで測定できます。これにより、テストがコードの全ての部分をカバーしているか確認できます。

```r
# カバレッジ測定スクリプト例 (coverage.R)
library(covr)

# カバレッジ計算
package_coverage <- covr::package_coverage()
print(package_coverage)

# カバレッジレポート生成
covr::report(package_coverage)

# 最小カバレッジのチェック
min_coverage <- 80
if (covr::percent_coverage(package_coverage) < min_coverage) {
  stop(paste("Test coverage is below", min_coverage, "%"))
}
```

## 2. デバッグとロギング例

### 2.1 futile.loggerの基本設定

ログ記録は`futile.logger`パッケージを使用します。以下は基本的な設定例です。

```r
# ロギングの設定
library(futile.logger)
library(cli)

# 解析開始時の設定
log_file <- paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M"), ".log")
flog.appender(appender.file(log_file))
flog.threshold(INFO)  # 通常実行時
# flog.threshold(DEBUG)  # デバッグ時
```

### 2.2 ログレベルの使い分け

適切なログレベルを使い分けることで、必要な情報だけを表示・記録できます。

```r
# 各ログレベルの例
flog.trace("Very detailed debugging: Iteration %d of inner loop", i)  # 最も詳細
flog.debug("Variable values: x=%f, y=%f", x, y)  # デバッグ情報
flog.info("Processing sample %s", sample_id)  # 一般的な情報
flog.warn("Missing values detected in %d genes", n_missing)  # 警告
flog.error("Failed to load file: %s", file_path)  # エラー
flog.fatal("Critical failure in core component: %s", component)  # 致命的エラー
```

### 2.3 cliを使ったユーザーインターフェース

ユーザーへの表示は`cli`パッケージを使用します。これにより、見やすく整形された出力ができます。

```r
# 解析スクリプト内での使用例
cli_h1("scRNA-seq Quality Control")
flog.info("Starting QC with %d cells and %d genes", ncells, ngenes)

# 長時間処理の進捗表示
cell_total <- ncol(sce)
cli_progress_bar("Processing cells", total = cell_total)
for (i in seq_len(cell_total)) {
  # セルごとの処理
  
  # 過剰なロギングを避ける（10%ごとに記録）
  if (i %% max(1, round(cell_total / 10)) == 0 || i == 1 || i == cell_total) {
    flog.debug("Processed %d/%d cells (%.1f%%)", i, cell_total, 100 * i / cell_total)
  }
  
  cli_progress_update()
}

# 結果の表示
qc_pass <- sum(sce$qc_pass)
flog.info("QC completed: %d/%d cells passed filters", qc_pass, cell_total)
cli_alert_success("{qc_pass} of {cell_total} cells passed QC filters")
```

### 2.4 エラーハンドリング

エラー処理は`tryCatch`を使用して行います。例外が発生した場合もログを残すことが重要です。

```r
# エラーハンドリング例
tryCatch({
  # 処理内容
  result <- perform_analysis(data)
  flog.info("Analysis completed successfully")
  return(result)
}, error = function(e) {
  flog.error("Error during analysis: %s", e$message)
  cli_alert_danger("Analysis failed: {e$message}")
  return(NULL)
}, warning = function(w) {
  flog.warn("Warning during analysis: %s", w$message)
  cli_alert_warning("Potential issue: {w$message}")
})
```

### 2.5 大規模データの効率的なロギング

大規模なオブジェクトをログに直接記録するのは避け、要約統計量のみを記録します。

```r
# 大規模データのロギング例
log_data_summary <- function(data, name = "data") {
  if (is.null(data)) {
    flog.warn("%s is NULL", name)
    return()
  }
  
  # データフレームやマトリックスの場合
  if (is.data.frame(data) || is.matrix(data)) {
    flog.info("%s dimensions: %d rows x %d columns", name, nrow(data), ncol(data))
    
    # 数値列の要約統計量
    if (is.data.frame(data)) {
      numeric_cols <- sapply(data, is.numeric)
      if (any(numeric_cols)) {
        for (col in names(data)[numeric_cols]) {
          stats <- summary(data[[col]])
          flog.debug("%s - Column %s summary: Min=%.2f, Median=%.2f, Max=%.2f, NA=%d", 
                     name, col, stats["Min."], stats["Median"], stats["Max."], 
                     sum(is.na(data[[col]])))
        }
      }
    }
    
    # 欠損値情報
    if (anyNA(data)) {
      na_count <- sum(is.na(data))
      na_percent <- 100 * na_count / (nrow(data) * ncol(data))
      flog.info("%s - Missing values: %d (%.2f%%)", name, na_count, na_percent)
    }
  }
  
  # Bioconductorオブジェクトの場合
  else if ("SummarizedExperiment" %in% class(data)) {
    flog.info("%s - SummarizedExperiment: %d features x %d samples", 
              name, nrow(data), ncol(data))
    
    # アッセイ情報
    assay_names <- assayNames(data)
    flog.debug("%s - Available assays: %s", name, paste(assay_names, collapse=", "))
    
    # メタデータ情報
    if (length(colData(data)) > 0) {
      flog.debug("%s - Sample metadata variables: %s", 
                name, paste(colnames(colData(data)), collapse=", "))
    }
  }
  
  # その他のオブジェクト
  else {
    flog.info("%s - Object of class: %s", name, paste(class(data), collapse=", "))
    flog.debug("%s - Object size: %.2f MB", name, object.size(data) / 1024^2)
  }
}

# 使用例
log_data_summary(expression_matrix, "expression_matrix")
log_data_summary(sample_metadata, "sample_metadata")
log_data_summary(sce_object, "single_cell_experiment")
```

### 2.6 ロギングとデバッグを活用した対話的な問題解決

対話的環境での問題解決例です。これにより、エラーの原因を迅速に特定できます。

```r
# 対話的デバッグの例
debug_differential_expression <- function(counts, metadata, design_formula) {
  flog.info("Starting interactive debugging of differential expression analysis")
  
  # 入力データのチェック
  flog.debug("Checking input data")
  log_data_summary(counts, "counts")
  log_data_summary(metadata, "metadata")
  
  # サンプル名の一貫性チェック
  count_samples <- colnames(counts)
  metadata_samples <- rownames(metadata)
  
  if (!identical(sort(count_samples), sort(metadata_samples))) {
    flog.error("Sample mismatch between counts and metadata")
    cli_alert_danger("Sample names in count matrix don't match metadata!")
    
    # 原因の詳細調査
    missing_in_counts <- setdiff(metadata_samples, count_samples)
    missing_in_metadata <- setdiff(count_samples, metadata_samples)
    
    if (length(missing_in_counts) > 0) {
      flog.debug("Samples in metadata but missing in counts: %s", 
                paste(missing_in_counts, collapse=", "))
      cli_alert_warning("Samples in metadata but missing in counts: {length(missing_in_counts)}")
    }
    
    if (length(missing_in_metadata) > 0) {
      flog.debug("Samples in counts but missing in metadata: %s", 
                paste(missing_in_metadata, collapse=", "))
      cli_alert_warning("Samples in counts but missing in metadata: {length(missing_in_metadata)}")
    }
    
    # 解決方法の提案
    cli_alert_info("Possible solutions:")
    cli_ul(c(
      "Check for typos in sample names",
      "Subset counts to only include samples in metadata",
      "Ensure metadata rownames match count column names"
    ))
    
    # サンプル名の比較を表示
    if (length(count_samples) < 20 && length(metadata_samples) < 20) {
      comparison <- data.frame(
        CountSamples = c(count_samples, rep(NA, max(0, length(metadata_samples) - length(count_samples)))),
        MetadataSamples = c(metadata_samples, rep(NA, max(0, length(count_samples) - length(metadata_samples))))
      )
      print(comparison)
    }
    
    return(NULL)
  }
  
  # デザイン行列の作成
  flog.debug("Creating design matrix from formula: %s", design_formula)
  tryCatch({
    design <- model.matrix(as.formula(design_formula), data = metadata)
    flog.debug("Design matrix dimensions: %d x %d", nrow(design), ncol(design))
    
    # デザイン行列の確認
    if (qr(design)$rank < ncol(design)) {
      flog.warn("Design matrix is not full rank - this may cause problems in model fitting")
      cli_alert_warning("Design matrix is not full rank!")
      
      # 原因の調査
      correlations <- cor(design)
      diag(correlations) <- 0
      high_cor <- which(abs(correlations) > 0.9, arr.ind = TRUE)
      
      if (nrow(high_cor) > 0) {
        for (i in 1:nrow(high_cor)) {
          var1 <- colnames(design)[high_cor[i, 1]]
          var2 <- colnames(design)[high_cor[i, 2]]
          cor_val <- correlations[high_cor[i, 1], high_cor[i, 2]]
          flog.debug("High correlation (%.3f) between variables: %s and %s", 
                    cor_val, var1, var2)
        }
        
        cli_alert_info("Detected high correlations between design variables.")
        cli_alert_info("Consider simplifying your design formula.")
      }
    }
    
  }, error = function(e) {
    flog.error("Error creating design matrix: %s", e$message)
    cli_alert_danger("Failed to create design matrix: {e$message}")
    
    # 一般的な問題の診断
    if (grepl("factor", e$message, ignore.case = TRUE)) {
      flog.debug("Design formula likely contains undefined factors")
      var_classes <- sapply(metadata, class)
      flog.debug("Metadata variable classes: %s", 
                paste(names(var_classes), var_classes, sep="=", collapse=", "))
      
      cli_alert_info("Check that all variables in your formula exist in metadata and have appropriate types")
    }
    
    return(NULL)
  })
  
  flog.info("Interactive debugging completed")
  cli_alert_success("Debugging completed - review log for details")
} 
# TCGA RNA-seqデータ解析パイプラインの開発例

#cursor #composer #プロンプト #開発指示 #R #tdd #git #実装例 #bioinformatics #TCGA

> [!NOTE]
> このファイルは[[Composer_プロンプトテンプレート]]のテンプレートを使用した具体的な開発例です。TCGAのRNA-seqデータを用いた解析パイプライン開発のプロセスを示しています。

## 概要

このドキュメントでは、TCGAからダウンロードしたRNA-seq STARカウントデータとクリニカルデータを用いて、以下の解析を行うパイプラインの開発プロセスを示します：
- PCA、ヒートマップによる探索的解析
- 特定遺伝子発現量でのグループ分け
- 差分発現遺伝子（DEG）解析
- パスウェイエンリッチメント解析
- 生存分析

開発は「最初にドキュメントで仕様をまとめて、それに従って開発し、適宜ドキュメントを修正する」という方針で、TDDアプローチとGitによるバージョン管理を適切に組み合わせて進めていきます。

このワークフローの作成は以下のプロンプトで行いました。
```
@{テンプレート}

具体例として、RNA-seq解析においてTCGAからダウンロードしたRNA-seq STAR カウントデータとクリニカルデータを基に、PCA, heatmapを作成し、特定の遺伝子の発現の大小でグループを分けてDEG, pathway enrichment analysis, survival analysisを行うパイプラインを作成するための一連のプロンプトを示して。
最初にドキュメントで仕様をまとめてそれに従って開発していき適宜ドキュメントを修正する方針としたい。
開発中に適宜バージョン管理もおこなうようにして。
```


## ステップ1: 仕様ドキュメントの作成


```
TCGA RNA-seqデータを用いた解析パイプラインの仕様ドキュメントを作成してください。
以下の解析ステップを含むパイプラインの詳細設計をR Markdownで記述してください：

1. 解析の背景と目的：
   - がん種別の遺伝子発現プロファイルの特徴付け
   - 特定の遺伝子（例：TP53）の発現量による患者層別化

2. データ取得と前処理：
   - TCGAデータ（STAR RNAカウント）のダウンロード方法
   - クリニカルデータの取得と結合
   - データの正規化とフィルタリング手順

3. 探索的データ解析：
   - PCA解析によるサンプル分布可視化
   - 遺伝子発現ヒートマップ作成
   - 特定遺伝子の発現分布と閾値設定

4. 差分発現解析：
   - DESeq2を用いた解析手法
   - 結果のフィルタリング基準（log2FC、p値）
   - 結果の可視化（Volcano plot、MA plot）

5. パスウェイ解析：
   - 使用するデータベースとパッケージ
   - エンリッチメント解析手法
   - 結果の可視化と解釈

6. 生存解析：
   - Kaplan-Meier曲線の作成手順
   - Cox比例ハザード解析の実施方法
   - 結果の解釈ガイドライン

7. 再現性と品質管理：
   - コードの構造と命名規則
   - ユニットテストの範囲
   - バージョン管理とドキュメント更新ポリシー

この仕様書は開発プロセスのガイドラインとして機能し、必要に応じて更新されることを想定してください。
TDDの原則に従い、各機能の実装前にテスト要件も明記してください。
```

## ステップ2: プロジェクト初期化とGit設定

```
TCGA_RNA_Analysis プロジェクトの初期化とGit設定を行うための手順を示してください。
以下の要素を含む初期設定スクリプトとコマンドを作成してください：

1. プロジェクトディレクトリ構造の構築：
   - R/（関数定義）
   - data/（入力データ保存用、gitignore対象）
   - results/（結果保存用、gitignore対象）
   - tests/（テストコード）
   - vignettes/（解析ドキュメント）
   - inst/extdata/（サンプルデータ）

2. 必要なパッケージのセットアップ：
   - DESeq2, TCGAbiolinks, survival, clusterProfiler, pheatmap, ggplot2など
   - renvを用いた環境管理の設定

3. Gitリポジトリの初期化：
   - .gitignoreファイルの設定（大きなデータファイル、中間結果などを除外）
   - READMEの作成
   - 初期コミット
   - develop, masterブランチの設定

4. CI/CD設定（オプション）：
   - GitHub Actionsの設定ファイル
   - テスト自動化スクリプト

すべてのコマンドを実行可能な形で提示し、各ステップの目的と意義を簡潔に説明してください。
```

## ステップ3: データ読み込みモジュールの実装

```
TDDルールに従って、TCGAデータの読み込みと前処理を行うモジュールを実装してください。
モジュール名は load_tcga_data とし、以下の機能を持つR関数群を開発してください：

1. TCGAからRNA-seqデータをダウンロードする関数：
   - 入力：がん種、バージョン、保存先ディレクトリ
   - TCGAbiolinksを使用
   - キャッシュ機能の実装（既存データの再利用）
   - エラーハンドリングと進捗報告

2. クリニカルデータを取得し加工する関数：
   - 入力：がん種、必要なフィールド一覧
   - 必要な前処理（日付形式の統一、カテゴリ変数の因子化など）
   - 欠損値の処理と報告

3. RNA-seqデータとクリニカルデータを結合する関数：
   - 適切なキーでのマージ処理
   - サンプルIDの一貫性確認
   - 結果のサマリー統計出力

4. データ正規化を行う関数：
   - DESeq2による正規化
   - TPM/FPKM変換オプション
   - 低発現遺伝子のフィルタリング

まず各関数のテストコードを testthat パッケージを使って記述し、次に実装を行ってください。
モックデータを適切に作成し、エッジケースを含むテストを行ってください。
関数はパイプライン全体での使いやすさを考慮した設計にしてください。
```

## ステップ4: 探索的データ解析モジュールの実装

```
TDDルールに従って、TCGA RNA-seqデータの探索的解析を行うモジュールを実装してください。
モジュール名は explore_expression_data とし、以下の機能を持つR関数群を開発してください：

1. PCA解析を実行する関数：
   - 入力：正規化済み発現量データ
   - prcomp関数を使用したPCA計算
   - ggplot2による可視化
   - サンプル属性によるポイントの色分け機能

2. 遺伝子発現ヒートマップを作成する関数：
   - 入力：正規化済み発現量データ、対象遺伝子リスト
   - pheatmapによる可視化
   - サンプルと遺伝子のクラスタリングオプション
   - カスタムアノテーションの追加機能

3. 特定遺伝子の発現分布を解析し、グループ分けを行う関数：
   - 入力：正規化済み発現量データ、対象遺伝子、閾値（任意）
   - 発現量分布のヒストグラム表示
   - 自動閾値決定アルゴリズム（オプション）
   - ハイ/ロー発現グループの割り当て

テストファーストの原則に従い、各関数のユニットテストを先に作成してください。
テストでは、様々な入力パターンとエッジケースを考慮してください。
関数は再利用性を高めるためmodular designを心がけ、適切なドキュメントを付けてください。
```

## ステップ5: 差分発現解析モジュールの実装

```
TDDルールに従って、遺伝子発現レベルによるグループ間の差分発現解析を行うモジュールを実装してください。
モジュール名は differential_expression とし、以下の機能を持つR関数群を開発してください：

1. DESeq2を用いた差分発現解析を実行する関数：
   - 入力：カウントデータ、グループ情報（因子型）、共変量（オプション）
   - DESeqDataSetオブジェクトの作成と前処理
   - 分散推定と差分発現検定
   - 結果のフィルタリングと整形

2. 差分発現結果を可視化する関数群：
   - Volcanoプロット作成関数
   - MAプロット作成関数
   - 上位DEGのヒートマップ作成関数
   - 特定遺伝子セットのboxplot作成関数

3. 結果を要約し、レポートを生成する関数：
   - 有意な遺伝子数のサマリー
   - 上位遺伝子リストの出力
   - 統計情報のテーブル作成
   - 結果ファイルのエクスポート機能

テストコードでは、モックデータセットを用意し、以下のケースをテストしてください：
- 単純な二群比較
- 共変量を含む複雑なデザイン
- エッジケース（ゼロ発現遺伝子、少数サンプル）
- 期待される出力フォーマットとデータタイプ

関数は再現性と柔軟性を重視し、パラメータのカスタマイズを可能にしてください。
すべての関数にはroxygen2形式の詳細なドキュメントを付けてください。
```

## ステップ6: パスウェイ解析モジュールの実装

```
TDDルールに従って、差分発現解析の結果を用いたパスウェイ解析モジュールを実装してください。
モジュール名は pathway_analysis とし、以下の機能を持つR関数群を開発してください：

1. 遺伝子セットエンリッチメント解析（GSEA）を実行する関数：
   - 入力：DESeq2の結果オブジェクト、ランキング方法
   - clusterProfilerを用いたGSEA実行
   - 複数のデータベース（KEGG, GO, Reactome等）のサポート
   - 結果のフィルタリングと整形

2. オーバーラップベースのエンリッチメント解析を実行する関数：
   - 入力：有意差のある遺伝子リスト、閾値パラメータ
   - エンリッチメント計算と統計的検定
   - 多重検定補正の複数オプション

3. パスウェイ解析結果を可視化する関数群：
   - エンリッチメントバープロット作成関数
   - エンリッチメントネットワーク図作成関数
   - ドットプロット作成関数
   - エンリッチメントマップ作成関数

4. 統合的な解釈を支援する関数：
   - 複数条件間のパスウェイ比較
   - 冗長性の除去と結果のクラスタリング
   - 生物学的解釈のためのサマリー生成

テストコードでは、以下のケースをカバーしてください：
- 様々なサイズの遺伝子リスト
- 異なるデータベースとパラメータ設定
- 結果の一貫性と再現性
- エッジケース（空のリスト、全ての遺伝子など）

各関数はパイプライン内での使用と独立した使用の両方に対応できるように設計し、
徹底したエラーチェックとユーザーフレンドリーなメッセージを実装してください。
```

## ステップ7: 生存解析モジュールの実装

```
TDDルールに従って、遺伝子発現データとクリニカルデータを用いた生存解析モジュールを実装してください。
モジュール名は survival_analysis とし、以下の機能を持つR関数群を開発してください：

1. 生存データを準備する関数：
   - 入力：クリニカルデータ、必要なイベント・時間変数
   - 生存データの妥当性チェック（欠損値、異常値）
   - Survオブジェクトの作成
   - 適切な打ち切りの処理

2. Kaplan-Meier解析を実行する関数：
   - 入力：Survオブジェクト、グループ変数
   - survminerを用いたKM曲線の作成
   - log-rank検定の実行と結果報告
   - カスタマイズ可能な視覚化オプション

3. Cox比例ハザード解析を実行する関数：
   - 入力：Survオブジェクト、共変量データフレーム
   - 単変量・多変量Coxモデルの構築
   - 比例ハザード仮定の検証
   - 結果の整形と解釈の補助

4. 予後予測モデルを構築・評価する関数：
   - 複数の遺伝子発現を組み合わせたスコアリング
   - モデルのパフォーマンス評価（C-index、AUC等）
   - 交差検証による堅牢性評価
   - ノモグラムの作成（オプション）

テストコードでは以下の点を検証してください：
- 様々なサンプルサイズと欠損パターン
- 異なるグループ数とバランス
- エッジケース（全生存、全死亡など）
- 結果の統計的妥当性

関数は再現性と明瞭な解釈を重視し、臨床的に意味のある出力を提供するよう設計してください。
すべての統計的検定には適切な多重検定補正を適用し、解釈の誤りを防止する機能を含めてください。
```

## ステップ8: パイプラインの統合とバージョン管理

```
これまでに実装した各モジュールを統合し、完全なパイプラインとして機能させるための手順を計画してください。
また、バージョン0.1.0としてリリースするためのGit管理手順も含めてください：

1. パイプライン統合：
   - メインスクリプトの作成（各モジュールを順番に呼び出す）
   - 設定ファイルの仕様設計（YAML/JSON形式）
   - モジュール間のデータ受け渡し方法の標準化
   - エラーハンドリングとロギングの一貫した実装

2. 統合テスト：
   - エンドツーエンドテストの実装方法
   - 小規模データセットでのテスト実行手順
   - パフォーマンステストと最適化

3. ドキュメント更新：
   - README.mdの完成
   - 使用方法ガイドの作成
   - パラメータリファレンスの作成
   - 実行例の追加

4. バージョン管理：
   - featureブランチからdevelopへのマージ手順
   - コードレビュー基準
   - バージョン0.1.0のタグ付け方法
   - CHANGELOG.mdの記載内容

5. リリース準備：
   - 最終テスト実行のチェックリスト
   - マスターブランチへのマージ手順
   - リリースノートの作成ガイドライン

各手順には具体的なGitコマンドとその目的を含め、TDD原則とGitフローに従った開発プロセスを維持してください。
パイプラインが再現可能で堅牢な結果を提供することを確認するための検証ステップも明記してください。
```

## ステップ9: ドキュメントの仕上げと次期開発計画

```
TCGA RNA-seq解析パイプラインのv0.1.0リリースに向けて、最終ドキュメントの仕上げと今後の開発計画を作成してください：

1. 最終ドキュメント：
   - 完成したパイプラインのR Markdownによる解説
   - TCGAデータを使った実際の解析例
   - パラメータ設定のベストプラクティス
   - 結果の解釈ガイドライン
   - トラブルシューティングセクション

2. リリースノート：
   - v0.1.0の主要機能リスト
   - 既知の制限事項と対処法
   - 依存パッケージとバージョン要件
   - インストール手順と初期設定

3. 次期開発計画（v0.2.0）：
   - 拡張予定の機能（他のがん種対応、追加の解析方法等）
   - パフォーマンス改善計画
   - ユーザーインターフェース改善案
   - コミュニティフィードバックの収集方法

4. 最終品質チェック：
   - コードスタイルの一貫性確認
   - テストカバレッジレポート
   - ドキュメントの完全性チェック
   - 再現性検証の最終手順

5. リリース後の計画：
   - バグ修正対応プロセス
   - ユーザーサポート体制
   - コミュニティ貢献ガイドライン
   - 次期リリースのタイムライン

ドキュメントはGitHubのマークダウン形式で準備し、コード例と出力例を十分に含めてください。
特に再現性と堅牢性に焦点を当て、パイプラインの限界と適切な使用方法を明確に伝えるようにしてください。
``` 
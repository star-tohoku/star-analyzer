# 計画・実装・比較記録の索引

整理日：2026-09-13。既存9文書を内容別に移動し、ファイル名・日付・解析結果は保持した。

## フォルダ構成

```text
mdfiles/
├── README.md                  この索引
├── kfparticle/
│   ├── plans/                 KFParticle導入・再構成の計画（3件）
│   ├── implementation/        実装と動作検証の記録（1件）
│   └── comparisons/           Λの手法・cut・候補数の比較（4件）
├── femto/
│   ├── plans/                 原子核別Femto解析の構成・移植計画
│   └── implementation/        実装・実イベント検証と実行手順
└── configuration/             設定の整理・保守記録（1件）
```

最新のΛ比較を確認する場合は、[2026-09-13：Imp5との候補数整合cut scan](kfparticle/comparisons/lambda_helix_kf_imp5_cut_scan_20260913.md) から読む。

## Femto解析の実装計画 — `femto/plans/`

| 文書 | 内容・位置付け |
|---|---|
| [plan_femto_lambda_nuclei_kfparticle_20260929.md](femto/plans/plan_femto_lambda_nuclei_kfparticle_20260929.md) | Λ–d/t/³He/⁴HeをStFemtoMakerへ移植し、保存済み高purity KFParticle選択を導入する実装前計画。旧histogram、全質量k*二次元分布、SE/ME、既存解析の保護と検証を定義 |

## Femto解析の実装・検証 — `femto/implementation/`

| 文書 | 内容・位置付け |
|---|---|
| [femto_lambda_nuclei_kfparticle_20260930.md](femto/implementation/femto_lambda_nuclei_kfparticle_20260930.md) | 4核種をStFemtoMaker＋保存済み高purity KFへ移植。比較基準は元のanaLambdaNuclearId。実装、10k検証、QA、設定と復元手順 |
| [femto_lambda_histograms_20260930.md](femto/implementation/femto_lambda_histograms_20260930.md) | Λ–原子核Femto解析のヒストグラム名・意味・軸・対象候補の日本語一覧。±3σの信号窓、全質量・サイドバンド、同一／混合イベントを区別 |
| [元解析のヒストグラム一覧（JSON）](femto/plans/femto_lambda_legacy_histogram_manifest_20260930.json) | 元解析のhistogramキー・型・軸・fill段階の機械可読一覧 |

## KFParticleの導入計画 — `kfparticle/plans/`

| 文書 | 内容・位置付け |
|---|---|
| [plan_kfparticle_pro_integration.md](kfparticle/plans/plan_kfparticle_pro_integration.md) | STAR pro / SL24yへの導入調査、参照パッケージ・他ユーザーの所在。初期の方針記録 |
| [plan_kfparticle_lambda.md](kfparticle/plans/plan_kfparticle_lambda.md) | 既存Helix解析を保護して、新しいKFParticle Λ Makerを追加する初期計画 |
| [plan_kfparticle_lambda_full_reconstruction.md](kfparticle/plans/plan_kfparticle_lambda_full_reconstruction.md) | Finder / TopoReconstructorを用いるfull再構成への拡張計画、各クラスの役割 |

計画書は当時の設計判断を残した履歴であり、すべてが現在の実装状態を表すわけではない。実装状態は次の検証記録と現行コード・設定を参照する。

## 実装・検証 — `kfparticle/implementation/`

| 文書 | 内容 |
|---|---|
| [kfparticle_full_implementation_20260907.md](kfparticle/implementation/kfparticle_full_implementation_20260907.md) | full再構成、ROOT5、mode別vertex条件、人工/実データ試験、既存Helixの非退行確認 |

## Λの比較・cut検討 — `kfparticle/comparisons/`

| 文書 | 内容・位置付け |
|---|---|
| [lambda_helix_kf_comparison_10000_20260907.md](kfparticle/comparisons/lambda_helix_kf_comparison_10000_20260907.md) | 標準KFと旧Helixの最初の10,000イベント質量比較 |
| [lambda_helix_kf_imp5_comparison_10000_20260908.md](kfparticle/comparisons/lambda_helix_kf_imp5_comparison_10000_20260908.md) | KFの共通選別を旧Imp5へ合わせた10,000イベント比較 |
| [kfparticle_lambda_qhu1_cut_comparison_20260908.md](kfparticle/comparisons/kfparticle_lambda_qhu1_cut_comparison_20260908.md) | qhu1コード・analysis note Fig.9とのcut比較と未確定事項。追加検討は一旦中断中 |
| [lambda_helix_kf_imp5_cut_scan_20260913.md](kfparticle/comparisons/lambda_helix_kf_imp5_cut_scan_20260913.md) | Λ近傍の候補数をそろえる距離/pointing cut scanとS/B比較。Notion Imp5表との照合も記録 |

## 設定の整理 — `configuration/`

| 文書 | 内容 |
|---|---|
| [unused_lambda_cut_cleanup_20260908.md](configuration/unused_lambda_cut_cleanup_20260908.md) | 未使用のPID/track/V0等の設定を除き、実効cutと既存解析の結果が変わらないことを確認した記録 |

## 関連する保存場所

- [analysisnote/](../analysisnote/)：日別の作業ログと解析レベルの詳細ノート。
- [KFParticle ROOT出力索引](../rootfile/auau13p5_anaLambda_KFParticle/mdfiles/README.md)：各ROOTファイルの生成条件と過去の移動対応。
- [2026-09-13 cut scan出力索引](../rootfile/auau13p5_anaLambda_KFParticle/mdfiles/cut_scan_20260913.md)：今回の候補数整合QAのROOT・CSV。
- [docs/REFERENCE.md](../docs/REFERENCE.md)：解析の実行・設定・QA手順。

## 今後の追加・過去のリンクについて

- 新しい文書は該当フォルダへ置き、この索引にも追加する。別の解析テーマが増えたら、`kfparticle/` と同じ階層にそのテーマのフォルダを作る。
- 日付付きの記録は残し、計画・実装状態・比較結果を混ぜない。日々の作業ログは従来どおり `analysisnote/YYYYMMDD/` に置く。
- 今回の9文書はすべて、旧 `mdfiles/<ファイル名>` から上表の位置へ移動した。ファイル名は変更していないため、この索引で旧ファイル名を検索できる。
- 直下には旧位置のコピーやsymlinkを残していない。過去のチャット・Notion・保存archive内の旧リンクは、この索引の新しいリンクで読み替える。
- リポジトリ内の現行Markdown参照を更新し、移動文書の相対リンクを補正した。以前のROOT整理で切れていた一部リンクも、実在する移動先へ修正した。
- 実行ログ・再現用archive・ROOT内部metadata・解析コード・YAMLは変更していない。解析の数値・結論も変更していない。

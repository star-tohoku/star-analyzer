# KFParticleを用いたΛ–原子核Femto解析の実装・検証記録

作成日：2026-09-30。実装の基準：[2026-09-29の計画](../plans/plan_femto_lambda_nuclei_kfparticle_20260929.md)。

実装・ローカル検証を完了した。Λ–d/t/³He/⁴Heは `StFemtoMaker` と保存済み高purity KF選別を使って動作する。比較基準は元の `anaLambdaNuclearId.C` であり、未デバッグだった旧原子核別マクロではない。両エネルギー各10,000入力イベントで、単独KF解析の候補・原子核QA・SE/ME出力との照合が成功した。低統計のCFやペアに含まれるΛのpurityの物理解釈は、別途検証が必要である。

ヒストグラム名と意味は、[Femto Lambdaのヒストグラム一覧](femto_lambda_histograms_20260930.md) にまとめた。

## 1. 対象範囲と比較基準

4つの入口 `anaFemtoLambda_{d,t,3He,4He}.C` は、full KFParticleによるΛ候補供給部を渡して、実際に **StFemtoMaker** を生成する構成へ変更した。未完成だった `StFemtoMakerLambdaNuclear` は新経路で使用しない。

**収量と原子核選別の比較基準は、変更していない `analysis/anaLambdaNuclearId.C` とした。** 旧 `anaFemtoLambda_4He.C` はデバッグ未完了のため、物理的な基準には使用していない。元のHelix解析、Lambda／NuclearId／Mix Maker、保存済み高purity KFカット、旧設定は残している。

今回はカットの再最適化、centrality較正の変更、ROOTの移行、farmへのjob投入、commit／pushは行っていない。

## 2. 実装した内容

| 構成要素 | 実装内容 |
|---|---|
| 共通Maker | `StMaker/StFemtoMaker/StFemtoMaker.{h,cxx}` にΛ専用経路を追加。既存Phi／Kaon経路は維持 |
| 任意に渡すインターフェース | `include/FemtoLambdaProvider.h`。候補と診断値だけを受け渡し、共通側にKF／SIMDのヘッダーを持ち込まない |
| full KF候補供給部 | `StMaker/kfparticle/FemtoLambdaKfProvider.{h,cxx}`。全PicoDst → Finder／Topo → 保存済み最終カット |
| 共通の最終選別 | `StMaker/kfparticle/KfLambdaSelector.h`。単独KF MakerとFemtoが同じ判定関数を使用 |
| 原子核選別・旧形式出力 | `StMaker/common/FemtoLambdaLegacy.{h,cxx}`。元のNuclearId選別とヒストを埋める段階を維持 |
| 入口・runner | 4組のマクロ／runnerと共通の `analysis/FemtoLambdaRunSupport.h`。mainconf、入力、出力、job ID、イベント上限を引数で指定 |
| 実行スクリプト | 4組の直接実行／コンテナ実行用スクリプト、`run_femtoLambda_common.sh`、`femtoLambda_environment.sh` |
| QA | `common/macro/checkHistAnaFemtoLambda.C`、runner、`script/singularity_checkHistAnaFemtoLambda.sh` |
| 設定 | エネルギー×原子核の8 mainconf、4種類の核種・channel定義、明示的なイベント／原子核／mixing設定、8つのヒストYAML |

`libStFemtoMaker.so` と `libStCommon.so` にKFParticleの共有ライブラリ依存を必須化していない。具体的なKF候補供給部は `libStKfParticleCommon.so` に置き、Λ用runnerだけが読み込む。Phi／Kaon用runnerには不要である。

### 維持した選別・物理量の定義

- 新mainconfは保存済み `kf_auau{13p5,3p85}_anaLambda_KFParticle_highpurity.yaml` を直接参照する。高purityカットの値は変更していない。
- イベント選別は `lambda_imp5_compat` と明記した。汎用Phi解析のイベント／trackカットをΛ娘粒子や原子核に追加適用しない。
- 両データセットはFXTである。13.5 GeVはビームエネルギー（2020年、P24iy、run 21033026）、3.85 GeVもビームエネルギー（2019年、P24iy、run 20160024）で、colliderの重心系エネルギーではない。
- 原子核PIDは、元のTPC選別、最適核種仮説の選択、rigidity上限を維持した。`nHitsDedx >= 15`、`gPt >= 0.1`、正のdE/dx、`maxPOverQ = 2.5`、`nSigma = 2` を使用する。原子核TOF選別は無効のままである。
- 不変質量軸には、質量制約をかけていないKF Λ質量を保存する。ペアのエネルギー計算には固定Λ質量1.115683 GeV/c²を使用し、Heの運動量にはZ=2補正を一度だけ適用する。
- 同一イベントのtrack共有判定には、確認済みのPico配列indexとイベント識別情報を使用する。別イベントの同じlocal indexは共有trackとして除外しない。
- Mixingは既存StFemtoMakerのpoolを使う。最新20イベント、両方向を独立処理、centrality 9 bin、Vz=198–202 cmを10 binとする。範囲外Vzは端のbinへ明示的に割り当てて数を記録し、黙ってイベントを除外しない。
- 元の0.5 T／1.4 mというmerging QA条件を残し、track IDと配列indexの混同を修正した。符号付き電荷・Heのrigidity・実イベント磁場を使うQAも別名で追加した。新しいclose-pair vetoは追加していない。

旧名の崩壊長・横方向崩壊長QAは、KF崩壊点とPVの幾何学距離を使う。別名の `hLambdaKF_DecayLength` はPV制約を使ったKF推定値である。`hDCA12`、`hDCAV0` もKF推定値となるため、Helixで得た値との数値一致を前提にはしない。

### Λのsignal選択：±3σ

全4核種・両エネルギーで、次の設定を共通に使う。

```yaml
lambdaSignalMean: 1.11596
lambdaSignalSigma: 0.00296
lambdaSignalNSigma: 3.0
lambdaSidebandOuterFactor: 4.0
```

**signalの条件は `abs(M(pπ) - 1.11596) <= 3 × 0.00296`、すなわち1.10708–1.12484 GeV/c²（両端を含む）である。** 今回の10,000イベント検証も既にこの±3σを使用しており、確認によってカット値や収量が変わったわけではない。

ここでのσはYAMLに設定した質量ピーク幅であり、各KF候補の `massError` ではない。今回、新しいKF質量分布を再fitしてσを求め直す処理は行っていない。

| 領域 | 質量範囲 [GeV/c²] | 用途 |
|---|---|---|
| 保存する全質量範囲 | 1.05–1.25 | inclusive質量QA、全質量のk*二次元分布 |
| signal（±3σ） | [1.10708, 1.12484] | signalのSE／ME、k*、q_lab、CF |
| 左sideband | [1.08044, 1.10708) | 低質量側の背景確認 |
| 右sideband | (1.12484, 1.15148] | 高質量側の背景確認 |
| 元のパラメータ | mean = 1.11596、sigma = 0.00296、nSigma = 3、外側係数 = 4 | sideband外端は中心から12σ |

`FemtoLambdaLegacy::MassRegion()` がこの条件でsignal／sidebandを分類する。各channelの上下限が同じ定義と一致することを初期化時に検査する。**全質量分布やsideband候補を保存する前に、全候補を±3σで削除することはしない。** これは元解析の `k* × invariant mass` と背景評価を残すためである。

この窓は以前のinclusive purity studyの [1.110, 1.122) とは異なる。また、比較図の表示・正規化範囲 [1.10, 1.14) も選別条件ではない。過去のinclusive purityを、そのまま原子核ペアや個別k*／centrality binのpurityとは扱わない。

### 旧解析と新解析で収量が一致しなくてよい理由

元の `anaLambdaNuclearId` mainconfはHelix Imp5カットを参照する。今回はΛ選別を**採用済みの高purity KF設定**へ置き換えた比較であり、旧解析と収量を一致させる再調整は行っていない。

| パラメータ | 旧Helix（両エネルギー） | KF 13p5 | KF 3p85 |
|---|---:|---:|---:|
| 陽子／πのnSigma | 3 / 3 | 3 / 3 | 3 / 3 |
| 陽子DCAの下限 [cm] | 0.7 | 0.7 | 1.0 |
| πのDCA下限 [cm] | 1.0 | 1.0 | 2.5 |
| nHitsFitの下限 | 15 | 20 | 25 |
| nHitsRatioの下限 | 0.52 | 0.52 | 0.52 |
| 娘粒子間距離の上限 [cm] | 0.5（Helix） | 0.5（KF） | 0.5（KF） |
| 親粒子とPVの距離上限 [cm] | 0.5（Helix） | 0.4（KF） | 0.5（KF） |
| cos(pointing)の下限 | 0.998 | 0.998 | 0.998 |
| Helix path lengthの上限 [cm] | 100 | 100 | 100 |
| KF崩壊長の下限 [cm] | 適用なし | 10 | 無効 |
| KF崩壊長有意度の下限 | 適用なし | 12 | 10 |

これ以外にFinder／interface／topologyのKF条件がある。全設定値は参照先のプリセットYAMLと出力metadataに保存する。この表はHelixとKFの距離推定量が同じという意味ではない。一致を要求したのは、原子核選別と同じKF設定を使う単独解析のΛ候補である。

## 3. ヒストグラム・出力の取り決め

[旧ヒストグラム定義の一覧（JSON）](../plans/femto_lambda_legacy_histogram_manifest_20260930.json) に、キー、directory、型、軸、埋める段階、任意だが未設定のキー、旧形式の空ヒストを記録した。日本語の説明は[ヒストグラム一覧](femto_lambda_histograms_20260930.md)を参照。

核種ごとのROOTファイルには次を保存する。

- トップ階層：Λ／共通QA、および共通形式の全質量・signal・左右sidebandのSE／MEヒスト。
- `true/`：元の共通NuclearId QAと、対象核種の同一イベントペア。
- `mix/`：対象核種の混合イベントペア。
- `true/hKstarMass_<S>_CentBin0..8` と `mix/hKstarMass_Mixed_<S>_CentBin0..8`：200 × 200 bin、k*=0–1 GeV/c、質量=1.05–1.25 GeV/c²。
- 全質量ペアはsignal／sideband分類の前に一度だけ埋める。MEの両方向を同じヒスト群へ反映する。
- ヒストとは別に `eventLedger`、`selectedLambda`、KF実効設定、mainconfの識別情報、イベント／候補数、完了／未完了状態を保存する。

既存出力pathへの上書きは拒否する。入力は、要求イベント上限より後のファイルも含めて事前検査し、covariance／PID branch欠落や入力件数不足では明示的に失敗させる。ROOT5のACLiCにはプロセス固有の作業directoryを使い、共有キャッシュの衝突を避ける。

## 4. 検証結果

根拠となるログ：`rootfile/femto_lambda_kf_validation_20260930/provenance/`。

| 確認項目 | 結果・根拠 |
|---|---|
| SL24y／ROOT 5.34/38／GCC 4.8.5による全ビルド | 成功。`build_retry.log`、`build_tests.log` |
| 共通最終選別・候補供給部 | 成功。4096通りの判定比較、境界・無効値・不正値・path・由来情報・固定質量の確認 |
| 入口・設定・QAの保護処理 | 6テストが成功。全8 mainconfで `lambdaSignalNSigma == 3.0`、正のσ、各channel境界の整合を確認する条件を追加し、再実行も成功 |
| 旧ヒストの定義・互換性 | 632項目の確認に成功 |
| 入力事前検査 | 正常fixtureと6種類の失敗ケースを確認 |
| Maker初期化・所有権・不正設定 | 25ケース成功。`maker_negative.log` |
| 4核種の旧互換helper合成試験 | 成功。`toy_legacy_*.log` |
| 元NuclearIdとの照合（13p5） | 10,000入力イベント、25,060原子核候補で一致。`nuclear_closure_13p5_10000.log` |
| 元NuclearIdとの照合（3p85） | 10,000入力イベント、47,699原子核候補で一致。`nuclear_closure_3p85_10000.log` |
| KFを読み込まないPhiDeuteron／Kaon | ACLiC・runner読込に成功。`phi_core_load.log`、`kaon_core_load.log` |
| 既存KF全再構成／Pico adapterの回帰試験 | 両方ともstatus=0。`kf_regression.log` |
| 共通ライブラリのKF非依存 | 確認済み。`core_library_dependencies.txt` |
| 保護対象の旧コード・設定 | 289ファイルがバイト単位で同一。`preservation.json`、`preservation_final.json` |
| KF候補・イベント記録の照合 | 両エネルギー・全4核種で成功。単独KFの候補情報が完全一致 |
| 元解析と新解析の最終出力監査 | 両エネルギー×4核種で成功。`audit_all.log`、`audit.csv` |
| QA PDF | v2の8ファイル×4ページと、旧Helix対KFの比較図を確認 |

NuclearIdの直接比較では、変更していないStNuclearIdMakerを基準に、全4核種のtrack ID、核種、Z補正済み運動量と、共通QA 27ヒストの全bin内容・誤差・underflow／overflowを照合した。同じイベントを二度評価したことによる乱数差を除くため、**試験内だけ**centrality乱数状態を再現した。解析本体でのcentrality判定はイベントごとに一度のままである。

設定307ファイルの中に既存PhiDeuteron／Kaon用mainconfはなかった。そのため、コンパイル・読込・依存関係の分離は確認したが、これらの解析について**実イベント収量が変わらないことまでは確認していない。**

## 5. 実データによる比較

入力は `config/picoDstList/auau{13p5,3p85}_femtoLambda_validation_20260930.list` の順序を固定したローカルリストである。data22／data24のファイルへ直接アクセスできたため、この実行は接続エラーが出ていたXRootD :1095には依存していない。各核種、元Helix解析、単独KFの参照解析に、それぞれ同じ先頭10,000入力イベントを指定した。

保存先：`rootfile/femto_lambda_kf_validation_20260930/{13p5,3p85}/`。

| ファイル | 内容 |
|---|---|
| `d.root`、`t.root`、`3He.root`、`4He.root` | 原子核別StFemtoMakerと保存済み高purity KFカットによる新解析 |
| `legacy.root` | 変更していない元のanaLambdaNuclearIdと、元のHelix設定 |
| `standalone.root` | 新Femtoと同じ保存済みカット・イベント方針による単独KF Λ解析 |

旧Helixと高purity KFのペア数の違いは、再構成・選別の変更による結果であり、一致を要求する値でも検出効率そのものでもない。Λの移植検証では、同設定の単独KF解析と選択候補が一致することを確認した。

### ペア収量の比較（各10,000入力イベント）

signal列には、前述の**±3σ**に入るペアだけを数えている。

| データセット | 原子核 | 旧Helix：signal SEペア | 新KF：signal SEペア | 新KF：全質量SEペア | 新KF：全質量MEペア |
|---|---|---:|---:|---:|---:|
| 13p5 FXT | d | 5460 | 2117 | 2483 | 106703 |
| 13p5 FXT | t | 724 | 265 | 326 | 15148 |
| 13p5 FXT | 3He | 270 | 83 | 107 | 5736 |
| 13p5 FXT | 4He | 37 | 18 | 20 | 792 |
| 3p85 FXT | d | 2764 | 733 | 819 | 34204 |
| 3p85 FXT | t | 338 | 93 | 110 | 4910 |
| 3p85 FXT | 3He | 348 | 103 | 113 | 3705 |
| 3p85 FXT | 4He | 63 | 24 | 28 | 788 |

上記はunderflow／overflowも含む候補・ペア数であり、背景差引き信号数やefficiencyではない。原子核との組合せに依存しないΛ候補数とイベント数は次のとおり。

| データセット | 各実行の入力イベント数 | 新解析でKF再構成を行ったイベント数 | 旧Helix：全質量のentries | 選択されたKF Λ候補数 |
|---|---:|---:|---:|---:|
| 13p5 FXT | 10000 | 9744 | 7216 | 873 |
| 3p85 FXT | 10000 | 9505 | 1983 | 144 |

同じエネルギーの新4核種実行では、イベント記録と選択Λ候補が一致した。

`audit.csv` と `provenance/audit_all.log` で、両エネルギー・4核種のイベント記録、単独KFの候補情報、元解析の原子核QA 27ヒスト、ヒスト型・軸・誤差、全質量の射影、centrality別の同内容ヒスト、ME方向別counterを確認した。

単独KFには全イベント記録が出力されない。読込／再構成件数と、選択候補のrun・event・娘track index・質量・運動量は一致したが、**除外イベントの識別情報を1件ずつ照合することはできない。** 新4核種同士のイベント記録は10,000行すべて一致している。

初回の読取専用監査では、名前空間内STL型のROOT5 dictionary生成に失敗した。型名を明示修飾して試験マクロを修正し、解析ROOTには変更を加えていない。13p5の先行KF専用監査 `audit_kf_13p5.csv` に残る `legacy_and_pair_audit = not_run` は、その時点の記録である。最終の全出力監査では8実行すべて成功した。

旧13p5 Helixの初回実行は、7,000イベントの進捗表示後、`legacy.root` 作成前に終了コード143でsignal終了した。成功扱いにはせず、ログを保存した。独立した再実行は10,000入力イベントをexit code 0で完了した（`13p5_legacy_retry1.log`／`.exitcode`、実時間1,377.06秒、CPU時間1,366.42秒）。単独KF参照もexit code 0で完了した。成功済みROOTは上書きしていない。

Λ質量比較図：[13p5](../../../share/figure/femto_lambda_kf_validation_20260930/legacy_vs_kf_mass_13p5.pdf)、[3p85](../../../share/figure/femto_lambda_kf_validation_20260930/legacy_vs_kf_mass_3p85.pdf)。

左は生の候補数、右は [1.10, 1.14) GeV/c²内のbin和がそれぞれ1になるように正規化した分布である。表示窓内の和は、13p5が2438／758、3p85が628／126（旧Helix／KF）。両PDFを画像化して目視確認した。これは再構成と最終カットが異なる設定の比較であり、KFアルゴリズムだけの効果を切り出した図ではない。

[核種別QA PDFの一覧](../../../share/figure/femto_lambda_kf_validation_20260930/README.md)では、各4ページの `*_v2.pdf` を使用する。全32ページの読取・非白紙確認と、質量・ペア分布、低統計⁴HeのCFの目視確認を行った。元PDFも保存している。v2はCFの一部の点が表示上限で切れていた問題を修正したもので、CFの計算値とROOTヒストは変更していない。

## 6. 実行方法

新しい出力pathを指定する。例：

```bash
bash script/singularity_run_anaFemtoLambda_4He.sh \
  config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml \
  config/picoDstList/auau3p85_femtoLambda_validation_20260930.list \
  rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_highpurity/local10000.root \
  local10000 10000
```

必要に応じて `4He` を `d`、`t`、`3He` に、`3p85` を `13p5` に置き換える。旧原子核別Maker用mainconfは新入口と互換ではないため、黙って別解析を実行せず、理由を表示して停止する。

```bash
bash script/singularity_checkHistAnaFemtoLambda.sh \
  rootfile/femto_lambda_kf_validation_20260930/3p85/4He.root \
  config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml \
  share/figure/femto_lambda_kf_validation_20260930/new_3p85_4He.pdf
```

QAにはinclusive／ペア中のΛ質量、SE／MEの質量対k*、centralityを合算したsignal／sidebandの生CFを表示する。MEまたは正規化領域が空の場合は「未定義」と表示し、測定値ゼロとは扱わない。これは診断用であり、fitしたpurityや最終的な背景差引きCFではない。核種・k*・centrality別のpurity fitは今回未実施で、特にHeのsidebandは低統計である。

## 7. 既存解析の保護と復元

実装前に、未commit／未追跡のコード・設定も含めた復元用archiveをリポジトリ外へ保存した。

`../implementation_backups/femto_lambda_20260930/before.tar.gz`

SHA-256：`39e4ad24999400b61637920a8337c52a89213542f28e73e918fd0abab6e1be15`。

新しい作業状態へarchive全体を上書き展開しないこと。別の一時directoryへ展開して比較し、その後の変更を保存したうえで、必要な入口・runner・設定だけを戻す。Git HEADだけでは元の未追跡Makerを復元できない。

実装時のソース・設定・試験（ignore対象の新YAMLを含む）も `rootfile/femto_lambda_kf_validation_20260930/provenance/implementation_sources.tar.gz` に保存した。SHA-256は `cf1cb2b1435ab668d288a7240a8c04b391ed8e8b1d43b25f60ecff2943dd6348`。これはその時点のソース・設定snapshotであり、ROOT／PDF成果物や、その後の文書更新を含むarchiveではない。

元の `anaLambdaNuclearId.C` はHelix参照として引き続き実行できる。既存高purity KF YAMLや過去のROOT出力は上書きしていない。新実行スクリプトは `.current_mainconf` を切り替えず、この移植によって既定の解析設定は変更されない。

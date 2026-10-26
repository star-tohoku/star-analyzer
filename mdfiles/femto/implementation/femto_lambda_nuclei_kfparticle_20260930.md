# KFParticleを用いたΛ–原子核Femto解析

`anaFemtoLambda_{d,t,3He,4He}` は `StFemtoMaker` にKFParticleのΛ候補供給部を渡す。マクロはStChainの構築・実行を担当し、再構成・選別・ヒスト充填はコンパイル済みMakerとhelperが行う。

## 構成と設定

| 構成要素 | 役割 |
|---|---|
| `include/FemtoLambdaProvider.h` | KF/SIMD型を持ち込まず候補・診断値を渡すインターフェース |
| `StMaker/kfparticle/FemtoLambdaKfProvider` | 全PicoDstをFinder/Topoへ渡し、保存済みKF最終カットでΛ候補を供給 |
| `StMaker/kfparticle/KfLambdaSelector.h` | 単独KF MakerとFemtoが共有する最終判定 |
| `StMaker/common/FemtoLambdaLegacy` | 原子核PID、単粒子QA、旧形式の核・pair出力 |
| `StMaker/common/FemtoNuclearDedxCalibration` | 明示したYAMLのTPC中心・幅。schema 2はlog10(R)多項式 |
| `StMaker/common/FemtoNuclearPairPidSparse` | 任意の8次元pair-PID同時分布 |
| `analysis/FemtoLambdaRunSupport.h` | 入力事前検査、明示的mainconf、完了判定、出力保護 |
| `common/macro/checkHistAnaFemtoLambda.C` | 保存ROOTからのQAとraw CF |

主入口 `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml` は2021 P24iyを指定する。13p5や各核種の基本設定は独立したmainconfを持つ。解析に必要な値はmainconfが参照するYAMLから取得し、欠落時に別の設定へ黙って切り替えない。

`libStFemtoMaker.so` と `libStCommon.so` はKFParticleを必須リンクしない。具体的な候補供給部は `libStKfParticleCommon.so` にあり、Λ用runnerが読み込む。

## 選別と運動量

`lambda_imp5_compat` は有限PV/maxNTrとcentralityを用いる。汎用Phi用のtrackカットをΛ娘・核へ追加適用しない。`kf` は保存済みhighpurityプリセットを直接参照し、同じ設定の単独KF解析と最終選別を共有する。

標準4He設定では4Heだけ校正多項式を使い、d/t/3Heは `NuclearIdDeDxVsMom.h` の係数を使用する。全核種でR<2はTOF不要、R>=2は有効TOFと核種別の厳密なm²窓を要求する。4Heは0.61<R<5 GeV/cで、dと同じm²窓を使う。TPC best-species判定後に領域・TOFで落ちた核を他核種へ再割当てしない。数値とvetoの正本はnuclearid YAMLである。

Heは元trackの剛性にZ=2を一度だけ掛けた物理運動量を使う。pair計算は固定Λ質量 `lambdaPairMass`、質量軸は拘束していないKF質量を使う。Pico配列indexとイベント識別子で共有trackを除き、別イベントの同じindexを誤って除かない。

## 質量領域と混合

標準maker YAMLのΛ信号条件は `abs(M-1.11596)<=3*0.00296` GeV/c²。σは設定されたピーク幅であり、候補ごとのmassErrorではない。信号両端を含み、隣接sidebandは信号側の端を除く。全質量・sidebandのQAを残すため、全候補を信号窓だけで削除しない。

同じpoolのcent9/Vz/EP区分、深さ、両方向混合をYAMLで指定する。範囲外Vzの明示的なclampは数を記録する。CFは足し合わせたSE/MEから後処理で計算する。

`lambdaHe3ControlEnabled` は物理3Heと、その同じ候補に4He質量を仮定したpair分布を追加する。後者は質量とk*/mTを変えるだけで、別のtrack・mixing channelではない。任意のpair-PID QAは物理poolに付随し、QAだけのイベントをpoolへ加えない。

## 出力と検証

全キー・型・軸・充填条件は[ヒストグラム一覧](femto_lambda_histograms_20260930.md)にまとめる。旧形式の軸契約は[回帰fixture](../../../tests/fixtures/femto_lambda_legacy_histogram_manifest.json)に保存する。標準4He設定は通常1,786ヒストと8次元Sparseを1個持つ。

入力のcovariance/PID branch、要求イベント数、出力先を事前検査する。既存ROOTを上書きしない。`eventLedger`、`selectedLambda`、実効設定、入力数・再構成数、`FemtoLambdaRunStatus`も保存する。schedulerの終了コードだけで解析完了と判断しない。

ビルド・回帰検証・実行方法は[REFERENCE](../../../docs/REFERENCE.md)、構成とQAの詳細は[Maker README](../../../StMaker/StFemtoMaker/README.md)を参照する。

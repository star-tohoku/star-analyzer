# KFParticleを使うΛ–原子核Femto解析のヒストグラム一覧

作成日: 2026-09-30。対象は `anaFemtoLambda_{d,t,3He,4He}` の新しいKF解析出力です。13p5/3p85で名称と軸は共通です。実際のYAML定義と `StFemtoMaker` / `FemtoLambdaLegacy` の記録処理を確認してまとめました。

## 表記と保存場所

| 表記 | 意味 |
|---|---|
| `{S}` | そのファイルの対象核種：`d`、`t`、`3He`、`4He` |
| `{C}` | 同じ順に `lambda_deuteron`、`lambda_triton`、`lambda_he3`、`lambda_he4` |
| `{i}` | `0`〜`8`。たとえば `CentBin2` はcent9=2のイベント |
| `{a,b}` | aまたはbへ置き換えた各名称。`{,_signal,...}` の先頭の空欄はsuffixなし |
| SE / ME | 同じイベントの組合せ / 異なるイベントを混合した組合せ |

cent9は大きいほど中心衝突です。詳細は[centralityの定義](../../../StRoot/StRefMultCorr/README.md)を参照してください。

| ROOT内の場所 | 内容 | 1ファイルのヒスト数 |
|---|---|---:|
| 最上位 `/` | Λ・イベントQA、共通Femto形式のSE/ME、混合QA | 67 |
| `true/` | 共通原子核QA、対象核種の旧形式SE、track-merging QA | 71 |
| `mix/` | 対象核種の旧形式ME | 42 |
| 合計 | TTreeやメタデータは含めない | 180 |

原子核QAには全4核種のPID仮説を残しますが、pairヒストはファイルごとの対象核種だけです。以下の1次元ヒストのY軸は、断りがなければ件数です。

## Λの質量窓：signalは±3σ

4核種とも、`lambdaSignalMean: 1.11596`、`lambdaSignalSigma: 0.00296`、`lambdaSignalNSigma: 3.0` です。したがってsignalの条件は **|M−1.11596| ≤ 3×0.00296 GeV/c²** です。

ここでσはYAMLに設定された質量ピーク幅であり、KF候補ごとの `massError` ではありません。各入力のピーク幅を自動フィットして窓を更新する処理でもありません。

| 区分 | 質量範囲 [GeV/c²] | 対応する名称 |
|---|---|---|
| 全質量 | 1.05〜1.25 | 共通形式のsuffixなし、`hKstarMass*`、inclusive Λ QA |
| signal | [1.10708, 1.12484]（±3σ、両端を含む） | `_signal`、旧形式のSB suffixなし |
| 左sideband | [1.08044, 1.10708)（−12σ〜−3σ） | `_leftSB`、`_SBNeg` |
| 右sideband | (1.12484, 1.15148]（+3σ〜+12σ） | `_rightSB`、`_SBPos` |

**signalのFemto解析は±3σですが、背景確認用の全質量・sidebandヒストは別に保存します。** `hLambda_InvMass` が±3σ外にも分布を持つのはこのためです。保存質量は質量制約をかけていないKF再構成値です。k*の計算では、別途固定Λ質量1.115683 GeV/c²を使用します。

## 最上位：Λ候補のQA

すべてKFの最終選別を通ったΛ候補について、原子核と組み合わせる前に1候補1回記録します。signal窓に限定せず、全保存質量範囲を含みます。娘PIDも「選ばれたΛの娘」であり、イベント内の全trackの分布ではありません。

| ヒストグラム名 | 軸と意味 |
|---|---|
| `hLambda_InvMass`、`hLambda_InvMass_CentBin{i}` | X：Λ不変質量。全centrality / 指定cent9 |
| `hLambda_{Pt,Eta,Phi}` | X：ΛのpT、擬ラピディティη、方位角φ |
| `hDCA12` | X：KFで求めたproton–pion間距離 [cm] |
| `hDCAV0` | X：KFで求めたΛとprimary vertex（PV）の距離 [cm] |
| `hCosPointing` | X：PV→崩壊点の方向とΛ運動量方向のなす角のcos |
| `hNSigma{Proton,Pion}` | X：娘のTPC proton/pion nσ |
| `hLambda_InvMass_vs_{Pt,DecayLength,TransDecayLength,Y}` | X：pT、PV→崩壊点の幾何学的距離、その横方向距離、再構成質量を用いたラピディティ。Y：不変質量 |
| `hDCAV0_vs_InvMass`、`hCosPointing_vs_InvMass` | X：不変質量。Y：PV距離 / pointing角のcos（上の行と質量軸の向きが異なる） |
| `hLambda_InvMass_vs_{Cent9,RefMultCorr}` | X：cent9 / 補正multiplicity。Y：不変質量 |
| `hLambdaKF_DecayLength` | X：PVを生成頂点として制約したKFコピーの崩壊長 [cm] |
| `hLambdaKF_DecayLengthSignificance` | X：上記崩壊長÷その誤差 |
| `hLambdaKF_Chi2Ndf` | X：質量制約なしKF候補のχ²/NDF |
| `hLambdaKF_TopoChi2Ndf` | X：PV制約を加えたKFコピーのχ²/NDF |

旧名の距離ヒストも、今回はKFから計算する量です。旧Helix値との数値一致を意味しません。`DecayLength` の幾何学的距離と `hLambdaKF_DecayLength` のPV制約付き推定量は区別してください。質量軸は200 bin、1.05〜1.25 GeV/c²です。

## 最上位：イベント・centralityのQA

| ヒストグラム名 | 軸と意味 | 記録段階 |
|---|---|---|
| `hVz`、`hRefMult` | X：PVのz [cm] / PicoEventのrefMult | bad-run等の除去前 |
| `hRefMultVsNTOFMatch` | X：nBTOFMatch、Y：raw multiplicity | 除去前 |
| `hRawMult` | X：centralityに使うraw multiplicity（今回のFXTではfxtMult） | bad-run・event条件通過後、pileup除去前 |
| `hCentralityRaw`、`hRefMultCorr`、`hRefMultWeight` | X：cent9、補正multiplicity、centrality weight | pileup等通過・centrality計算後、許容centrality選択前 |
| `hCentralityVsVz` | X：PVのz、Y：cent9 | 同上 |
| `hRefMultVsNTOFMatchAfter` | X：nBTOFMatch、Y：raw multiplicity | 同上 |
| `hCentrality`、`hCentrality16` | X：選択後cent9 / cent16 | 許容centrality通過後。`useWeight` が有効ならweight付き |
| `hRawMult_vs_RefMultCorr` | X：補正multiplicity、Y：raw multiplicity | KF処理成功後 |
| `h{RawMult,RefMultCorr,NTracks,TofMatchMult}_vs_Cent9` | X：cent9、Y：raw/補正multiplicity、入力track数、nBTOFMatch | 同上 |
| `hN{Proton,Pion}Cand_vs_Cent9` | X：cent9、Y：KFへ投入したproton/pionのPID仮説数（絶対PDGで集計） | 同上。最終Λ娘数ではない |
| `hNLambdaPairs_vs_Cent9` | X：cent9、Y：Finder/TopoのΛ粒子数（絶対PDG=3122） | 同上。adapter検査・最終cutの前の数であり、採用Λ数ではない |
| `hN` | 1 binのイベントカウンタ | KF・核種・pair処理が正常終了したイベント。Λが0候補でも加算 |

`hRawMult` とKF処理後のmultiplicity対cent9 QAは `fillCentralityQA` が有効なときに記録します。`hN` は全入力数ではありません。全入力数は後述の `inputEvents` を参照してください。物理候補・pairヒストへcentrality weightを自動適用しているわけではありません。

## true/：共通の原子核QA

KF処理まで通過したイベントのtrackについて記録します。基本条件はnHitsDedx ≥ 15、gPt ≥ 0.1、正のdE/dxです。以下のTPC仮説QAは、最終best-species選択・p/q上限・任意TOF選択の**前**を含み、最終pair用核種の分布とは異なります。

| ヒストグラム名（すべて `true/`） | 軸と意味・記録対象 |
|---|---|
| `hDedxP` | X：TPCのraw p/q相当運動量、Y：dE/dx。基本条件通過track |
| `hDedxP_{e,pi,K,p}` | 同じ軸。各通常粒子のTPC nσ条件を通るtrack（仮説間の重複あり） |
| `hDedxP_else` | 同じ軸。e/π/K/pのいずれからも離れたtrack |
| `hDedxP_{d,t,3He,4He}` | 同じ軸。各原子核のTPC仮説を通るtrack。全4仮説を保存 |
| `hPvs{Eta,Y,Pt}_{d,t,3He,4He}` | X：仮説ごとの物理運動量p、Y：η / y / pT。HeではpとpTをZ=2補正。TPC仮説通過時 |
| `hM2P` | X：raw p/q、Y：TOF βからraw運動量で計算したm²。利用可能なTOFを持つtrack |
| `hDedxP_{d,t,3He}_m2` | X：raw p/q、Y：dE/dx。各TOF m²窓に入るtrack。最終核種採用とは独立したQA |
| `hDedxP_cut` | X：raw p/q、Y：dE/dx。TOF未match/無効、または広いm²・運動量のQA条件を通らないtrack |

`hDedxP_cut` は「最終採用核種」や「TOFなしだけ」の分布ではありません。また核種の `m2_selection: false` でもTOF QAは記録されます。`hM2P` のHeのm²は物理質量の二乗へZ²補正した値ではありません。

## SE/MEのpairヒスト：旧形式と共通Femto形式

同じtrackをΛ娘と原子核に重ねるSE pairは除きます。MEは異なるイベントのみで、current Λ×過去の核種と過去のΛ×current核種の両方向を含みます。k*はpair静止系での片粒子運動量、q_labは実験室系の3運動量差の大きさです。

| 保存場所 | ヒストグラム名 | 軸・対象 |
|---|---|---|
| `true/` | `hKstar_{S}`、`hKstar_{S}_CentBin{i}` | X：k*。signal（±3σ）のSE、全cent / cent別 |
| `mix/` | `hKstar_Mixed_{S}`、`hKstar_Mixed_{S}_CentBin{i}` | 同上、ME |
| `true/` | `hKstar_{S}_{SBPos,SBNeg}`、`hKstar_{S}_{SBPos,SBNeg}_CentBin{i}` | X：k*。右/左sidebandのSE |
| `mix/` | `hKstar_Mixed_{S}_{SBPos,SBNeg}`、`hKstar_Mixed_{S}_{SBPos,SBNeg}_CentBin{i}` | 同上、ME |
| `true/` | `hQlab_{S}{,_SBPos,_SBNeg}` | X：q_lab。signal / 右SB / 左SBのSE（cent別版はない） |
| `mix/` | `hQlab_Mixed_{S}{,_SBPos,_SBNeg}` | 同上、ME |
| `true/` | `hKstarMass_{S}_CentBin{i}` | X：k*、Y：Λ質量。**signal/SB分類前の全質量SE** |
| `mix/` | `hKstarMass_Mixed_{S}_CentBin{i}` | 同上、全質量ME |
| 最上位 | `hKstar{SE,ME}_{C}{,_signal,_leftSB,_rightSB}` | X：k*。suffixなし＝全質量、他は各質量窓 |
| 最上位 | `hKstar{SE,ME}VsCent_{C}{,_signal,_leftSB,_rightSB}` | X：k*、Y：cent9。同じ分類をcent別に保存 |

k*軸は200 bin・0〜1 GeV/c、q_lab軸は200 bin・0〜2 GeV/c、全質量TH2は200×200 binです。軸外はROOTのunderflow/overflowに入ります。

旧形式と共通形式は**同じ採用pairを別名で保存したもの**です。たとえば `true/hKstar_d` と最上位の `hKstarSE_lambda_deuteron_signal` は同じ母集団なので、足し合わせてはいけません。全質量TH2も1 pairにつき1回で、signal/SBごとに再加算しません。これらはSE/ME分布であり、Maker出力に完成済みCFやpurityヒストがあるわけではありません。

## true/：track-mergingの確認

| ヒストグラム名 | 軸・意味 |
|---|---|
| `hDphiDeta_proton_{S}` | X：Δφ*、Y：Δη。signal SEのΛのproton娘と核種。旧固定磁場0.5 T・評価半径1.4 mの規約（ID/index参照は修正済み） |
| `hDphiDeta_proton_{S}_EventField` | 同じ組合せ。イベント磁場、符号付き電荷、Heの物理pTを整合させた追加QA |

各軸200 bin、−0.25〜0.25です。これらを保存するだけで、close-pair vetoは追加していません。

## 最上位：混合処理の確認

`hMixSamplerQA` はXが下表の状態番号0〜14、Yがchannel番号です。今回の設定ではbase channel（Y=0）だけを混合し、同じpairを各質量窓へ振り分けるため、他のY binを個別に埋める必要はありません。bin内容はイベント数ではなく、各状態のpair数またはbuffer方向数の累積です。

| Xの値 | 意味 |
|---|---|
| 0 / 1 / 2 / 3 | 試行pair数 / 混合可能pair数 / 記録pair数 / 未記録pair数 |
| 4 / 5 | sampling上限による未処理数 / pair条件による除外数 |
| 6 / 7 / 8 | 利用可能なbuffer方向数 / 空のbuffer方向数 / 誤って選ばれた空方向数 |
| 9 / 10 | forward / reverseで記録したpair数 |
| 11 / 12 | forward / reverseの混合可能pair数 |
| 13 / 14 | track共有等による除外数 / 質量窓による除外数 |

今回の `bufferAll` と全質量base channelでは、sampling上限やsignal窓でbase pairを落としません。X=4、8、14は通常0です。forwardはcurrent Λ×buffer核種、reverseはbuffer Λ×current核種です。

## ヒストグラムではない保存項目

| ROOT内の名称 | 種類・意味 |
|---|---|
| `eventLedger` | TTree。入力番号、run/event、処理段階、cent9、採用Λ数 |
| `selectedLambda` | TTree。採用Λの元イベント、娘のPico配列index、質量、運動量 |
| `KFParticleEffectiveConfiguration` | 実際に使用したKF設定を保存した文字列 |
| `FemtoLambdaMainconf`、`FemtoLambdaRunStatus` | mainconfの識別情報 / `completed`・`incomplete` |
| `inputEvents`、`reconstructedEvents`、`selectedLambdaCandidates` | 全入力数、KF再構成成功イベント数、採用Λ総数 |
| `mixingVzClampedEvents` | 混合Vz範囲外で端のbinへ割り当てたイベント数 |

現在のYAMLには `hDedxP_4He_m2`、原子核QAの `*_CentBin{i}`、`hVz_{S}`、`hMult_{S}` 等の旧コード上だけのoptionalヒストはありません。未作成のものを一覧の実出力に数えていません。`mix/` のヒストを空の複製として `true/` へ書く旧来の構成も引き継いでいません。

参照：[実装・検証まとめ](femto_lambda_nuclei_kfparticle_20260930.md)、[全キー・軸の機械可読一覧](../plans/femto_lambda_legacy_histogram_manifest_20260930.json)、[代表のΛ/common YAML](../../../config/hist/hist_femtoLambda_deuteron_kf.yaml)、[代表の核種/pair YAML](../../../config/hist/hist_femtoLambda_deuteron_nuclear.yaml)、[記録処理](../../../StMaker/common/FemtoLambdaLegacy.cxx)。


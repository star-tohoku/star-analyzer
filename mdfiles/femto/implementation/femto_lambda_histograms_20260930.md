# KFParticleを使うΛ–原子核Femto解析のヒストグラム一覧

対象は `anaFemtoLambda_{d,t,3He,4He}` のKF解析出力。以下は3.85 GeV・2021 P24iyの標準4He設定の名称・軸・充填条件をまとめる。


## 現行anaFemtoLambda_4He_KFParticleの全出力

対象入口：[main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml](../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml)。実際に選別する核は⁴Heと³Heで、³He→⁴He質量仮定は同じ³He候補・pairから派生する比較用表現。d/tのTPC仮説QAも保存するが、この設定にΛ–d/tのCFチャネルはない。

| ROOT保存場所 | 通常TH1/TH2 | 内容 |
|---|---:|---|
| `/`（トップ階層） | 145 | イベント・Λ/KF QA、Λ/娘/核の単粒子アクセプタンス、3表現の共通形式SE/ME、mixing QA |
| `true/` | 840 | 原子核QAと3表現のSE・mT・merging（300）、selected/prePIDペアQA（540） |
| `mix/` | 801 | 3表現のME・mT（261）、selected/prePIDペアQA（540） |
| 合計 | **1,786** | TTree・TNamed・TParameterは含めない |

さらにトップ階層に8次元THnSparseD **`hNuclearPidPairJoint`を1個**保存し、ヒストグラム総数は**1,787個**。この数は本節の標準4He設定に対応し、他のprofileとは異なる。

**全ヒストグラムの実名・型・軸・簡単な説明は、本ページ末尾の「現行設定の全ヒストグラムキー一覧」に収録した。** 1,786キーとSparseを3つのヒストYAMLから列挙する。

### 全体の読み方

| 分類 | 対象・主な用途 | 詳細説明 |
|---|---|---|
| イベント・centrality QA | event選択、multiplicity、cent9/16、正常終了数 | 本文「最上位：イベント・centralityのQA」 |
| Λ/KF QA | 全保存質量範囲の採用Λ・娘PID・距離・pointing・KF品質 | 本文「最上位：Λ候補のQA」 |
| 単粒子アクセプタンス | Λ信号、娘、⁴He、³He、³He→⁴He質量仮定のpT–Y、追加β分布 | 本節の下表と本文pT–Y節 |
| 共通原子核QA | 全4核種・通常粒子のdE/dx仮説、TOF m²、運動学 | 本文「true/：共通の原子核QA」 |
| ⁠⁴He専用PID/TOF QA | PID前の対照試料、TPC best後/TOF後、段階別数 | 本文高剛性TOF節と下記現行条件 |
| Λ–核pair | ⁴He・³He・³HeAs4HeのSE/ME、信号/左右SB、k*、q_lab、全Λ質量、cent9、mT | 本文pair・k*–mT節 |
| 近接pair QA | Λ娘陽子と各核表現のΔφ*–Δη | 本文track-merging節 |
| 混合QA | 8チャネルのmixing状態別カウント。base indexは⁴He=0、³He=4 | 本文混合処理節 |
| ペアPID QA | selected/prePID、k*対nσ/Rの1,080 TH2Dと同時8次元Sparse | 本文「原子核PIDとk*の相関QA」 |

### 現在のPID・TOF条件

- d/t/³HeのdE/dx中心・幅は`NuclearIdDeDxVsMom.h`の関数・係数・有効域を使用。⁴Heだけ確定したlog10(R)の中心3次・幅2次係数へ更新し、補間校正は使わない。
- 全4核種でR<2 GeV/cはTOF不要、R≥2は有効TOFと核種別m²窓を要求。⁴Heはdと同じ`3.268773 < (m/q)² < 3.693147`。tは`7.2564275～8.2816925`、³Heは`1.78285～2.06485`、いずれも両端除外。
- ⁠⁴Heは厳密にR>0.61、R<5 GeV/c、旧⁴He vetoなし。d/t/³Heの既存最小R・vetoは保持（特にR≥2.5のveto）。TOF要件だけで全域が採用されるわけではない。
- TPC best-species判定を先に行い、後段のR/TOF条件で落ちたwinnerを他核種へ再割当てしない。³He物理質量版と⁴He質量仮定版は同じ採用候補を使う。
- `hHe4PidStage`の0～4は最終R下限/上限等の前を含み、5は最終候補。`0=1+2`だが、現在のR制限では`5=1+4`を一般には仮定しない。


### 現在のβ・アクセプタンス追加分

既存のpT–YはΛ信号・陽子娘・π娘、⁴He、³He、³HeAs4Heについて保持。核は最終選別済み候補、Λと娘はΛ信号候補につき1回で、pair数による重みなし。

| TH2D名 | X軸 | Y軸 | 対象 |
|---|---|---|---|
| `hLambda_BetaVsEta_signal` | η_lab | β=物理P/E | Λ信号候補 |
| `hLambda_BetaVsYLab_signal` | Y_lab | β=物理P/E | 同上 |
| `hNucleus_BetaVsEta_he4` | η_lab | β=物理P/E | 最終⁴He候補、全R |
| `hNucleus_BetaVsYLab_he4` | Y_lab | β=物理P/E | 同上 |

βはTOF測定値ではない。⁴HeのPとpTはZ=2補正済み物理運動量。η/Yは実験室系、CMシフトなし。η/Y軸200 bin・−5～5、β軸200 bin・0～1。

さらに`hNucleus_{BetaVsEta,BetaVsYLab,PtVsYLab}_he4`の各基底名へ次の末尾を付けた**12個のTH2D**を保存する。対応する全体分布と同じ軸である。

| 末尾 | 条件 |
|---|---|
| `_lowR_noTof` | R<2、TOF不要。TOFがある候補も含む |
| `_lowR_tofMatch` | R<2、有効TOFあり、m²窓なし |
| `_lowR_tofM2` | R<2、有効TOF＋既存⁴He m²窓 |
| `_highR_tofM2` | R≥2、有効TOF＋既存⁴He m²窓 |

低R3区分は`tofM2 ⊆ tofMatch ⊆ noTof`で重複するため足し合わせない。各変数で`lowR_noTof + highR_tofM2 = 全体分布`がflowを含めて成立する。境界判定は元trackのdouble R・m²を使用。Λのβ分布は全体signal候補のまま。ROOTの`FemtoHe4AcceptanceSamples`へ区分定義・R境界・m²窓を記録する。


### 現在の定義・実装

- [KF・イベント・アクセプタンス・共通pair YAML](../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_control_acceptanceSamples20261022_kf.yaml)
- [原子核QA・旧形式pair YAML](../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_control_nuclear.yaml)
- [ペアPID QAとSparse YAML](../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_pid_pair_qa.yaml)
- [単粒子・原子核・pair充填とROOT保存](../../../StMaker/common/FemtoLambdaLegacy.cxx)、[event/mixingとledger保存](../../../StMaker/StFemtoMaker/StFemtoMaker.cxx)

`FemtoNuclearDedxCalibration`には消費した⁴He係数YAML全文と出典、`FemtoNuclearPidDomain`には核種別R・veto・TOF条件、`FemtoHe4AcceptanceSamples`には追加アクセプタンス試料の定義が入る。これらはTNamedで、上のヒストグラム数には含めない。その他の非ヒストグラム出力は本文末尾の保存項目とPID QA節を参照。

以下に分類別の充填条件と全キーを示す。


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

## 最上位：採用粒子のpT–Y_Lab分布

全4核種・両エネルギーに追加しました。4つともTH2Dで、**横軸は実験室系ラピディティ `Y_Lab`、縦軸はpT [GeV/c]**です。
`Y_Lab = 0.5 × ln((E+pz)/(E−pz))` とし、擬ラピディティηではありません。重心系への移動やビームラピディティの差引きは行いません。
Y_Lab軸は200 bin・−5〜5、pT軸は200 bin・0〜10 GeV/cです。軸の定義はヒストYAMLにあり、範囲外はunderflow／overflowに保存します。

| ヒストグラム名 | 記録対象・回数 | 運動量・質量の定義 |
|---|---|---|
| `hLambda_PtVsYLab_signal` | 高purity KF選別と±3σ質量窓を通ったΛ。原子核との組合せ前に1候補1回 | KFのΛ運動量と、ペア計算にも使う固定Λ質量 `lambdaPairMass` |
| `hLambdaProton_PtVsYLab_signal` | 上記Λを構成するp。Λ候補ごとに1回 | Pico配列indexで対応づけた娘trackの `gMom()` と `lambdaDaughterProtonMass` |
| `hLambdaPion_PtVsYLab_signal` | 上記Λを構成するπ⁻。Λ候補ごとに1回 | 同じく娘trackの `gMom()` と `lambdaDaughterPionMass` |
| `hNucleus_PtVsYLab_{deuteron,triton,he3,he4}` | そのファイルの対象核種で、best-species・p/q上限・設定されたTOF条件まで通過した候補。1track1回 | 原子核ペアに渡す四元運動量。HeはZ=2補正を一度だけ適用し、質量は核種YAMLの値 |

娘粒子にはsecondary trackのglobal運動量を使い、PVに結び付けたprimary運動量 `pMom()` やKF崩壊点へ輸送した娘運動量は使いません。
娘質量は各maker YAMLに `lambdaDaughterProtonMass: 0.9382720813`、`lambdaDaughterPionMass: 0.13957039` [GeV/c²] と明記しています。これらは同梱KFParticleの粒子質量定義と同じ値で、QA計算にだけ使用します。既存のKF再構成やk*計算の質量は変更しません。

原子核の分布は、同じイベントにΛがない場合も記録します。後で混合イベントにも使用する、最終選別後の候補全体を示すためです。
Λ・原子核とも、**ペア数や混合回数による重み付けはしません**。ただし娘trackが複数のΛ候補で共有された場合、娘の分布には各Λ候補について1回ずつ入ります。

これらは「実データで選択された粒子が占める運動学的範囲」を示します。生成粒子数で割ったアクセプタンス・再構成効率ではなく、検出器の効率そのものではありません。
また、Λの分布は±3σだけが対象なので、全質量の `hLambda_Pt` 等と母集団が異なります。

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


### 従来の核QA

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

## 最上位：信号窓を低質量側・高質量側に分けたk*

従来の±3σ信号、全質量、左右サイドバンドをすべて維持し、**同じ採用済みsignal pair**を次の2群へ追加記録します。
実際の判定は従来の `MassRegion(mass)==1` を先に確認し、質量中心より低いかどうかで分けます。

| suffix | Λ質量の条件 [GeV/c²] | 意味 |
|---|---|---|
| `_signalLow` | [1.10708, 1.11596) | −3σ〜0の低質量側 |
| `_signalHigh` | [1.11596, 1.12484] | 0〜+3σの高質量側。中心ちょうどはこちらだけに入る |

| ヒストグラム名（最上位） | 軸・対象 |
|---|---|
| `hKstarSE_{C}_signalLow` | X：k*。低質量側の同一イベントペア |
| `hKstarME_{C}_signalLow` | X：k*。低質量側の混合イベントペア |
| `hKstarSE_{C}_signalHigh` | X：k*。高質量側の同一イベントペア |
| `hKstarME_{C}_signalHigh` | X：k*。高質量側の混合イベントペア |
| `hKstar{SE,ME}VsCent_{C}_{signalLow,signalHigh}` | X：k*、Y：cent9。上の4分布をcentrality別に保存 |

`{C}` は `lambda_deuteron`、`lambda_triton`、`lambda_he3`、`lambda_he4` です。
k*軸は従来と同じ200 bin・0〜1 GeV/cです。質量判定にはbin中心ではなく候補の再構成質量を使います。保存質量はFloat_tであり、数値境界は従来と同じ浮動小数点比較に従います。

SE／MEの各bin（underflow／overflowを含む）で、**`signalLow + signalHigh = signal`**になります。
両群は重複しませんが、元の `signal` の部分集合なので、3つを足し合わせて収量にしないでください。
追加のpair loopやmixing channelは作らず、元の4 channelとmixing pool・両方向処理を維持します。新しいcutではありません。

各群のCFは、それぞれのSE／MEを独立に正規化して作ります。
`checkHistAnaFemtoLambda.C` は、設定されたcentrality範囲を合算・rebinした後、
従来どおり0.5〜1.0 GeV/c（YAMLの `normQMin/Max`）でMEをSEへ正規化し、低側／高側それぞれのCFを描画します。
CFを足し合わせて元のsignal CFにしてはいけません。合算する場合はSEとMEの段階で合算し、改めて正規化します。
正規化区間の統計がゼロならCFは未定義と表示し、MEがゼロのbinも測定値ゼロとして描きません。

アクセプタンス／質量窓二分割追加時点のQA PDFは、従来4ページに3ページを加えた7ページです。さらに本ページ末尾のk*–m_T追加により、最新出力では計12ページになります。
既存ROOTには新しい12ヒストは自動追加されません。**Makerを再ビルドして再解析する必要があります**。
旧ROOTのQAは警告を出して従来4ページを表示します。新ヒストが一部だけ存在する不完全な出力はエラーにします。

## true/：track-mergingの確認

| ヒストグラム名 | 軸・意味 |
|---|---|
| `hDphiDeta_proton_{S}` | X：Δφ*、Y：Δη。signal SEのΛのproton娘と核種。旧固定磁場0.5 T・評価半径1.4 mの規約（ID/index参照は修正済み） |
| `hDphiDeta_proton_{S}_EventField` | 同じ組合せ。イベント磁場、符号付き電荷、Heの物理pTを整合させた追加QA |

各軸200 bin、−0.25〜0.25です。これらを保存するだけで、close-pair vetoは追加していません。

## 最上位：混合処理の確認

`hMixSamplerQA` はXが下表の状態番号0〜14、Yがchannel番号です。基本設定ではbase channel（Y=0）を混合する。現行³He control付き設定では⁴HeのY=0と³HeのY=4を混合し、同じpairを各質量窓へ振り分けるため、他のY binを個別に埋める必要はありません。bin内容はイベント数ではなく、各状態のpair数またはbuffer方向数の累積です。

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

参照：[実装・検証まとめ](femto_lambda_nuclei_kfparticle_20260930.md)、[全キー・軸の機械可読一覧](../../../tests/fixtures/femto_lambda_legacy_histogram_manifest.json)、[代表のΛ/common YAML](../../../config/hist/hist_femtoLambda_deuteron_kf.yaml)、[代表の核種/pair YAML](../../../config/hist/hist_femtoLambda_deuteron_nuclear.yaml)、[記録処理](../../../StMaker/common/FemtoLambdaLegacy.cxx)。


## centrality別のk*–m_T

後段でcentralityとm_Tを指定してCFを作るため、**SE／MEを別々に、cent9=0〜8の各binにTH2Dを保存**します。
既存192ヒストを保持し、90ヒストを追加して計282ヒストです。全4核種・両エネルギーの設定に追加しています。

### m_Tの定義

ユーザー指定の平均質量・ペア横運動量の半分を使います（以下c=1）。

```text
m       = (m_Lambda + m_nucleus) / 2
pT_pair = |pT_vector_Lambda + pT_vector_nucleus|
pT      = pT_pair / 2
m_T     = sqrt(m^2 + pT^2)
```

これは **sqrt(M_pair^2 + pT_pair^2)というペア全体の横質量ではありません**。
横運動量はベクトル和であり、個々のpTの大きさの和でもありません。再構成Λの候補質量や、Λ／原子核各々の横質量の平均も使いません。

質量はk*と共通の固定仮説を使います。Λはmaker YAMLの `lambdaPairMass: 1.115683`、原子核はnuclearid YAMLの `nuclearMassDeuteron: 1.87561`、`nuclearMassTriton: 2.80892`、`nuclearMassHe3: 2.80839`、`nuclearMassHe4: 3.72742` [GeV/c²] です。
運動量は既存のペア候補に保存した物理運動量で、HeのZ=2補正を再適用しません。
signal、左右sideband、低側／高側、SE／MEすべてに同じ式と質量仮説を使用します。
出力ROOTの `FemtoLambdaPairMtDefinition`（TNamed）に式・使用質量・単位を記録します。

### 名前・軸・選択領域

`{S}` は `d,t,3He,4He`、`{i}` はcent9の0〜8です。

| 保存場所 | ヒストグラム名 | Λの質量領域 |
|---|---|---|
| `true/` | `hKstarMt_{S}_CentBin{i}` | signal：±3σ |
| `mix/` | `hKstarMt_Mixed_{S}_CentBin{i}` | signal：±3σ |
| `true/` | `hKstarMt_{S}_SBPos_CentBin{i}` | 右sideband |
| `mix/` | `hKstarMt_Mixed_{S}_SBPos_CentBin{i}` | 右sideband |
| `true/` | `hKstarMt_{S}_SBNeg_CentBin{i}` | 左sideband |
| `mix/` | `hKstarMt_Mixed_{S}_SBNeg_CentBin{i}` | 左sideband |
| `true/` | `hKstarMt_{S}_signalLow_CentBin{i}` | signalの低質量半分 |
| `mix/` | `hKstarMt_Mixed_{S}_signalLow_CentBin{i}` | signalの低質量半分 |
| `true/` | `hKstarMt_{S}_signalHigh_CentBin{i}` | signalの高質量半分 |
| `mix/` | `hKstarMt_Mixed_{S}_signalHigh_CentBin{i}` | signalの高質量半分 |

- X：k* [GeV/c]、既存と同じ200 bin・0〜1（bin幅0.005）。
- Y：m_T [GeV/c²]、200 bin・0〜10（bin幅0.05）。各nuclearヒストYAMLの `PairMt` 軸で変更可能。
- 合計：2（SE／ME）×5（質量領域）×9（cent9）＝90ヒスト／核種別ROOT。
- 各採用pairを、従来のk*と同じ選別・同じcent9で1回ずつ記録する。追加のcut、mixing pool、pair loopはない。
- MEの両方向も同じ保存処理を通り、同じcent9のイベント同士で構成する。
- m_Tの範囲外はunderflow／overflowへ残る。m_Tを全域（flow込み）で射影すれば元のcentrality別k*を再現する。
- 各k*–m_T binでも `signalLow + signalHigh = signal`。分割版とsignalを重ねて加算しない。
- 全質量のk*–massは従来どおり別に残す。今回のk*–m_Tには全質量版はなく、5つの質量領域を保存する。

### 後段でのm_T選択とCF

1. 希望のcentrality binのTH2を、SEとMEそれぞれで加算する。
2. 両者で**同じm_T範囲**を指定し、X軸（k*）へ射影する。
3. k*のrebinを揃え、そのcentrality・m_T範囲内のSE／MEで改めて正規化してCFを作る。
4. sideband CFも同じcentrality・m_T選択で作る。完成したCF同士を平均してcentralityを統合しない。

たとえばTH2ポインタを `hSE`／`hME` とし、m_Tのbin 31〜40（[1.50,2.00) GeV/c²）を選ぶ場合：

```cpp
TH1D* se = hSE->ProjectionX("se_mt150_200", 31, 40, "e");
TH1D* me = hME->ProjectionX("me_mt150_200", 31, 40, "e");
// 選択後のSE／MEを同じk* binningへ揃え、設定した正規化域でCFを作る。
```

これは射影例であり、全核種でこのm_T範囲に統計があるという意味ではありません。
**任意範囲の選択精度は保存bin幅0.05 GeV/c²まで**です。bin途中を境界にしても、その内部を正確に分割できません。
境界はbin edgeに合わせ、隣接領域を[下限,上限)として重複させないでください。もっと細かい選択には `PairMt` のbin数を変えて再解析が必要です。

CFの正規化は `C = (SE/ME) × (ME_norm/SE_norm)`。既存設定は0.5〜1.0 GeV/cで、必ず選択したm_T・centralityの分布から求めます。
MEがゼロのbinや正規化統計がゼロの範囲はCF未定義として扱います。今回Makerへ完成済みCFを保存する処理は追加していません。

### QA・旧ROOTとの互換性

QAの末尾に、signal／SBPos／SBNeg／signalLow／signalHighのk*–m_T（SEとME）を各1ページ追加し、最新ROOTは12ページです。
図は `cfCent9Min..cfCent9Max` を合算したものですが、ROOTには9 bin個別のTH2が残っています。
旧ROOTに90個の新ヒストが全くなければ警告して従来ページだけを表示します。一部欠損・異次元・軸不整合はPDF作成前に停止します。
既存ROOTにm_T情報は後付けできないため、再ビルド・再解析した新ROOTを使用してください。


## 原子核PIDとk*の相関QA

原子核のpurityをk*・centrality・Λ質量領域別に調べるための、**追加の入力分布**である。
Makerがpurity fitやCFの背景差引きを実行する機能ではない。
mainconfの任意キー `nuclearPidPairQaHist:` で専用ヒストYAMLを指定した場合だけ有効にし、
キーがない従来profileでは追加ヒストの生成・QA候補の保持を行わない。
専用定義は [原子核PIDペアQA YAML](../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_pid_pair_qa.yaml) に保存する。
標準3p85⁴He mainconfがこのQAを明示的に有効にする。
旧統合ROOTにはこれらの新分布はないため、再ビルド・再解析後のROOTが必要になる。


### 保存する2段階と運動量の定義

| 段階 | 対象 | 用途と制限 |
|---|---|---|
| `selected`（名称に追加suffixなし） | 元のTPC nσ窓・best-species判定・核種ごとの剛性範囲・TOF条件まで通った最終核候補 | 従来CFと同じ採用pairを追跡する。nσ窓で切断済みなので、これだけから窓外の背景やピーク形を復元しない |
| `prePID`（名称に `PrePid`） | 対象核種自身のTPC nσ窓とbest-species競合だけを外したQA仮説候補。track品質、剛性の許可／除外領域、適用されるTOF条件は維持 | 背景を含むnσ全形状のfit用入力。無条件の全track分布ではない。同一trackが³He仮説と⁴He仮説の両方に入ることがあり、仮説同士を独立粒子数として足さない |

同一核種・pair条件の `selected` は `prePID` に含まれる。
両stageを独立標本として尤度へ二重に加算したり、ヒストを足して統計を増やしたりしない。
³He実質量と `3HeAs4He` も同じtrack／pairの派生表現なので統計的に独立ではない。

⁴He表現は `nSigmaHe4`、³Heと `3HeAs4He` は**選別元の `nSigmaHe3`**を使う。
`3HeAs4He` は同じ³He trackの運動量を保って⁴He質量でk*を計算する比較用表現であり、
⁴Heのpullへ取り替えたり、別の核種／mixing channelとして追加したりしない。
剛性 `R=p/|q|` は元trackの物理運動量の大きさをHeの電荷数2で割った値で、
運動量のZ=2補正を二重に適用しない。追加pair QAでは保存済みFloat_t P4からRを求めるため、
元trackのdouble剛性の完全な再現ではなく、丸めを含む診断値である。
TOF量も未補正剛性から計算した `(m/q)^2` である。

QA候補は物理CF用の候補storeとは別に保持する。
MEでは**従来の物理候補によるpool参加条件・深さ・centrality/Vz区分・両方向を維持した同じpool**にQA情報を添付し、
QA候補だけのイベントを新たにpoolへ参加させない。
元イベントとtrack indexによる共有track vetoを維持する。
最終選択とQA選択の判定には元のdouble値を使い、保存済みFloat_tから境界判定をやり直さない。

### TH2Dの名前・軸・保存場所

次の表の `{S}` は `4He,3He,3HeAs4He`、`{i}` はcent9の `0..8`、
`{M}` は空欄（±3σsignal）、`_SBPos`、`_SBNeg`、`_signalLow`、`_signalHigh` のいずれかである。

| 保存場所 | ヒストグラム名 | X軸 | Y軸 |
|---|---|---|---|
| `true/` | `hKstarNSigma_{S}{M}_CentBin{i}` | k* | 選別元核種のdE/dx nσ、最終選択後 |
| `mix/` | `hKstarNSigma_Mixed_{S}{M}_CentBin{i}` | k* | 同上、ME |
| `true/` | `hKstarPOverQ_{S}{M}_CentBin{i}` | k* | 元trackのR、最終選択後 |
| `mix/` | `hKstarPOverQ_Mixed_{S}{M}_CentBin{i}` | k* | 同上、ME |
| `true/` | `hKstarNSigmaPrePid_{S}{M}_CentBin{i}` | k* | prePID仮説のdE/dx nσ |
| `mix/` | `hKstarNSigmaPrePid_Mixed_{S}{M}_CentBin{i}` | k* | 同上、ME |
| `true/` | `hKstarPOverQPrePid_{S}{M}_CentBin{i}` | k* | prePID仮説のR |
| `mix/` | `hKstarPOverQPrePid_Mixed_{S}{M}_CentBin{i}` | k* | 同上、ME |

全て明示YAML軸を使う。k*は200 bin・0〜1 GeV/c（5 MeV/c幅）、
Rは120 bin・0〜6 GeV/c（0.05幅）、nσは120 bin・−12〜12（0.2幅）。
軸範囲外の有限値はflowへ残し、通常のY軸可視範囲だけを射影するとその分は落ちる。
各段階で `2変数×3表現×2(SE/ME)×5質量区分×9cent9=540` 件、2段階で1080 TH2Dである。
signalとその低側・高側は部分集合として重複保存されるため、5区分を全て足して収量にしない。
最終選択後TH2のY全域・flow込みProjectionXは対応する従来のcentrality別k*を再現し、
各段階・各二次元binで `signalLow+signalHigh=signal` となる設計である。

2種類のTH2はそれぞれ周辺分布であり、両者だけから `(k*,R,nσ)` の同時分布は復元できない。
R依存の応答・purityをk*へ移すときにRとnσの独立性を暗黙に仮定しないよう、次の疎同時分布も残す。

### THnSparseDの同時分布と射影

最上位の `hNuclearPidPairJoint` は、同じQA pairを次の8軸で保存する1個の `THnSparseD` である。
疎形式なので全binを密に確保しないが、埋まったbinに応じてメモリ・ROOT容量・hadd時間が増える。

| 軸番号 | 軸名 | binning／カテゴリ |
|---:|---|---|
| 0 | `Kstar` | 200 bin、0〜1 GeV/c |
| 1 | `POverQ` | 120 bin、0〜6 GeV/c |
| 2 | `NSigma` | 120 bin、−12〜12 |
| 3 | `Cent9` | 9 bin、値0〜8 |
| 4 | `Representation` | 0=`he4`、1=`he3`、2=`he3as4` |
| 5 | `MixingMode` | 0=`SE`、1=`ME` |
| 6 | `MassClass` | 0=`SBNeg`、1=`signalLow`、2=`signalHigh`、3=`SBPos` |
| 7 | `PidStage` | 0=`selected`、1=`prePID` |

categorical軸もbin数・範囲をYAMLに明記し、上記意味との不一致は初期化時に拒否する。
signalを独立カテゴリとして二重に埋めず、signalLow／signalHighを**排他的に1回だけ**保存する。
signalの射影はMassClassの1と2を加算する。左右SBとは重ならない。
カテゴリを選び、他の物理軸のflowを含めて加算した射影は上のTH2と対応する。
各pairは単位weightで記録し、疎分布は `Sumw2` を有効にする。
非有限のk*・R・nσまたは不正カテゴリはQAだけで拒否し、通常CFのpair判定へ戻さない。
同時分布と二次元射影の閉包のため、新QAは有限な同一tupleを記録する。

ROOT5でflowを含む範囲を選んだ軸へ射影するときは、元の軸を維持する `"EA"` を使う。
`"E"` だけでは、出力先の物理軸がflowを含む幅へ拡張され、元TH2の軸と一致しなくなる。
たとえばcent9=0、⁴He、SE、±3σsignal、prePIDを選ぶ場合：

```cpp
for (int axis = 0; axis < 3; ++axis)
  joint->GetAxis(axis)->SetRange(0, joint->GetAxis(axis)->GetNbins() + 1);
joint->GetAxis(3)->SetRange(1, 1); // cent9=0
joint->GetAxis(4)->SetRange(1, 1); // he4
joint->GetAxis(5)->SetRange(1, 1); // SE
joint->GetAxis(6)->SetRange(2, 3); // signalLow + signalHigh
joint->GetAxis(7)->SetRange(2, 2); // prePID
TH2D* nsigmaVsKstar = joint->Projection(2, 0, "EA"); // X=k*, Y=nSigma
TH2D* rigidityVsKstar = joint->Projection(1, 0, "EA"); // X=k*, Y=R
// 使用後、次の選択へ影響させないため全8軸のSetRange(0, 0)を行う。
```

`A` は射影先軸のrange指定を無視して元の全軸を保つoptionである。
射影先のk*などをさらに制限したい場合は、このTH2を得た後で目的のbinを選ぶ。
独立toyでは1080個のカテゴリ／質量区分別射影が、別に組み立てたTH2の全flow・内容・誤差²と一致した。

このR幅0.05・nσ幅0.2は保存後の分割精度の上限であり、
以前の低R studyのnative R幅0.01や元dE/dxスペクトルを復元するものではない。
有限のflowは保持するが、flowの内部をさらにR・nσで分けられない。
解析時には選択範囲、flow方針、TPC winner条件、TOF条件、応答モデルを明示し、
³He coreのCFをそのまま⁴He窓へ漏れる³He tailのCFと同一視しない。

| 最上位のメタデータ名 | 内容 |
|---|---|
| `FemtoLambdaNuclearPairPidQaDefinition` | TH2の母集団・stage・元核種pull・剛性・質量表現の定義 |
| `NuclearPidPairJointMetadata` | 疎8軸のカテゴリ対応、exact QA YAML、flow／weight方針、helperが受け取ったtuple内の正常fill・非有限値拒否・不正カテゴリ拒否の件数。caller側共通finite gateより前の除外数は含まない |

TH2Dの内容配列は2段階で約203 MiB（1080×202×122×8 byte、flow込み）となり、
全TH2にSumw2配列も確保すると合計約406 MiBになる。疎分布とpool payloadの占有量は別に増える。
これは配列の概算であり、実RSSや本番ROOTサイズではない。新規設定では短時間smokeと高占有toyの容量・時間を確認し、
本番投入前に利用可能メモリとhadd負荷を評価する。offなら追加分を確保しない。


## 現行設定の全ヒストグラムキー一覧

現在のmainconfが参照する3つのYAMLを展開した**全1,786通常ヒストグラムと1 Sparse**の一覧。キーはROOT内のフルパス。各区分を開くと全名称を検索できる。型は`HistManager.cpp`の生成規則（無指定の1DはTH1F、2DはTH2F等）を適用し、実ROOTのクラス名でも照合した。軸は「軸ラベル / bin数 [下端, 上端]」。上端は軸範囲の端で、値の選別境界そのものではない。候補cutとflowの扱いは本文を参照。

同じpairを別名で保存する旧形式/共通形式、同じ信号の全体/上下半分、同じ低R候補のTOF比較群は互いに重複し、独立な収量として合計しない。ペアPID QAは完成したPurityやCFではなく後処理の入力分布である。

<details>
<summary>トップ階層（145個）：全キー・型・軸・説明</summary>

| ROOTキー | 型 | 軸 | 説明 |
|---|---|---|---|
| `hCentrality` | TH1F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5] | 選択後cent9 |
| `hCentrality16` | TH1F | X: cent16 / 16 bin [-0.5, 15.5] | 選択後cent16 |
| `hCentralityRaw` | TH1F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5] | centrality計算後・許容範囲選択前のcent9 |
| `hCentralityVsVz` | TH2F | X: V_{z} [cm] / 400 bin [-300, 300]<br>Y: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5] | centrality計算後のPV z–cent9 |
| `hCosPointing` | TH1F | X: cos(#theta) / 100 bin [0.99, 1.01] | PV→崩壊点とΛ運動量のpointing cos |
| `hCosPointing_vs_InvMass` | TH2F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25]<br>Y: cos(#theta) / 100 bin [0.99, 1.01] | Λ不変質量–pointing cos |
| `hDCA12` | TH1F | X: DCA [cm] / 200 bin [0, 5] | 採用Λの娘陽子–π間KF距離 |
| `hDCAV0` | TH1F | X: DCA [cm] / 200 bin [0, 5] | 採用Λ–PV間KF距離 |
| `hDCAV0_vs_InvMass` | TH2F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25]<br>Y: DCA [cm] / 200 bin [0, 5] | Λ不変質量–PV距離 |
| `hKstarMEVsCent_lambda_he3` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・ME・Λ全保存質量・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_as_he4` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_as_he4_leftSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_as_he4_rightSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_as_he4_signal` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_as_he4_signalHigh` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_as_he4_signalLow` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_leftSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・ME・Λ低側sideband・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_rightSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・ME・Λ高側sideband・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_signal` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・ME・Λ信号窓・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_signalHigh` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・ME・Λ信号高半分・最終核選別後 |
| `hKstarMEVsCent_lambda_he3_signalLow` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・ME・Λ信号低半分・最終核選別後 |
| `hKstarMEVsCent_lambda_he4` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・ME・Λ全保存質量・最終核選別後 |
| `hKstarMEVsCent_lambda_he4_leftSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・ME・Λ低側sideband・最終核選別後 |
| `hKstarMEVsCent_lambda_he4_rightSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・ME・Λ高側sideband・最終核選別後 |
| `hKstarMEVsCent_lambda_he4_signal` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・ME・Λ信号窓・最終核選別後 |
| `hKstarMEVsCent_lambda_he4_signalHigh` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・ME・Λ信号高半分・最終核選別後 |
| `hKstarMEVsCent_lambda_he4_signalLow` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・ME・Λ信号低半分・最終核選別後 |
| `hKstarME_lambda_he3` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・ME・Λ全保存質量・最終核選別後 |
| `hKstarME_lambda_he3_as_he4` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後 |
| `hKstarME_lambda_he3_as_he4_leftSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後 |
| `hKstarME_lambda_he3_as_he4_rightSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後 |
| `hKstarME_lambda_he3_as_he4_signal` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後 |
| `hKstarME_lambda_he3_as_he4_signalHigh` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後 |
| `hKstarME_lambda_he3_as_he4_signalLow` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後 |
| `hKstarME_lambda_he3_leftSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後 |
| `hKstarME_lambda_he3_rightSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後 |
| `hKstarME_lambda_he3_signal` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後 |
| `hKstarME_lambda_he3_signalHigh` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号高半分・最終核選別後 |
| `hKstarME_lambda_he3_signalLow` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号低半分・最終核選別後 |
| `hKstarME_lambda_he4` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ全保存質量・最終核選別後 |
| `hKstarME_lambda_he4_leftSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後 |
| `hKstarME_lambda_he4_rightSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後 |
| `hKstarME_lambda_he4_signal` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後 |
| `hKstarME_lambda_he4_signalHigh` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号高半分・最終核選別後 |
| `hKstarME_lambda_he4_signalLow` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号低半分・最終核選別後 |
| `hKstarSEVsCent_lambda_he3` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・SE・Λ全保存質量・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_as_he4` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_as_he4_leftSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_as_he4_rightSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_as_he4_signal` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_as_he4_signalHigh` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_as_he4_signalLow` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_leftSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・SE・Λ低側sideband・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_rightSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・SE・Λ高側sideband・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_signal` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・SE・Λ信号窓・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_signalHigh` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・SE・Λ信号高半分・最終核選別後 |
| `hKstarSEVsCent_lambda_he3_signalLow` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–³He・SE・Λ信号低半分・最終核選別後 |
| `hKstarSEVsCent_lambda_he4` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・SE・Λ全保存質量・最終核選別後 |
| `hKstarSEVsCent_lambda_he4_leftSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・SE・Λ低側sideband・最終核選別後 |
| `hKstarSEVsCent_lambda_he4_rightSB` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・SE・Λ高側sideband・最終核選別後 |
| `hKstarSEVsCent_lambda_he4_signal` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・SE・Λ信号窓・最終核選別後 |
| `hKstarSEVsCent_lambda_he4_signalHigh` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・SE・Λ信号高半分・最終核選別後 |
| `hKstarSEVsCent_lambda_he4_signalLow` | TH2F | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: cent9 / 9 bin [-0.5, 8.5] | k*–cent9・Λ–⁴He・SE・Λ信号低半分・最終核選別後 |
| `hKstarSE_lambda_he3` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・SE・Λ全保存質量・最終核選別後 |
| `hKstarSE_lambda_he3_as_he4` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後 |
| `hKstarSE_lambda_he3_as_he4_leftSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後 |
| `hKstarSE_lambda_he3_as_he4_rightSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後 |
| `hKstarSE_lambda_he3_as_he4_signal` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後 |
| `hKstarSE_lambda_he3_as_he4_signalHigh` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後 |
| `hKstarSE_lambda_he3_as_he4_signalLow` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後 |
| `hKstarSE_lambda_he3_leftSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後 |
| `hKstarSE_lambda_he3_rightSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後 |
| `hKstarSE_lambda_he3_signal` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後 |
| `hKstarSE_lambda_he3_signalHigh` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号高半分・最終核選別後 |
| `hKstarSE_lambda_he3_signalLow` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号低半分・最終核選別後 |
| `hKstarSE_lambda_he4` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ全保存質量・最終核選別後 |
| `hKstarSE_lambda_he4_leftSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後 |
| `hKstarSE_lambda_he4_rightSB` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後 |
| `hKstarSE_lambda_he4_signal` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後 |
| `hKstarSE_lambda_he4_signalHigh` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号高半分・最終核選別後 |
| `hKstarSE_lambda_he4_signalLow` | TH1F | X: k* [GeV/c] / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号低半分・最終核選別後 |
| `hLambdaKF_Chi2Ndf` | TH1F | X: Chi2Ndf / 200 bin [0, 100] | 採用ΛのKF品質QA（詳細は本文） |
| `hLambdaKF_DecayLength` | TH1F | X: DecayLength / 200 bin [0, 100] | 採用ΛのKF品質QA（詳細は本文） |
| `hLambdaKF_DecayLengthSignificance` | TH1F | X: DecayLengthSignificance / 200 bin [0, 100] | 採用ΛのKF品質QA（詳細は本文） |
| `hLambdaKF_TopoChi2Ndf` | TH1F | X: TopoChi2Ndf / 200 bin [0, 100] | 採用ΛのKF品質QA（詳細は本文） |
| `hLambdaPion_PtVsYLab_signal` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: p_{T} [GeV/c] / 200 bin [0, 10] | Λ信号候補のπ娘アクセプタンス・1候補1回 |
| `hLambdaProton_PtVsYLab_signal` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: p_{T} [GeV/c] / 200 bin [0, 10] | Λ信号候補の陽子娘アクセプタンス・1候補1回 |
| `hLambda_BetaVsEta_signal` | TH2D | X: #eta_{lab} / 200 bin [-5, 5]<br>Y: #beta = &#124;p&#124;/E / 200 bin [0, 1] | Λ信号候補の運動学的βアクセプタンス・1候補1回 |
| `hLambda_BetaVsYLab_signal` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: #beta = &#124;p&#124;/E / 200 bin [0, 1] | Λ信号候補の運動学的βアクセプタンス・1候補1回 |
| `hLambda_Eta` | TH1F | X: #eta / 100 bin [-2.5, 2.5] | 採用Λの単粒子運動学・全保存質量範囲 |
| `hLambda_InvMass` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin0` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin1` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin2` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin3` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin4` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin5` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin6` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin7` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_CentBin8` | TH1F | X: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_vs_DecayLength` | TH2F | X: Decay Length [cm] / 200 bin [0, 100]<br>Y: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_vs_Pt` | TH2F | X: p_{T} [GeV/c] / 200 bin [0, 10]<br>Y: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_vs_RefMultCorr` | TH2F | X: refMult_{corr} / 1000 bin [0, 1000]<br>Y: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_vs_TransDecayLength` | TH2F | X: L_{T} [cm] / 200 bin [0, 50]<br>Y: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_InvMass_vs_Y` | TH2F | X: y / 100 bin [-2, 2]<br>Y: M_{p#pi^{-}} [GeV/c^{2}] / 200 bin [1.05, 1.25] | 採用Λの全保存質量QA（詳細は本文・軸） |
| `hLambda_Phi` | TH1F | X: #phi [rad] / 100 bin [-3.15, 3.15] | 採用Λの単粒子運動学・全保存質量範囲 |
| `hLambda_Pt` | TH1F | X: p_{T} [GeV/c] / 200 bin [0, 10] | 採用Λの単粒子運動学・全保存質量範囲 |
| `hLambda_PtVsYLab_signal` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: p_{T} [GeV/c] / 200 bin [0, 10] | Λ信号候補のpT–Y_labアクセプタンス・1候補1回 |
| `hMixSamplerQA` | TH2D | X: status (FemtoMixingSampler enum) / 15 bin [-0.5, 14.5]<br>Y: channel index / 8 bin [-0.5, 7.5] | mixing処理の状態別カウント×8チャネル |
| `hN` | TH1I | X: dummy / 1 bin [0, 1] | 正常終了event数。全入力数とは別 |
| `hNLambdaPairs_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: N_{pairs} / 200 bin [0, 200] | KF Finder/Topo内Λ数。最終cut前 |
| `hNPionCand_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: N_{#pi^{-}} / 200 bin [0, 200] | KFへ投入したπ仮説数。最終娘数ではない |
| `hNProtonCand_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: N_{p} / 200 bin [0, 200] | KFへ投入した陽子仮説数。最終娘数ではない |
| `hNSigmaPion` | TH1F | X: n#sigma_{#pi} / 100 bin [-5, 5] | 採用Λのπ娘TPC nσ |
| `hNSigmaProton` | TH1F | X: n#sigma_{p} / 100 bin [-5, 5] | 採用Λの陽子娘TPC nσ |
| `hNTracks_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: N_{tracks} / 500 bin [0, 2000] | KF処理成功後のcent9–入力track数 |
| `hNucleus_BetaVsEta_he4` | TH2D | X: #eta_{lab} / 200 bin [-5, 5]<br>Y: #beta = &#124;p&#124;/E / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–η・全R・1候補1回 |
| `hNucleus_BetaVsEta_he4_highR_tofM2` | TH2D | X:  valid TOF and He4 mass2 window) / 200 bin [-5, 5]<br>Y: #eta_{lab} / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–η・高R・TOF＋m²窓・1候補1回 |
| `hNucleus_BetaVsEta_he4_lowR_noTof` | TH2D | X:  no TOF requirement) / 200 bin [-5, 5]<br>Y: #eta_{lab} / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–η・低R・TOF不要・1候補1回 |
| `hNucleus_BetaVsEta_he4_lowR_tofM2` | TH2D | X:  valid TOF and He4 mass2 window) / 200 bin [-5, 5]<br>Y: #eta_{lab} / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–η・低R・TOF＋m²窓・1候補1回 |
| `hNucleus_BetaVsEta_he4_lowR_tofMatch` | TH2D | X:  valid TOF / 200 bin [-5, 5]<br>Y:  no mass2 cut) / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–η・低R・有効TOFのみ・1候補1回 |
| `hNucleus_BetaVsYLab_he4` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: #beta = &#124;p&#124;/E / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–Y_lab・全R・1候補1回 |
| `hNucleus_BetaVsYLab_he4_highR_tofM2` | TH2D | X:  valid TOF and He4 mass2 window) / 200 bin [-5, 5]<br>Y: y_{lab} / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–Y_lab・高R・TOF＋m²窓・1候補1回 |
| `hNucleus_BetaVsYLab_he4_lowR_noTof` | TH2D | X:  no TOF requirement) / 200 bin [-5, 5]<br>Y: y_{lab} / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–Y_lab・低R・TOF不要・1候補1回 |
| `hNucleus_BetaVsYLab_he4_lowR_tofM2` | TH2D | X:  valid TOF and He4 mass2 window) / 200 bin [-5, 5]<br>Y: y_{lab} / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–Y_lab・低R・TOF＋m²窓・1候補1回 |
| `hNucleus_BetaVsYLab_he4_lowR_tofMatch` | TH2D | X:  valid TOF / 200 bin [-5, 5]<br>Y:  no mass2 cut) / 200 bin [0, 1] | ⁴He候補アクセプタンス・β–Y_lab・低R・有効TOFのみ・1候補1回 |
| `hNucleus_PtVsYLab_he3` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: p_{T} [GeV/c] / 200 bin [0, 10] | ³He候補アクセプタンス・pT–Y_lab・全R・1候補1回 |
| `hNucleus_PtVsYLab_he3_as_he4` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: p_{T} [GeV/c] / 200 bin [0, 10] | ³He→⁴He質量仮定候補アクセプタンス・pT–Y_lab・全R・1候補1回 |
| `hNucleus_PtVsYLab_he4` | TH2D | X: y_{lab} / 200 bin [-5, 5]<br>Y: p_{T} [GeV/c] / 200 bin [0, 10] | ⁴He候補アクセプタンス・pT–Y_lab・全R・1候補1回 |
| `hNucleus_PtVsYLab_he4_highR_tofM2` | TH2D | X:  valid TOF and He4 mass2 window) / 200 bin [-5, 5]<br>Y: y_{lab} / 200 bin [0, 10] | ⁴He候補アクセプタンス・pT–Y_lab・高R・TOF＋m²窓・1候補1回 |
| `hNucleus_PtVsYLab_he4_lowR_noTof` | TH2D | X:  no TOF requirement) / 200 bin [-5, 5]<br>Y: y_{lab} / 200 bin [0, 10] | ⁴He候補アクセプタンス・pT–Y_lab・低R・TOF不要・1候補1回 |
| `hNucleus_PtVsYLab_he4_lowR_tofM2` | TH2D | X:  valid TOF and He4 mass2 window) / 200 bin [-5, 5]<br>Y: y_{lab} / 200 bin [0, 10] | ⁴He候補アクセプタンス・pT–Y_lab・低R・TOF＋m²窓・1候補1回 |
| `hNucleus_PtVsYLab_he4_lowR_tofMatch` | TH2D | X:  valid TOF / 200 bin [-5, 5]<br>Y:  no mass2 cut) / 200 bin [0, 10] | ⁴He候補アクセプタンス・pT–Y_lab・低R・有効TOFのみ・1候補1回 |
| `hRawMult` | TH1F | X: mult / 1000 bin [0, 1000] | centrality入力raw multiplicity（FXT: fxtMult） |
| `hRawMult_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: rawMult / 1000 bin [0, 1000] | KF処理成功後のcent9–raw multiplicity |
| `hRawMult_vs_RefMultCorr` | TH2F | X: refMult_{corr} / 1000 bin [0, 1000]<br>Y: rawMult / 1000 bin [0, 1000] | KF処理成功後の補正–raw multiplicity |
| `hRefMult` | TH1F | X: RefMult / 1000 bin [0, 1000] | PicoEventのrefMult。event除去前 |
| `hRefMultCorr` | TH1F | X: refMult_{corr} / 1000 bin [0, 1000] | 補正multiplicity |
| `hRefMultCorr_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: refMult_{corr} / 1000 bin [0, 1000] | KF処理成功後のcent9–補正multiplicity |
| `hRefMultVsNTOFMatch` | TH2F | X: nBTOFMatch / 500 bin [0, 500]<br>Y: refMult / 1000 bin [0, 1000] | event除去前のnBTOFMatch–raw multiplicity |
| `hRefMultVsNTOFMatchAfter` | TH2F | X: nBTOFMatch / 500 bin [0, 500]<br>Y: refMult / 1000 bin [0, 1000] | event/pileup選択後のnBTOFMatch–raw multiplicity |
| `hRefMultWeight` | TH1F | X: weight / 125 bin [0.5, 3] | centrality補正weightの分布 |
| `hTofMatchMult_vs_Cent9` | TH2F | X: cent9 (0=70-80%%, 8=0-5%%) / 9 bin [-0.5, 8.5]<br>Y: nBTOFMatch / 500 bin [0, 500] | KF処理成功後のcent9–nBTOFMatch |
| `hVz` | TH1F | X: V_{z} [cm] / 400 bin [-300, 300] | PVのz。event除去前 |

</details>

<details>
<summary>true/（840個）：全キー・型・軸・説明</summary>

| ROOTキー | 型 | 軸 | 説明 |
|---|---|---|---|
| `true/hDedxP` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | 品質通過trackのR–dE/dx、最終PID前 |
| `true/hDedxP_3He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・3He・最終核種選別とは独立 |
| `true/hDedxP_3He_m2` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TOF m²仮説QA・R–dE/dx・3He・最終核種選別とは独立 |
| `true/hDedxP_4He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・4He・最終核種選別とは独立 |
| `true/hDedxP_4He_TOFControl` | TH2F | X: R=&#124;p/q&#124; [GeV/c] / 600 bin [0, 6]<br>Y: dE/dx [keV/cm] / 1200 bin [0, 120] | TPC PID前の条件付きTOF対照試料のR–dE/dx |
| `true/hDedxP_4He_selected` | TH2F | X: R=&#124;p/q&#124; [GeV/c] / 600 bin [0, 6]<br>Y: dE/dx [keV/cm] / 1200 bin [0, 120] | 最終⁴He候補のR–dE/dx |
| `true/hDedxP_K` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・K・最終核種選別とは独立 |
| `true/hDedxP_cut` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TOF未match/無効または広いTOF QA条件を通らないtrack |
| `true/hDedxP_d` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・d・最終核種選別とは独立 |
| `true/hDedxP_d_m2` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TOF m²仮説QA・R–dE/dx・d・最終核種選別とは独立 |
| `true/hDedxP_e` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・e・最終核種選別とは独立 |
| `true/hDedxP_else` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | 通常e/π/K/p仮説から離れたtrackのR–dE/dx |
| `true/hDedxP_p` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・p・最終核種選別とは独立 |
| `true/hDedxP_pi` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・pi・最終核種選別とは独立 |
| `true/hDedxP_t` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TPC仮説QA・R–dE/dx・t・最終核種選別とは独立 |
| `true/hDedxP_t_m2` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  dE/dx (keV/cm) / 1200 bin [0, 120] | TOF m²仮説QA・R–dE/dx・t・最終核種選別とは独立 |
| `true/hDphiDeta_proton_3He` | TH2F | X: #Delta#phi* / 200 bin [-0.25, 0.25]<br>Y: #Delta#eta / 200 bin [-0.25, 0.25] | Λ娘陽子–³Heの近接pair QA・signal SE・旧固定磁場規約 |
| `true/hDphiDeta_proton_3HeAs4He` | TH2F | X: #Delta#phi* / 200 bin [-0.25, 0.25]<br>Y: #Delta#eta / 200 bin [-0.25, 0.25] | Λ娘陽子–³He→⁴He質量仮定の近接pair QA・signal SE・旧固定磁場規約 |
| `true/hDphiDeta_proton_3HeAs4He_EventField` | TH2F | X: Delta phi* / 200 bin [-0.25, 0.25]<br>Y: Delta eta / 200 bin [-0.25, 0.25] | Λ娘陽子–³He→⁴He質量仮定の近接pair QA・signal SE・イベント磁場 |
| `true/hDphiDeta_proton_3He_EventField` | TH2F | X: Delta phi* / 200 bin [-0.25, 0.25]<br>Y: Delta eta / 200 bin [-0.25, 0.25] | Λ娘陽子–³Heの近接pair QA・signal SE・イベント磁場 |
| `true/hDphiDeta_proton_4He` | TH2F | X: #Delta#phi* / 200 bin [-0.25, 0.25]<br>Y: #Delta#eta / 200 bin [-0.25, 0.25] | Λ娘陽子–⁴Heの近接pair QA・signal SE・旧固定磁場規約 |
| `true/hDphiDeta_proton_4He_EventField` | TH2F | X: Delta phi* / 200 bin [-0.25, 0.25]<br>Y: Delta eta / 200 bin [-0.25, 0.25] | Λ娘陽子–⁴Heの近接pair QA・signal SE・イベント磁場 |
| `true/hHe4PidStage` | TH1F | X: Stage / 6 bin [-0.5, 5.5] | TPC/低R/高R/有効TOF/m²窓/最終候補の段階別track数 |
| `true/hKstarMass_3HeAs4He_CentBin0` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=0 |
| `true/hKstarMass_3HeAs4He_CentBin1` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=1 |
| `true/hKstarMass_3HeAs4He_CentBin2` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=2 |
| `true/hKstarMass_3HeAs4He_CentBin3` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=3 |
| `true/hKstarMass_3HeAs4He_CentBin4` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=4 |
| `true/hKstarMass_3HeAs4He_CentBin5` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=5 |
| `true/hKstarMass_3HeAs4He_CentBin6` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=6 |
| `true/hKstarMass_3HeAs4He_CentBin7` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=7 |
| `true/hKstarMass_3HeAs4He_CentBin8` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・SE・Λ全保存質量・最終核選別後・cent9=8 |
| `true/hKstarMass_3He_CentBin0` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=0 |
| `true/hKstarMass_3He_CentBin1` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=1 |
| `true/hKstarMass_3He_CentBin2` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=2 |
| `true/hKstarMass_3He_CentBin3` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=3 |
| `true/hKstarMass_3He_CentBin4` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=4 |
| `true/hKstarMass_3He_CentBin5` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=5 |
| `true/hKstarMass_3He_CentBin6` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=6 |
| `true/hKstarMass_3He_CentBin7` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=7 |
| `true/hKstarMass_3He_CentBin8` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・SE・Λ全保存質量・最終核選別後・cent9=8 |
| `true/hKstarMass_4He_CentBin0` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=0 |
| `true/hKstarMass_4He_CentBin1` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=1 |
| `true/hKstarMass_4He_CentBin2` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=2 |
| `true/hKstarMass_4He_CentBin3` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=3 |
| `true/hKstarMass_4He_CentBin4` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=4 |
| `true/hKstarMass_4He_CentBin5` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=5 |
| `true/hKstarMass_4He_CentBin6` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=6 |
| `true/hKstarMass_4He_CentBin7` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=7 |
| `true/hKstarMass_4He_CentBin8` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・SE・Λ全保存質量・最終核選別後・cent9=8 |
| `true/hKstarMt_3HeAs4He_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarMt_3HeAs4He_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarMt_3HeAs4He_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarMt_3HeAs4He_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarMt_3HeAs4He_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarMt_3HeAs4He_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarMt_3HeAs4He_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarMt_3HeAs4He_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarMt_3HeAs4He_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarMt_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarMt_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarMt_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarMt_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarMt_3He_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarMt_3He_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarMt_3He_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarMt_3He_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarMt_3He_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarMt_3He_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarMt_3He_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarMt_3He_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarMt_3He_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarMt_3He_SBNeg_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarMt_3He_SBNeg_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarMt_3He_SBNeg_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarMt_3He_SBNeg_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarMt_3He_SBNeg_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarMt_3He_SBNeg_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarMt_3He_SBNeg_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarMt_3He_SBNeg_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarMt_3He_SBNeg_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarMt_3He_SBPos_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarMt_3He_SBPos_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarMt_3He_SBPos_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarMt_3He_SBPos_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarMt_3He_SBPos_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarMt_3He_SBPos_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarMt_3He_SBPos_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarMt_3He_SBPos_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarMt_3He_SBPos_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarMt_3He_signalHigh_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarMt_3He_signalHigh_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarMt_3He_signalHigh_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarMt_3He_signalHigh_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarMt_3He_signalHigh_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarMt_3He_signalHigh_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarMt_3He_signalHigh_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarMt_3He_signalHigh_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarMt_3He_signalHigh_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarMt_3He_signalLow_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarMt_3He_signalLow_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarMt_3He_signalLow_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarMt_3He_signalLow_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarMt_3He_signalLow_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarMt_3He_signalLow_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarMt_3He_signalLow_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarMt_3He_signalLow_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarMt_3He_signalLow_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarMt_4He_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarMt_4He_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarMt_4He_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarMt_4He_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarMt_4He_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarMt_4He_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarMt_4He_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarMt_4He_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarMt_4He_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarMt_4He_SBNeg_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarMt_4He_SBNeg_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarMt_4He_SBNeg_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarMt_4He_SBNeg_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarMt_4He_SBNeg_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarMt_4He_SBNeg_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarMt_4He_SBNeg_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarMt_4He_SBNeg_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarMt_4He_SBNeg_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarMt_4He_SBPos_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarMt_4He_SBPos_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarMt_4He_SBPos_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarMt_4He_SBPos_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarMt_4He_SBPos_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarMt_4He_SBPos_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarMt_4He_SBPos_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarMt_4He_SBPos_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarMt_4He_SBPos_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarMt_4He_signalHigh_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarMt_4He_signalHigh_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarMt_4He_signalHigh_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarMt_4He_signalHigh_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarMt_4He_signalHigh_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarMt_4He_signalHigh_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarMt_4He_signalHigh_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarMt_4He_signalHigh_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarMt_4He_signalHigh_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarMt_4He_signalLow_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarMt_4He_signalLow_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarMt_4He_signalLow_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarMt_4He_signalLow_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarMt_4He_signalLow_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarMt_4He_signalLow_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarMt_4He_signalLow_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarMt_4He_signalLow_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarMt_4He_signalLow_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=8 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=0 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=1 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=2 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=3 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=4 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=5 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=6 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=7 |
| `true/hKstarNSigmaPrePid_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=8 |
| `true/hKstarNSigma_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarNSigma_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarNSigma_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarNSigma_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarNSigma_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarNSigma_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarNSigma_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarNSigma_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarNSigma_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarNSigma_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarNSigma_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarNSigma_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarNSigma_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarNSigma_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarNSigma_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarNSigma_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarNSigma_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarNSigma_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarNSigma_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarNSigma_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarNSigma_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarNSigma_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarNSigma_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarNSigma_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarNSigma_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarNSigma_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarNSigma_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarNSigma_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarNSigma_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarNSigma_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarNSigma_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarNSigma_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarNSigma_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarNSigma_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarNSigma_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarNSigma_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarNSigma_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarNSigma_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarNSigma_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarNSigma_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarNSigma_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarNSigma_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarNSigma_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarNSigma_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarNSigma_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarNSigma_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarNSigma_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarNSigma_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarNSigma_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarNSigma_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarNSigma_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarNSigma_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarNSigma_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarNSigma_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarNSigma_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・prePID・cent9=8 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=0 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=1 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=2 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=3 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=4 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=5 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=6 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=7 |
| `true/hKstarPOverQPrePid_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・prePID・cent9=8 |
| `true/hKstarPOverQ_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号高半分・最終核選別後・cent9=8 |
| `true/hKstarPOverQ_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=0 |
| `true/hKstarPOverQ_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=1 |
| `true/hKstarPOverQ_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=2 |
| `true/hKstarPOverQ_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=3 |
| `true/hKstarPOverQ_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=4 |
| `true/hKstarPOverQ_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=5 |
| `true/hKstarPOverQ_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=6 |
| `true/hKstarPOverQ_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=7 |
| `true/hKstarPOverQ_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・SE・Λ信号低半分・最終核選別後・cent9=8 |
| `true/hKstar_3He` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後 |
| `true/hKstar_3HeAs4He` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後 |
| `true/hKstar_3HeAs4He_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstar_3HeAs4He_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstar_3HeAs4He_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstar_3HeAs4He_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstar_3HeAs4He_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstar_3HeAs4He_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstar_3HeAs4He_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstar_3HeAs4He_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstar_3HeAs4He_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstar_3HeAs4He_SBNeg` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstar_3HeAs4He_SBNeg_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstar_3HeAs4He_SBPos` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後 |
| `true/hKstar_3HeAs4He_SBPos_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstar_3HeAs4He_SBPos_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstar_3HeAs4He_SBPos_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstar_3HeAs4He_SBPos_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstar_3HeAs4He_SBPos_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstar_3HeAs4He_SBPos_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstar_3HeAs4He_SBPos_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstar_3HeAs4He_SBPos_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstar_3HeAs4He_SBPos_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstar_3He_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstar_3He_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstar_3He_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstar_3He_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstar_3He_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstar_3He_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstar_3He_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstar_3He_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstar_3He_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstar_3He_SBNeg` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後 |
| `true/hKstar_3He_SBNeg_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstar_3He_SBNeg_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstar_3He_SBNeg_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstar_3He_SBNeg_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstar_3He_SBNeg_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstar_3He_SBNeg_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstar_3He_SBNeg_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstar_3He_SBNeg_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstar_3He_SBNeg_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstar_3He_SBPos` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後 |
| `true/hKstar_3He_SBPos_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstar_3He_SBPos_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstar_3He_SBPos_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstar_3He_SBPos_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstar_3He_SBPos_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstar_3He_SBPos_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstar_3He_SBPos_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstar_3He_SBPos_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstar_3He_SBPos_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hKstar_4He` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後 |
| `true/hKstar_4He_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=0 |
| `true/hKstar_4He_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=1 |
| `true/hKstar_4He_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=2 |
| `true/hKstar_4He_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=3 |
| `true/hKstar_4He_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=4 |
| `true/hKstar_4He_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=5 |
| `true/hKstar_4He_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=6 |
| `true/hKstar_4He_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=7 |
| `true/hKstar_4He_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ信号窓・最終核選別後・cent9=8 |
| `true/hKstar_4He_SBNeg` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後 |
| `true/hKstar_4He_SBNeg_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=0 |
| `true/hKstar_4He_SBNeg_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=1 |
| `true/hKstar_4He_SBNeg_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=2 |
| `true/hKstar_4He_SBNeg_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=3 |
| `true/hKstar_4He_SBNeg_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=4 |
| `true/hKstar_4He_SBNeg_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=5 |
| `true/hKstar_4He_SBNeg_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=6 |
| `true/hKstar_4He_SBNeg_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=7 |
| `true/hKstar_4He_SBNeg_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ低側sideband・最終核選別後・cent9=8 |
| `true/hKstar_4He_SBPos` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後 |
| `true/hKstar_4He_SBPos_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=0 |
| `true/hKstar_4He_SBPos_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=1 |
| `true/hKstar_4He_SBPos_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=2 |
| `true/hKstar_4He_SBPos_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=3 |
| `true/hKstar_4He_SBPos_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=4 |
| `true/hKstar_4He_SBPos_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=5 |
| `true/hKstar_4He_SBPos_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=6 |
| `true/hKstar_4He_SBPos_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=7 |
| `true/hKstar_4He_SBPos_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・SE・Λ高側sideband・最終核選別後・cent9=8 |
| `true/hM2P` | TH2F | X: p (GeV/c) / 600 bin [0, 6]<br>Y: m^{2} (GeV^{2}/c^{4}) / 1000 bin [-0.5, 25] | 有効TOFを持つ品質通過trackのR–(m/q)² |
| `true/hM2P_4He_TPC` | TH2F | X: R=&#124;p/q&#124; [GeV/c] / 600 bin [0, 6]<br>Y: (m/q)^{2} [(GeV/c^{2})^{2}] / 1000 bin [-0.5, 25] | TPC best=⁴He、高R・有効TOF、m²窓前 |
| `true/hM2P_4He_selected` | TH2F | X: R=&#124;p/q&#124; [GeV/c] / 600 bin [0, 6]<br>Y: (m/q)^{2} [(GeV/c^{2})^{2}] / 1000 bin [-0.5, 25] | 最終⁴He候補のうち有効TOFありのR–(m/q)² |
| `true/hNSigmaHe4VsP_TOFControl` | TH2F | X: R=&#124;p/q&#124; [GeV/c] / 600 bin [0, 6]<br>Y: n#sigma_{He4} / 400 bin [-10, 10] | 同対照試料のR–⁴He nσ |
| `true/hPvsEta_3He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  #eta / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsEta_4He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  #eta / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsEta_d` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  #eta / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsEta_t` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  #eta / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsPt_3He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  p_{T} (GeV/c) / 600 bin [0, 6] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsPt_4He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  p_{T} (GeV/c) / 600 bin [0, 6] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsPt_d` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  p_{T} (GeV/c) / 600 bin [0, 6] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsPt_t` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  p_{T} (GeV/c) / 600 bin [0, 6] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsY_3He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  y / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsY_4He` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  y / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsY_d` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  y / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hPvsY_t` | TH2F | X:  p (GeV/c) / 600 bin [0, 6]<br>Y:  y / 100 bin [-3, 3] | 各TPC核仮説通過trackの物理P対η/Y/pT・HeはZ=2補正 |
| `true/hQlab_3He` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He・SE・Λ信号窓・最終核選別後 |
| `true/hQlab_3HeAs4He` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He→⁴He質量仮定・SE・Λ信号窓・最終核選別後 |
| `true/hQlab_3HeAs4He_SBNeg` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He→⁴He質量仮定・SE・Λ低側sideband・最終核選別後 |
| `true/hQlab_3HeAs4He_SBPos` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He→⁴He質量仮定・SE・Λ高側sideband・最終核選別後 |
| `true/hQlab_3He_SBNeg` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He・SE・Λ低側sideband・最終核選別後 |
| `true/hQlab_3He_SBPos` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He・SE・Λ高側sideband・最終核選別後 |
| `true/hQlab_4He` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–⁴He・SE・Λ信号窓・最終核選別後 |
| `true/hQlab_4He_SBNeg` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–⁴He・SE・Λ低側sideband・最終核選別後 |
| `true/hQlab_4He_SBPos` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–⁴He・SE・Λ高側sideband・最終核選別後 |

</details>

<details>
<summary>mix/（801個）：全キー・型・軸・説明</summary>

| ROOTキー | 型 | 軸 | 説明 |
|---|---|---|---|
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin0` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=0 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin1` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=1 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin2` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=2 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin3` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=3 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin4` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=4 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin5` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=5 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin6` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=6 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin7` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=7 |
| `mix/hKstarMass_Mixed_3HeAs4He_CentBin8` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He→⁴He質量仮定・ME・Λ全保存質量・最終核選別後・cent9=8 |
| `mix/hKstarMass_Mixed_3He_CentBin0` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=0 |
| `mix/hKstarMass_Mixed_3He_CentBin1` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=1 |
| `mix/hKstarMass_Mixed_3He_CentBin2` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=2 |
| `mix/hKstarMass_Mixed_3He_CentBin3` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=3 |
| `mix/hKstarMass_Mixed_3He_CentBin4` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=4 |
| `mix/hKstarMass_Mixed_3He_CentBin5` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=5 |
| `mix/hKstarMass_Mixed_3He_CentBin6` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=6 |
| `mix/hKstarMass_Mixed_3He_CentBin7` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=7 |
| `mix/hKstarMass_Mixed_3He_CentBin8` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–³He・ME・Λ全保存質量・最終核選別後・cent9=8 |
| `mix/hKstarMass_Mixed_4He_CentBin0` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=0 |
| `mix/hKstarMass_Mixed_4He_CentBin1` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=1 |
| `mix/hKstarMass_Mixed_4He_CentBin2` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=2 |
| `mix/hKstarMass_Mixed_4He_CentBin3` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=3 |
| `mix/hKstarMass_Mixed_4He_CentBin4` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=4 |
| `mix/hKstarMass_Mixed_4He_CentBin5` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=5 |
| `mix/hKstarMass_Mixed_4He_CentBin6` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=6 |
| `mix/hKstarMass_Mixed_4He_CentBin7` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=7 |
| `mix/hKstarMass_Mixed_4He_CentBin8` | TH2F | X: X / 200 bin [0, 1]<br>Y: 件数軸 / 200 bin [1.05, 1.25] | k*–Λ不変質量（全質量）・Λ–⁴He・ME・Λ全保存質量・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3HeAs4He_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3He_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3He_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3He_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3He_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3He_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3He_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3He_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3He_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3He_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3He_SBNeg_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3He_SBPos_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3He_signalHigh_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_3He_signalLow_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_4He_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_4He_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_4He_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_4He_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_4He_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_4He_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_4He_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_4He_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_4He_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_4He_SBNeg_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_4He_SBPos_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_4He_signalHigh_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin0` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin1` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin2` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin3` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin4` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin5` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin6` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin7` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarMt_Mixed_4He_signalLow_CentBin8` | TH2D | X: k* (GeV/c) / 200 bin [0, 1]<br>Y: m_{T} (GeV/c^{2}) / 200 bin [0, 10] | k*–mT・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=8 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=0 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=1 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=2 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=3 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=4 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=5 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=6 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=7 |
| `mix/hKstarNSigmaPrePid_Mixed_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=8 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{3He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarNSigma_Mixed_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: n#sigma_{4He} / 120 bin [-12, 12] | k*–TPC nσ・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・prePID・cent9=8 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=0 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=1 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=2 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=3 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=4 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=5 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=6 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=7 |
| `mix/hKstarPOverQPrePid_Mixed_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・prePID・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3HeAs4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He→⁴He質量仮定・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_3He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–³He・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_4He_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_4He_SBNeg_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_4He_SBPos_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_4He_signalHigh_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号高半分・最終核選別後・cent9=8 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin0` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=0 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin1` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=1 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin2` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=2 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin3` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=3 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin4` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=4 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin5` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=5 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin6` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=6 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin7` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=7 |
| `mix/hKstarPOverQ_Mixed_4He_signalLow_CentBin8` | TH2D | X: k* [GeV/c] / 200 bin [0, 1]<br>Y: raw &#124;p/q&#124; [GeV/c] / 120 bin [0, 6] | k*–R・Λ–⁴He・ME・Λ信号低半分・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_3He` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後 |
| `mix/hKstar_Mixed_3HeAs4He` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_3HeAs4He_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_3HeAs4He_SBNeg_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_3HeAs4He_SBPos_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_3He_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_3He_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_3He_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_3He_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_3He_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_3He_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_3He_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_3He_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_3He_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_3He_SBNeg` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_3He_SBNeg_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_3He_SBPos` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_3He_SBPos_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–³He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_4He` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後 |
| `mix/hKstar_Mixed_4He_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_4He_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_4He_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_4He_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_4He_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_4He_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_4He_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_4He_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_4He_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ信号窓・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_4He_SBNeg` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_4He_SBNeg_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ低側sideband・最終核選別後・cent9=8 |
| `mix/hKstar_Mixed_4He_SBPos` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin0` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=0 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin1` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=1 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin2` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=2 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin3` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=3 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin4` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=4 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin5` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=5 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin6` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=6 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin7` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=7 |
| `mix/hKstar_Mixed_4He_SBPos_CentBin8` | TH1F | X: k* (GeV/c) / 200 bin [0, 1] | k*・Λ–⁴He・ME・Λ高側sideband・最終核選別後・cent9=8 |
| `mix/hQlab_Mixed_3He` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He・ME・Λ信号窓・最終核選別後 |
| `mix/hQlab_Mixed_3HeAs4He` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He→⁴He質量仮定・ME・Λ信号窓・最終核選別後 |
| `mix/hQlab_Mixed_3HeAs4He_SBNeg` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He→⁴He質量仮定・ME・Λ低側sideband・最終核選別後 |
| `mix/hQlab_Mixed_3HeAs4He_SBPos` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He→⁴He質量仮定・ME・Λ高側sideband・最終核選別後 |
| `mix/hQlab_Mixed_3He_SBNeg` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He・ME・Λ低側sideband・最終核選別後 |
| `mix/hQlab_Mixed_3He_SBPos` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–³He・ME・Λ高側sideband・最終核選別後 |
| `mix/hQlab_Mixed_4He` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–⁴He・ME・Λ信号窓・最終核選別後 |
| `mix/hQlab_Mixed_4He_SBNeg` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–⁴He・ME・Λ低側sideband・最終核選別後 |
| `mix/hQlab_Mixed_4He_SBPos` | TH1F | X: q_lab (GeV/c) / 200 bin [0, 2] | q_lab・Λ–⁴He・ME・Λ高側sideband・最終核選別後 |

</details>

### トップ階層：hNuclearPidPairJoint（THnSparseD、1個）

全8軸を同時保存。カテゴリの意味と投影コードは本文PID QA節を参照。

| 軸番号 | 軸名 | bin数・範囲 | 意味 |
|---|---|---|---|
| 0 | Kstar | 200 [0, 1] | k* [GeV/c] |
| 1 | POverQ | 120 [0, 6] | raw &#124;p/q&#124; [GeV/c] |
| 2 | NSigma | 120 [-12, 12] | source-nucleus nSigma |
| 3 | Cent9 | 9 [-0.5, 8.5] | cent9 |
| 4 | Representation | 3 [-0.5, 2.5] | 0=he4, 1=he3, 2=he3as4 |
| 5 | MixingMode | 2 [-0.5, 1.5] | 0=SE, 1=ME |
| 6 | MassClass | 4 [-0.5, 3.5] | 0=SBNeg, 1=signalLow, 2=signalHigh, 3=SBPos |
| 7 | PidStage | 2 [-0.5, 1.5] | 0=selected, 1=prePID |


## 大規模統合時のSparse保存

`hNuclearPidPairJoint` は占有bin数に応じて増大し、ROOTの約1 GB/単一オブジェクトのシリアライズ上限を超えることがある。通常1,786ヒストを主ROOTに統合し、Sparseは同一軸・Sumw2のまま入力群ごとの別ROOTに保存する。この場合は各Sparseへ同じ範囲条件を適用し、低次元へ射影してから加算する。巨大Sparseを先に1個へ加算・保存しない。各入力群を重複なく網羅し、Entries、座標、内容・誤差、underflow/overflow、軸定義を保存後に確認する。

`FemtoNuclearDedxCalibration` の出典フィールドは係数を同定するための記録であり、fitスクリプトや出典ROOTを実行時に開かない。本番が必要とする係数は校正YAMLに全て含む。

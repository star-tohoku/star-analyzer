# KFParticle Λ再構成の構成と設定

## 再構成の構成とROOT環境

`anaLambda_KFParticle.C → StLambdaKFParticleMaker → StPicoKFParticleInterface → KFParticleTopoReconstructor → KFParticleFinder`。

macroはStChainの構築・イベントループ・mainconfパス引渡しだけを行う。
track covariance、PID仮説、PV分類、Finderによる±3122探索、raw massとPV拘束コピーの量はコンパイル済みKFライブラリ内で処理する。
scalar手動pair fitへのfallbackはない。既存StLambda／StLambdaMakerのAPI・配置・出力は変更しない。

- STAR SL24y、sl73_x8664_gcc485、ROOT 5.34/38、GCC 4.8.5、STAR batch相当Singularity。
- .DEV2の39ファイル／11実装を同じsnapshotから取り込み、`star_analyzer_kfp`へ隔離した。
  snapshot ID: `033c38b7128ea7cb351a2785a3a65ff0ec2e19e391ca3680d8ffd4d17c27ceb5`。
- raw parent massは親質量拘束なし。ただし採用TopoのPV χ²/NDF < 3等の内部選別後である。

## analysis.modeとvertex cut

- 既存 `setup_config_from_analysisinfo.py` は、analysis-infoの `analysis.mode` から**設定生成時にVzのみ**を切り替える。
  `fxtmult` は198–202 cm、`refmult` は−100–100 cm。
- 既存 `ConfigManager` はanalysis-infoからanaNameを読むが、実行時のmode別vertex変更はしない。
- `EventCutConfig` のXY中心は指定がなければ(0,0)。

KF専用 `KfEventSelection` がmainconfの `analysis:` と `event:` を読み、
`analysis.mode` に対応する `event.vertexByMode.<mode>` を実行時に選択する。

| mode | Vz [cm] | XY中心 [cm] | 半径 [cm] |
| --- | --- | --- | --- |
| refmult（collider） | [-100,100] | (0,0) | 2 |
| fxtmult（fixed target） | [198,202] | (-0.4,-2.0) | 2 |

これらの数値はC++に固定せず、KF専用event YAML内に置いた。

設定例（mainconfが参照するevent YAML）:

```yaml
vertexByMode:
  refmult:
    vzRange: [-100.0, 100.0]
    center: [0.0, 0.0]
    radius: 2.0
  fxtmult:
    vzRange: [198.0, 202.0]
    center: [-0.4, -2.0]
    radius: 2.0
```

優先順位・安全性:

- KFではprofileのVz／XY／半径がflatなminVz／maxVz／vtxCenterX/Y／maxVrより優先。
- RefMult、valid VPD、maxNTrは既存event設定を引き続き適用。
- 共有EventCutConfig singletonを変更せず、KF内へコピーする。
- unknown／missing mode、必要profile欠落、重複、不正型・次元・非有限値・範囲違反は初期化エラー。
- enabledなcentrality.modeがanalysis.modeと不一致なら拒否。誤設定後に以前のprofileを使い続けない。
- 既存parserは階層を平坦化するため、profile内ではflatキーと異なるvzRange／center／radiusを用いた。
  KF helperはyaml-cppで階層を明示的に読む。
- AuAu19もKF専用eventコピーを配置し、旧Helix設定の共有参照から分離。
- モードだけを変えても別データ用のPID、centrality校正、histogram軸まで切り替わるわけではない。
- 汎用設定生成scriptは変更していない。KFへのmode適用にこのscriptを再実行しないこと。
  generic macro生成は独自KF macroを上書きする可能性がある。


## ビルドと検証

```bash
bash script/singularity_make.sh config/mainconf/main_auau13p5_anaLambda_KFParticle.yaml all -j2
bash script/singularity_make.sh config/mainconf/main_auau13p5_anaLambda_KFParticle.yaml --no-clean test-kfparticle-full-chain test-kfparticle-pico-adapter KF_TEST_CUTS=config/cuts/kf/kf_auau13p5_anaLambda_KFParticle.yaml -j2
```

full-chainテストはΛ/anti-Λ、磁場符号、ID/Sort、複数PID、raw/拘束コピー、resetを確認する。Pico adapterは6 parameter/21 covariance、TOFのindex 0と不正参照、dE/dx pull、strict YAML、mode別vertexと境界・不正値を検査する。

出力の実効設定・backend fingerprintと完了状態を読み戻して確認する。小統計の動作検証は物理purity・検出効率の評価とは別である。

- [KF core](../../../StRoot/KFParticle/README.md)
- [Snapshot provenance](../../../StRoot/KFParticle/PROVENANCE.md)
- [Pico adapter](../../../StMaker/kfparticle/README.md)
- [Lambda–nucleus profile](../../femto/implementation/femto_lambda_nuclei_kfparticle_20260930.md)

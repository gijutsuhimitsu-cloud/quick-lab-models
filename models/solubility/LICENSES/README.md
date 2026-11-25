# このモデルに関わるライセンス

*English: this folder holds no license texts of its own. Which terms apply to what is described in the repository-root LICENSE (scope of the MIT License for program code) and LICENSE-DATA.md (weights, tables, reference values and example data). RDKit.js and the third-party components and fonts bundled in its WebAssembly have their own notices in vendor/rdkit/. Upstream data sources and their terms (CC BY 4.0 or CC0 1.0) are listed in ../ATTRIBUTION.md.*

このフォルダには、原文を置いていません。どの物にどの条件がかかるかは、次の場所に書いてあります。

| 何について | 見る所 |
|---|---|
| プログラム（`js/`・`train/` の道具・`tests/` のプログラム・`examples/`）に MIT がかかる範囲 | リポジトリの直下の [LICENSE](../../../LICENSE) |
| 重み（`weights/`）、学習の土台の表（`train/base_data/`）、照合用の値（`tests/reference_*.json`）、例の CSV（`train/examples/`）、除外の key の一覧（`train/rebuild/bradley_excluded_keys.txt`）の条件 | リポジトリの直下の [LICENSE-DATA.md](../../../LICENSE-DATA.md) |
| 元のデータごとの作者・版・ライセンス・変えたこと | [../ATTRIBUTION.md](../ATTRIBUTION.md) |
| RDKit.js（このフォルダには写さず、共有の物を使う）と、その中に入っている第三者の部品・フォントの表示 | リポジトリの [vendor/rdkit/](../../../vendor/rdkit/) |

元のデータのライセンス（詳しくは [../ATTRIBUTION.md](../ATTRIBUTION.md)）：

| データ | ライセンス | 義務 |
|---|---|---|
| BigSolDB 2.2、MixtureSolDB、Compounds with known Abraham descriptors | CC BY 4.0 | 出典の表示と変更の説明（ATTRIBUTION.md を残す） |
| AqSolDB、Abraham Solvent Master Dataset（AbraLlama）、Jean-Claude Bradley Open Melting Point Dataset | CC0 1.0 | なし（CC0 は登録した人が手放せる権利の範囲。元の資料ごとの第三者の権利は ATTRIBUTION.md の注意を見る） |

原文を同梱する義務のあるデータ（Apache 2.0 など）は、このモデルでは使っていません。

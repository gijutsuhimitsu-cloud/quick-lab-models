# 出典（溶けるか・混ざるか 組み合わせ版 0.2）

*English: data sources used to train and build the weights and tables in `weights/` and the reference values in `tests/`. Our own contributions (trained weights, the selection, processing and added parts of the tables, reference values) are offered under CC BY 4.0; upstream material keeps its original terms (CC BY 4.0 or CC0 1.0). When you redistribute the weights or tables, keep this attribution. All sources allow commercial use without share-alike terms.*

こちらが許諾できる物（学習済みの重み、表の選び方・加工・追加の部分、照合用の値）は CC BY 4.0 で提供します。元のデータには、それぞれの配布元の条件（CC BY 4.0 か CC0 1.0）が引き続きかかります。重みや表を配る時は、このページを一緒に付けてください。

版・作者のくわしい記録と、RDKit.js・式の出典は、同じデータを使う [../solubility/ATTRIBUTION.md](../solubility/ATTRIBUTION.md) にもあります。**溶質の名前は、配布元の表から写していません**（どの表も、溶質は構造〈SMILES〉で持ちます。名前があるのは、例の固体の日本語の名前〈こちらで付けた物〉と溶媒の名前だけです）。

| 用途 | 出典（作者、配布元、DOI） | ライセンス | こちらで変えたこと |
|---|---|---|---|
| 溶解度の実測（25 °C 付近と全温度）、融点、融解熱（DSC）、モデルとネットワークの学習、式の係数の当てはめ | Lev Krasnov, Dmitry Malikov, Marina Kiseleva, Sergei Tatarin, Sergey Sosnin, Stanislav Bezzubov. **BigSolDB 2.2**. Zenodo, doi:[10.5281/zenodo.22648301](https://doi.org/10.5281/zenodo.22648301)（論文 doi:[10.1038/s41597-025-05559-8](https://doi.org/10.1038/s41597-025-05559-8)） | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | 溶質 × 溶媒の平均、構造の正規化、両性・分子量 500 以上の有機溶媒の行と誤りの疑いのある行の除外、有機溶媒の値の溶液 1 L あたりへの換算 |
| 水の溶解度の実測（BigSolDB に値がない分子）、水のモデルとネットワークの教師の学習 | Murat Cihan Sorkun, Abhishek Khetan, Süleyman Er. **AqSolDB**. Harvard Dataverse, doi:[10.7910/DVN/OVHAW8](https://doi.org/10.7910/DVN/OVHAW8)（論文 doi:[10.1038/s41597-019-0151-1](https://doi.org/10.1038/s41597-019-0151-1)） | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | 塩・混合物・電荷を持つ物などの除外、同じ構造の平均 |
| 融点の照合の表、融点の推算モデル、ネットワークの教師の学習 | Jean-Claude Bradley, Antony Williams, Andrew Lang. **Jean-Claude Bradley Open Melting Point Dataset**. figshare, doi:[10.6084/m9.figshare.1031637](https://doi.org/10.6084/m9.figshare.1031637) | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | 出所が DrugBank・Aldrich Catalog-Handbook・Sigma-Aldrich の行を除外、同じ構造の中央値 |
| 融点（上の 2 つのどちらにも値がない学習の溶質 24 種） | Antony Williams, Daniel Lowe, Igor Tetko. **Melting Point and Pyrolysis Point Data for Tens of Thousands of Chemicals**. figshare, doi:[10.6084/m9.figshare.2007426](https://doi.org/10.6084/m9.figshare.2007426)（論文 doi:[10.1186/s13321-016-0113-y](https://doi.org/10.1186/s13321-016-0113-y)） | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | 学習の溶質と一致した値だけを使用 |
| 誤りの疑いのある文献の行の選び出し | V. Ramani ほか. **SC3: A Multi-Solvent Solubility Challenge and Benchmark**. Zenodo, doi:[10.5281/zenodo.23999757](https://doi.org/10.5281/zenodo.23999757) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | 文献どうしの一致の度合いを参考に、学習から外す行を選んだ（値そのものは重みに入れていない） |
| 混合溶媒の補正の係数、溶媒の沸点 | Dmitry Malikov, Lev Krasnov, Marina Kiseleva, Elizaveta Meshcheriakova, Fedor Kuznetsov, Vladimir Elistratov, Matvei Vasiyarov, Sergei Tatarin, Stanislav Bezzubov. **MixtureSolDB**. Zenodo, doi:[10.5281/zenodo.18660057](https://doi.org/10.5281/zenodo.18660057)（論文 doi:[10.1038/s41597-026-07047-z](https://doi.org/10.1038/s41597-026-07047-z)） | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | 25 °C 付近の系列から、溶媒の組ごとの係数を当てはめた |
| 溶媒の Abraham の係数 | Andrew Lang, Youngmin Lee. **Abraham Solvent Master Dataset（AbraLlama）**. figshare, doi:[10.6084/m9.figshare.25898476](https://doi.org/10.6084/m9.figshare.25898476) | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | 14 溶媒の行を使用 |
| 溶質の Abraham 記述子の推算式 | Jean-Claude Bradley, William Acree, Andrew Lang. **Compounds with known Abraham descriptors**. figshare, doi:[10.6084/m9.figshare.1176994](https://doi.org/10.6084/m9.figshare.1176994) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | 値の欠けた行などを除き、構造ごとに平均して一次式を当てはめた |
| 構造の読み込み、記述子と指紋の計算 | RDKit.js（@rdkit/rdkit 2025.3.4）。リポジトリの `vendor/rdkit/` | BSD 3-Clause（同梱の部品の表示は `vendor/rdkit/THIRD_PARTY_NOTICES.md`） | — |

このモデルは、非営利限定（NC）・改変禁止（ND）・継承（SA）のデータと、HSP（ハンセン溶解度パラメータ）の表を使っていません。元の作者・権利者がこのモデルを推奨していることを示すものではありません。

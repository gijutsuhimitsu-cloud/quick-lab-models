# 学習の土台の表（train/base_data/）

*English: processed public training tables (gzip CSV). Licensed CC BY 4.0 like the weights; see ../../ATTRIBUTION.md for sources.*

学習し直す時に、元のデータ（数十 MB）を取らなくても済むように、加工した表を同梱しています。ライセンスは重みと同じ **CC BY 4.0**（出典は [../../ATTRIBUTION.md](../../ATTRIBUTION.md)）。

| ファイル | 中身 | 出所 |
|---|---|---|
| `bigsol_25C.csv.gz` | 25 °C 付近（298.15 ± 5 K）の実測を、溶質 × 溶媒ごとに平均した表（11,488 組、溶質 1,400）。列：構造・溶媒・log10 mol/L・点の数・出典の論文の DOI | BigSolDB 2.2（CC BY 4.0） |
| `bigsol_solutes.csv.gz` | 溶質ごとの分子量・Crippen の log P・融点（`tmSrc`：sol = BigSolDB、mp = 融点の公開データ）・推算した Abraham 記述子（kE〜kB は表にある記述子） | BigSolDB 2.2、Abraham 記述子の表（CC BY 4.0）、融点の公開データ（CC0 1.0） |
| `aqsol_rows.csv.gz` | 水への溶解度（8,718 構造。同じ構造は平均。id は AqSolDB の番号） | AqSolDB（CC0 1.0） |
| `aqsol_solutes.csv.gz` | BigSolDB に無い溶質の入力（融点は融点の公開データから、3,267 溶質） | AqSolDB、融点の公開データ（CC0 1.0） |
| `tm_train.csv.gz` | 融点のモデルの学習用の表（in_bigsol・融点 °C・構造の記述子 97 列、19,287 構造） | 融点の公開データ（CC0 1.0。非営利の薬のデータベースと試薬会社のカタログの行を除く） |
| `mix_series_25C.csv.gz` | 混合溶媒の系列の途中の点（29,609 点）と、対数線形からのずれ | MixtureSolDB（CC BY 4.0）、密度は BigSolDB 2.2 |
| `abraham_solutes_clean.csv.gz` | Abraham 記述子の表を構造ごとに平均した物（1,835 構造） | Compounds with known Abraham descriptors（CC BY 4.0） |

- 構造は RDKit.js 2025.3.4 で正規化した SMILES（明示の水素を外した形）。
- 溶質の名前は、配布元の表から写していません（構造だけ）。
- 元のデータから作り直す道具は `../rebuild/`。

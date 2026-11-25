# 元のデータの置き場（train/raw/）

*English: put the original public datasets here only if you want to rebuild everything from scratch. This folder is git-ignored except for this file.*

ふだんの学習し直し（自社のデータを足すなど）には要りません。`train/base_data/` の表から学習できます。
全部を元のデータから作り直したい時だけ、下の形で置きます（場所は環境変数 `SOL_RAW` でも変えられます）。取り方とライセンスは `../README.md` の「元のデータから全部を作り直す」。

```
raw/
  bigsoldb/            BigSolDBv2.2.csv・BigSolDBv2.2_gamma.csv・BigSolDBv2.2_densities.csv
  mixturesoldb/        MixtureSolDB.csv・solvent_boiling_point.csv
  aqsoldb/             curated-solubility-dataset.csv（Harvard Dataverse から手で取る）
  bradley_mp/          BradleyMeltingPointDataset.xlsx
  abraham_descriptors/ 20140922AbrahamDescriptorsBackUp.xls
  abraham_solvents/    20240524 Abraham Solvent Master Dataset.xlsx
```

このフォルダの中身は `.gitignore` でリポジトリに入りません（この README を除く）。

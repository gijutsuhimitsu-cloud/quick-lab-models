# 自社のデータの例（train/examples/）

*English: `own_data_example.csv` shows the own-data CSV format. It is a real slice of public data (pure-solvent end points from MixtureSolDB, CC BY 4.0) for 11 solutes that are not in the public training tables, relabeled as "例-01" … "例-11" to look like in-house data.*

- `own_data_example.csv`：50 行・溶質 11 種・溶媒 19 種。形は [../OWN_DATA.md](../OWN_DATA.md)。
- 中身は**本物の公開データ**です：MixtureSolDB（Malikov ほか、doi:10.5281/zenodo.18660057、CC BY 4.0）の、純粋な溶媒での点（混合溶媒の系列の端）。学習の土台の表（BigSolDB 2.2 の 25 °C 付近・AqSolDB）に入っていない溶質だけを選びました。
- こちらで変えたこと：溶質の名前を「例-01」〜「例-11」に置き換えた（社内の番号に見立てた）。濃度（g/100 g 溶媒）を、純粋な溶媒の密度（BigSolDB 2.2）で mg/mL に換算し、一部の行は mol/L で書いた（単位の換算の例）。各行の `note` に、元の論文の DOI を残した。
- 47 行が 20〜30 °C、3 行が 40 °C（温度の外で飛ばされる例）。融点の列は空（書けば使われる）。
- 試し方：`node tests/check_own_data.cjs`、`node train/run_all.cjs --own train/examples/own_data_example.csv`。
- このファイルのライセンスは CC BY 4.0（元のデータに従う）。

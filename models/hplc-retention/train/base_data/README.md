# 基礎データ（公開の重みの学習に使った表）

*English summary: processed rows of the six open datasets used to train the published weights (no CC BY-SA, no non-commercial data). Each part keeps the license of its source (SMRT: Apache 2.0; HSM2, Weber, MCaRT, Wiczling group: CC BY 4.0; Zhang et al.: CC0 1.0). Keep `../../ATTRIBUTION.md` and `../../LICENSES/Apache-2.0.txt` when you redistribute these files. Formats: CSV / gzip CSV / JSON only (no pickle).*

公開の重み（`weights/`）は、この表だけから作りました。元のファイル（約 400 MB）を取り直さなくても学習し直せるように、学習に使う行を小さな表にして同梱しています。作り直しの手順は `../README.md` の「元のデータから作り直す」。

## ファイル

| ファイル | 中身 | 大きさ |
|---|---|---|
| `methods.csv` | 1 行 = 1 つの分析法（1,130）。出所・カラム（一般名と `weights/columns.json` の id）・カラムの性質（HSM の H〜D・logkEB）・固定相の種類・寸法・温度・流量・死時間・pH・アセトニトリルの割合・添加剤の印・勾配の表（JSON）か等組成の割合・補正を固定するか・注記 | 0.5 MB |
| `rows.csv.gz` | 1 行 = 1 つの測定（199,712）。`method_id`・`ck`（化合物の鍵＝InChIKey の前半 14 文字）・`rt_min`（勾配の保持時間、分）か `logk_iso`（等組成の log10 k）。**並びは学習の時の並び**（変えると木がわずかに変わる） | 1.9 MB |
| `compounds.csv.gz` | 1 行 = 1 つの化合物（80,988）。標準の SMILES と、`js/features.js` で計算した特徴量：記述子 41（`d_`、小数 4 桁）・官能基の数 56（`g_`）・イオン化の規則ごとの数 30（`ion_`）・512 ビットの指紋（`fp`、16 進。予測の幅の分割に使う） | 13.9 MB |
| `stock.csv` | ページのストックの化合物 102 種（こちらで作った名前 → 構造の表と、その特徴量） | 0.04 MB |
| `names.csv.gz` | 名前の索引の材料（出所の表の名前と構造のまま 845 行。こちらの表・Weber ら・Wiczling グループ・Zhang ら） | 0.01 MB |
| `names_manual.tsv`・`names_weber.tsv` | こちらで書いた名前 → SMILES の表（ストックと HSM2 の試験溶質、Weber らの化合物） | 小 |
| `names_exclude.tsv` | 名前と構造が合わないので、名前の索引に載せない名前（Wiczling グループの元の表の 4 件） | 小 |
| `retention_form.json` | 保持の式の定数（アセトニトリル: Weber らの等組成の実測に当てはめ、メタノール: Wiczling グループの勾配時間を変えた測定で決めた） | 小 |

## 出所ごとの行の数

| 出所 | ライセンス | 分析法 | 測定 | 化合物 | 中身 |
|---|---|---|---|---|---|
| METLIN SMRT（`smrt_figshare`） | Apache 2.0 | 1 | 79,955 | 79,888 | 勾配 1 条件の保持時間 |
| HSM2 database（`hsm2_mendeley`） | CC BY 4.0 | 551 | 8,816 | 15 | カラム 551 本 × 試験溶質 16 の log k（パラメータからの再構成。独立の実測ではない） |
| Weber ら（`weber_isocratic_mendeley`） | CC BY 4.0 | 126 | 5,996 | 105 | 等組成の log k（割合 × 温度。繰り返しは平均） |
| MCaRT（`mcart_msd_zenodo`） | CC BY 4.0 | 8 | 3,614 | 474 | 4 カラム × 2 緩衝液の勾配の保持時間 |
| Wiczling グループ（`wiczling_osf`） | CC BY 4.0 | 418 | 92,593 | 283 | 5 カラム × pH・溶媒・温度・勾配時間を変えた保持時間 |
| Zhang ら（`zhang2024_scidb`） | CC0 1.0 | 26 | 8,738 | 341 | 26 の勾配の分析法の保持時間（Dataset 10 は入れていない＝評価専用） |
| 合計 | | 1,130 | 199,712 | 80,988 | |

## こちらで変えたこと（CC BY 4.0・Apache 2.0 の「変更の表示」）

- 単位の変換（SMRT の秒 → 分、保持係数 → log10）、印の付いた値・候補のピークがあいまいな組・正でない保持係数の除外、同じ条件の繰り返しの平均（Weber）、カラム外の時間の補正（Wiczling）。
- 構造から特徴量を計算した（RDKit.js 2025.3.4 と `js/features.js`）。化合物は InChIKey の前半 14 文字でまとめた（立体異性体は同じ化合物）。
- 測定条件の整理：勾配の表の再平衡の区間を落とした、添加剤からの pH の見積もり、死時間と装置の遅れの見積もり（SMRT）。見積もりは `methods.csv` の `note` に書いた。
- **カラムの製品名とメーカー名は載せていません。** カラムは「固定相の種類＋HSM2 データベースの番号」の一般名（例 `C18 No.5`）と、`weights/columns.json` の id で表します。HSM2 の表に無いカラム（MCaRT の 4 本）は、同じ固定相の中央値の値です。
- SMRT の測定条件（カラムの種類・寸法・勾配・流量・ギ酸）は、SMRT の論文（doi:10.1038/s41467-019-13680-7）に書かれた条件を、RepoRT（CC BY-SA 4.0、https://github.com/michaelwitting/RepoRT ）の記録 0186 の記載を経由して読み、こちらの形で書いた個別の事実の値です。論文の本文との 1 項目ずつの照合は、まだしていません。死時間・装置の遅れ・pH はこちらの見積もりです。RepoRT の保持値・構造は、どの行にも入っていません（SMRT の保持値と構造は SMRT v1 から）。

## ライセンス

こちらで整えた部分（表の組み方・特徴量・こちらで作った表）は、重みと同じ CC BY 4.0 で配ります（リポジトリの直下の `LICENSE-DATA.md`）。元のデータには、元のライセンスがそのまま及びます（上の表。SMRT は Apache 2.0、CC BY の物は CC BY 4.0、CC0 の物は CC0）。表をまとめて、または抜き出して配る時は：

- 出典の表示（`../../ATTRIBUTION.md`）を付ける。
- `rows.csv.gz` の `method_id=smrt` の保持時間は、SMRT v1 の秒を分に直し、有効な構造を持つ行を選んだ物です。`compounds.csv.gz` は複数の出所の構造をまとめ、特徴量を足した物です。行の出所は、`rows.csv.gz` の `ck` と `methods.csv` の `method_id`・`source` をつなぐと分かります（同じ構造が複数の出所に出ることがあります）。この 2 つのファイル（や、その抜き出し）を配る時は、Apache License 2.0 の文（`../../LICENSES/Apache-2.0.txt`）と、上の「こちらで変えたこと」を一緒に付けてください（SMRT の配布物に NOTICE のファイルはありません）。

これは事実の整理で、法律の意見ではありません。

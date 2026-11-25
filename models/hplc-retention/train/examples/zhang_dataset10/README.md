# 例：自社のデータの形（Zhang らの Dataset 10 を書き直した物）

*English summary: an example "own data" set in the format of `../../OWN_DATA.md`, made from a real slice of public CC0 data (Zhang et al., Science Data Bank doi:10.57760/sciencedb.15823, Dataset 10: one methanol gradient method, 335 measurements). It is NOT anyone's in-house data. Dataset 10 is not in the base training data (it was kept for testing), so it mimics "adding one new method of your own".*

**本当の自社のデータではありません。** 公開データ（CC0 1.0）の一部を、自社のデータの形（`methods.csv`・`measurements.csv`）に書き直した例です。

- 出所：Zhang, Y. ほか, “Retention time dataset for heterogeneous molecules in reversed-phase liquid chromatography”, Science Data Bank, doi:10.57760/sciencedb.15823（論文 doi:10.1038/s41597-024-03780-5）。CC0 1.0（表示の義務はありませんが、出所として書いています）。
- 中身：Dataset 10 の 1 つの勾配の分析法（C18 No.93 の値、100 × 2.1 mm、2.2 µm、30 °C、0.2 mL/min、メタノール／水、ギ酸＋ギ酸アンモニウム）と、335 化合物の保持時間。
- こちらで変えたこと：移動相 A に入っている 10% のメタノールを勾配の表に入れた（`0:0.1;2:0.1;24:1;27:1`）、pH 3.7 は添加剤からの見積もり、カラムは製品名ではなく `weights/columns.json` の一般名（同じ固定相の粒径違いの値）。
- Dataset 10 は、公開の重みの学習に入れていない（評価専用にした）ファイルです。なので、この例は「自社の分析法を 1 つ足す」状況に近くなります。化合物は 333 個すべてが、Zhang らのほかの分析法で基礎データに入っています。

## 動かし方

`models/hplc-retention/` で：

```bash
python train/scripts/eval_external.py train/examples/zhang_dataset10            # 学習し直さずに、公開の重みで評価
python train/run.py --own train/examples/zhang_dataset10 --name example         # 基礎データ＋この例で学習（約 10 分）
python train/run.py --own train/examples/zhang_dataset10 --name example_fast --trees 60 --qtrees 30   # 形を確かめるだけの速い試し
```

## この道具を作った PC での結果（木 900 本）

保持時間の誤差（平均の絶対誤差。かっこの中は分析時間 27 分に対する割合）：

| | 公開の重み（補正なし） | 例を足して学習した重み（補正あり） |
|---|---|---|
| 評価用の化合物 32 個（どちらの学習にも使っていない） | 2.08 分（7.7%）、順位の相関 0.945 | 1.83 分（6.8%）、0.948 |
| 全部の 335 化合物（学習に使った化合物を含む。良く見えすぎる） | 1.22 分（4.5%）、0.990 | 0.90 分（3.3%）、0.989 |

- 良くなった分の多くは、分析法ごとの補正（死時間の倍率・装置の遅れ）を、この分析法の測定から決められたことによります。この例では補正が格子の端（死時間 × 0.6、遅れ 2 分）になりました。装置の遅れが分かっている時は、`methods.csv` の `dwell_min` に書くと、その値で固定できます。
- 学習し直さずに、同じカラムの基礎データの分析法の補正の中央値（死時間 × 0.8、遅れ 1 分）を使うだけでも、全部の化合物で 4.5% → 3.8% になります（`eval_external.py --cal train/work/free/calibration.json`）。
- 評価用の化合物が 32 個と少ないので、数字のぶれは大きい（目安）です。

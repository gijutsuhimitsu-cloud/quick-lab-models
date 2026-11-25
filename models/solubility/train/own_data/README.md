# 自社のデータの置き場（train/own_data/）

*English: put your own solubility CSV files here. Everything in this folder except this README is git-ignored, so your data never leaves your machine through this repository.*

- ここに置いた CSV は、リポジトリの `.gitignore` でコミットされません（この README だけが入ります）。
- 形は [../OWN_DATA.md](../OWN_DATA.md)。例は [../examples/own_data_example.csv](../examples/own_data_example.csv)。
- 学習は手元の計算機の中だけで行います。データをどこかへ送る処理はありません。
- 学習の途中の物（`train/work/`）とできた重み（`train/out/`）にも、自社のデータから作った値が入ります。これらも `.gitignore` でコミットされません。社外に出す時は、中身を確かめてから出してください。

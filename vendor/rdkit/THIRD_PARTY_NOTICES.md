# RDKit.js と、その中の部品のライセンス

*English: `RDKit_minimal.js` / `RDKit_minimal.wasm` (@rdkit/rdkit 2025.3.4-1.0.0) are the RDKit core (BSD 3-Clause, `LICENSE`) compiled with Emscripten. The WebAssembly file also contains third-party code and fonts. Their notices and license texts are listed below and kept in `LICENSES/`. If you ship these two files (for example inside a product), ship this file and `LICENSES/` with them.*

`RDKit_minimal.js` と `RDKit_minimal.wasm` は、RDKit の本体を Emscripten で WebAssembly にした物です。WebAssembly の中には、RDKit のほかに、フォントや第三者のライブラリも入っています。この 2 つのファイルを配る時（製品に組み込む時も）は、このファイルと `LICENSES/` を一緒に付けてください。

## 部品の一覧

| 部品 | 中での役目 | 著作権の表示 | ライセンス | 原文 |
|---|---|---|---|---|
| RDKit（本体） | 構造の読み込み・記述子の計算 | Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti and others | BSD 3-Clause | [LICENSE](LICENSE) |
| Telex-Regular（フォント） | 構造の絵を描く時の既定の文字 | Copyright 2011 The Telex Project Authors (info@andrestorresi.com.ar), with Reserved Font Name Telex. | SIL Open Font License 1.1 | [LICENSES/Telex-OFL-1.1.txt](LICENSES/Telex-OFL-1.1.txt) |
| Roboto-Regular（フォント） | 構造の絵を描く時の文字（選べる物） | Copyright 2011 Google Inc. All Rights Reserved. | Apache License 2.0 | [LICENSES/Roboto-NOTICE.txt](LICENSES/Roboto-NOTICE.txt)・[LICENSES/Apache-2.0.txt](LICENSES/Apache-2.0.txt) |
| InChI API Library 1.07.3 | InChI・InChIKey の計算 | Copyright (c) 2024 InChI Project | MIT | [LICENSES/InChI-MIT.txt](LICENSES/InChI-MIT.txt) |
| CoordGen（coordgenlibs 3.0.2） | 2 次元の座標の生成 | Copyright (c) 2017, Schrödinger, Inc. | BSD 3-Clause | [LICENSES/CoordGen-BSD-3-Clause.txt](LICENSES/CoordGen-BSD-3-Clause.txt) |
| FreeType | フォントの読み込み（構造の絵） | Portions of this software are copyright © The FreeType Project (https://freetype.org). All rights reserved. | FreeType License（FTL。FreeType は FTL と GPL v2 の選択で、ここでは FTL） | [LICENSES/FreeType-FTL.txt](LICENSES/FreeType-FTL.txt) |
| Emscripten の実行時のコード | `RDKit_minimal.js` の読み込みの部分と、WebAssembly の中の実行時の部品 | Copyright (c) 2010-2014 Emscripten authors | MIT と University of Illinois/NCSA Open Source License の 2 つ | [LICENSES/Emscripten-MIT-NCSA.txt](LICENSES/Emscripten-MIT-NCSA.txt) |
| musl libc（Emscripten の C の標準ライブラリ） | メモリ・文字列などの基本の関数 | Copyright © 2005-2020 Rich Felker, et al. | MIT | [LICENSES/musl-MIT.txt](LICENSES/musl-MIT.txt) |
| LLVM libc++・libc++abi（Emscripten の C++ の標準ライブラリ） | C++ の基本の部品・例外の処理 | LLVM Project | Apache License 2.0 with LLVM Exceptions | https://llvm.org/LICENSE.txt |
| Boost（ヘッダーのライブラリ。property_tree など） | JSON の読み書きなど | 各ファイルの作者 | Boost Software License 1.0 | https://www.boost.org/LICENSE_1_0.txt |

- **フォント**：OFL は、フォントそのものに付く条件です。フォントを使うプログラムを OFL にする条件ではありません。フォントだけを取り出して売ることはできません。
- **LLVM の libc++・libc++abi**：LLVM Exceptions により、コンパイルされてバイナリに埋め込まれた部分は、Apache 2.0 の第 4 条 (a)(b)(d)（表示の添付など）を満たさずに配れます。ここではリンクだけを書いています。
- **Boost**：Boost Software License は、コンパイルされた実行形式だけを配る時は、ライセンスの文を付けなくてよいとしています。ここではリンクだけを書いています。

## 原文の取り先（2026-11-25 に取得）

| ファイル | 取り先（版） |
|---|---|
| `LICENSE` | rdkit/rdkit の `license.txt`（Release_2025_03_4）と同じ文 |
| `LICENSES/Telex-OFL-1.1.txt` | rdkit/rdkit の `Data/Fonts/telex_font_license.txt`（Release_2025_03_4） |
| `LICENSES/Roboto-NOTICE.txt` | rdkit/rdkit の `Data/Fonts/roboto_regular_license.txt`（Release_2025_03_4）。中身は「Apache 2.0 で使える」という 3 行のメモ |
| `LICENSES/Apache-2.0.txt` | Apache License 2.0 の全文（https://www.apache.org/licenses/LICENSE-2.0.txt と同じ文） |
| `LICENSES/InChI-MIT.txt` | IUPAC-InChI/InChI の `LICENSE`（v1.07.3） |
| `LICENSES/CoordGen-BSD-3-Clause.txt` | schrodinger/coordgenlibs の `LICENSE`（v3.0.2。RDKit の Release_2025_03_4 の `External/CoordGen/CMakeLists.txt` が取りに行く版） |
| `LICENSES/FreeType-FTL.txt` | freetype/freetype の `docs/FTL.TXT`（master） |
| `LICENSES/Emscripten-MIT-NCSA.txt` | emscripten-core/emscripten の `LICENSE`（main） |
| `LICENSES/musl-MIT.txt` | emscripten-core/emscripten の `system/lib/libc/musl/COPYRIGHT`（main） |

## 確かめたことと、確かめられなかったこと

- 中に何が入っているかは、`RDKit_minimal.wasm` の中の文字（フォントの名前の表、エラーの文、ビルドの時のファイルの場所など）から確かめました。Telex と Roboto のフォント、InChI、CoordGen（sketcherMinimizer）、FreeType、Boost、Emscripten の libc++abi の印があります。maeparser など、ほかの部品の印は見つかりませんでした。
- 確かめたファイル（SHA-256）：`RDKit_minimal.js` は `58d3c996ade7e0b0137d4f9363ff6c204b545689428fa0898eeaed579ef788d9`、`RDKit_minimal.wasm` は `e0967d44fed59e44a2d07bfc08f3a86821c54a14c95bb7e6e392fe867333b8c8`。
- この RDKit.js のビルドに使われた FreeType と Emscripten の正確な版は、ファイルの中からは分かりませんでした。そのため、この 2 つの原文は配布元の今の版から取りました（どちらもライセンスの文は長く変わっていません）。
- RDKit.js を新しい版に入れ替えた時は、この表を見直してください。

これは事実の整理で、法律の意見ではありません。

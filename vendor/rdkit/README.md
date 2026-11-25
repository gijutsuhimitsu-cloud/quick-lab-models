# RDKit.js（同梱）

- 版：@rdkit/rdkit 2025.3.4-1.0.0（`RDKit_minimal.js`・`RDKit_minimal.wasm`）
- ライセンス：RDKit の本体は BSD 3-Clause（`LICENSE`）。RDKit の著作権表示はそのまま残してください。
- 中の部品：WebAssembly の中には、フォント（Telex・Roboto）、InChI、CoordGen、FreeType、Boost、Emscripten の実行時の部品も入っています。それぞれの表示とライセンスの原文は `THIRD_PARTY_NOTICES.md` と `LICENSES/` にあります。
- 配る時（製品に組み込む時も）：`RDKit_minimal.js`・`RDKit_minimal.wasm` と一緒に、`LICENSE`・`THIRD_PARTY_NOTICES.md`・`LICENSES/` を付けてください。
- 配布元：https://github.com/rdkit/rdkit-js ・ https://www.npmjs.com/package/@rdkit/rdkit
- このリポジトリでは、構造（SMILES）を読んで記述子を計算するためだけに使います。どのモデルもこの 1 つを共有します。

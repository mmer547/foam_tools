# foam_cylindrical_slice_unfold

OpenFOAM の**体積メッシュ**と半径 `r` の**円筒面**（`vtkCylinder`）を **`vtkCutter`** で交差させ、得られた**表面（PolyData）**を**周方向に展開**し、平面 `(u, v)` 上に表示します。

- **u** = `r * θ`（周方向の弧長。θ はラジアン）
- **v** = 軸方向座標（軸が `z` なら `z - cz`。`--axis-origin` で軸上の基準点を指定）

θ は主に `[-π, π]` なので、全周展開では **θ = ±π で接縫**が出ることがあります。

## 依存

```bash
pip install pyvista numpy matplotlib
```

（`vtk` は `pyvista` に同梱されます。）

## 使い方

```bash
python foam_cylindrical_slice_unfold.py path/to/case --radius 0.05
python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --scalars p --vtk unfolded.vtp
python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --triangles
```

- `--radius` / `-r`: 円筒半径（カッターに使用。展開の `u = r*θ` にも使用）
- `--triangles`: カッター出力を三角形化
- `--axis`: 円筒軸 `x` / `y` / `z`（既定: `z`）
- `--scalars`: コンタ用の場（未指定時は `U` から `|U|`、なければ先頭のスカラー）
- `--vtk`: 展開後メッシュを保存（**`.vtp`** など PolyData 向け。拡張子が合わない場合は自動で置換）

`--help` でオプション一覧を表示します。

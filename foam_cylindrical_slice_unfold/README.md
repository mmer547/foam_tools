# foam_cylindrical_slice_unfold

OpenFOAM の**体積メッシュ**と半径 `r` の**円筒面**（`vtkCylinder`）を **`vtkCutter`** で交差させ、得られた**表面（PolyData）**を**周方向に展開**し、平面 `(u, v)` 上に表示します。

- **u** = `r * θ`（周方向の弧長。θ はラジアン）
- **v** = 軸方向座標（軸が `z` なら `z - cz`。`--axis-origin` で軸上の基準点を指定）

θ は主に `[-π, π]` です。全周展開では **θ = ±π の接縫**をまたぐセルをそのまま写すと、画面を横断する引き伸ばし面になりコンターが非連続に見えます。本ツールは該当セルを**除去**してから展開します（接縫付近に細い隙間が残ることがあります）。コンター用スカラーが cell データのときは point へ平均化し、滑らかに塗ります。

### SurfaceLIC 用の速度ベクトル

ケースに速度 **`U`** があるとき、展開**前**の円筒接線基底で

- 第1成分 = `U` の**周方向**（θ 増加方向）への射影  
- 第2成分 = `U` の**軸方向**への射影  
- 第3成分 = `0`  

となる **point ベクトル**を既定名 **`U_unfold`** で書き込みます。メッシュは展開後 **XY 平面**（`z=0`）上にあります。

- **`--surface-lic`**: VTK の `vtkSurfaceLICMapper` で **SurfaceLIC 画像**を `-o` に保存（PyVista はカメラ・スクリーンショット補助。要 `U_unfold`）
- **`--no-u-unfold`**: `U_unfold` を追加しない（`--surface-lic` と併用不可）
- **`--u-unfold-name NAME`**: 配列名を変更（既定: `U_unfold`）

`--vtk` でメッシュを保存し、ParaView 側で SurfaceLIC を掛ける使い方もできます。

`U` が **cell データ**だけの場合は、内部で **point へ平均化**（`cell_data_to_point_data`）してから変換します。

## 依存

```bash
pip install pyvista numpy matplotlib
```

（`vtk` は `pyvista` に同梱されます。`--surface-lic` も同梱の VTK を使います。）

## 使い方

```bash
python foam_cylindrical_slice_unfold.py path/to/case --radius 0.05
python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --scalars p --vtk unfolded.vtp
python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --triangles
python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --surface-lic -o lic.png
```

- `--radius` / `-r`: 円筒半径（カッターに使用。展開の `u = r*θ` にも使用）
- `--triangles`: カッター出力を三角形化
- `--axis`: 円筒軸 `x` / `y` / `z`（既定: `z`）
- `--scalars`: コンタ用の場（未指定時は `U` から `|U|`、なければ先頭のスカラー）
- `--vtk`: 展開後メッシュを保存（**`.vtp`** など PolyData 向け。拡張子が合わない場合は自動で置換）
- `--surface-lic`: VTK SurfaceLIC で画像出力
- `--no-u-unfold` / `--u-unfold-name`: SurfaceLIC 用ベクトル（既定 `U_unfold`）

`--help` でオプション一覧を表示します。

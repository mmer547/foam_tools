# foam_velocity_slice_position

OpenFOAM の結果を読み込み、**任意の点を通る座標平面**でスライスした流速（`|U|`）のコンタ図を画像として保存します。平面の法線は `--normal`、通る点は `-p` / `--origin` で指定します。

## 依存

```bash
pip install pyvista numpy matplotlib
```

## 使い方

スクリプトがあるディレクトリで実行する例です。

```bash
python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z
python foam_velocity_slice_position.py path/to/case --origin 1.0 0 0 --normal x -o slice.png
python foam_velocity_slice_position.py path/to/case -p 0 2 0 --normal y --time 500
python foam_velocity_slice_position.py path/to/case --origin "(1.0, 0, 0)" --normal x
python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z --scalars Ux
python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z --surface-lic -o lic.png
```

- `-p` / `--origin`: `OX OY OZ` の3語、または `"(OX, OY, OZ)"` / `OX,OY,OZ` の1語。未指定時は原点 `(0,0,0)`。
- `--scalars`: コンタスカラー。`U_mag`（既定）/ `Ux` / `Uy` / `Uz`（別名: `mag`, `U_x`, `U_y`, `U_z`）
- `--surface-lic`: VTK の `vtkSurfaceLICMapper` で SurfaceLIC 画像を `-o` に保存（面内速度と `--scalars` のコンタ）
- ケース直下に空の `case.foam` があると読みやすいです（無くてもディレクトリを渡せる場合があります）。

`--help` でオプション一覧を表示します。

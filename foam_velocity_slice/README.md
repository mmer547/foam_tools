# foam_velocity_slice

OpenFOAM の結果を読み込み、**座標平面**（既定では原点を通る面）でスライスした流速（`|U|`）のコンタ図を画像として保存します。法線方向は `--normal`（例: `z` なら z=0 付近の XY 的な断面）で選べます。

## 依存

```bash
pip install pyvista numpy matplotlib
```

## 使い方

```bash
python foam_velocity_slice.py path/to/case
python foam_velocity_slice.py path/to/case -o u_mag_z0.png --normal z
python foam_velocity_slice.py path/to/case --time 500 --normal x
python foam_velocity_slice.py path/to/case --origin "(0, 0, 0.5)" --normal z
```

- `--origin`: 面が通る点（既定: `0 0 0`）。3語または `"(OX, OY, OZ)"` 形式でも指定可能です。

任意位置でのスライスを主に使う場合は、リポジトリ内の `foam_velocity_slice_position` を参照してください。

`--help` でオプション一覧を表示します。

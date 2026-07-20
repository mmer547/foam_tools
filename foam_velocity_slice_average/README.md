# foam_velocity_slice_average

OpenFOAM の結果を読み込み、**領域の下端**から指定ピッチでスライス断面を取得し、各断面内の流速 3 成分（Ux, Uy, Uz）に加え、`--normal` で指定した軸まわりの**半径方向流速（Ur）**と**回転方向流速（Utheta）**の面積加重平均を CSV に出力します。

スライス面は `--normal` で選んだ軸に垂直な平面です（既定: `z` → XY 面を z 方向に等間隔で取得）。下端はメッシュの最小座標（`--start` で上書き可）です。

## 依存

```bash
pip install pyvista numpy matplotlib
```

## 使い方

```bash
python foam_velocity_slice_average.py path/to/case --pitch 0.1
python foam_velocity_slice_average.py path/to/case --pitch 0.05 --normal z -o u_avg.csv
python foam_velocity_slice_average.py path/to/case --pitch 0.1 --time 500
python foam_velocity_slice_average.py path/to/case --pitch 0.1 --offset 0.05
python foam_velocity_slice_average.py path/to/case --pitch 0.1 --start 0 --end 2.0
python foam_velocity_slice_average.py path/to/case --pitch 0.1 --plot
python foam_velocity_slice_average.py path/to/case --pitch 0.1 --plot profile.png
```

`--plot` を付けると、取得した断面平均流速をプロファイル図として PNG 等で保存します（縦軸=座標、横軸=Ux/Uy/Uz/Ur/Utheta）。

Ur / Utheta はメッシュ中心を通る `--normal` 軸を基準に、各セル位置で速度を分解してから面積加重平均した値です。回転方向は右ねじの法線が +軸方向となる向き（例: `--normal z` なら +z から見て反時計回りが正）です。

## 出力 CSV

1 行目はヘッダ、2 行目以降が各断面の結果です。

| 列 | 内容 |
|----|------|
| `x` / `y` / `z` | 断面位置（`--normal` に応じた軸座標） |
| `Ux` | 断面内 Ux の面積加重平均 |
| `Uy` | 断面内 Uy の面積加重平均 |
| `Uz` | 断面内 Uz の面積加重平均 |
| `Ur` | 断面内の半径方向流速の面積加重平均 |
| `Utheta` | 断面内の回転方向流速の面積加重平均 |

断面がメッシュと交差しない場合は `NaN` を出力します。

## 主なオプション

| オプション | 説明 |
|------------|------|
| `--pitch` | 断面間隔（必須） |
| `--normal` | 法線方向 `x` / `y` / `z`（既定: `z`） |
| `--offset` | 下端から最初の断面までのオフセット（既定: `0`） |
| `--start` / `--end` | スライス範囲（未指定時はメッシュ境界） |
| `--time` | 時刻（数値 / `first` / `latest`、既定: `latest`） |
| `-o` / `--output` | 出力 CSV パス |
| `--plot [PLOT]` | 流速プロファイル図を出力（パス省略時: `foam_velocity_slice_average.png`） |
| `--dpi` | `--plot` 出力の DPI（既定: `150`） |

`--help` でオプション一覧を表示します。

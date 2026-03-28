# foam_log_plot

OpenFOAM の **foamLog** が生成した `logs/` 内の残差データを読み、1 枚の画像（PNG など）として残差グラフを出力します。

描画対象は **Ux, Uy, Uz, k, epsilon, omega** のみです（該当ファイルが無ければその系列は省略）。`*FinalRes_*` のような Final 残差ファイルは描画しません。

## 依存

```bash
pip install matplotlib numpy
```

## 使い方

```bash
python foam_log_plot.py ./logs
python foam_log_plot.py . -o residuals.png
```

foamLog をケースのカレントで実行した場合は `logs` ができるので、そのパスを渡します。

`--help` でオプション一覧を表示します。

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenFOAM の foamLog が生成した logs/ ディレクトリ内の残差データを読み込み、
1 枚の画像（PNG など）として残差グラフを出力する。
描画するのは Ux, Uy, Uz, k, epsilon, omega のみ（該当ファイルが無ければその系列は省略）。
Final 残差（foamLog の *FinalRes_* ファイル、例: pFinalRes_0）は描画しない。

依存: pip install matplotlib numpy

使用例:
  python foam_log_plot.py ./logs
  python foam_log_plot.py . -o residuals.png
  foamLog をカレントで実行した場合は logs ができるので:
  python foam_log_plot.py logs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# foamLog のファイル名は「変数名_連番」（例: Ux_0, epsilon_0）
PLOT_VARIABLES: frozenset[str] = frozenset(
    {"Ux", "Uy", "Uz", "k", "epsilon", "omega"}
)
PLOT_VARIABLE_ORDER: tuple[str, ...] = (
    "Ux",
    "Uy",
    "Uz",
    "k",
    "epsilon",
    "omega",
)


def _variable_from_stem(stem: str) -> str | None:
    """ファイル名（拡張子なし）から foamLog の変数名を取り出す。"""
    if "_" not in stem:
        return None
    base, suffix = stem.rsplit("_", 1)
    if not suffix.isdigit():
        return None
    return base


def _parse_columns(path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """2 列以上の数値行を (x, residual) として読む。コメント行・空行をスキップ。"""
    xs: list[float] = []
    ys: list[float] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        print(f"警告: 読み込み失敗 {path}: {e}", file=sys.stderr)
        return None

    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            x = float(parts[0])
            y = float(parts[1])
        except ValueError:
            continue
        if y <= 0:
            # 半対数プロットのため正の残差のみ
            continue
        xs.append(x)
        ys.append(y)

    if not xs:
        return None
    return np.asarray(xs), np.asarray(ys)


def _is_final_residual_file(name: str) -> bool:
    """foamLog が付ける Final 残差ファイル名（変数名 + FinalRes + _番号）を除外する。"""
    return "finalres" in name.lower()


def _collect_series(logs_dir: Path) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """logs_dir 内から Ux/Uy/Uz/k/epsilon/omega の系列のみ収集。ラベルはファイル名（拡張子なし）。"""
    if not logs_dir.is_dir():
        raise NotADirectoryError(f"ディレクトリではありません: {logs_dir}")

    series: list[tuple[str, np.ndarray, np.ndarray]] = []
    for p in sorted(logs_dir.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_file():
            continue
        if _is_final_residual_file(p.name):
            continue
        var = _variable_from_stem(p.stem)
        if var is None or var not in PLOT_VARIABLES:
            continue
        parsed = _parse_columns(p)
        if parsed is None:
            continue
        x, y = parsed
        label = p.stem
        series.append((label, x, y))

    order_rank = {name: i for i, name in enumerate(PLOT_VARIABLE_ORDER)}

    def _sort_key(item: tuple[str, np.ndarray, np.ndarray]) -> tuple[int, str]:
        stem = item[0]
        v = _variable_from_stem(stem) or ""
        return (order_rank.get(v, 99), stem.lower())

    series.sort(key=_sort_key)
    return series


def _resolve_logs_dir(path: Path) -> Path:
    """path が logs ならそのまま、親で logs があればそれを優先。"""
    if path.is_dir() and path.name == "logs":
        return path
    if path.is_dir() and (path / "logs").is_dir():
        return path / "logs"
    if path.is_dir():
        return path
    raise FileNotFoundError(f"logs ディレクトリが見つかりません: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="foamLog の logs/ から残差を読み、1 枚のグラフ画像を保存します。"
    )
    parser.add_argument(
        "path",
        type=Path,
        help="logs ディレクトリ、またはその親ディレクトリ（./ で foamLog 実行した場所など）",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("foam_residuals.png"),
        help="出力画像パス（既定: foam_residuals.png）",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="画像の DPI（既定: 150）",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="OpenFOAM residuals (foamLog)",
        help="グラフタイトル",
    )
    parser.add_argument(
        "--no-grid",
        action="store_true",
        help="グリッドを表示しない",
    )
    args = parser.parse_args()

    try:
        logs_dir = _resolve_logs_dir(args.path.resolve())
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        return 1

    series = _collect_series(logs_dir)
    if not series:
        print(
            f"エラー: {logs_dir} から有効な残差データを読み込めませんでした。",
            file=sys.stderr,
        )
        return 1

    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")

    for label, x, y in series:
        ax.semilogy(x, y, label=label, linewidth=1.0)

    ax.set_xlabel("Time / iteration")
    ax.set_ylabel("Residual")
    ax.set_title(args.title)
    if not args.no_grid:
        ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(loc="best", fontsize=8, ncol=2)

    out = args.output.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    print(f"保存しました: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

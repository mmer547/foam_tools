#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenFOAM の結果を読み込み、領域の下端から指定ピッチでスライス断面を取得し、
各断面内の流速 3 成分（Ux, Uy, Uz）に加え、--normal で指定した軸まわりの
半径方向流速（Ur）・回転方向流速（Utheta）の面積加重平均を CSV に出力する。

スライス面は --normal で選んだ軸に垂直な平面とし、下端（その軸の最小座標）から
--pitch 間隔で上端まで断面を取る。

依存: pip install pyvista numpy matplotlib

使用例:
  python foam_velocity_slice_average.py path/to/case --pitch 0.1
  python foam_velocity_slice_average.py path/to/case --pitch 0.05 --normal z -o u_avg.csv
  python foam_velocity_slice_average.py path/to/case --pitch 0.1 --time 500 --offset 0.05
  python foam_velocity_slice_average.py path/to/case --pitch 0.1 --start 0 --end 2.0
  python foam_velocity_slice_average.py path/to/case --pitch 0.1 --plot
  python foam_velocity_slice_average.py path/to/case --pitch 0.1 --plot profile.png
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np


def _resolve_openfoam_path(case_path: Path) -> Path:
    """OpenFOAMReader に渡すパス（.foam ファイル優先、無ければケースディレクトリ）。"""
    p = case_path.resolve()
    if p.is_file():
        return p
    if not p.is_dir():
        raise FileNotFoundError(f"パスが見つかりません: {case_path}")
    foams = sorted(p.glob("*.foam"))
    if foams:
        return foams[0]
    ctrl = p / "system" / "controlDict"
    if not ctrl.is_file():
        raise FileNotFoundError(
            f"OpenFOAM ケースではなさそうです（system/controlDict がありません）: {p}"
        )
    return p


def _normal_axis(name: str) -> tuple[float, float, float]:
    m = name.strip().lower()
    if m == "x":
        return (1.0, 0.0, 0.0)
    if m == "y":
        return (0.0, 1.0, 0.0)
    if m == "z":
        return (0.0, 0.0, 1.0)
    raise ValueError(f"normal は x, y, z のいずれか: {name!r}")


def _axis_index(name: str) -> int:
    m = name.strip().lower()
    if m == "x":
        return 0
    if m == "y":
        return 1
    if m == "z":
        return 2
    raise ValueError(f"normal は x, y, z のいずれか: {name!r}")


def _extract_internal_mesh(multiblock) -> "object":
    import pyvista as pv

    if not isinstance(multiblock, pv.MultiBlock):
        return multiblock
    if "internalMesh" in multiblock.keys():
        m = multiblock["internalMesh"]
        if m is not None:
            return m
    for i in range(multiblock.n_blocks):
        b = multiblock[i]
        if b is not None and b.n_cells > 0:
            return b
    raise RuntimeError("internalMesh を取得できませんでした。")


def _pick_time_index(reader, time_arg: str | None) -> float:
    times = list(reader.time_values)
    if not times:
        return 0.0
    if time_arg is None or time_arg.lower() == "latest":
        return float(times[-1])
    if time_arg.lower() == "first":
        return float(times[0])
    try:
        t = float(time_arg)
    except ValueError:
        raise SystemExit(f"無効な --time: {time_arg!r}（数値 / first / latest）")
    nearest = min(times, key=lambda x: abs(float(x) - t))
    return float(nearest)


def _axis_bounds(mesh, axis_idx: int) -> tuple[float, float]:
    bounds = mesh.bounds
    return float(bounds[2 * axis_idx]), float(bounds[2 * axis_idx + 1])


def _slice_positions(
    lo: float,
    hi: float,
    pitch: float,
    offset: float,
) -> list[float]:
    if pitch <= 0:
        raise ValueError(f"pitch は正の数である必要があります: {pitch}")
    start = lo + offset
    if start > hi:
        return []
    tol = max(1e-9, abs(hi) * 1e-9)
    count = int(np.floor((hi - start) / pitch + tol)) + 1
    positions = [start + i * pitch for i in range(count)]
    return [p for p in positions if p <= hi + tol]


def _slice_origin(mesh, axis_idx: int, position: float) -> tuple[float, float, float]:
    origin = list(mesh.center)
    origin[axis_idx] = position
    return (float(origin[0]), float(origin[1]), float(origin[2]))


def _axis_center(mesh) -> tuple[float, float, float]:
    c = mesh.center
    return (float(c[0]), float(c[1]), float(c[2]))


def _cell_cylindrical_velocity(
    positions: np.ndarray,
    velocities: np.ndarray,
    axis_idx: int,
    axis_center: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """
    軸に垂直な平面内で、各点の速度を半径方向 Ur と回転方向 Utheta に分解する。
    回転方向は右ねじの法線が +軸方向となる向き（+z 軸まわりなら反時計回り）。
    """
    perp = [i for i in range(3) if i != axis_idx]
    i0, i1 = perp[0], perp[1]
    centers = axis_center

    d0 = positions[:, i0] - centers[i0]
    d1 = positions[:, i1] - centers[i1]
    r = np.hypot(d0, d1)

    u0 = velocities[:, i0]
    u1 = velocities[:, i1]

    eps = 1e-12
    safe_r = np.where(r > eps, r, 1.0)
    ur = (u0 * d0 + u1 * d1) / safe_r
    utheta = (-u0 * d1 + u1 * d0) / safe_r
    on_axis = r <= eps
    ur = np.where(on_axis, 0.0, ur)
    utheta = np.where(on_axis, 0.0, utheta)
    return ur, utheta


def _slice_area_weighted_velocity_mean(
    slice_mesh,
    axis_idx: int,
    axis_center: tuple[float, float, float],
) -> tuple[float, float, float, float, float]:
    """スライス面上の Ux, Uy, Uz, Ur, Utheta の面積加重平均を返す。"""
    nan5 = (float("nan"),) * 5
    if slice_mesh.n_cells == 0:
        return nan5

    work = slice_mesh
    if "U" in work.point_data and "U" not in work.cell_data:
        work = work.point_data_to_cell_data()
    elif "U" not in work.cell_data:
        keys = list(work.point_data.keys()) + list(work.cell_data.keys())
        raise RuntimeError(
            "速度場 U がスライス上に見つかりません。"
            f"利用可能なデータ: {keys[:20]}{'...' if len(keys) > 20 else ''}"
        )

    work = work.compute_cell_sizes(area=True, length=False, volume=False)
    areas = np.asarray(work.cell_data["Area"], dtype=float)
    u = np.asarray(work.cell_data["U"], dtype=float)[:, :3]
    centers = np.asarray(work.cell_centers().points, dtype=float)

    valid = areas > 0
    if not np.any(valid):
        return nan5

    areas = areas[valid]
    u = u[valid]
    centers = centers[valid]
    ur, utheta = _cell_cylindrical_velocity(centers, u, axis_idx, axis_center)

    total_area = float(areas.sum())
    return (
        float(np.dot(u[:, 0], areas) / total_area),
        float(np.dot(u[:, 1], areas) / total_area),
        float(np.dot(u[:, 2], areas) / total_area),
        float(np.dot(ur, areas) / total_area),
        float(np.dot(utheta, areas) / total_area),
    )


def _write_csv(
    path: Path,
    axis_name: str,
    rows: list[tuple[float, float, float, float, float, float]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([axis_name, "Ux", "Uy", "Uz", "Ur", "Utheta"])
        for position, ux, uy, uz, ur, utheta in rows:
            writer.writerow([position, ux, uy, uz, ur, utheta])


def _save_profile_plot(
    path: Path,
    axis_name: str,
    rows: list[tuple[float, float, float, float, float, float]],
    *,
    time: float,
    dpi: int,
) -> None:
    """縦軸=座標、横軸=流速成分のプロファイル図を保存する。"""
    import matplotlib.pyplot as plt

    valid = [
        row
        for row in rows
        if not any(np.isnan(v) for v in row[1:])
    ]
    if not valid:
        raise RuntimeError("プロット可能なデータがありません（すべて NaN）。")

    positions = np.array([r[0] for r in valid], dtype=float)
    ux = np.array([r[1] for r in valid], dtype=float)
    uy = np.array([r[2] for r in valid], dtype=float)
    uz = np.array([r[3] for r in valid], dtype=float)
    ur = np.array([r[4] for r in valid], dtype=float)
    utheta = np.array([r[5] for r in valid], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 8))
    ax.plot(ux, positions, label="Ux", marker="o", markersize=3)
    ax.plot(uy, positions, label="Uy", marker="o", markersize=3)
    ax.plot(uz, positions, label="Uz", marker="o", markersize=3)
    ax.plot(ur, positions, label="Ur", marker="o", markersize=3)
    ax.plot(utheta, positions, label="Utheta", marker="o", markersize=3)
    ax.set_xlabel("Velocity")
    ax.set_ylabel(axis_name)
    ax.set_title(f"Slice-averaged velocity profile  (time = {time})")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "OpenFOAM 結果を下端からピッチ指定でスライスし、"
            "各断面の Ux/Uy/Uz および Ur/Utheta の面積加重平均を CSV に出力します。"
        )
    )
    parser.add_argument(
        "case",
        type=Path,
        help="ケースディレクトリ、または .foam ファイル",
    )
    parser.add_argument(
        "--pitch",
        type=float,
        required=True,
        help="スライス断面の間隔（法線軸方向の長さ）",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("foam_velocity_slice_average.csv"),
        help="出力 CSV（既定: foam_velocity_slice_average.csv）",
    )
    parser.add_argument(
        "--normal",
        choices=("x", "y", "z"),
        default="z",
        help="スライス面の法線方向（既定: z）",
    )
    parser.add_argument(
        "--offset",
        type=float,
        default=0.0,
        help="下端から最初の断面までのオフセット（既定: 0）",
    )
    parser.add_argument(
        "--start",
        type=float,
        default=None,
        help="スライス開始位置（未指定時はメッシュ下端）",
    )
    parser.add_argument(
        "--end",
        type=float,
        default=None,
        help="スライス終了位置（未指定時はメッシュ上端）",
    )
    parser.add_argument(
        "--time",
        type=str,
        default="latest",
        help="時刻（数値）、または first / latest（既定: latest）",
    )
    parser.add_argument(
        "--plot",
        nargs="?",
        const="",
        metavar="PLOT",
        help=(
            "流速プロファイル図を出力する（縦軸=座標、横軸=Ux/Uy/Uz/Ur/Utheta）。"
            " パス省略時は foam_velocity_slice_average.png"
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="--plot 出力の DPI（既定: 150）",
    )
    args = parser.parse_args()

    if args.pitch <= 0:
        print("エラー: --pitch は正の数を指定してください。", file=sys.stderr)
        return 1
    if args.start is not None and args.end is not None and args.start > args.end:
        print("エラー: --start は --end 以下である必要があります。", file=sys.stderr)
        return 1

    try:
        import pyvista as pv
    except ImportError:
        print(
            "エラー: pyvista が必要です。  pip install pyvista",
            file=sys.stderr,
        )
        return 1

    try:
        case_file = _resolve_openfoam_path(args.case)
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        return 1

    reader = pv.OpenFOAMReader(str(case_file))
    t = _pick_time_index(reader, args.time)
    reader.set_active_time_value(t)
    multiblock = reader.read()
    mesh = _extract_internal_mesh(multiblock)

    axis_idx = _axis_index(args.normal)
    mesh_lo, mesh_hi = _axis_bounds(mesh, axis_idx)
    lo = mesh_lo if args.start is None else args.start
    hi = mesh_hi if args.end is None else args.end

    if lo > hi:
        print("エラー: スライス範囲が無効です（開始 > 終了）。", file=sys.stderr)
        return 1

    try:
        positions = _slice_positions(lo, hi, args.pitch, args.offset)
    except ValueError as e:
        print(f"エラー: {e}", file=sys.stderr)
        return 1

    if not positions:
        print(
            "警告: 生成されたスライス位置がありません。"
            f"（範囲 [{lo}, {hi}], pitch={args.pitch}, offset={args.offset}）",
            file=sys.stderr,
        )

    normal = _normal_axis(args.normal)
    axis_center = _axis_center(mesh)
    rows: list[tuple[float, float, float, float, float, float]] = []
    empty_count = 0

    for position in positions:
        origin = _slice_origin(mesh, axis_idx, position)
        sl = mesh.slice(normal=normal, origin=origin)
        ux, uy, uz, ur, utheta = _slice_area_weighted_velocity_mean(
            sl, axis_idx, axis_center
        )
        if np.isnan(ux):
            empty_count += 1
        rows.append((position, ux, uy, uz, ur, utheta))

    out = args.output.resolve()
    _write_csv(out, args.normal, rows)

    print(
        f"保存しました: {out}  "
        f"(time={t}, normal={args.normal}, slices={len(rows)}, empty={empty_count})"
    )
    if empty_count:
        print(
            f"警告: {empty_count} 断面で有効なセルがありませんでした（NaN を出力）。",
            file=sys.stderr,
        )

    if args.plot is not None:
        plot_path = (
            Path(args.plot)
            if args.plot
            else Path("foam_velocity_slice_average.png")
        )
        plot_out = plot_path.resolve()
        try:
            _save_profile_plot(
                plot_out,
                args.normal,
                rows,
                time=t,
                dpi=args.dpi,
            )
        except RuntimeError as e:
            print(f"警告: グラフ出力をスキップしました: {e}", file=sys.stderr)
        else:
            print(f"保存しました: {plot_out}  (time={t}, normal={args.normal})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

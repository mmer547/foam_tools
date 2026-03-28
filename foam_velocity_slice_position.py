#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenFOAM の結果を読み込み、任意の点を通る座標平面でスライスした
流速（速度ベクトル U の大きさ）のコンタ図を画像として保存する派生版。

平面の法線を --normal で、平面が通る点を -p / --origin で指定する。
例: 法線 z、原点 (0,0,1.2) なら「z=1.2 の XY 面」と同じ向きの断面。

依存: pip install pyvista numpy matplotlib

使用例:
  python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z
  python foam_velocity_slice_position.py path/to/case --origin 1.0 0 0 --normal x -o slice.png
  python foam_velocity_slice_position.py path/to/case -p 0 2 0 --normal y --time 500
  python foam_velocity_slice_position.py path/to/case --origin "(1.0, 0, 0)" --normal x
  # -p / --origin は OX OY OZ の3語、または "(OX, OY, OZ)" / OX,OY,OZ の1語でも可
  # ケース直下に空の case.foam があると読みやすい（無くてもディレクトリを渡せば可な場合あり）
"""

from __future__ import annotations

import argparse
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


def _parse_origin_token(s: str) -> tuple[float, float, float]:
    """括弧付き・カンマ区切り・空白区切りのいずれかで3成分を解釈する。"""
    t = s.strip()
    if t.startswith("(") and t.endswith(")"):
        t = t[1:-1].strip()
    if not t:
        raise ValueError("空の座標です")
    if "," in t:
        parts = [p.strip() for p in t.split(",") if p.strip() != ""]
    else:
        parts = t.split()
    if len(parts) != 3:
        raise ValueError(f"成分は3つである必要があります（得られた成分数: {len(parts)}）")
    try:
        return (float(parts[0]), float(parts[1]), float(parts[2]))
    except ValueError as e:
        raise ValueError(f"数値として解釈できません: {s!r}") from e


class _OriginTripletAction(argparse.Action):
    """-p/--origin に 3 数値、または (OX,OY,OZ) 形式の 1 引数を渡す。"""

    def __init__(self, option_strings, dest, nargs="+", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 1:
            try:
                triplet = _parse_origin_token(values[0])
            except ValueError as e:
                parser.error(f"{option_string}: {e}")
        elif len(values) == 3:
            try:
                triplet = (
                    float(values[0]),
                    float(values[1]),
                    float(values[2]),
                )
            except ValueError:
                parser.error(f"{option_string}: 数値として解釈できません {values!r}")
        else:
            parser.error(
                f"{option_string}: 3 つの数値、または括弧・カンマの1引数で指定してください "
                f'(例: 0 0 1 または "(0, 0, 1)")。受け取った引数: {values!r}'
            )
        setattr(namespace, self.dest, triplet)


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


def _ensure_velocity_magnitude(mesh) -> tuple[str, str]:
    """point_data / cell_data の U から U_mag を作り、(場の所在, スカラー名) を返す。"""
    u = None
    loc = ""
    if "U" in mesh.point_data:
        u = np.asarray(mesh.point_data["U"])
        loc = "point"
    elif "U" in mesh.cell_data:
        u = np.asarray(mesh.cell_data["U"])
        loc = "cell"
    else:
        keys = list(mesh.point_data.keys()) + list(mesh.cell_data.keys())
        raise RuntimeError(
            "速度場 U が見つかりません。"
            f"利用可能なデータ: {keys[:20]}{'...' if len(keys) > 20 else ''}"
        )

    if u.ndim != 2 or u.shape[1] < 3:
        raise RuntimeError(f"U の形状が想定外です: {u.shape}（(N,3) を想定）")

    mag = np.linalg.norm(u[:, :3], axis=1)
    name = "U_mag"
    if loc == "point":
        mesh.point_data[name] = mag
    else:
        mesh.cell_data[name] = mag
    return loc, name


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "OpenFOAM 結果を読み、指定した点を通る平面で流速コンタを保存します。"
            "平面の向きは --normal、通る点は -p（長い形 --origin）で指定します。"
        )
    )
    parser.add_argument(
        "case",
        type=Path,
        help="ケースディレクトリ、または .foam ファイル",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("foam_velocity_slice.png"),
        help="出力画像（既定: foam_velocity_slice.png）",
    )
    parser.add_argument(
        "--normal",
        choices=("x", "y", "z"),
        default="z",
        help=(
            "スライス面の法線（x→面の法線が±X、実際の面は -p/--origin の x 座標を通る）"
        ),
    )
    parser.add_argument(
        "-p",
        "--origin",
        action=_OriginTripletAction,
        nargs="+",
        metavar="ORIGIN",
        default=None,
        help=(
            "断面が通る点: OX OY OZ の3語、または "
            '"(OX, OY, OZ)" / OX,OY,OZ の1語。未指定時は 0 0 0'
        ),
    )
    parser.add_argument(
        "--time",
        type=str,
        default="latest",
        help="時刻（数値）、または first / latest（既定: latest）",
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="turbo",
        help="カラーマップ名（matplotlib）",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="出力 DPI",
    )
    args = parser.parse_args()
    origin = args.origin if args.origin is not None else (0.0, 0.0, 0.0)

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

    _ensure_velocity_magnitude(mesh)
    normal = _normal_axis(args.normal)
    origin_t = tuple(origin)

    sl = mesh.slice(normal=normal, origin=origin_t)
    if sl.n_cells == 0:
        print(
            "警告: スライスが空です。指定位置付近にメッシュが無い、"
            "または面が領域と交差していません。",
            file=sys.stderr,
        )

    # コンタ図（2D スライスを off-screen で保存）
    plotter = pv.Plotter(off_screen=True, window_size=(900, 700))
    plotter.add_mesh(
        sl,
        scalars="U_mag",
        cmap=args.cmap,
        show_edges=False,
        scalar_bar_args={"title": "|U|", "vertical": True},
    )
    plotter.view_vector(normal)
    plotter.camera.parallel_projection = True
    plotter.add_axes()
    plotter.add_text(f"time = {t}", font_size=10)

    out = args.output.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(out), transparent_background=False)
    plotter.close()
    ox, oy, oz = origin_t
    print(
        f"保存しました: {out}  (time={t}, normal={args.normal}, origin=({ox}, {oy}, {oz}))"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

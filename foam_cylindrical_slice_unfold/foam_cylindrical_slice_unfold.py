#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenFOAM の体積メッシュを読み、vtkCylinder + vtkCutter で半径 r の円筒面と交差させ、
得られた表面を周方向角 theta と軸方向座標で展開し平面 (u, v) に写像する。
速度 U があるときは展開面内の接線成分 (U·e_theta, U·e_axial, 0) を point ベクトルとして追加する
（ParaView の SurfaceLIC など用。既定の名前は U_unfold）。

写像（既定: 軸が +Z、軸上の基準点は --axis-origin）:
  u = r * theta   （theta = atan2(y - cy, x - cx)、ラジアン、展開後の周方向弧長）
  v = z - cz  （軸方向。軸が x や y のときは対応する軸方向座標）

依存: pip install pyvista numpy matplotlib（vtk は pyvista に同梱）

使用例:
  python foam_cylindrical_slice_unfold.py path/to/case --radius 0.05
  python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --triangles -o unfold.png
  python foam_cylindrical_slice_unfold.py path/to/case -r 0.05 --axis y --axis-origin 0 0 0 --time latest
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


def _parse_origin_token(s: str) -> tuple[float, float, float]:
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
    """--axis-origin に 3 数値、または (OX,OY,OZ) 形式の 1 引数を渡す。"""

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
                f'(例: 0 0 0 または "(0, 0, 0)")。受け取った引数: {values!r}'
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


def _pick_scalar_array_name(mesh, explicit: str | None) -> str:
    if explicit:
        if explicit in mesh.point_data:
            return explicit
        if explicit in mesh.cell_data:
            return explicit
        raise RuntimeError(
            f"--scalars で指定した {explicit!r} が point_data / cell_data にありません。"
        )
    if "U" in mesh.point_data or "U" in mesh.cell_data:
        _ensure_velocity_magnitude(mesh)
        return "U_mag"
    for name in mesh.point_data.keys():
        arr = mesh.point_data[name]
        if np.asarray(arr).ndim == 1 and name not in ("vtkOriginalPointIds",):
            return name
    for name in mesh.cell_data.keys():
        arr = mesh.cell_data[name]
        if np.asarray(arr).ndim == 1:
            return name
    raise RuntimeError("表示用スカラーが見つかりません。--scalars で名前を指定してください。")


def _cylindrical_rt_axial(
    points: np.ndarray,
    axis: str,
    origin: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """各点について円筒座標 r, theta, axial（軸に沿った座標）を返す。"""
    ox, oy, oz = origin
    x = points[:, 0] - ox
    y = points[:, 1] - oy
    z = points[:, 2] - oz
    axis_l = axis.strip().lower()
    if axis_l == "z":
        r = np.sqrt(x * x + y * y)
        theta = np.arctan2(y, x)
        axial = z
    elif axis_l == "y":
        r = np.sqrt(x * x + z * z)
        theta = np.arctan2(z, x)
        axial = y
    elif axis_l == "x":
        r = np.sqrt(y * y + z * z)
        theta = np.arctan2(z, y)
        axial = x
    else:
        raise ValueError(f"axis は x, y, z のいずれか: {axis!r}")
    return r, theta, axial


def _unfold_points(
    points: np.ndarray,
    axis: str,
    axis_origin: tuple[float, float, float],
    radius_nominal: float,
) -> np.ndarray:
    _, theta, axial = _cylindrical_rt_axial(points, axis, axis_origin)
    u = radius_nominal * theta
    v = axial
    return np.column_stack([u, v, np.zeros(points.shape[0], dtype=points.dtype)])


def _cylinder_uv_basis(
    points: np.ndarray,
    axis: str,
    origin: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """
    円筒面上の接線基底（単位ベクトル）。
    e_theta: 周方向（theta 増加方向）。e_axial: 円筒軸方向（+X / +Y / +Z）。
    いずれも (N, 3)。
    """
    ox, oy, oz = origin
    dx = points[:, 0] - ox
    dy = points[:, 1] - oy
    dz = points[:, 2] - oz
    n = points.shape[0]
    ax = axis.strip().lower()
    if ax == "z":
        r = np.sqrt(dx * dx + dy * dy)
        e_theta = np.column_stack([-dy, dx, np.zeros(n, dtype=points.dtype)])
        e_axial = np.tile(np.array([0.0, 0.0, 1.0], dtype=points.dtype), (n, 1))
    elif ax == "y":
        r = np.sqrt(dx * dx + dz * dz)
        e_theta = np.column_stack([-dz, np.zeros(n, dtype=points.dtype), dx])
        e_axial = np.tile(np.array([0.0, 1.0, 0.0], dtype=points.dtype), (n, 1))
    elif ax == "x":
        r = np.sqrt(dy * dy + dz * dz)
        e_theta = np.column_stack([np.zeros(n, dtype=points.dtype), -dz, dy])
        e_axial = np.tile(np.array([1.0, 0.0, 0.0], dtype=points.dtype), (n, 1))
    else:
        raise ValueError(f"axis は x, y, z のいずれか: {axis!r}")
    r_safe = np.maximum(r, 1e-30)
    e_theta = e_theta / r_safe[:, np.newaxis]
    return e_theta, e_axial


def _velocity_unfold_xyz(
    points: np.ndarray,
    U: np.ndarray,
    axis: str,
    axis_origin: tuple[float, float, float],
) -> np.ndarray:
    """世界座標の U を展開平面上の (u方向, v方向, 0) 成分に射影（各点で内積）。"""
    if U.ndim != 2 or U.shape[1] < 3:
        raise RuntimeError(f"U の形状が想定外です: {U.shape}（(N,3) を想定）")
    e_th, e_ax = _cylinder_uv_basis(points, axis, axis_origin)
    uu = np.einsum("ij,ij->i", U[:, :3], e_th)
    vv = np.einsum("ij,ij->i", U[:, :3], e_ax)
    z0 = np.zeros(points.shape[0], dtype=U.dtype)
    return np.column_stack([uu, vv, z0])


def _ensure_point_velocity(mesh) -> np.ndarray | None:
    """スライスメッシュの各頂点での U (N,3)。無ければ None。cell のみのときは point に複写する。"""
    if "U" in mesh.point_data:
        return np.asarray(mesh.point_data["U"], dtype=float)
    if "U" not in mesh.cell_data:
        return None
    c2p = mesh.cell_data_to_point_data()
    if "U" not in c2p.point_data:
        return None
    u_p = np.asarray(c2p.point_data["U"], dtype=float)
    mesh.point_data["U"] = u_p
    return u_p


def _vtk_cylinder_implicit(
    radius: float,
    axis: str,
    axis_origin: tuple[float, float, float],
):
    """
    無限円筒面（半径 radius）の vtkImplicitFunction。
    VTK の vtkCylinder は既定で軸が +Y。--axis に合わせて回転し、軸が --axis-origin を通るように平行移動する。
    """
    import vtk

    cyl = vtk.vtkCylinder()
    cyl.SetRadius(float(radius))
    cyl.SetCenter(0.0, 0.0, 0.0)
    cx, cy, cz = (float(axis_origin[0]), float(axis_origin[1]), float(axis_origin[2]))
    tf = vtk.vtkTransform()
    tf.PostMultiply()
    tf.Identity()
    ax = axis.strip().lower()
    if ax == "z":
        tf.RotateX(90.0)
    elif ax == "x":
        tf.RotateZ(-90.0)
    elif ax == "y":
        pass
    else:
        raise ValueError(f"axis は x, y, z のいずれか: {axis!r}")
    tf.Translate(cx, cy, cz)
    cyl.SetTransform(tf)
    return cyl


def _resolve_mesh_save_path(path: Path, mesh) -> Path:
    """PolyData 等、mesh.save で許可される拡張子に合わせる。"""
    writers = getattr(mesh, "_WRITERS", None) or {}
    allowed = frozenset(writers.keys())
    path = path.expanduser().resolve()
    ext = path.suffix.lower()
    if ext == "":
        path = path.with_suffix(".vtp")
        ext = ".vtp"
    if ext in allowed:
        return path
    order = (".vtp", ".vtu", ".vtk", ".vtkhdf")
    fallback = next((e for e in order if e in allowed), None)
    if fallback is None and allowed:
        fallback = sorted(allowed)[0]
    if fallback is None:
        raise RuntimeError(
            "このメッシュ型に対応する save 拡張子がありません。"
            f"（利用可能: {sorted(allowed)}）"
        )
    new_path = path.with_suffix(fallback)
    print(
        f"注意: {ext or '(拡張子なし)'} はこのメッシュ型では保存できません。"
        f" {new_path.name}（{fallback}）に保存します。",
        file=sys.stderr,
    )
    return new_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "OpenFOAM 結果を円筒面でスライス（vtkCutter）し、周方向・軸方向に展開した平面表示を保存します。"
        )
    )
    parser.add_argument(
        "case",
        type=Path,
        help="ケースディレクトリ、または .foam ファイル",
    )
    parser.add_argument(
        "-r",
        "--radius",
        type=float,
        required=True,
        help="円筒半径（展開の周方向スケール r*theta に使う目標半径）",
    )
    parser.add_argument(
        "--triangles",
        action="store_true",
        help="カッター出力を三角形に分割する",
    )
    parser.add_argument(
        "--axis",
        choices=("x", "y", "z"),
        default="z",
        help="円筒軸の向き（+X / +Y / +Z）。既定: z",
    )
    parser.add_argument(
        "--axis-origin",
        action=_OriginTripletAction,
        nargs="+",
        metavar="XYZ",
        default=(0.0, 0.0, 0.0),
        help='軸が通る基準点 (cx cy cz)、または "(cx, cy, cz)"（既定: 0 0 0）',
    )
    parser.add_argument(
        "--time",
        type=str,
        default="latest",
        help="時刻（数値）、または first / latest（既定: latest）",
    )
    parser.add_argument(
        "--scalars",
        type=str,
        default=None,
        help="塗り分けに使うスカラー名（未指定時は U から |U|、なければ先頭の1成分スカラー）",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("foam_cylindrical_unfold.png"),
        help="出力画像（既定: foam_cylindrical_unfold.png）",
    )
    parser.add_argument(
        "--vtk",
        type=Path,
        default=None,
        help="展開後メッシュを保存するパス（省略時は書き出さない）。PolyData のため .vtp 等",
    )
    parser.add_argument(
        "--no-u-unfold",
        action="store_true",
        help="速度 U の展開面ベクトル（SurfaceLIC 用）を追加しない",
    )
    parser.add_argument(
        "--u-unfold-name",
        type=str,
        default="U_unfold",
        metavar="NAME",
        help="展開面内速度の point ベクトル名（既定: U_unfold）。第1成分周方向第2成分軸方向",
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

    r0 = float(args.radius)
    if r0 <= 0:
        print("エラー: --radius は正の数である必要があります。", file=sys.stderr)
        return 1

    try:
        import pyvista as pv
    except ImportError:
        print("エラー: pyvista が必要です。  pip install pyvista", file=sys.stderr)
        return 1

    try:
        case_file = _resolve_openfoam_path(args.case)
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        return 1

    axis_origin = tuple(args.axis_origin)

    reader = pv.OpenFOAMReader(str(case_file))
    t = _pick_time_index(reader, args.time)
    reader.set_active_time_value(t)
    multiblock = reader.read()
    mesh = _extract_internal_mesh(multiblock)

    try:
        impl = _vtk_cylinder_implicit(r0, args.axis, axis_origin)
    except ImportError:
        print(
            "エラー: vtk が必要です（通常は pip install pyvista で入ります）。",
            file=sys.stderr,
        )
        return 1

    cut = mesh.slice_implicit(
        impl, generate_triangles=args.triangles, progress_bar=False
    )
    if cut.n_cells == 0:
        print(
            "エラー: 円筒スライスが空です。--radius・--axis・--axis-origin と計算領域の交差を確認してください。",
            file=sys.stderr,
        )
        return 1

    scalar_name = _pick_scalar_array_name(cut, args.scalars)
    pts = np.asarray(cut.points, dtype=float)

    if not args.no_u_unfold:
        u_pts = _ensure_point_velocity(cut)
        if u_pts is not None:
            if u_pts.shape[0] != pts.shape[0]:
                print(
                    "警告: U の点数が頂点数と一致しません。U_unfold は追加しません。",
                    file=sys.stderr,
                )
            else:
                name_u = (args.u_unfold_name.strip() or "U_unfold")
                cut.point_data[name_u] = _velocity_unfold_xyz(
                    pts, u_pts, args.axis, axis_origin
                )
                print(
                    f"SurfaceLIC 用ベクトル {name_u!r} を追加しました（周方向, 軸方向, 0）。"
                )

    cut.points = _unfold_points(pts, args.axis, axis_origin, r0)

    if args.vtk is not None:
        out_mesh = _resolve_mesh_save_path(args.vtk, cut)
        out_mesh.parent.mkdir(parents=True, exist_ok=True)
        cut.save(str(out_mesh))
        print(f"展開メッシュを保存しました: {out_mesh}")

    plotter = pv.Plotter(off_screen=True, window_size=(1000, 700))
    plotter.add_mesh(
        cut,
        scalars=scalar_name,
        cmap=args.cmap,
        show_edges=False,
        scalar_bar_args={"title": scalar_name, "vertical": True},
    )
    plotter.camera_position = "xy"
    plotter.camera.parallel_projection = True
    plotter.add_axes()
    u_span = 2.0 * np.pi * r0
    plotter.add_text(
        f"time = {t}  |  u = r*theta (幅 ~ {u_span:.4g}), v = axial  |  r = {r0:g} (cyl slice)",
        font_size=9,
    )

    out = args.output.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(out), transparent_background=False)
    plotter.close()
    print(f"保存しました: {out}  (time={t}, axis={args.axis}, r={r0:g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

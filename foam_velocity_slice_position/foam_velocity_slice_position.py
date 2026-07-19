#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenFOAM の結果を読み込み、任意の点を通る座標平面でスライスした
流速のコンタ図を画像として保存する派生版。

平面の法線を --normal で、平面が通る点を -p / --origin で指定する。
例: 法線 z、原点 (0,0,1.2) なら「z=1.2 の XY 面」と同じ向きの断面。
--scalars で |U| または U の各成分（Ux / Uy / Uz）を選べる。

依存: pip install pyvista numpy matplotlib

使用例:
  python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z
  python foam_velocity_slice_position.py path/to/case --origin 1.0 0 0 --normal x -o slice.png
  python foam_velocity_slice_position.py path/to/case -p 0 2 0 --normal y --time 500
  python foam_velocity_slice_position.py path/to/case --origin "(1.0, 0, 0)" --normal x
  python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z --scalars Ux
  python foam_velocity_slice_position.py path/to/case -p 0 0 0.5 --normal z --surface-lic -o lic.png
  # -p / --origin は OX OY OZ の3語、または "(OX, OY, OZ)" / OX,OY,OZ の1語でも可
  # --surface-lic 指定時は VTK の vtkSurfaceLICMapper で SurfaceLIC 画像を出力する
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


# --scalars の別名 → (配列名, カラーバータイトル)
_SCALAR_CHOICES: dict[str, tuple[str, str]] = {
    "u_mag": ("U_mag", "|U|"),
    "mag": ("U_mag", "|U|"),
    "|u|": ("U_mag", "|U|"),
    "ux": ("Ux", "Ux"),
    "u_x": ("Ux", "Ux"),
    "uy": ("Uy", "Uy"),
    "u_y": ("Uy", "Uy"),
    "uz": ("Uz", "Uz"),
    "u_z": ("Uz", "Uz"),
}


def _resolve_scalar_choice(name: str) -> tuple[str, str]:
    """--scalars の指定を (配列名, カラーバータイトル) に正規化する。"""
    key = name.strip().lower()
    if key in _SCALAR_CHOICES:
        return _SCALAR_CHOICES[key]
    raise ValueError(
        f"不明な --scalars: {name!r}。"
        " U_mag / Ux / Uy / Uz（別名: mag, U_x, U_y, U_z）を指定してください。"
    )


def _ensure_velocity_scalars(mesh) -> str:
    """
    point_data / cell_data の U から U_mag / Ux / Uy / Uz を作る。
    場の所在（"point" / "cell"）を返す。
    """
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

    u3 = u[:, :3]
    fields = {
        "U_mag": np.linalg.norm(u3, axis=1),
        "Ux": u3[:, 0].copy(),
        "Uy": u3[:, 1].copy(),
        "Uz": u3[:, 2].copy(),
    }
    target = mesh.point_data if loc == "point" else mesh.cell_data
    for name, values in fields.items():
        target[name] = values
    return loc


def _prepare_slice_for_surface_lic(sl, normal: tuple[float, float, float]):
    """
    SurfaceLIC 用にスライスへ point_data の U / スカラー / U_plane を整える。
    U_plane は法線成分を除いた面内速度。
    """
    if "U" not in sl.point_data:
        if "U" in sl.cell_data:
            sl = sl.cell_data_to_point_data()
        else:
            keys = list(sl.point_data.keys()) + list(sl.cell_data.keys())
            raise RuntimeError(
                "SurfaceLIC 用の速度場 U が見つかりません。"
                f"利用可能なデータ: {keys[:20]}{'...' if len(keys) > 20 else ''}"
            )

    u = np.asarray(sl.point_data["U"], dtype=float)
    if u.ndim != 2 or u.shape[1] < 3:
        raise RuntimeError(f"U の形状が想定外です: {u.shape}（(N,3) を想定）")

    n = np.asarray(normal, dtype=float)
    n_norm = float(np.linalg.norm(n))
    if n_norm <= 0.0:
        raise RuntimeError("法線ベクトルがゼロです。")
    n = n / n_norm
    u3 = u[:, :3]
    sl.point_data["U_plane"] = u3 - np.outer(u3 @ n, n)
    _ensure_velocity_scalars(sl)
    return sl


def _matplotlib_lookup_table(cmap_name: str, smin: float, smax: float, n: int = 256):
    """matplotlib カラーマップから vtkLookupTable を作る。"""
    import vtk
    from matplotlib import colormaps

    try:
        cmap = colormaps[cmap_name]
    except KeyError as e:
        raise RuntimeError(f"不明なカラーマップです: {cmap_name!r}") from e

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(int(n))
    lut.SetRange(float(smin), float(smax))
    lut.Build()
    denom = max(n - 1, 1)
    for i in range(n):
        r, g, b, a = cmap(i / denom)
        lut.SetTableValue(i, float(r), float(g), float(b), float(a))
    return lut


def _in_plane_size(
    mesh, view_normal: tuple[float, float, float]
) -> tuple[float, float]:
    """スライス面内でのメッシュ幅・高さ（ワールド単位）。"""
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    ax = abs(view_normal[0])
    ay = abs(view_normal[1])
    az = abs(view_normal[2])
    if az >= ax and az >= ay:
        return float(xmax - xmin), float(ymax - ymin)
    if ay >= ax and ay >= az:
        return float(xmax - xmin), float(zmax - zmin)
    return float(ymax - ymin), float(zmax - zmin)


# 右端のカラーバー＋数値ラベル用に確保する画面幅の割合
_COLORBAR_RIGHT_RESERVE = 0.24


def _window_size_for_slice(
    mesh,
    view_normal: tuple[float, float, float],
    *,
    base: int = 1000,
    right_reserve: float = _COLORBAR_RIGHT_RESERVE,
    max_side: int = 1600,
) -> tuple[int, int]:
    """
    メッシュの面内アスペクト比に合わせたウィンドウサイズ。
    右端のカラーバー分だけ横幅を余分に取る。
    """
    dx, dy = _in_plane_size(mesh, view_normal)
    mesh_aspect = max(dx, 1e-30) / max(dy, 1e-30)
    usable = max(1.0 - right_reserve, 0.5)
    # ウィンドウ全体の幅/高さ（メッシュ表示領域が mesh_aspect になるよう）
    win_aspect = mesh_aspect / usable
    if win_aspect >= 1.0:
        h = base
        w = int(round(base * win_aspect))
    else:
        w = base
        h = int(round(base / win_aspect))
    w = int(min(max(w, 640), max_side))
    h = int(min(max(h, 640), max_side))
    return (w, h)


def _fit_parallel_camera(
    plotter,
    mesh,
    view_normal: tuple[float, float, float],
    *,
    margin: float = 0.06,
    right_reserve: float = _COLORBAR_RIGHT_RESERVE,
) -> None:
    """
    平行投影カメラをメッシュにフィットさせる。
    右端 right_reserve をカラーバー用に空け、メッシュは左側領域に収める。
    """
    plotter.view_vector(view_normal)
    plotter.camera.parallel_projection = True
    plotter.reset_camera()

    dx, dy = _in_plane_size(mesh, view_normal)
    win_w, win_h = plotter.window_size
    aspect = float(win_w) / max(float(win_h), 1.0)
    usable = max(1.0 - right_reserve, 0.5)
    pad = 1.0 + margin
    # parallel_scale = ビューポート高さの半分（ワールド単位）
    scale_h = 0.5 * dy * pad
    scale_w = 0.5 * dx * pad / max(aspect * usable, 1e-30)
    cam = plotter.camera
    cam.parallel_scale = max(scale_h, scale_w, 1e-30)

    # WindowCenter: 投影中心をずらしてメッシュを左寄せ（右端を空ける）
    # 座標は [-1, 1]。正の x で内容が左へ寄る。
    cam.SetWindowCenter(float(right_reserve), 0.0)


def _screenshot_surface_lic(
    mesh,
    output: Path,
    *,
    vector_name: str,
    scalar_name: str,
    scalar_title: str,
    cmap: str,
    title: str,
    view_normal: tuple[float, float, float],
    window_size: tuple[int, int] | None = None,
) -> None:
    """
    VTK の vtkSurfaceLICMapper で SurfaceLIC を描画し PNG 保存する。
    PyVista は Plotter（カメラ・オフスクリーン）補助に使う。
    """
    import pyvista as pv
    import vtk

    if vector_name not in mesh.point_data:
        raise RuntimeError(
            f"SurfaceLIC 用ベクトル {vector_name!r} が point_data にありません。"
        )
    if scalar_name not in mesh.point_data:
        raise RuntimeError(
            f"コンター用スカラー {scalar_name!r} が point_data にありません。"
        )

    vectors = np.asarray(mesh.point_data[vector_name], dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] < 3:
        raise RuntimeError(
            f"{vector_name!r} の形状が想定外です: {vectors.shape}（(N,3) を想定）"
        )
    if not np.isfinite(vectors[:, :3]).any() or np.allclose(vectors[:, :3], 0.0):
        raise RuntimeError(
            f"{vector_name!r} がゼロまたは非有限です。SurfaceLIC を描画できません。"
        )

    # LIC は三角形サーフェス向け
    surf = mesh
    if not getattr(surf, "is_all_triangles", False):
        try:
            surf = surf.triangulate()
        except Exception:
            pass

    scalars = np.asarray(surf.point_data[scalar_name], dtype=float).reshape(-1)
    smin = float(np.nanmin(scalars))
    smax = float(np.nanmax(scalars))
    if not np.isfinite(smin) or not np.isfinite(smax) or smin == smax:
        smin, smax = 0.0, 1.0

    mapper = vtk.vtkSurfaceLICMapper()
    mapper.SetInputData(surf)
    pd = surf.GetPointData()
    pd.SetActiveVectors(vector_name)
    pd.SetActiveScalars(scalar_name)
    mapper.SetScalarVisibility(True)
    mapper.SetScalarModeToUsePointData()
    mapper.SetColorModeToMapScalars()
    lut = _matplotlib_lookup_table(cmap, smin, smax)
    mapper.SetLookupTable(lut)
    mapper.SetUseLookupTableScalarRange(True)

    lic = mapper.GetLICInterface()
    lic.SetNoiseTextureSize(384)
    lic.SetNoiseGrainSize(1)
    lic.SetNumberOfSteps(44)
    lic.SetStepSize(0.22)
    lic.SetLICIntensity(0.6)
    lic.SetEnhancedLIC(1)
    lic.SetAntiAlias(1)
    lic.SetEnhanceContrast(lic.ENHANCE_CONTRAST_LIC)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle(scalar_title)
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.UnconstrainedFontSizeOn()
    scalar_bar.SetWidth(0.05)
    scalar_bar.SetHeight(0.55)
    # 右端予約領域内に配置。数値はバーの右側へ（メッシュ側に食い込まない）
    scalar_bar.SetPosition(0.86, 0.22)
    scalar_bar.SetTextPositionToSucceedScalarBar()
    title_prop = scalar_bar.GetTitleTextProperty()
    title_prop.SetFontSize(12)
    title_prop.BoldOff()
    title_prop.ItalicOff()
    title_prop.ShadowOff()
    title_prop.SetColor(0.1, 0.1, 0.1)
    label_prop = scalar_bar.GetLabelTextProperty()
    label_prop.SetFontSize(10)
    label_prop.BoldOff()
    label_prop.ItalicOff()
    label_prop.ShadowOff()
    label_prop.SetColor(0.1, 0.1, 0.1)

    if window_size is None:
        window_size = _window_size_for_slice(surf, view_normal)

    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.set_background("white")
    plotter.renderer.AddActor(actor)
    plotter.renderer.AddActor2D(scalar_bar)
    _fit_parallel_camera(plotter, surf, view_normal)
    plotter.add_axes()
    if title:
        plotter.add_text(title, font_size=10)

    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(output), transparent_background=False)
    plotter.close()


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
        "--scalars",
        type=str,
        default="U_mag",
        help=(
            "コンタに使うスカラー: U_mag（既定）/ Ux / Uy / Uz"
            "（別名: mag, U_x, U_y, U_z）"
        ),
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
    parser.add_argument(
        "--surface-lic",
        action="store_true",
        help=(
            "出力画像を VTK vtkSurfaceLICMapper による SurfaceLIC にする"
            "（面内速度 U_plane と --scalars のコンタ）"
        ),
    )
    args = parser.parse_args()
    origin = args.origin if args.origin is not None else (0.0, 0.0, 0.0)
    try:
        scalar_name, scalar_title = _resolve_scalar_choice(args.scalars)
    except ValueError as e:
        print(f"エラー: {e}", file=sys.stderr)
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

    _ensure_velocity_scalars(mesh)
    normal = _normal_axis(args.normal)
    origin_t = tuple(origin)

    sl = mesh.slice(normal=normal, origin=origin_t)
    if sl.n_cells == 0:
        print(
            "警告: スライスが空です。指定位置付近にメッシュが無い、"
            "または面が領域と交差していません。",
            file=sys.stderr,
        )

    out = args.output.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    ox, oy, oz = origin_t
    title = f"time = {t}  |  {scalar_title}"
    meta = (
        f"(time={t}, normal={args.normal}, origin=({ox}, {oy}, {oz}), "
        f"scalars={scalar_name})"
    )

    if args.surface_lic:
        try:
            sl = _prepare_slice_for_surface_lic(sl, normal)
            _screenshot_surface_lic(
                sl,
                out,
                vector_name="U_plane",
                scalar_name=scalar_name,
                scalar_title=scalar_title,
                cmap=args.cmap,
                title=title + " | SurfaceLIC",
                view_normal=normal,
            )
        except Exception as e:
            print(f"エラー: SurfaceLIC 画像の出力に失敗しました: {e}", file=sys.stderr)
            return 1
        print(f"SurfaceLIC を保存しました: {out}  {meta}")
        return 0

    # コンタ図（2D スライスを off-screen で保存）
    win = _window_size_for_slice(sl, normal)
    plotter = pv.Plotter(off_screen=True, window_size=win)
    plotter.add_mesh(
        sl,
        scalars=scalar_name,
        cmap=args.cmap,
        show_edges=False,
        scalar_bar_args={
            "title": scalar_title,
            "vertical": True,
            "width": 0.05,
            "height": 0.55,
            "position_x": 0.86,
            "position_y": 0.22,
            "title_font_size": 12,
            "label_font_size": 10,
            "n_labels": 5,
        },
    )
    _fit_parallel_camera(plotter, sl, normal)
    # 数値ラベルをバー右側へ（メッシュ側へ食い込まない）
    for bar in plotter.scalar_bars.values():
        try:
            bar.SetTextPositionToSucceedScalarBar()
        except Exception:
            pass
    plotter.add_axes()
    plotter.add_text(title, font_size=10)
    plotter.screenshot(str(out), transparent_background=False)
    plotter.close()
    print(f"保存しました: {out}  {meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

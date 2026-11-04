#!/usr/bin/env python3
"""BGMN/Profex の .lst から Rietveld 結果を直接集計して Excel に保存する。

実行例:
  全件実行:
    uv run --no-project --with pandas --with openpyxl --with pymatgen python 2_Program/compile_rietveld_from_lst.py

  1ファイルだけ確認:
    uv run --no-project --with pandas --with openpyxl --with pymatgen python 2_Program/compile_rietveld_from_lst.py --sample 0_Data/1_Rietveld/SPring-8/Mafic/SP-MTN-43/SP-MTN-43.lst --output-dir /private/tmp/rietveld_check
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd
from pymatgen.core import Element


# BGMN の .lst テキストを読むための正規表現。
# 例:
#   A=0.4913667+-0.0000021
#   Rp=3.43%  Rpb=6.31%  R=4.46%  Rwp=4.52% Rexp=1.92%
#   Local parameters and GOALs for phase SiO2p3221
NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
VALUE_RE = re.compile(rf"^(?P<name>[^=\s]+)\s*=\s*(?P<value>{NUMBER_RE}|UNDEF)(?:\+\-(?P<esd>{NUMBER_RE}))?")
FIT_RE = re.compile(rf"\b(?P<name>Rp|Rpb|R|Rwp|Rexp)=(?P<value>{NUMBER_RE})%")
PHASE_RE = re.compile(r"^Local parameters and GOALs for phase\s+(.+?)\s*$")

# Atomic positions の行は、先頭に multiplicity、後半に E=(...) の元素情報を持つ。
# BGMN は SI, AL, FE のように元素記号を大文字で出すため、ここでは大文字1-2文字を
# 受け取り、後で canonical_element() で Si/Al/Fe の表記へ正規化する。
ATOM_RE = re.compile(r"^\s*(?P<multiplicity>\d+)\s+")
EL_RE = re.compile(r"(?P<element>[A-Z]{1,2})(?:[+-]?\d*)?\((?P<occupancy>[-+]?\d+(?:\.\d*)?|\.\d+)\)")

# Param 出力で抽出する格子定数・粒径パラメータ。
# Chem の酸化物列は、新規元素への汎用性を保つため固定せず、実データから動的に作る。
PARAM_COLUMNS = ["A", "B", "C", "ALPHA", "BETA", "GAMMA", "GrainSize(1,1,1)"]

# 陽イオン量を、従来の酸化物表記の質量へ変換する係数。
# 結晶構造中の O の数をそのまま使うのではなく、旧 Chem.CSV と同じ酸化物種
# (SiO2, Al2O3, Fe2O3 など) として報告する。
# ここにない元素は酸化物名を推定せず、元素名そのものの wt% として出力する。
OXIDE_FORMULAS = {
    "H": ("H2O", {"H": 2, "O": 1}, 2),
    "C": ("CO2", {"C": 1, "O": 2}, 1),
    "F": ("F", {"F": 1}, 1),
    "Na": ("Na2O", {"Na": 2, "O": 1}, 2),
    "Mg": ("MgO", {"Mg": 1, "O": 1}, 1),
    "Al": ("Al2O3", {"Al": 2, "O": 3}, 2),
    "Si": ("SiO2", {"Si": 1, "O": 2}, 1),
    "P": ("P2O5", {"P": 2, "O": 5}, 2),
    "Cl": ("Cl", {"Cl": 1}, 1),
    "K": ("K2O", {"K": 2, "O": 1}, 2),
    "Ca": ("CaO", {"Ca": 1, "O": 1}, 1),
    "Ti": ("TiO2", {"Ti": 1, "O": 2}, 1),
    "Mn": ("MnO", {"Mn": 1, "O": 1}, 1),
    "Fe": ("Fe2O3", {"Fe": 2, "O": 3}, 2),
}

# GOAL 名は phase block 名と完全一致しないことがある。
# よく出る短縮名だけをこの対応表で解決する。
# Qquartz は、同じ .lst 内に QSiO2p3221 と SiO2p3221-2 が併存する場合があるため、
# 固定 alias ではなく assign_goal_phase() で未割当の quartz 系 phase へ対応させる。
GOAL_PHASE_ALIASES = {
    "musc2m1": "MUSCOVITE_2M1",
    "chlorite2b": "CHLORITE_IIb-2",
    "fluorapa": "Fluorapatite",
    "albite": "Plagioclase_Albite",
    "calcite": "Calcite",
    "titanite": "Titanite",
    "epidote": "Epidote",
    "actinolite": "Actinolite",
    "glaucophane": "Glaucophane",
    "ilmenite": "Ilmenite",
    "tio2r": "TiO2r",
    "hematite": "Hematite",
    "graphite2h": "Graphite-2H",
}


@dataclass
class PhaseBlock:
    """`.lst` 内の 'Local parameters and GOALs for phase ...' 1ブロック分。"""

    name: str
    params: dict[str, tuple[float | None, float | None]] = field(default_factory=dict)
    composition: dict[str, float] = field(default_factory=lambda: defaultdict(float))


@dataclass
class ParsedLst:
    """1つの `.lst` から Chem/Param/Goal 表を再構築するための抽出結果。"""

    path: Path
    fit: dict[str, float]
    goals: dict[str, tuple[float, float | None]]
    phases: list[PhaseBlock]


def parse_args() -> argparse.Namespace:
    """全ファイル一括処理または1ファイル確認用の CLI オプションを読む。"""

    parser = argparse.ArgumentParser(
        description="BGMN/Profex の .lst から集計 Excel を作成します。",
        epilog=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-root", type=Path, default=Path("0_Data/1_Rietveld"))
    parser.add_argument("--output-dir", type=Path, default=Path("0_Data/2_Compiled_Data"))
    parser.add_argument("--sample", type=Path, help="Compile one .lst file instead of scanning input-root.")
    parser.add_argument("--include-old", action="store_true", help="Include paths containing OLD/old.")
    parser.add_argument("--save-intermediate-csv", action="store_true")
    return parser.parse_args()


def as_float(value: str | None) -> float | None:
    """BGMN の数値文字列を float にし、UNDEF や空欄は欠損値として扱う。"""

    if value is None or value == "" or value == "UNDEF":
        return None
    return float(value)


def parse_value_line(line: str) -> tuple[str, float | None, float | None] | None:
    """`A=0.491+-0.001` のような行を name, value, ESD に分解する。"""

    match = VALUE_RE.match(line.strip())
    if not match:
        return None
    return match.group("name"), as_float(match.group("value")), as_float(match.group("esd"))


def normalize_name(value: str) -> str:
    """`Qmusc2m1` と `MUSCOVITE_2M1` のような名前を比較しやすくする。"""

    return re.sub(r"[^a-z0-9]", "", value.lower())


def sample_metadata(path: Path, input_root: Path) -> dict[str, str]:
    """旧 CSV export に含まれていたメタデータ列をパスから作る。"""

    rel_parts = path.relative_to(input_root).parts
    if len(rel_parts) < 3:
        raise ValueError(f"Cannot infer instrument/rock type from path: {path}")
    instrument, rock_type = rel_parts[0], rel_parts[1]
    return {
        "File": str(path),
        "Sample": path.stem,
        "Sample ID": f"{path.stem}.txt Data 1",
        "Instrument": instrument,
        "Rock_Type": rock_type,
    }


def find_lst_files(input_root: Path, include_old: bool) -> list[Path]:
    """現行の Laboratory/SPring-8 かつ Pelitic/Mafic 配下の `.lst` を探す。"""

    paths: list[Path] = []
    for path in input_root.glob("*/*/**/*.lst"):
        if path.name.startswith("._"):
            continue
        rel_parts = path.relative_to(input_root).parts
        if len(rel_parts) < 3:
            continue
        if rel_parts[0] not in {"Laboratory", "SPring-8"}:
            continue
        if rel_parts[1] not in {"Pelitic", "Mafic"}:
            continue
        if not include_old and any(part.lower() == "old" for part in rel_parts):
            continue
        paths.append(path)
    return sorted(paths)


def parse_lst(path: Path) -> ParsedLst:
    """1つの `.lst` を fit 統計量、GOAL 値、phase block に分解する。"""

    lines = path.read_text(errors="replace").splitlines()
    fit: dict[str, float] = {}
    goals: dict[str, tuple[float, float | None]] = {}
    phases: list[PhaseBlock] = []
    current: PhaseBlock | None = None
    in_global = False
    in_atomic = False

    for line in lines:
        # Rwp/Rexp などの fit 統計量は冒頭の1行にまとまっているため、見つかったものを回収する。
        fit.update({match.group("name"): float(match.group("value")) for match in FIT_RE.finditer(line)})

        # Global block には、鉱物量に相当する Q... 値が並んでいる。
        if line.startswith("Global parameters and GOALs"):
            in_global = True
            current = None
            in_atomic = False
            continue

        # Local phase block には、格子定数と atomic positions が含まれる。
        phase_match = PHASE_RE.match(line)
        if phase_match:
            in_global = False
            in_atomic = False
            current = PhaseBlock(phase_match.group(1).strip())
            phases.append(current)
            continue

        # この marker 以降では、E=(...) を元素組成として読む。
        if current and line.startswith("Atomic positions for phase"):
            in_atomic = True
            continue

        if in_global:
            # ここでは Q... 行だけを残す。EPS2 など鉱物モードではない値は、
            # 旧 Goal.CSV 型の鉱物量テーブルでは使わない。
            parsed = parse_value_line(line)
            if parsed:
                name, value, esd = parsed
                if value is not None and name.startswith("Q"):
                    goals[name] = (value, esd)
            continue

        if not current:
            continue

        parsed = parse_value_line(line)
        if parsed and not in_atomic:
            # Param_compile_df に必要な格子定数・粒径パラメータだけを保存する。
            name, value, esd = parsed
            if name in PARAM_COLUMNS:
                current.params[name] = (value, esd)
            continue

        if in_atomic:
            # 各元素について multiplicity * occupancy を足し、Refined Composition に
            # 相当する化学式量を作る。
            atom_match = ATOM_RE.match(line)
            element_text = extract_element_text(line)
            if not atom_match or element_text is None:
                continue
            multiplicity = float(atom_match.group("multiplicity"))
            for element_match in EL_RE.finditer(element_text):
                element = canonical_element(element_match.group("element"))
                occupancy = float(element_match.group("occupancy"))
                current.composition[element] += multiplicity * occupancy

    return ParsedLst(path=path, fit=fit, goals=goals, phases=phases)


def extract_element_text(line: str) -> str | None:
    """occupancy の括弧を含む E=(...) から、中身の元素情報だけを取り出す。"""

    marker = "E=("
    start = line.find(marker)
    if start == -1:
        return None
    text = line[start + len(marker) :].strip()
    if text.endswith(")"):
        text = text[:-1]
    return text


def canonical_element(element: str) -> str:
    """BGMN の大文字元素記号、例 SI/FE を、通常の Si/Fe 表記へ直す。"""

    element = element.strip()
    if len(element) == 1:
        return element.upper()
    return element[0].upper() + element[1:].lower()


def format_amount(value: float) -> str:
    """旧 Param.CSV に近い見た目で化学式中の量を整形する。"""

    if abs(value - round(value)) < 1e-8:
        return str(int(round(value)))
    return f"{value:.4f}"


def refined_composition(composition: dict[str, float]) -> str:
    """集計した元素数から `Refined Composition` 文字列を作る。"""

    if not composition:
        return ""
    items = sorted(composition.items(), key=lambda item: (item[0] == "O", item[0]))
    return " ".join(f"{element.upper()}{format_amount(amount)}" for element, amount in items if amount > 1e-10)


def atomic_weight(element: str) -> float:
    """pymatgen から元素の原子量を取得する。"""

    return float(Element(element).atomic_mass)


def oxide_factor(element: str) -> tuple[str, float]:
    """元素量を出力用の酸化物または元素 wt% へ変換する係数を返す。"""

    if element not in OXIDE_FORMULAS:
        return element, 1.0

    oxide_name, formula, element_count = OXIDE_FORMULAS[element]
    oxide_mass = sum(atomic_weight(formula_element) * count for formula_element, count in formula.items())
    element_mass = atomic_weight(element) * element_count
    return oxide_name, oxide_mass / element_mass


def oxide_wt_percent(composition: dict[str, float]) -> dict[str, float]:
    """元素組成を酸化物 wt% に変換する。

    対応表にない新規元素は、酸化物名を推定せず元素名の wt% として残す。
    """

    oxide_masses: dict[str, float] = defaultdict(float)
    for element, amount in composition.items():
        if element == "O" or amount <= 0:
            continue
        try:
            oxide, factor = oxide_factor(element)
            oxide_masses[oxide] += amount * atomic_weight(element) * factor
        except ValueError:
            continue

    total = sum(oxide_masses.values())
    if total <= 0:
        return {}
    return {oxide: mass / total * 100.0 for oxide, mass in oxide_masses.items()}


def goal_phase_name(goal_name: str, phase_names: Iterable[str]) -> str | None:
    """Q... の GOAL 名を、安全に対応できる local phase 名へ変換する。

    完全一致または明示 alias だけを採用する。Qquartz のように割当済み phase を
    見て判断する必要がある GOAL は、assign_goal_phase() 側で処理する。
    """

    raw = goal_name[1:] if goal_name.startswith("Q") else goal_name
    phase_names = list(phase_names)
    raw_norm = normalize_name(raw)

    # まず Q を除いた名前と phase block 名が完全一致する場合を採用する。
    for phase in phase_names:
        if raw_norm == normalize_name(phase):
            return phase

    # 完全一致しない短縮名だけを、明示した対応表で解決する。
    alias = GOAL_PHASE_ALIASES.get(normalize_name(raw))
    if alias:
        for phase in phase_names:
            if normalize_name(phase) == normalize_name(alias):
                return phase
    return None


def is_quartz_phase_name(phase_name: str) -> bool:
    """phase 名が quartz 系 SiO2 phase かを判定する。"""

    normalized = normalize_name(phase_name)
    return normalized.startswith("sio2p3221")


def assign_goal_phase(goal_name: str, phase_names: Iterable[str], phase_quantities: dict[str, float]) -> str | None:
    """Chem の phase quantity に使う phase 名を決める。

    通常は完全一致または明示 alias で対応付ける。Qquartz は、すでに
    QSiO2p3221 などで使われた quartz phase を避け、未割当の quartz 系 phase へ入れる。
    """

    phase_names = list(phase_names)
    raw = goal_name[1:] if goal_name.startswith("Q") else goal_name
    if normalize_name(raw) == "quartz":
        for phase in phase_names:
            if is_quartz_phase_name(phase) and phase_quantities.get(phase, 0.0) == 0.0:
                return phase
        for phase in phase_names:
            if is_quartz_phase_name(phase):
                return phase
        return None
    return goal_phase_name(goal_name, phase_names)


def build_rows(parsed_files: list[ParsedLst], input_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """parse 済み `.lst` から Chem_df, Param_df, Goal_df 相当の表を作る。"""

    chem_rows: list[dict[str, object]] = []
    param_rows: list[dict[str, object]] = []
    goal_rows: list[dict[str, object]] = []

    for parsed in parsed_files:
        meta = sample_metadata(parsed.path, input_root)
        phase_by_name = {phase.name: phase for phase in parsed.phases}
        phase_names = list(phase_by_name)

        # Q... 値は Goal_df には全て保存し、Chem の重み付けにも phase quantity として使う。
        # Qquartz は割当済みの SiO2p3221 を上書きせず、未割当の quartz 系 phase へ入れる。
        phase_quantities = {phase: 0.0 for phase in phase_names}

        # Goal_df: 鉱物量 Q 値と refinement quality の統計量を格納する。
        for goal, (value, esd) in parsed.goals.items():
            phase = assign_goal_phase(goal, phase_names, phase_quantities)
            goal_rows.append({**meta, "Parameter, Goal": goal, "Value": value, "ESD": esd})
            if phase is None:
                print(f"Skipped Chem phase quantity assignment: {meta['Sample']} {goal}")
                continue
            if phase_quantities[phase] != 0.0:
                print(f"Skipped duplicate Chem phase quantity assignment: {meta['Sample']} {goal} -> {phase}")
                continue
            phase_quantities[phase] = value

        if "Rwp" in parsed.fit:
            goal_rows.append({**meta, "Parameter, Goal": "Rwp", "Value": parsed.fit["Rwp"], "ESD": None})
        if "Rexp" in parsed.fit:
            goal_rows.append({**meta, "Parameter, Goal": "Rexp", "Value": parsed.fit["Rexp"], "ESD": None})
        if parsed.fit.get("Rwp") is not None and parsed.fit.get("Rexp"):
            gof = parsed.fit["Rwp"] / parsed.fit["Rexp"]
            goal_rows.append({**meta, "Parameter, Goal": "Chi2", "Value": gof * gof, "ESD": None})
            goal_rows.append({**meta, "Parameter, Goal": "GOF", "Value": gof, "ESD": None})

        # Param_df と、phase ごとの Chem_df 行を作る。
        for phase in parsed.phases:
            param_row: dict[str, object] = {
                **meta,
                "Phase": phase.name,
                "Refined Composition": refined_composition(phase.composition),
            }
            for column in PARAM_COLUMNS:
                value, esd = phase.params.get(column, (None, None))
                param_row[column] = value
                param_row[f"ESD({column})"] = esd
            param_rows.append(param_row)

            # Chem_df は旧 export と同じく、phase x oxide ごとに1行を持つ。
            oxides = oxide_wt_percent(phase.composition)
            phase_quantity = phase_quantities.get(phase.name, 0.0)
            for oxide, quantity in oxides.items():
                chem_rows.append(
                    {
                        **meta,
                        "Phase": phase.name,
                        "Phase Quantity": phase_quantity,
                        "Oxide": oxide,
                        "Quantity (wt-%)": quantity,
                    }
                )

        # Weighted total は、鉱物量で重み付けして再構築したバルク組成。
        weighted = weighted_bulk_oxides(parsed.phases, phase_quantities)
        weighted_total = sum(phase_quantities.values())
        for oxide, quantity in weighted.items():
            chem_rows.append(
                {
                    **meta,
                    "Phase": "Weighted total",
                    "Phase Quantity": weighted_total,
                    "Oxide": oxide,
                    "Quantity (wt-%)": quantity,
                }
            )

    return pd.DataFrame(chem_rows), pd.DataFrame(param_rows), pd.DataFrame(goal_rows)


def weighted_bulk_oxides(phases: list[PhaseBlock], phase_quantities: dict[str, float]) -> dict[str, float]:
    """phase quantity で重み付けした合計として、バルク酸化物 wt% を計算する。

    旧 Chem.CSV の Weighted total は、対応できた phase quantity の合計で再正規化せず、
    各 phase の酸化物 wt% に phase quantity を掛けた合計として出力されている。
    """

    totals: dict[str, float] = defaultdict(float)
    total_quantity = sum(phase_quantities.get(phase.name, 0.0) for phase in phases)
    if total_quantity <= 0:
        return {}
    for phase in phases:
        quantity = phase_quantities.get(phase.name, 0.0)
        oxides = oxide_wt_percent(phase.composition)
        for oxide, value in oxides.items():
            totals[oxide] += value * quantity
    return dict(totals)


def write_outputs(chem_df: pd.DataFrame, param_df: pd.DataFrame, goal_df: pd.DataFrame, output_dir: Path) -> None:
    """元の notebook と同じ5つの Excel ファイルを書き出す。"""

    output_dir.mkdir(parents=True, exist_ok=True)

    # 鉱物化学組成テーブルと、そのうち Weighted total だけを抜き出した表。
    pivot_chem_df = chem_df.pivot_table(
        index=["Sample", "Instrument", "Rock_Type", "Phase", "Phase Quantity"],
        columns="Oxide",
        values="Quantity (wt-%)",
    )
    pivot_bulk_comp_df = pivot_chem_df.loc[pivot_chem_df.index.get_level_values("Phase") == "Weighted total"]

    param_compile_df = param_df[["Sample", "Instrument", "Rock_Type", "Phase", *PARAM_COLUMNS]]

    # 鉱物モード表。Rwp/Rexp/Chi2/GOF は下で pivot_error_df として分ける。
    pivot_bulk_mode_df = goal_df.pivot_table(
        index=["Sample", "Instrument", "Rock_Type"],
        columns="Parameter, Goal",
        values="Value",
    ).fillna(0)

    error_col = ["Rwp", "Rexp", "Chi2", "GOF"]
    for column in error_col:
        if column not in pivot_bulk_mode_df.columns:
            pivot_bulk_mode_df[column] = 0.0
    pivot_error_df = pivot_bulk_mode_df[error_col]
    pivot_bulk_mode_df = pivot_bulk_mode_df.drop(error_col, axis=1).fillna(0)

    pivot_chem_df.to_excel(output_dir / "pivot_chem_df.xlsx")
    pivot_bulk_comp_df.to_excel(output_dir / "pivot_bulk_comp_df.xlsx")
    param_compile_df.to_excel(output_dir / "Param_compile_df.xlsx", index=False)
    pivot_error_df.to_excel(output_dir / "pivot_error_df.xlsx")
    pivot_bulk_mode_df.to_excel(output_dir / "pivot_bulk_mode_df.xlsx")


def main() -> None:
    """CLI の実行入口。"""

    args = parse_args()
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()

    if args.sample:
        # 1ファイル指定モード。新しく refine したサンプルの確認に使う。
        lst_paths = [args.sample.resolve()]
        if not lst_paths[0].is_relative_to(input_root):
            input_root = common_parent(input_root, lst_paths[0])
    else:
        # 一括処理モード。旧 notebook と同じく、現行の Rietveld 結果フォルダを走査する。
        lst_paths = [path.resolve() for path in find_lst_files(input_root, args.include_old)]

    if not lst_paths:
        raise SystemExit("No .lst files found.")

    parsed_files = [parse_lst(path) for path in lst_paths]
    chem_df, param_df, goal_df = build_rows(parsed_files, input_root)
    write_outputs(chem_df, param_df, goal_df, output_dir)

    if args.save_intermediate_csv:
        # デバッグ用。再構築した3つの中間テーブルを確認したい場合だけ CSV 保存する。
        chem_df.to_csv(output_dir / "Chem_from_lst.csv", index=False)
        param_df.to_csv(output_dir / "Param_from_lst.csv", index=False)
        goal_df.to_csv(output_dir / "Goal_from_lst.csv", index=False)

    print(f"Compiled {len(lst_paths)} .lst file(s) to {output_dir}")


def common_parent(root: Path, child: Path) -> Path:
    """--sample が --input-root の外にある場合の fallback root を求める。"""

    root_parts = root.parts
    child_parts = child.parts
    common: list[str] = []
    for left, right in zip(root_parts, child_parts):
        if left != right:
            break
        common.append(left)
    return Path(*common) if common else Path("/")


if __name__ == "__main__":
    main()

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from compile_rietveld_str_parameters import (
    PARAMETERS,
    SUMMARY_COLUMNS,
    find_structure_files,
    output_path_for,
    parse_structure,
    prefer_generic_structures,
    summarize_records,
)


class CompileRietveldStrParametersTests(unittest.TestCase):
    def test_refinement_flags_and_gewicht(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "sample" / "Albite.str"
            path.parent.mkdir()
            path.write_text(
                "PHASE=Plagioclase_Albite //\n"
                "PARAM=A=8.14E-01_0.812^0.815 PARAM=B=1.27_1.2^1.3 "
                "PARAM=C=0.716_0.7^0.8 //\n"
                "PARAM=ALPHA=94_93^95 PARAM=BETA=116_115^117 "
                "PARAM=GAMMA=87_86^88 //\n"
                "k1=0 PARAM=k2=0_0^0.0001 PARAM=B1=0_0^0.01 "
                "GEWICHT=SPHAR4 //\n",
                encoding="utf-8",
            )

            parsed = parse_structure(path, root)

            self.assertEqual(parsed.phase, "Plagioclase_Albite")
            self.assertEqual(parsed.gewicht, 4)
            self.assertEqual(parsed.flags["k1"], 0)
            for name in set(PARAMETERS) - {"k1"}:
                self.assertEqual(parsed.flags[name], 1)

    def test_zero_value_with_param_is_refined(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "sample" / "phase.str"
            path.parent.mkdir()
            path.write_text(
                "PHASE=Test\nPARAM=k1=0_0^1 k2=ANISOSQR^0.0001 "
                "B1=ANISO^0.01 GEWICHT=SPHAR10\n",
                encoding="utf-8",
            )

            parsed = parse_structure(path, root)

            self.assertEqual(parsed.flags["k1"], 1)
            self.assertEqual(parsed.flags["k2"], 0)
            self.assertEqual(parsed.flags["B1"], 0)
            self.assertEqual(parsed.gewicht, 10)

    def test_generic_file_wins_over_sample_suffixed_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = root / "sample"
            folder.mkdir()
            generic = folder / "Albite.str"
            suffixed = folder / "Albite-sample.str"
            content = "PHASE=Albite //\nGEWICHT=SPHAR4\n"
            generic.write_text(content, encoding="utf-8")
            suffixed.write_text(content, encoding="utf-8")

            selected = prefer_generic_structures(
                [parse_structure(generic, root), parse_structure(suffixed, root)]
            )

            self.assertEqual([record.path.name for record in selected], ["Albite.str"])

    def test_uppercase_extension_and_old_are_handled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            current = root / "sample" / "PHASE.STR"
            old = root / "OLD_2024" / "old.str"
            sidecar = root / "sample" / "._PHASE.str"
            current.parent.mkdir()
            old.parent.mkdir()
            current.write_text("PHASE=Test\nGEWICHT=SPHAR0\n", encoding="utf-8")
            old.write_text("PHASE=Old\nGEWICHT=SPHAR0\n", encoding="utf-8")
            sidecar.write_text("PHASE=Sidecar\nGEWICHT=SPHAR0\n", encoding="utf-8")

            self.assertEqual(find_structure_files(root), [current])

    def test_output_path_identifies_instrument_and_rock_type(self) -> None:
        root = Path("0_Data/1_Rietveld/SPring-8/Mafic")
        output_dir = Path("0_Data/2_Compiled_Data/Refined_parameter")
        self.assertEqual(
            output_path_for(root, output_dir),
            output_dir / "SPring-8_Mafic_Rietveld_refinement_parameters.xlsx",
        )

    def test_phase_summary_uses_zero_one_and_two(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = []
            for folder, k1, k2, gewicht in (
                ("sample1", False, True, 4),
                ("sample2", False, False, 6),
                ("sample3", False, False, 4),
            ):
                path = root / folder / "phase.str"
                path.parent.mkdir()
                tokens = ["PHASE=TestPhase", f"GEWICHT=SPHAR{gewicht}"]
                if k1:
                    tokens.append("PARAM=k1=0_0^1")
                if k2:
                    tokens.append("PARAM=k2=0_0^0.001")
                path.write_text(" ".join(tokens), encoding="utf-8")
                records.append(parse_structure(path, root))

            summary = summarize_records(records)

            self.assertEqual(len(summary), 1)
            row = dict(zip(SUMMARY_COLUMNS, summary[0]))
            self.assertEqual(row["Phase"], "TestPhase")
            self.assertEqual(row["Record_Count"], 3)
            self.assertEqual(row["k1"], 0)
            self.assertEqual(row["k2"], 2)
            self.assertEqual(row["A"], 0)
            self.assertEqual(row["GEWICHT_Mode"], 4)
            self.assertEqual(row["GEWICHT_Min"], 4)
            self.assertEqual(row["GEWICHT_Max"], 6)

    def test_gewicht_reports_all_tied_modes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = []
            for folder, gewicht in (("sample1", 2), ("sample2", 4)):
                path = root / folder / "phase.str"
                path.parent.mkdir()
                path.write_text(
                    f"PHASE=TestPhase //\nGEWICHT=SPHAR{gewicht}\n", encoding="utf-8"
                )
                records.append(parse_structure(path, root))

            row = dict(zip(SUMMARY_COLUMNS, summarize_records(records)[0]))

            self.assertEqual(row["GEWICHT_Mode"], "2, 4")
            self.assertEqual(row["GEWICHT_Min"], 2)
            self.assertEqual(row["GEWICHT_Max"], 4)


if __name__ == "__main__":
    unittest.main()

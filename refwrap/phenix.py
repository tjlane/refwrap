from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gemmi
import reciprocalspaceship as rs
from reciprocalspaceship.utils import add_rfree


@dataclass
class PhenixParams:
    number_of_macrocycles: int = 3

    def format(self) -> list[str]:
        formatted_parameters: list[str] = [
            f"main.number_of_macro_cycles={self.number_of_macrocycles}",
        ]
        return formatted_parameters


@dataclass
class RefinementStats:
    r_work: float
    r_free: float

    @staticmethod
    def extract_r_values(text: str) -> tuple[float, float]:
        match = re.search(r"R-work\s*=\s*([0-9.]+).*?R-free\s*=\s*([0-9.]+)", text)
        if match:
            r_work = float(match.group(1))
            r_free = float(match.group(2))
            return r_work, r_free
        msg = "R values not found in the text."
        raise ValueError(msg)

    @classmethod
    def read_from_log(cls, logfile: Path) -> RefinementStats:
        logtext = logfile.read_text()
        r_work, r_free = cls.extract_r_values(logtext)
        return RefinementStats(
            r_work=r_work,
            r_free=r_free,
        )


def phenixrefine(
    *,
    structure: gemmi.Structure,
    amplitudes: rs.DataSeries,
    uncertainties: rs.DataSeries,
    rflags: rs.DataSeries,
    cell: Any,
    spacegroup: Any,
    phenix_params: PhenixParams,
) -> tuple[gemmi.Structure, rs.DataSet, RefinementStats]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        refinement_data = rs.DataSet(
            {"F-obs": amplitudes, "SIGF-obs": uncertainties, "R-free-flags": rflags}
        )
        refinement_data.cell = cell
        refinement_data.spacegroup = spacegroup
        refinement_data.write_mtz(str(tmpdir_path / "input.mtz"))

        structure.write_pdb(str(tmpdir_path / "input.pdb"))

        cmd = ["phenix.refine", "input.pdb", "input.mtz", *phenix_params.format()]

        process_return = subprocess.run(cmd, cwd=tmpdir, check=False, capture_output=True)
        print(process_return)

        if process_return.returncode != 0:
            raise subprocess.CalledProcessError(process_return.returncode, cmd=cmd)

        refined_pdb = tmpdir_path / "input_refine_001.pdb"
        refined_mtz = tmpdir_path / "input_refine_001.mtz"
        refinment_log = tmpdir_path / "input_refine_001.log"

        shutil.copy(tmpdir_path / "input_refine_001.log", "/Users/tjlane/Desktop")

        assert refined_pdb.exists()
        assert refined_mtz.exists()
        assert refinment_log.exists()

        structure = gemmi.read_structure(str(refined_pdb))
        mtz_dataset = rs.read_mtz(str(refined_mtz))
        stats = RefinementStats.read_from_log(refinment_log)

        return structure, mtz_dataset, stats


if __name__ == "__main__":
    p = PhenixParams(number_of_macrocycles=1)
    mtz = Path("/Users/tjlane/Desktop/ocp-xfel-bench/models/ech/10ps/light-F.mtz")
    pdb = Path(
        "/Users/tjlane/Desktop/ocp-xfel-bench/models/ech/10ps/10ps_extr_realspace-occ20_1.pdb"
    )

    ds = rs.read_mtz(str(mtz))
    add_rfree(ds, fraction=0.05, ccp4_convention=False, inplace=True)
    init_structure = gemmi.read_structure(str(pdb))

    structure, mtz_dataset, stats = phenixrefine(
        structure=init_structure,
        amplitudes=ds["FW-F"],
        uncertainties=ds["FW-SIGF"],
        rflags=ds["R-free-flags"],
        cell=ds.cell,
        spacegroup=ds.spacegroup,
        phenix_params=p,
    )

    print(structure)
    print(mtz_dataset)
    print(stats)

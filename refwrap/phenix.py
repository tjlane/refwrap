from __future__ import annotations

import subprocess
import tempfile
import gemmi
import os
import re

from pathlib import Path
from dataclasses import dataclass
import reciprocalspaceship as rs
import shutil


@dataclass
class MyPhenixParams:
    data_labels: str | None = None
    number_of_macrocycles: int = 3

    def format(self) -> list[str]:
        formatted_parameters: list[str] = [
            f"main.number_of_macro_cycles={self.number_of_macrocycles}",
        ]

        if self.data_labels:
            formatted_parameters.append(f"refinement.input.xray_data.labels={self.data_labels}")

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
        else:
            raise ValueError("R values not found in the text.")

    @classmethod
    def read_from_log(cls, logfile: Path) -> RefinementStats:
        with open(logfile, "r") as f:
            logtext = logfile.read_text()

        r_work, r_free = cls.extract_r_values(logtext)

        return RefinementStats(
            r_work=r_work,
            r_free=r_free,
        )


def phenixrefine(
        *,
        pdb_file: Path,
        mtz_file: Path,
        phenix_params: MyPhenixParams,
) -> tuple[gemmi.Structure, rs.DataSet, RefinementStats]:

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        os.symlink(pdb_file, tmpdir_path / 'input.pdb')
        os.symlink(mtz_file, tmpdir_path / 'input.mtz')

        cmd = [
            "phenix.refine",
            "input.pdb",
            "input.mtz",
            *phenix_params.format()
        ]
        
        process_return = subprocess.run(cmd, cwd=tmpdir, check=False, capture_output=True)
        print(process_return)

        if process_return.returncode != 0:
            raise subprocess.CalledProcessError(process_return.returncode, cmd=cmd)

        refined_pdb = tmpdir_path / "input_refine_001.pdb"
        refined_mtz = tmpdir_path / "input_refine_001.mtz"
        refinment_log = tmpdir_path / "input_refine_001.log"

        shutil.copy(tmpdir_path / 'input_refine_001.log', "/Users/tjlane/Desktop")

        assert refined_pdb.exists()
        assert refined_mtz.exists()
        assert refinment_log.exists()

        structure = gemmi.read_structure(str(refined_pdb))
        mtz_dataset = rs.read_mtz(str(refined_mtz))
        stats = RefinementStats.read_from_log(refinment_log)

        return structure, mtz_dataset, stats
    

if __name__ == "__main__":
    p = MyPhenixParams(data_labels="FW-F,FW-SIGF", number_of_macrocycles=1)
    mtz = Path("/Users/tjlane/Desktop/ocp-xfel-bench/models/ech/10ps/light-F.mtz")
    pdb = Path("/Users/tjlane/Desktop/ocp-xfel-bench/models/ech/10ps/10ps_extr_realspace-occ20_1.pdb")
    
    structure, mtz_dataset, stats = phenixrefine(pdb_file=pdb, mtz_file=mtz, phenix_params=p)

    print(structure)
    print(mtz_dataset)
    print(stats)
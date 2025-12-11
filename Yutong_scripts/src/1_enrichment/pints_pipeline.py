import subprocess
import pandas as pd

from typing import Iterable, Optional
from utils import dir_exists

class PintsPipeline:
    """
    A class to run the PINTS peak calling pipeline
    including pints_caller and pints_boundary_extender
    """
    def __init__(self, bw_pl: str, bw_mn: str, peak_dir: str,
                 job_name: str, alpha: float = 0.3, 
                 disable_sticks: bool = True, thread: int = 16,
                 promoter_file: Optional[str] = None):
        """
        Parameters
        ----------
        bw_pl : str
            Path to the plus strand bigwig file.
        bw_mn : str
            Path to the minus strand bigwig file.
        peak_dir : str
            Directory to save output files.
        job_name : str
            Prefix for output files.
        alpha : float, optional
            The stringency for PINTS to cluster nearby TSSs into a peak. 
            0 is the least stringent; 1 is the most stringent. (default: 0.3)
        disable_sticks : bool, optional
            Set this switch to remove stick-like peaks (signal on a single position) 
            (default: True)
        thread : int, optional
            Max number of threads PINTS can create (default: 16)
        promoter_file : Optional[str], optional
            Path to the promoter BED file for boundary extension, by default None.
        """
        self.bw_pl = bw_pl
        self.bw_mn = bw_mn
        self.peak_dir = peak_dir
        self.file_prefix = job_name
        self.alpha = alpha
        self.disable_sticks = disable_sticks
        self.thread = thread
        self.promoter_file = promoter_file

    def run_pints_caller(self):
        """
        Call peaks from the provided plus and minus strand bigwig files.
        """
        argv = [
            "pints_caller",
            "--bw-pl", self.bw_pl,
            "--bw-mn", self.bw_mn,
            "--save-to", self.peak_dir,
            "--file-prefix", self.file_prefix,
            "--alpha", str(self.alpha),
            "--thread", str(self.thread),
        ]
        if self.disable_sticks:
            argv.append("--remove-sticks")
            
        subprocess.run(argv)

    def run_pints_boundary_extender(self,
                                    peak_dir: Optional[str] = None,
                                    element_dir: Optional[str] = None,
                                    div_ext_left: int = 60,
                                    div_ext_right: int = 60,
                                    single_ext_left: int = 60,
                                    single_ext_right: int = 60,
                                    promoter_file: Optional[str] = None):
        """
        Extend peak boundaries from the called peaks. Need to run after 
        the peak calling step. This function find the called peaks directly in the 
        output directory. The bigwigs are still needed for finding summits 
        and prominent bases for extension.

        Parameters
        ----------
        peak_dir : Optional[str], optional
            Directory where the peak files are located. If None, use the peak_dir
            by default None.
        element_dir : Optional[str], optional
            Directory to save the extended peak files (elements). If None,
            use the peak_dir by default None.
        div_ext_left : int, optional
            divergent extension left (default: (60,))
        div_ext_right : int, optional
            divergent extension right (default: (60,))
        single_ext_left : int, optional
            unidirectional extension left (default: (60,))
        single_ext_right : int, optional
            unidirectional extension right (default: (60,))
        promoter_file : Optional[str], optional
            Promoter bed, if specified, will create a separate bed file for
            distal elements (default: None)
        """
        if peak_dir is None:
            peak_dir = self.peak_dir
        if element_dir is None:
            element_dir = self.peak_dir
        dir_exists(element_dir)

        promoter = self.promoter_file if promoter_file is None else promoter_file
        argv = [
            "pints_boundary_extender",
            "--bw-pl", self.bw_pl,
            "--bw-mn", self.bw_mn,
            "--divergent-files", f"{peak_dir}/{self.file_prefix}_1_divergent_peaks.bed",
            "--bidirectional-files", f"{peak_dir}/{self.file_prefix}_1_bidirectional_peaks.bed",
            "--unidirectional-files", f"{peak_dir}/{self.file_prefix}_1_unidirectional_peaks.bed",
            "--save-to", element_dir,
            "--div-ext-left", str(div_ext_left),
            "--div-ext-right", str(div_ext_right),
            "--single-ext-left", str(single_ext_left),
            "--single-ext-right", str(single_ext_right),
            "--promoter-bed", promoter
        ]

        subprocess.run(argv)

    def _bed_path(self, file_prefix: str, category: str, 
                  tag: str, distal: bool) -> str:
        """
        Helper function to get the bed file path.
        """
        suffix = "_e" if distal else ""
        return f"{file_prefix}_1_{category}_peaks_element_{tag}bp{suffix}.bed"

    def _read_bed3(self, path: str) -> pd.DataFrame:
        """
        Read a bed3 file into a DataFrame.
        """
        df = pd.read_csv(path, sep="\t", header=None, usecols=[0,1,2])
        df.columns = ["chr", "start", "end"]
        df["start"] = df["start"].astype(int)
        df["end"]   = df["end"].astype(int)
        return df

    def _union_beds(self, paths: Iterable[str]) -> pd.DataFrame:
        """
        Merge multiple bed3 files into one DataFrame and drop duplicates.
        """
        parts = [ self._read_bed3(p) for p in paths ]
        if not parts:
            raise ValueError("No bed files provided for merging.")
        return pd.concat(parts, ignore_index=True).drop_duplicates()

    def _write_bed(self, df: pd.DataFrame, out: str) -> None:
        """
        Write a DataFrame to a bed3 file without header.
        """
        df.to_csv(out, sep="\t", header=False, index=False)

    def merge_and_save_elements(self, element_dir: Optional[str] = None,
                             div_ext_left: int = 60, div_ext_right: int = 60,
                             single_ext_left: int = 60, single_ext_right: int = 60,
                             enhancer: bool = False, promoter: bool = False,
                             merged_prefix: str = "all_peaks_element",
                             include_unidirectional: bool = True):
        """
        Concatenate peaks and save as a final output.

        Parameters
        ----------
        element_dir : Optional[str], optional
            Directory where the extended peak files (elements) are located. 
            If None, use the peak_dir by default None.
        div_ext_left : int, optional
            divergent extension left (default: (60,))
        div_ext_right : int, optional
            divergent extension right (default: (60,))
        single_ext_left : int, optional
            unidirectional extension left (default: (60,))
        single_ext_right : int, optional
            unidirectional extension right (default: (60,))
        enhancer : bool, optional
            If True, will create a separate merged bed file for 
            distal elements (default: False)
        promoter : bool, optional
            If True, will create a separate merged bed file for 
            proximal elements (default: False)
        merged_prefix : str, optional
            Prefix for the merged output files (default: "all_peaks_element")
        include_unidirectional : bool, optional
            If True, include unidirectional peaks in the merging (default: True)
        """
        categories = ["bidirectional", "divergent"]
        if include_unidirectional:
            categories.append("unidirectional")

        if element_dir is None:
            element_dir = self.peak_dir
        dir_exists(element_dir)

        if div_ext_left != div_ext_right or single_ext_left != single_ext_right or \
            div_ext_left != single_ext_left:
            print("Warning: Right now only support symmetric extension for merging peaks.")
            exit(1)

        cat_paths = [
            self._bed_path(f"{element_dir}/{self.file_prefix}", 
                           cat, div_ext_left, distal=False)
            for cat in categories
        ]
        all_df = self._union_beds(cat_paths)
        out_all = f"{element_dir}/{merged_prefix}_{div_ext_left}bp.bed"
        self._write_bed(all_df, out_all)

        if enhancer or promoter:
            distal_paths = [
                self._bed_path(f"{element_dir}/{self.file_prefix}", 
                               cat, div_ext_left, distal=True)
                for cat in categories
            ]
            enh_df = self._union_beds(distal_paths)

        if enhancer:
            out_enh = f"{element_dir}/{merged_prefix}_{div_ext_left}bp_e.bed"
            self._write_bed(enh_df, out_enh)
        if promoter:
            prom_df = all_df.merge(enh_df, how='outer', indicator=True)
            prom_df = prom_df[prom_df['_merge'] == 'left_only'].drop(columns=['_merge'])
            out_prom = f"{element_dir}/{merged_prefix}_{div_ext_left}bp_p.bed"
            self._write_bed(prom_df, out_prom)
    
    def clear_intermediate_files(self, element_dir: Optional[str] = None,
                                 peaks: bool = True, elements: bool = True):
        """
        Remove intermediate files generated during peak calling and boundary extension.

        Parameters
        ----------
        element_dir : Optional[str], optional
            Directory where the extended peak files (elements) are located. 
            If None, use the peak_dir by default None.
        peaks : bool, optional
            If True, remove peak calling intermediate files (default: True).
        elements : bool, optional
            If True, remove boundary extension intermediate files (default: True).
        """
        if element_dir is None:
            element_dir = self.peak_dir

        patterns = []
        if peaks:
            patterns.append(f"{element_dir}/{self.file_prefix}_1_*_peaks.bed")

        if elements:
            patterns.append(f"{element_dir}/{self.file_prefix}_1_*_peaks_element_*bp*.bed")

        for pattern in patterns:
            subprocess.run(f"rm -f {pattern}", shell=True)
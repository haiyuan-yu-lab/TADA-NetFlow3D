import os
import pysam
import tempfile
import pybedtools
import pandas as pd
import numpy as np

from scipy.stats import fisher_exact
from utils import ensure_parent_dir, run_command
from typing import Dict, Iterable, Tuple, Union, Set, List, Optional


class VCFHandler:
    """
    A class to handle VCF files and perform various operations.
    """
    def __init__(self):
        pass

    def set_vcf(self, vcf_path: str):
        """
        Set the VCF file.
        """
        self.vcf_path = vcf_path

    def get_vcf(self) -> str:
        """
        Get the path of the VCF file.
        """
        return self.vcf_path
    
    def ensure_index(self, vcf_path: str = None):
        """
        Ensure that the VCF file is indexed. If not, create an index.
        """
        if vcf_path is None:
            vcf_path = self.vcf_path

        index_suffix = ".tbi" if vcf_path.endswith(".vcf.gz") else ".csi"
        index_path = f"{vcf_path}{index_suffix}"
        if not os.path.exists(index_path):
            command = ["bcftools", "index", "-f", vcf_path]
            run_command(command)

    def _info_keys_and_csq_fields(self, vcf_path: str
                                  ) -> Tuple[Set[str], Optional[Set[str]]]:
        """
        Return (info_keys, csq_fields) from VCF header.
        If no CSQ field, csq_fields is None.
        """
        with pysam.VariantFile(vcf_path) as vcf:
            info_keys = set(vcf.header.info.keys())
            csq_fields = None
            if "CSQ" in info_keys:
                desc = vcf.header.info["CSQ"].description
                if "Format: " in desc:
                    csq_part = desc.split("Format: ")[1]
                    csq_fields = {field.strip() for field in csq_part.split("|")}
        return info_keys, csq_fields

    def _norm_type_in_preds(self, preds: Iterable[Tuple[str, str, object]],
                           info_tags: Set[str]) -> Iterable[Tuple[str, str, object]]:
        """
        Normalize TYPE predicates to use the virtual TYPE field.
        """
        aliases = {"SNV": "snp", "SNP": "snp", "snv": "snp", 
                   "INDEL": "indel", "Indel": "indel", "InDel": "indel"}
        out = []
        for _, (field, op, value) in enumerate(preds):
            if field == "TYPE": # use virtual field instead of INFO field
                out.append((field, op, aliases.get(value)))
            else:
                out.append((field, op, value))
        return out

    def _norm_type_in_groups(self, groups: Iterable[Iterable[Tuple[str, str, object]]],
                            info_tags: Set[str]
                            ) -> Iterable[Iterable[Tuple[str, str, object]]]:
        return [
            self._norm_type_in_preds(grp, info_tags) for grp in (groups or [])
        ]
    
    def _fmt_pred(self, field: str, op: str, value: object) -> str:
        """
        Format a single predicate into a bcftools expression.
        """
        if isinstance(value, str):
            value = f'"{value}"'
        return f"({field} {op} {value})"

    def _build_and(self, preds: Iterable[Tuple[str, str, object]]) -> str:
        """
        Build an AND expression from predicates.
        """
        parts = []
        for p in preds or ():
            if not (isinstance(p, tuple) and len(p) == 3):
                raise ValueError("Each predicate must be a (field, op, value) tuple.")
            parts.append(self._fmt_pred(*p))
        return " && ".join(parts) if parts else ""

    def _build_or_of_ands(self, groups: Iterable[Iterable[Tuple[str, str, object]]]
                          ) -> str:
        """
        Build an OR expression of AND groups.
        """
        alts = []
        for grp in groups or ():
            s = self._build_and(grp)
            if s:
                alts.append(f"({s})")
        return " || ".join(alts) if alts else ""

    def _build_expr(self, include: Iterable[Tuple[str, str, object]] = (), 
                    exclude: Iterable[Tuple[str, str, object]] = (),
                    include_any: Iterable[Iterable[Tuple[str, str, object]]] = (),
                    exclude_any: Iterable[Iterable[Tuple[str, str, object]]] = ()
                    ) -> Tuple[str, str]:
        """
        Build a bcftools expression from include/exclude predicates.
        Returns (include_expr, exclude_expr).
        """
        inc_all = self._build_and(include)
        inc_any = self._build_or_of_ands(include_any)

        exc_all = self._build_and(exclude)
        exc_any = self._build_or_of_ands(exclude_any)

        inc_parts = " && ".join(
            part for part in (inc_all, f"({inc_any})" if inc_any else "")
            if part
        )
        exc_parts = " || ".join(
            part for part in (exc_all, f"({exc_any})" if exc_any else "")
            if part
        )
        return inc_parts, exc_parts

    def _split_vep_info(self, vcf_path: str, workdir: str, csq_tags: Set[str]) -> str:
        """
        Use the bcftools +split-vep plugin to split VEP annotations.
        """
        fd, tmp_vcf = tempfile.mkstemp(prefix="splitvep_", suffix=".vcf.gz", dir=workdir)
        os.close(fd)

        c_arg = ",".join(f"{t}:float" for t in sorted(csq_tags))

        run_command([
            "bcftools", "+split-vep",
            "-s", "worst", # annotation for variants will be the same
            "-c", c_arg,
            "-Oz", "-o", tmp_vcf,
            vcf_path
        ])
        run_command(["tabix", "-f", "-p", "vcf", tmp_vcf])
        return tmp_vcf
    
    def pass_filters(self, output_path: str, vcf_path: str = None, 
                     include: Iterable[Tuple[str, str, object]] = (),
                     exclude: Iterable[Tuple[str, str, object]] = (),
                     include_any: Iterable[Iterable[Tuple[str, str, object]]] = (),
                     exclude_any: Iterable[Iterable[Tuple[str, str, object]]] = (),
                     threads: int = 4, bgzip_index: bool = True,
                     tmpdir: str = ".",
                     keep_intermediate: bool = False):
        """
        Filter variants based on selected criteria.
        Only support filtering the INFO fields for now.
        
        Parameters
        ----------
        output_path : str
            Path to save the filtered VCF file.
        vcf_path : str, optional
            Path to the input VCF file. If None, use the VCF set in the instance.
        include : Iterable[Tuple[str, str, object]], optional
            Conditions to include variants, combined with AND.
        exclude : Iterable[Tuple[str, str, object]], optional
            Conditions to exclude variants, combined with AND.
        include_any : Iterable[Iterable[Tuple[str, str, object]]], optional
            Groups of conditions to include variants, combined with OR of ANDs.
        exclude_any : Iterable[Iterable[Tuple[str, str, object]]], optional
            Groups of conditions to exclude variants, combined with OR of ANDs.
        threads : int, optional
            Number of threads to use (default is 4).
        bgzip_index : bool, optional
            Whether to create a bgzip compressed and indexed VCF (default is True).
        tmpdir : str, optional
            Directory for temporary files (default is current directory).
        keep_intermediate : bool, optional
            Whether to keep intermediate files (default is False).
        """
        if vcf_path is None:
            vcf_path = self.vcf_path

        out_parent = ensure_parent_dir(output_path)
        intermediate = None
        vcf_in = vcf_path

        info_tags, csq_fields = self._info_keys_and_csq_fields(vcf_path)
        all_tags = info_tags.union(csq_fields or set())
        input_tags = {pred[0] for pred in include} | {pred[0] for pred in exclude} | \
                     {tag for grp in include_any for pred in grp for tag in (pred[0],)} | \
                     {tag for grp in exclude_any for pred in grp for tag in (pred[0],)}

        to_validate = {t for t in input_tags if t != "TYPE"}  # allow virtual TYPE
        if not to_validate.issubset(all_tags):
            missing = to_validate - all_tags
            raise ValueError(f"Input tags not found in VCF INFO fields: {missing}\n"
                             f"Available tags: {all_tags}")

        if csq_fields and any(tag in csq_fields for tag in input_tags):
            workdir = tmpdir if tmpdir else out_parent
            csq_tags = {tag for tag in input_tags if tag in csq_fields}
            intermediate = self._split_vep_info(vcf_path, workdir, csq_tags)
            vcf_in = intermediate

        inc_expr, exc_expr = self._build_expr(
            self._norm_type_in_preds(include, info_tags), 
            self._norm_type_in_preds(exclude, info_tags), 
            self._norm_type_in_groups(include_any, info_tags), 
            self._norm_type_in_groups(exclude_any, info_tags)
        )

        argv = ["bcftools", "view", "--threads", str(threads)]
        if inc_expr:
            argv += ["-i", inc_expr]
        if exc_expr:
            argv += ["-e", exc_expr]
        argv += [vcf_in, "-Oz", "-o", output_path]
        run_command(argv)

        if bgzip_index:
            run_command(["tabix", "-f", "-p", "vcf", output_path])
        if intermediate and not keep_intermediate:
            for p in (intermediate, f"{intermediate}.tbi"):
                try: os.remove(p)
                except FileNotFoundError: pass

    def bcftools_count(self, vcf_path: str) -> str:
        """
        Return the total number of variant record in vcf_path.
        """
        self.ensure_index(vcf_path)
        return run_command(["bcftools", "index", "--nrecords", vcf_path], capture=True).strip()

    def isec_unique(self, a_vcf_path: str, b_vcf_path: str, output_path: str):
        """
        Extract variants unique to a_vcf_path compared to b_vcf_path.
        """
        run_command(["bcftools", "isec", "-C",  "-w1", "-Oz", 
                     "-o", output_path, a_vcf_path, b_vcf_path])
        run_command(["tabix", "-f", "-p", "vcf", output_path]) 


class VCFBEDOverlapper:
    """
    A class to handle overlap between VCF and BED files.
    """
    def __init__(self, config: Dict[str, List[str]],
                 output_dir: str = ".", 
                 enhancer_dir: str = "data/enhancer/",
                 enhancer_prefix: str = "all_peaks_element_200bp_e",
                 ref: str = "/fs/cbsuhy01/storage/yz2676/data/DetailedProCapNet/wdir/reference/hg38.22XY.fasta"):
        """
        Args
        ----
        config : Dict[str, List[str]]
            A dictionary mapping labels to lists of BED directory names.
        output_dir : str
            Directory to save output files.
        enhancer_dir : str
            Directory containing enhancer BED files.
        enhancer_prefix : str
            Saved BED file name.
        ref : str
            Path to the reference genome file, keeps only the major contigs (22+XY)
        """
        self.config = config
        self.output_dir = output_dir
        self.enhancer_dir = enhancer_dir
        self.enhancer_prefix = enhancer_prefix
        self.ref = ref

    def _write_nonoverlap(self, bed: pd.DataFrame, label: str) -> str:
        """
        Merge the BED file to ensure non-overlapping regions.

        Args
        ----
        bed : pd.DataFrame
            Concatenated BED regions as a DataFrame.
        label : str
            Label for the BED file.
        """
        merged_bed = os.path.join(self.output_dir, label, f"{self.enhancer_prefix}.merged.bed")
        ensure_parent_dir(merged_bed)
        faidx_path = f"{self.ref}.fai"
        if not os.path.exists(faidx_path):
            run_command(["samtools", "faidx", self.ref])
        
        pybedtools.BedTool.from_dataframe(bed
                ).filter(
                    lambda x: x.chrom in {f"chr{str(i)}" for i in range(1, 23)} # autosomes only (stat model assumption)
                ).sort(faidx=faidx_path).merge().saveas(merged_bed)
        return merged_bed

    def _harmonize_chrom_names(self, vcf_path: str, bed_path: str) -> pd.DataFrame:
        """
        Harmonize chromosome names in the BED file to match those in the VCF.
        E.g., add 'chr' prefix if needed.
        """
        vcf_chroms = set()
        with pysam.VariantFile(vcf_path) as vcf:
            for rec in vcf.header.contigs.values():
                vcf_chroms.add(rec.name)

        def harmonize(chrom: str) -> str:
            if chrom in vcf_chroms:
                return chrom
            elif f"chr{chrom}" in vcf_chroms:
                return f"chr{chrom}"
            elif chrom.startswith("chr") and chrom[3:] in vcf_chroms:
                return chrom[3:]
            else:
                return chrom  # leave unchanged if no match found
        bed = pd.read_csv(bed_path, sep="\t", header=None)
        bed[0] = bed[0].apply(harmonize)
        return bed

    def _is_autosome(self, chrom: str) -> bool:
        """
        Check if the chromosome is an autosome (1-22).
        """
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        return chrom.isdigit() and 1 <= int(chrom) <= 22
    
    def map_variants_to_bed(self, vcf_path: str, names: List[str],
                            file: str, write_overlap_vcf: bool = False) -> None:
        """
        For each enhancer label,
            - Concatenate corresponding BED files.
            - Merge to ensure non-overlapping regions.
            - Find overlaps between VCF variants and BED regions.
            (pybedtools handled 1-based vs 0-based coordinates)
            - Save the overlapping BED regions to output files.
        
        Args
        ----
        vcf_path : str
            Path to the VCF file.
        names : List[str]
            List of enhancer names to analyze.
        file : str
            Output file name for saving overlapping bed regions.
        write_overlap_vcf : bool
            Whether to write the full overlapping VCF records (default is False).
        """
        for name in names:
            bed_names = self.config.get(name, [])
            if not bed_names:
                raise ValueError(f"Enhancer name '{name}' not found in config.")
            bed_region = pd.concat([pd.read_csv(
                os.path.join(self.enhancer_dir, bed_name, 
                             f"{self.enhancer_prefix}.bed"), sep="\t", header=None
            ) for bed_name in bed_names], ignore_index=True)
            sorted_bed_path = self._write_nonoverlap(bed_region, name) # keep only autosomes

            # Keep only autosomal chromosomes (statistical model assumption)
            bed_df = self._harmonize_chrom_names(vcf_path, sorted_bed_path)
            bed_bt = pybedtools.BedTool.from_dataframe(bed_df)

            vcf_bt = pybedtools.BedTool(vcf_path)
            if write_overlap_vcf:
                overlap_bed = bed_bt.intersect(vcf_bt, wa=True, wb=True)
            else:
                overlap_bed = bed_bt.intersect(vcf_bt, wa=True)
            overlap_bed.saveas(os.path.join(self.output_dir, name, file))

class StatisticalTests:
    """
    A class for statistical tests on variants.
    """
    def __init__(self):
        pass

    def fisher_exact_onesided(self, a: int, b: int, c: int, d: int
                              ) -> Tuple[float, float, float, float]:
        """
        Test whether the association between two categorical variables (A, B) is significant
        using one-sided Fisher's exact test. 

        Args
        ----
        a : int
            Count in group A with condition.
        b : int
            Total count in group A.
        c : int
            Count in group B with condition.
        d : int
            Total count in group B.
        """
        table = np.array([[a, b-a], [c, d-c]])
        odds_ratio, p_value = fisher_exact(table, alternative='greater')
        SE_log_OR = np.sqrt(1/a + 1/(b-a) + 1/c + 1/(d-c))

        return (odds_ratio, np.log(odds_ratio), SE_log_OR, p_value)

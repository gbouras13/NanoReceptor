# read and set the version

from .input_commands import get_input, instantiate_dirs, validate_fastq
from .processes import run_minimap2, keep_primary_supplementary_mappings_convert_sam, keep_primary_supplementary_mappings_convert_sam_direct, run_flagstat
from .post_processing import get_total_read_count, parse_bam, pivot_df
from .main import run
from .version import __version__

__all__ = ['get_input',
           'instantiate_dirs', 'validate_fastq',
           'run_minimap2', 'keep_primary_supplementary_mappings_convert_sam', 'keep_primary_supplementary_mappings_convert_sam_direct', 'run_flagstat',
           'get_total_read_count', 'parse_bam', 'pivot_df',
           'run',
            '__version__'
           ]
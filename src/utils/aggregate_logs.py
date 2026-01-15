#!/usr/bin/env python3
"""
Universal log aggregator for OP3_v2 pipelines.

Modes:
  1) Recursive: scan ./logs (or --log_dir) and aggregate per dataset/pipeline folder.
  2) Targeted: aggregate a specific directory with a given pipeline type.
"""

import os
import re
import argparse
from typing import List, Tuple, Optional, Dict, Set
from dataclasses import dataclass
from pathlib import Path

# Import logger from parsing_utils
from src.utils.parsing_utils import logger


# Constants
SENTINEL_JOB_ID = 999999999
OUT_EXTENSION = '.out'
ERR_EXTENSION = '.err'
LOG_EXTENSIONS = (OUT_EXTENSION, ERR_EXTENSION)  # Ordered: .out before .err

# Priority offsets (.out files come before .err files)
PRIORITY_OUT_FILE = 0.0
PRIORITY_ERR_FILE = 0.1

# Valid pipeline parameters
VALID_PARAMS = {'group_all_replicates', 'separate_replicates'}
FILTER_PREFIX = 'filter_min_cells_'


@dataclass
class PipelineConfig:
    """Configuration for pipeline-specific log priorities."""
    name: str
    main_pattern_template: str  # e.g., "^{name}_{dataset}(\.PID\d+)?\.(out|err)$"
    main_pattern_generic: str   # e.g., "^{name}_\w+(\.PID\d+)?\.(out|err)$"
    keyword_priorities: Dict[str, int]  # keywords -> priority level


# Pipeline configurations
PIPELINE_CONFIGS = {
    'pseudobulk': PipelineConfig(
        name='pseudobulk',
        main_pattern_template=r'^pseudobulk_{dataset}(\.pid\d+|\.\d+)?\.(out|err)$',
        main_pattern_generic=r'^pseudobulk_\w+(\.pid\d+|\.\d+)?\.(out|err)$',
        keyword_priorities={
            'download': 1,
            'process': 2,
            'unite': 3,
            'combine': 3,
        }
    ),
    'deg': PipelineConfig(
        name='deg',
        main_pattern_template=r'^deg_{dataset}(\.pid\d+|\.\d+)?\.(out|err)$',
        main_pattern_generic=r'^deg_\w+(\.pid\d+|\.\d+)?\.(out|err)$',
        keyword_priorities={
            'preprocessing': 1,
            'processing_pseudobulk': 1,
            'deg_analysis': 2,
            'analysis': 2,
            'aggregat': 3,
        }
    ),
}


def extract_job_id_and_index(filename: str) -> Tuple[int, int]:
    """
    Extract job ID and array index from filename for sorting.
    
    Args:
        filename: Log filename to parse
        
    Returns:
        Tuple of (job_id, array_index), or sentinel values if not found
    """
    # Try array job pattern first: name_JOBID_INDEX.ext
    match = re.search(r'[._](\d+)_(\d+)\.(out|err)$', filename)
    if match:
        return (int(match.group(1)), int(match.group(2)))
    
    # Try single job pattern: name_JOBID.ext
    match = re.search(r'[._](\d+)\.(out|err)$', filename)
    if match:
        return (int(match.group(1)), 0)
    
    return (SENTINEL_JOB_ID, SENTINEL_JOB_ID)


def get_file_priority(filename: str, config: PipelineConfig, dataset: Optional[str] = None) -> float:
    """
    Determine priority order for pipeline log files based on configuration.
    
    Priority order within each level:
    1. Files with PID (e.g., pseudobulk_dataset.PID123.out) - highest priority
    2. Files without any ID (e.g., pseudobulk_dataset.out)
    3. Files with job IDs (e.g., pseudobulk_dataset.123456.out) - lowest priority
    
    Args:
        filename: Log filename to prioritize
        config: Pipeline configuration
        dataset: Optional dataset name for more specific pattern matching
        
    Returns:
        Priority value (lower = higher priority)
    """
    basename = os.path.basename(filename).lower()
    is_out = basename.endswith(OUT_EXTENSION)
    offset = PRIORITY_OUT_FILE if is_out else PRIORITY_ERR_FILE
    
    # Check if this is the main pipeline file
    if dataset:
        pattern = config.main_pattern_template.format(dataset=re.escape(dataset))
    else:
        pattern = config.main_pattern_generic
    
    if re.match(pattern, basename):
        # Sub-priority within main pattern: PID files > no ID > job ID
        base_priority = 0.0
        
        # Check if file has PID pattern (e.g., .PID123456.)
        if re.search(r'\.pid\d+\.', basename):
            sub_priority = 0.0  # Highest priority
        # Check if file has job ID (numbers before the extension)
        elif re.search(r'\.\d+\.(out|err)$', basename):
            sub_priority = 0.4  # Lowest priority (job ID files)
        # No ID at all
        else:
            sub_priority = 0.2  # Second priority
        
        return base_priority + sub_priority + offset
    
    # Check keyword-based priorities
    for keyword, priority_level in config.keyword_priorities.items():
        if keyword in basename:
            return float(priority_level) + offset
    
    # Default priority for unrecognized files
    return 10.0 + offset


def extract_dataset_from_path(log_dir: str) -> Optional[str]:
    """
    Extract dataset name from log directory path.
    
    Args:
        log_dir: Path to log directory
        
    Returns:
        Dataset name if found, None otherwise
    """
    if not log_dir:
        return None
    
    normalized_path = os.path.normpath(log_dir)
    path_parts = normalized_path.split(os.sep)
    
    try:
        logs_idx = path_parts.index('logs')
        if logs_idx + 1 < len(path_parts):
            return path_parts[logs_idx + 1]
    except ValueError:
        pass
    
    return None


def sort_log_files(files: List[str], pipeline_type: str, log_dir: Optional[str] = None) -> List[str]:
    """
    Sort log files by priority, then by array index and job ID.
    
    Args:
        files: List of log file paths
        pipeline_type: Type of pipeline ('deg' or 'pseudobulk')
        log_dir: Optional log directory for dataset extraction
        
    Returns:
        Sorted list of file paths
        
    Raises:
        ValueError: If pipeline_type is invalid
    """
    if pipeline_type not in PIPELINE_CONFIGS:
        valid_types = ', '.join(PIPELINE_CONFIGS.keys())
        raise ValueError(f"Invalid pipeline_type: {pipeline_type}. Must be one of: {valid_types}")
    
    config = PIPELINE_CONFIGS[pipeline_type]
    dataset = extract_dataset_from_path(log_dir) if log_dir else None
    
    def sort_key(filepath: str) -> Tuple[float, int, int, str]:
        filename = os.path.basename(filepath)
        priority = get_file_priority(filename, config, dataset)
        job_id, array_idx = extract_job_id_and_index(filename)
        # Sort by: priority, array_idx, job_id, filename
        # This groups all _0 files together, then all _1 files, etc.
        return (priority, array_idx, job_id, filename)
    
    return sorted(files, key=sort_key)


def collect_logs(log_dir: str, pipeline_type: str) -> List[str]:
    """
    Collect .out/.err files under log_dir filtered by pipeline_type in filename.
    
    Args:
        log_dir: Directory to search for log files
        pipeline_type: Pipeline type to filter by
        
    Returns:
        List of matching log file paths
    """
    results = []
    pipeline_lower = pipeline_type.lower()
    
    for root, _, files in os.walk(log_dir):
        for filename in files:
            filename_lower = filename.lower()
            
            # Check file extension
            if not any(filename.endswith(ext) for ext in LOG_EXTENSIONS):
                continue
            
            # Check if pipeline type is in filename
            if pipeline_lower not in filename_lower:
                continue
            
            # Exclude files containing other pipeline types
            # e.g., exclude "deg_processing_pseudobulk" when collecting "pseudobulk" logs
            if pipeline_lower == 'pseudobulk' and 'deg' in filename_lower:
                continue
            if pipeline_lower == 'deg' and filename_lower.startswith('pseudobulk'):
                continue
            
            results.append(os.path.join(root, filename))
    
    return results


def pair_log_files(sorted_files: List[str]) -> List[str]:
    """
    Pair .out and .err files together, maintaining sort order.
    
    Args:
        sorted_files: Pre-sorted list of log files
        
    Returns:
        List with .out/.err pairs grouped together
    """
    processed: Set[str] = set()
    all_files: Set[str] = set(sorted_files)
    file_pairs: List[str] = []
    
    for filepath in sorted_files:
        if filepath in processed:
            continue
        
        # Get base path without extension
        base = filepath.rsplit('.', 1)[0]
        
        # Add both .out and .err if they exist
        for ext in LOG_EXTENSIONS:
            candidate = base + ext
            if candidate in all_files and candidate not in processed:
                file_pairs.append(candidate)
                processed.add(candidate)
    
    return file_pairs


def write_combined_log(
    files: List[str],
    output_path: str,
    pipeline_type: str,
    log_dir: str
) -> Optional[str]:
    """
    Write a combined log from the provided files.
    
    Args:
        files: List of log files to combine
        output_path: Path for the output combined log
        pipeline_type: Type of pipeline
        log_dir: Directory containing the logs
        
    Returns:
        Output path if successful, None if no files processed
    """
    if not files:
        return None
    
    try:
        sorted_files = sort_log_files(files, pipeline_type=pipeline_type, log_dir=log_dir)
        file_pairs = pair_log_files(sorted_files)
        
        total_lines = 0
        files_written = 0
        
        with open(output_path, 'w', encoding='utf-8') as outfile:
            for filepath in file_pairs:
                # Skip empty files
                if os.path.getsize(filepath) == 0:
                    continue
                
                filename = os.path.basename(filepath)
                
                try:
                    with open(filepath, 'r', encoding='utf-8', errors='replace') as infile:
                        for line in infile:
                            outfile.write(f"{filename}:{line.rstrip()}\n")
                            total_lines += 1
                    files_written += 1
                    
                except IOError as e:
                    logger.warning(f"Failed to read {filepath}: {e}")
                    continue
        
        logger.info(f"Created {output_path} ({total_lines} lines from {files_written} files)")
        return output_path
        
    except IOError as e:
        logger.error(f"Failed to write combined log to {output_path}: {e}")
        return None


def is_valid_deg_subdir(root: str, deg_root: str) -> bool:
    """
    Check if a directory is a valid DEG subdirectory for aggregation.
    
    Args:
        root: Directory path to check
        deg_root: Root DEG directory
        
    Returns:
        True if directory should be processed
    """
    rel = os.path.relpath(root, deg_root)
    parts = rel.split(os.sep)
    
    # Must have exactly 3 levels; expected layout: param/qc/filter
    if len(parts) != 3:
        return False
    
    param, _, filter_dir = parts
    
    # Validate parameter name
    if param not in VALID_PARAMS:
        return False
    
    # Validate filter directory
    if not filter_dir.startswith(FILTER_PREFIX):
        return False
    
    return True


def process_deg_logs(deg_root: str, dataset_name: str) -> List[str]:
    """
    Process DEG logs in the specified directory structure.
    
    Args:
        deg_root: Root directory for DEG logs
        dataset_name: Name of the dataset for filename
        
    Returns:
        List of created combined log paths
    """
    created = []
    
    for root, _, _ in os.walk(deg_root):
        if not is_valid_deg_subdir(root, deg_root):
            continue
        
        deg_files = collect_logs(root, 'deg')
        if deg_files:
            output_path = os.path.join(root, f'deg_{dataset_name}_combined.log')
            result = write_combined_log(deg_files, output_path, 'deg', root)
            if result:
                created.append(result)
    
    return created


def process_version_directory(version_dir: Path, dataset_name: str) -> List[str]:
    """
    Process logs for a single version directory (full or subsample).
    
    Args:
        version_dir: Path to version directory
        dataset_name: Name of the dataset
        
    Returns:
        List of created combined log paths
    """
    created = []
    version_dir_str = str(version_dir)
    
    # Process pseudobulk logs
    pb_files = collect_logs(version_dir_str, 'pseudobulk')
    if pb_files:
        output_path = version_dir / f"pseudobulk_{dataset_name}_combined.log"
        result = write_combined_log(pb_files, str(output_path), 'pseudobulk', version_dir_str)
        if result:
            created.append(result)
    
    # Process DEG logs
    deg_root = version_dir / 'deg'
    if deg_root.is_dir():
        created.extend(process_deg_logs(str(deg_root), dataset_name))
    
    return created


def aggregate_recursive(base_dir: str) -> List[str]:
    """
    Aggregate logs across datasets and pipelines from a logs root directory.
    
    Args:
        base_dir: Root directory containing logs
        
    Returns:
        List of created combined log paths
        
    Raises:
        ValueError: If base_dir doesn't exist
    """
    base_path = Path(base_dir).resolve()
    
    if not base_path.is_dir():
        raise ValueError(f"Directory not found: {base_dir}")
    
    created = []
    
    # Iterate through dataset directories
    dataset_dirs = sorted([d for d in base_path.iterdir() if d.is_dir()])
    
    for dataset_dir in dataset_dirs:
        dataset_name = dataset_dir.name
        
        # Process both full and subsample versions
        for version in ('full', 'subsample'):
            version_dir = dataset_dir / version
            if version_dir.is_dir():
                created.extend(process_version_directory(version_dir, dataset_name))
    
    return created


def aggregate_target(
    log_dir: str,
    pipeline_type: str,
    output_name: Optional[str] = None
) -> Optional[str]:
    """
    Aggregate logs for a specific directory and pipeline type.
    
    Args:
        log_dir: Directory containing logs to aggregate
        pipeline_type: Type of pipeline ('deg' or 'pseudobulk')
        output_name: Optional custom output filename
        
    Returns:
        Path to combined log if successful, None otherwise
        
    Raises:
        ValueError: If log_dir doesn't exist
    """
    log_path = Path(log_dir).resolve()
    
    if not log_path.is_dir():
        raise ValueError(f"Directory not found: {log_dir}")
    
    if output_name is None:
        output_name = 'combined.log'
    
    output_path = log_path / output_name
    files = collect_logs(str(log_path), pipeline_type)
    
    if not files:
        logger.warning(f"No log files found in {log_dir} for pipeline {pipeline_type}")
        return None
    
    return write_combined_log(files, str(output_path), pipeline_type, str(log_path))


def main() -> None:
    """Main entry point for the log aggregator."""
    parser = argparse.ArgumentParser(
        description='Log aggregator for OP3_v2 pipelines',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Recursive mode: aggregate all logs under ./logs
  %(prog)s --recursive
  
  # Recursive mode with custom directory
  %(prog)s --recursive --log_dir /path/to/logs
  
  # Target mode: aggregate specific directory
  %(prog)s --log_dir /path/to/logs/dataset1/full --pipeline_type pseudobulk
  
  # Target mode with custom output name
  %(prog)s --log_dir ./logs/deg/param1 --pipeline_type deg --output_name my_log.txt
        """
    )
    parser.add_argument(
        '--recursive',
        action='store_true',
        help='Aggregate logs under ./logs (or --log_dir) for all datasets/pipelines'
    )
    parser.add_argument(
        '--log_dir',
        type=str,
        default='./logs',
        help='Logs root (recursive mode) or target directory (non-recursive mode)'
    )
    parser.add_argument(
        '--pipeline_type',
        type=str,
        choices=list(PIPELINE_CONFIGS.keys()),
        help='Pipeline type for non-recursive mode'
    )
    parser.add_argument(
        '--output_name',
        type=str,
        default=None,
        help='Output filename for non-recursive mode'
    )
    
    args = parser.parse_args()
    
    try:
        if args.recursive:
            created = aggregate_recursive(args.log_dir)
            logger.info(f"\nCreated {len(created)} combined log files")
            for path in created:
                print(f"  - {path}")
        else:
            if not args.pipeline_type:
                parser.error("--pipeline_type is required when not using --recursive")
            
            result = aggregate_target(args.log_dir, args.pipeline_type, args.output_name)
            if result:
                print(f"Created: {result}")
            else:
                raise SystemExit(1)
                
    except (ValueError, IOError) as e:
        logger.error(str(e))
        raise SystemExit(1)


if __name__ == '__main__':
    main()

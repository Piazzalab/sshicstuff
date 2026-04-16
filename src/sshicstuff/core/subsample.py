"""

"""
import logging
import os
import re
import subprocess
import sys

logger = logging.getLogger(__name__)


def _check_gzip():
    """
    Check if gzip is installed and retrieve its version.
    """
    try:
        # Execute gzip with the --version argument to capture the version information
        result = subprocess.run(["gzip", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)
        # Parse the output to find the version number
        version_match = re.search(r"gzip (\d+\.\d+)", result.stdout)
        if version_match:
            version = version_match.group(1)
            logger.info("gzip version %s is installed.", version)
            return version
        else:
            logger.error("Unable to determine gzip version from the output.")
            return None
    except subprocess.CalledProcessError:
        logger.error("gzip is not installed or not functioning correctly. "
                      "Please install or fix gzip before running this function.")
        return None
    except FileNotFoundError as e:
        logger.error("gzip not found: %s", e)
        return None
    except subprocess.SubprocessError as e:
        logger.error("Subprocess error when checking gzip version: %s", e)
        return None


def _check_seqtk():
    """
    Check if seqtk is installed and retrieve its version.
    """
    try:
        # Execute seqtk without arguments to capture the usage information which contains the version
        result = subprocess.run("seqtk", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False, text=True)
        # Parse the output to find the version number
        version_match = re.search(r"Version: (\S+)", result.stdout)
        if version_match:
            version = version_match.group(1)
            logger.info("seqtk version %s is installed.", version)
            return version
        else:
            logger.error("Unable to determine seqtk version from the output.")
            sys.exit(1)
    except subprocess.CalledProcessError:
        logger.error("seqtk is not installed or not functioning correctly. "
                      "Please install or fix seqtk before running this function.")
        sys.exit(1)


def subsample_fastq(
    input_path: str,
    seed: int = 100,
    size: int = 4000000,
    compress: bool = True,
    force: bool = False,
):
    """
    Subsample and compress FASTQ file using seqtk.

    Parameters
    ----------
    input_path : str
        Path to the input FASTQ file.
    seed : int
        Seed for random number generator.
    size : int
        Number of reads to subsample randomly.
    compress : bool
        Whether to compress the output file (with gzip).
    force : bool
        Whether to overwrite the output file if it already exists.
    """

    _check_seqtk()

    # Determine the appropriate suffix for output file based on size
    if size >= 1000000000:
        suffix = f"sub{size // 1000000000}G"
    elif size >= 1000000:
        suffix = f"sub{size // 1000000}M"
    elif size >= 1000:
        suffix = f"sub{size // 1000}K"
    else:
        suffix = f"sub{size}"

    # Construct the output path
    pattern = r"(.+)(\.end[12]\.fastq)\.gz"
    match = re.match(pattern, input_path)
    if match:
        prefix = match.group(1)
        end_part = match.group(2)
        output_path = f"{prefix}_{suffix}{end_part}"
    else:
        logger.error("Input file does not match expected pattern")
        raise ValueError("Input file does not match expected pattern")

    if os.path.exists(output_path) and not force:
        logger.warning(
            "Output file %s already exists. skipping subsampling.", output_path
        )
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    # Log the start of subsampling
    logger.info(
        "Starting subsampling of %s to %s.gz with seed %d and size %d",
        input_path,
        output_path,
        seed,
        size,
    )

    # Run seqtk to subsample the FASTQ file
    seqtk_command = f"seqtk sample -s {seed} {input_path} {size} > {output_path}"
    try:
        subprocess.run(seqtk_command, shell=True, check=True)
        logger.info("Subsampling completed successfully for %s", output_path)
    except subprocess.CalledProcessError as e:
        logger.error("Error in subsampling %s: %s", input_path, e)
        raise

    # Optionally compress the output file
    if compress:
        _check_gzip()
        if os.path.exists(output_path + ".gz"):
            logger.warning(
                "Output file %s.gz already exists. removing it.", output_path
            )
            os.remove(output_path + ".gz")

        gzip_command = f"gzip {output_path}"
        logger.info("Compressing %s", output_path)
        try:
            subprocess.run(gzip_command, shell=True, check=True)
            output_path += ".gz"
            logger.info("Compression completed successfully for %s", output_path)
        except subprocess.CalledProcessError as e:
            logger.error("Error in compressing %s: %s", output_path, e)
            raise
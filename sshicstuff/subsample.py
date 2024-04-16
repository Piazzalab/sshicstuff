import os.path
import subprocess
import logging
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')



def check_seqtk():
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
            logging.info(f"seqtk version {version} is installed.")
            return version
        else:
            logging.error("Unable to determine seqtk version from the output.")
            return None
    except subprocess.CalledProcessError:
        logging.error("seqtk is not installed or not functioning correctly. "
                      "Please install or fix seqtk before running this function.")
        return None


def check_gzip():
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
            logging.info(f"gzip version {version} is installed.")
            return version
        else:
            logging.error("Unable to determine gzip version from the output.")
            return None
    except subprocess.CalledProcessError:
        logging.error("gzip is not installed or not functioning correctly. "
                      "Please install or fix gzip before running this function.")
        return None
    except Exception as e:
        logging.error(f"Unexpected error when checking gzip version: {e}")
        return None


def subsample(
        input_path: str,
        seed: int = 100,
        size: int = 4000000,
        compress: bool = True,
        force: bool = False
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

    check_seqtk()

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
        logging.error("Input file does not match expected pattern")
        raise ValueError("Input file does not match expected pattern")

    if os.path.exists(output_path) and not force:
        logging.warning(f"Output file {output_path} already exists. skipping subsampling.")
        logging.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    # Log the start of subsampling
    logging.info(f"Starting subsampling of {input_path} to {output_path}.gz with seed {seed} and size {size}")

    # Run seqtk to subsample the FASTQ file
    seqtk_command = f"seqtk sample -s {seed} {input_path} {size} > {output_path}"
    try:
        subprocess.run(seqtk_command, shell=True, check=True)
        logging.info(f"Subsampling completed successfully for {output_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in subsampling {input_path}: {e}")
        raise

    # Optionally compress the output file
    if compress:
        check_gzip()
        if os.path.exists(output_path + ".gz"):
            logging.warning(f"Output file {output_path}.gz already exists. removing it.")
            os.remove(output_path + ".gz")

        gzip_command = f"gzip {output_path}"
        logging.info(f"Compressing {output_path}")
        try:
            subprocess.run(gzip_command, shell=True, check=True)
            output_path += ".gz"
            logging.info(f"Compression completed successfully for {output_path}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in compressing {output_path}: {e}")
            raise


if __name__ == "__main__":

    """
    Example usage
    
    python3 ./main.py subsample -s 100 -z 1000000 -c -F ./data/test.fastq
    """

    pass



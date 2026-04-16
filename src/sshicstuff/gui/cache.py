"""
Session-aware file cache utilities for the Dash GUI.
"""

from __future__ import annotations

import base64
import os
from uuid import uuid4

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

#: Root cache directory.  Override via the ``SSHICSTUFF_CACHE_DIR`` env var.
__CACHE_DIR__: str = os.environ.get("SSHICSTUFF_CACHE_DIR", "/tmp/sshicstuff_cache")

#: Unique identifier for this process instance.  Separates parallel server
#: processes from one another under the same base cache directory.
APP_INSTANCE_ID: str = os.environ.get("SSHICSTUFF_APP_INSTANCE_ID", str(uuid4()))

os.makedirs(__CACHE_DIR__, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def uploaded_files_cache(cache_dir: str) -> list[str]:
    """Return the list of file names present in *cache_dir*."""
    return [
        f for f in os.listdir(cache_dir)
        if os.path.isfile(os.path.join(cache_dir, f))
    ]


def save_file_cache(name: str, content: str, cache_dir: str) -> None:
    """Decode a Dash upload payload and write it to *cache_dir*/*name*.

    Parameters
    ----------
    name:
        Original file name.
    content:
        Base-64 data URI string supplied by ``dcc.Upload``.
    cache_dir:
        Destination directory (must already exist).
    """
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(cache_dir, name), "wb") as fh:
        fh.write(base64.decodebytes(data))
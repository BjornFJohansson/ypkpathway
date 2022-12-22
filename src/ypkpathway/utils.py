#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""
import shutil
from pathlib import Path
from typing import Iterable


def _copy2(src, dst, *, follow_symlinks=True):
    """docstring."""
    try:
        dst = shutil.copy2(src, dst)
    except shutil.SameFileError:
        pass
    return Path(dst)


def _find_path_to_file(folders: Iterable[Path],
                       fn: str):
    """docstring."""
    foundpath = None
    for folder in folders:
        path = folder/fn
        if path.is_file():
            foundpath = path
            break
    return foundpath

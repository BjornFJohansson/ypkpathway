#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""
import os
import shutil
from pathlib import Path
from typing import Iterable
import copy
from itertools import chain

class Pth(type(Path())):
    """docstring."""

    red = "\u001b[31m"
    blue = "\u001b[34m"
    green = "\u001b[32m"
    plain = "\x1b[0m"

    def __new__(cls,
                path: str | Path,
                datadirs: list[Path] = None,
                *args, **kwargs):
        """docstring."""
        self = super().__new__(cls, path)
        dd = datadirs or []
        self.datadirs = [Path(d) for d in dd]
        self.origin = None
        self.others = []
        self.status = chr(32)
        return self
    
    def __truediv__(self, other):
        result = Pth(super().__truediv__(other), self.datadirs)
        return result

    def __rtruediv__(self, other):
        result = Pth(super().__rtruediv__(other), self.datadirs)
        return result

    def with_suffix(self, suffix):
        result = Pth(super().with_suffix(suffix), self.datadirs)
        return result

    @property
    def parent(self):
        result = Pth(super().parent, self.datadirs)
        return result
    
    def find_in_dirs(self):
        """docstring."""
        self.origin, *self.others = list(
            chain(*[d.rglob(self.name) for d in self.datadirs])) or [None]
        self.status = "f" if self.origin else "-"
        return bool(self.origin)

    def copy_to_dir(self):
        """docstring."""
        if not self.origin:
            self.find_in_dirs()
        try:
            dst = shutil.copy2(self.origin, self.parent)
        except (shutil.SameFileError, TypeError):
            dst = None
        self.status = "+" if dst else self.status
        return dst

    def link_in_dir(self):
        """docstring."""
        src = self.find_in_dirs()
        try:
            os.symlink(src, self)
        except TypeError:
            pass
        return None

    def __repr__(self):
        """docstring."""
        color = {"-": self.red,
                 "+": self.green,
                 "f": self.blue,
                 chr(32): self.plain}[self.status]
        return f"{color}{super().__repr__()}({self.status}){self.plain}"

if __name__ == "__main__":

    datadirs = ("/home/bjorn/Desktop/mec@githb/YeastPathwayKit/sequences",
                "/home/bjorn/Desktop/mec@githb/YeastPathwayKitPrivate/sequences",)
    
    self = Pth("xyzzy.py", datadirs)
    
    a = self/Path("qwerty")
    b = "self"/self
    
    c = a.with_suffix(".gb")



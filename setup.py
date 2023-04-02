# -*- coding: utf-8 -*-

import sys

from setuptools import setup


if __name__ == "__main__":
    # clean --all does not remove extensions automatically
    if "clean" in sys.argv and "--all" in sys.argv:
        import pathlib
        import shutil

        # Delete complete build directory
        path = pathlib.Path.cwd() / "build"
        shutil.rmtree(str(path), ignore_errors=True)
    else:
        setup()

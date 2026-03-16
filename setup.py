import setuptools

import os
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Extension
from setuptools.command.build_ext import build_ext

class RustPyo3Extension(Extension):
    def __init__(self, name: str, manifest_path: str, crate_name: str) -> None:
        super().__init__(name, sources=[])
        self.manifest_path = os.fspath(Path(manifest_path).resolve())
        self.crate_name = crate_name


class RustBuild(build_ext):
    def build_extension(self, ext: Extension) -> None:
        if not isinstance(ext, RustPyo3Extension):
            raise RuntimeError(f"Unsupported extension type: {type(ext)!r}")
        self._build_rust_extension(ext)

    def _build_rust_extension(self, ext: RustPyo3Extension) -> None:
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        ext_fullpath.parent.mkdir(parents=True, exist_ok=True)

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug

        cargo_target_dir = (Path(self.build_temp) / "rust_cargo_target").resolve()
        cargo_target_dir.mkdir(parents=True, exist_ok=True)
        cargo_tmp_dir = cargo_target_dir / "tmp"
        cargo_tmp_dir.mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        env["PYO3_PYTHON"] = sys.executable
        env["CARGO_TARGET_DIR"] = str(cargo_target_dir)
        env["TMPDIR"] = str(cargo_tmp_dir)

        cargo_cmd = [
            "cargo",
            "build",
            "--manifest-path",
            ext.manifest_path,
            "--features",
            "python",
        ]
        if not debug:
            cargo_cmd.append("--release")

        subprocess.run(cargo_cmd, env=env, check=True)

        profile = "debug" if debug else "release"
        if sys.platform.startswith("linux"):
            lib_name = f"lib{ext.crate_name}.so"
        elif sys.platform.startswith("darwin"):
            lib_name = f"lib{ext.crate_name}.dylib"
        elif sys.platform == "win32":
            lib_name = f"{ext.crate_name}.dll"
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}")

        built_lib = cargo_target_dir / profile / lib_name
        if not built_lib.exists():
            raise RuntimeError(f"Rust extension build output not found: {built_lib}")

        shutil.copyfile(built_lib, ext_fullpath)

setuptools.setup(
    name="pythonknot",
    version="0.1.16",
    author="yjianzhu",
    author_email="yjianzhu@mail.ustc.edu.cn",
    description="pythonknot for knot theory calculation",
    long_description="Python Knot is a Python package for creating and manipulating knots and links, calculating knot invariants, and more.",
    url="",
    ext_modules=[
        RustPyo3Extension(
            name="pythonknot.alexander_poly",
            manifest_path="vendor/rust_knot/Cargo.toml",
            crate_name="rust_knot",
        ),
    ],
    cmdclass={"build_ext": RustBuild},

    packages=setuptools.find_packages(where='src'), 
    package_dir={"": "src"},
    package_data={
        "pythonknot": ["data/table_knot_Alexander_polynomial.txt"],
    },
    zip_safe=False,
    install_requires=[
        'numpy',
        'scipy',
        'rmsd'
    ],
    include_package_data=True,
)

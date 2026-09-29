#!/usr/bin/env python3
"""Test the core Makefile HTSLIB=bundled|external selection with a tiny HTSlib fixture."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
MAKEFILE = REPO / "core/legacy/source/Makefile"


class HtslibBuildTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="star-htslib-build-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.prefix = self.root / "external"
        headers = self.prefix / "include/htslib"
        headers.mkdir(parents=True)
        (headers / "khash.h").write_text("#define TEST_EXTERNAL_HTSLIB 1\n")
        (headers / "sam.h").write_text(
            'extern "C" const char *hts_version(void);\n'
        )
        src = self.root / "hts.c"
        src.write_text('const char *hts_version(void) { return "fixture"; }\n')
        libdir = self.prefix / "lib"
        libdir.mkdir()
        obj = self.root / "hts.o"
        subprocess.run(["cc", "-c", str(src), "-o", str(obj)], check=True)
        self.lib = libdir / "libhts.a"
        subprocess.run(["ar", "rcs", str(self.lib), str(obj)], check=True)
        pcdir = libdir / "pkgconfig"
        pcdir.mkdir()
        (pcdir / "htslib.pc").write_text(
            "Name: htslib\nDescription: build discovery fixture\nVersion: 1.0\n"
            f"Cflags: -I{self.prefix}/include\nLibs: {self.lib}\n"
        )
        (self.root / "parametersDefault").write_text("fixture default\n")
        (self.root / "test.cpp").write_text(
            "#include <htslib/khash.h>\n"
            "#ifndef TEST_EXTERNAL_HTSLIB\n#error wrong headers\n#endif\n"
        )
        self.env = dict(os.environ, PKG_CONFIG_LIBDIR=str(pcdir), PKG_CONFIG_PATH="")
        for name in ("CPATH", "CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH", "MAKEFLAGS"):
            self.env.pop(name, None)

    def make(self, *args, ok=True, target="Depend.list", selection="external"):
        selection_args = [f"HTSLIB={selection}"] if selection else []
        result = subprocess.run(
            ["make", "--no-print-directory", "-f", str(MAKEFILE),
             f"ROOT_DIR={REPO}", "SOURCES=test.cpp", *selection_args,
             target, *args],
            cwd=self.root, env=self.env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30,
        )
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def test_pkg_config_prefix_and_cached_scan(self):
        self.make()
        dep = self.root / "Depend.list"
        self.assertIn(str(self.prefix / "include/htslib/khash.h"), dep.read_text())
        before = dep.stat().st_mtime_ns
        self.make()
        self.assertEqual(dep.stat().st_mtime_ns, before)
        self.assertFalse((self.root / "htslib").exists())

    def test_cppflags_and_extra_flags_reach_scan(self):
        for name in ("CPPFLAGS", "CXXFLAGSextra"):
            with self.subTest(name=name):
                self.make("HTSLIB_CFLAGS=", f"HTSLIB_LIBS={self.lib}",
                          f"{name}=-I{self.prefix}/include", "-W", "test.cpp")

    def test_custom_prefix_compiles_after_dependency_scan(self):
        # Match Carl's indirect include, with khash.h only below include/htslib.
        (self.root / "ParametersSolo.h").write_text('#include "htslib/khash.h"\n')
        (self.root / "test.cpp").write_text(
            '#include "ParametersSolo.h"\n'
            '#ifndef TEST_EXTERNAL_HTSLIB\n#error wrong headers\n#endif\n'
            'int fixture() { return TEST_EXTERNAL_HTSLIB; }\n'
        )
        self.make(target="test.o")
        self.assertGreater((self.root / "test.o").stat().st_size, 0)
        self.assertIn(str(self.prefix / "include/htslib/khash.h"),
                      (self.root / "Depend.list").read_text())

    def test_missing_headers_are_actionable(self):
        result = self.make("HTSLIB_CFLAGS=-nostdinc", ok=False)
        self.assertIn("HTSLIB=external requires installed HTSlib headers and library", result.stdout)
        self.assertIn("libhts-dev", result.stdout)
        self.assertIn("HTSLIB=bundled", result.stdout)
        self.assertFalse((self.root / "Depend.list").exists())

    def test_missing_library_is_actionable(self):
        result = self.make(f"HTSLIB_LIBS={self.root}/absent.so", ok=False)
        self.assertIn("HTSLIB=external requires installed HTSlib headers and library", result.stdout)
        self.assertFalse((self.root / "Depend.list").exists())

    def test_unknown_selection_is_rejected(self):
        result = self.make(ok=False, selection="system")
        self.assertIn("HTSLIB must be 'bundled' or 'external'", result.stdout)
        self.assertFalse((self.root / "Depend.list").exists())

    def test_switching_selection_rescans_dependencies(self):
        bundled = self.root / "htslib/htslib"
        bundled.mkdir(parents=True)
        (bundled / "khash.h").write_text("#define TEST_EXTERNAL_HTSLIB 1\n")
        (bundled.parent / "libhts.a").touch()
        self.make()
        dep = self.root / "Depend.list"
        self.assertIn(str(self.prefix / "include/htslib/khash.h"), dep.read_text())
        self.make("PKG_CONFIG=false", "CPPFLAGS=-nostdinc", selection="bundled")
        self.assertIn("htslib/htslib/khash.h", dep.read_text())
        self.assertNotIn(str(self.prefix), dep.read_text())

    def test_scan_failure_keeps_previous_dependencies(self):
        self.make()
        dep = self.root / "Depend.list"
        good = dep.read_bytes()
        (self.root / "bad.cpp").write_text('#include "nonexistent-fixture.h"\n')
        self.make("SOURCES=test.cpp bad.cpp", "-W", "test.cpp", ok=False)
        self.assertEqual(dep.read_bytes(), good)
        self.assertFalse(list(self.root.glob("Depend.list.*")))

    def test_default_uses_bundled_headers_without_external_dependencies(self):
        bundled = self.root / "htslib/htslib"
        bundled.mkdir(parents=True)
        (bundled / "khash.h").write_text("#define TEST_EXTERNAL_HTSLIB 1\n")
        (bundled.parent / "libhts.a").touch()
        self.make("PKG_CONFIG=false", "CPPFLAGS=-nostdinc", selection=None)
        self.assertIn("htslib/htslib/khash.h", (self.root / "Depend.list").read_text())


if __name__ == "__main__":
    unittest.main()

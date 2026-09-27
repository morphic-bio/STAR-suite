#!/usr/bin/env python3
"""Test real core Makefile rules using a tiny external HTSlib fixture."""
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
        self.chromap = self.root / "Chromap-suite"
        (self.chromap / "src").mkdir(parents=True)
        (self.chromap / "Makefile").write_text("# fixture, no build needed\n")
        (self.chromap / "src/libchromap.h").write_text("// fixture\n")
        (self.root / "parametersDefault").write_text("fixture default\n")
        (self.root / "test.cpp").write_text(
            "#include <htslib/khash.h>\n"
            "#ifndef TEST_EXTERNAL_HTSLIB\n#error wrong headers\n#endif\n"
        )
        self.env = dict(os.environ, PKG_CONFIG_LIBDIR=str(pcdir), PKG_CONFIG_PATH="")
        for name in ("CPATH", "CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH", "MAKEFLAGS"):
            self.env.pop(name, None)

    def make(self, *args, ok=True, target="Depend.list"):
        result = subprocess.run(
            ["make", "--no-print-directory", "-f", str(MAKEFILE),
             f"ROOT_DIR={REPO}", f"CHROMAP_SUITE_DIR={self.chromap}",
             "SOURCES=test.cpp", target, *args],
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
                self.make("CHROMAP_HTSLIB_CFLAGS=", f"CHROMAP_SYS_HTS={self.lib}",
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
        result = self.make("CHROMAP_HTSLIB_CFLAGS=-nostdinc", ok=False)
        self.assertIn("requires external HTSlib headers and library", result.stdout)
        self.assertIn("libhts-dev", result.stdout)
        self.assertIn("core-portable", result.stdout)
        self.assertFalse((self.root / "Depend.list").exists())

    def test_missing_library_is_actionable(self):
        result = self.make(f"CHROMAP_SYS_HTS={self.root}/absent.so", ok=False)
        self.assertIn("requires external HTSlib headers and library", result.stdout)
        self.assertFalse((self.root / "Depend.list").exists())

    def test_missing_chromap_is_actionable(self):
        result = self.make(f"CHROMAP_SUITE_DIR={self.root}/absent",
                           "CHROMAP_SOURCE_REQUIRED=1", ok=False)
        self.assertIn("Chromap-suite source not found", result.stdout)
        self.assertFalse((self.root / "Depend.list").exists())

    def test_non_chromap_partial_target_does_not_require_chromap_checkout(self):
        self.make(f"CHROMAP_SUITE_DIR={self.root}/absent")

    def test_scan_failure_keeps_previous_dependencies(self):
        self.make()
        dep = self.root / "Depend.list"
        good = dep.read_bytes()
        (self.root / "bad.cpp").write_text('#include "nonexistent-fixture.h"\n')
        self.make("SOURCES=test.cpp bad.cpp", "-W", "test.cpp", ok=False)
        self.assertEqual(dep.read_bytes(), good)
        self.assertFalse(list(self.root.glob("Depend.list.*")))

    def test_portable_uses_bundled_headers_without_external_dependencies(self):
        bundled = self.root / "htslib/htslib"
        bundled.mkdir(parents=True)
        (bundled / "khash.h").write_text("#define TEST_EXTERNAL_HTSLIB 1\n")
        (bundled.parent / "libhts.a").touch()
        self.make("WITH_CHROMAP=0", "PKG_CONFIG=false",
                  f"CHROMAP_SUITE_DIR={self.root}/absent",
                  "CPPFLAGS=-nostdinc")
        self.assertIn("htslib/htslib/khash.h", (self.root / "Depend.list").read_text())


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Guard STAR's library builds against compile-time and host-CPU inputs."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ReproducibleLibraryBuildTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="star-library-build-test-")
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.env = dict(os.environ, LC_ALL="C", TZ="UTC")
        for name in ("CXXFLAGS", "CPPFLAGS", "CPATH", "CPLUS_INCLUDE_PATH",
                     "MAKEFLAGS", "MFLAGS", "MAKEOVERRIDES"):
            self.env.pop(name, None)

    def run_command(self, argv, **kwargs):
        result = subprocess.run(argv, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=60, **kwargs)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def test_pcg_strict_build_and_seeded_sequence_across_paths_and_epochs(self):
        # This sequence is frozen from the unmodified vendored PCG header.
        # STAR's callers pass explicit seeds; the arbitrary-seed helper must
        # remain independent of build location/time without altering that RNG.
        program = r'''
#include "pcg_random.hpp"
#include <cstdint>
#include <iostream>
int main() {
    constexpr auto a = pcg_extras::static_arbitrary_seed<uint32_t>::value;
    constexpr auto b = pcg_extras::static_arbitrary_seed<uint64_t>::value;
    std::cout << a << " " << b << "\n";
    pcg32 rng(42u, 54u);
    for (int i = 0; i < 6; ++i) std::cout << std::hex << rng() << "\n";
}
'''
        outputs = []
        for directory, epoch in (("first", "0"), ("other/location", "946684800")):
            build = self.work / directory
            build.mkdir(parents=True)
            for name in ("pcg_random.hpp", "pcg_extras.hpp", "pcg_uint128.hpp"):
                shutil.copy2(ROOT / "core/features/libscrna/include" / name, build)
            source = build / "probe.cpp"
            source.write_text(program)
            executable = build / "probe"
            env = dict(self.env, SOURCE_DATE_EPOCH=epoch)
            self.run_command(["g++", "-std=c++11", "-O2", "-Werror=date-time",
                              str(source), "-o", str(executable)], env=env)
            output = self.run_command([str(executable)], env=env)
            self.assertEqual(output.splitlines()[1:], [
                "a15c02b7", "7b47f409", "ba1d3330", "83d2f293", "bfa4784b", "cbed606e"])
            outputs.append(output)
        self.assertEqual(outputs[0], outputs[1])

    def test_libem_compile_commands_do_not_detect_host_cpu(self):
        makefile = ROOT / "core/features/vbem/source/libem/Makefile"
        # libem's own sources normally resolve from its working directory;
        # copying them keeps this dry run independent of existing objects.
        for source in makefile.parent.glob("*.cpp"):
            shutil.copy2(source, self.work)
        for flags in (None, "-pipe -fno-omit-frame-pointer"):
            with self.subTest(flags=flags):
                command = ["make", "--no-print-directory", "-n", "-B",
                           "-f", str(makefile), f"ROOT_DIR={ROOT}", "CXX=g++", "libem.a"]
                if flags is not None:
                    command.append("CXXFLAGS=" + flags)
                output = self.run_command(command, cwd=self.work, env=self.env)
                compilations = [shlex.split(line) for line in output.splitlines()
                                if line.startswith("g++ ")]
                self.assertGreater(len(compilations), 0, output)
                for argv in compilations:
                    for flag in argv:
                        self.assertFalse(flag.startswith(("-march=", "-mtune=", "-mcpu="))
                                         and flag.split("=", 1)[1] == "native", argv)
                    self.assertIn("-std=c++11", argv)
                    self.assertIn("-fopenmp", argv)
                    if flags is not None:
                        self.assertIn("-fno-omit-frame-pointer", argv)


if __name__ == "__main__":
    unittest.main()

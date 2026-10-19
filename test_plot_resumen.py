"""Pruebas del agrupamiento y del promedio de presiones."""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from plot_resumen import average_summaries, find_summaries, main, plot_averages, read_summary


class SummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def write(self, name, rows):
        path = self.directory / name
        path.write_text("# box_length pressure\n" + rows)
        return path

    def test_finds_only_matching_size_and_coordination_in_numeric_run_order(self):
        names = [
            "network_4096_4.300_0010_summary.dat",
            "network_4096_4.300_0002_summary.dat",
            "network_1024_4.300_0000_summary.dat",
            "network_4096_4.310_0000_summary.dat",
            "network_4096_4.300_0000_step_0001.dat",
            "network_4096_4.300_mean.dat",
        ]
        for name in names:
            (self.directory / name).touch()
        found = find_summaries(self.directory, 4096, 4.3)
        self.assertEqual([path.name for path in found], [names[1], names[0]])

    def test_accepts_legacy_names_and_prefers_summary_for_the_same_run(self):
        for name in ["network_4096_4.300_0000_resumen.dat",
                     "network_4096_4.300_0000_summary.dat",
                     "network_4096_4.300_0001_resumen.dat"]:
            (self.directory / name).touch()
        found = find_summaries(self.directory, 4096, 4.3)
        self.assertEqual([path.name for path in found], [
            "network_4096_4.300_0000_summary.dat",
            "network_4096_4.300_0001_resumen.dat",
        ])

    def test_reads_single_row_as_strain_and_preserves_pressure_sign(self):
        path = self.write("one.dat", "1.2 -0.5\n")
        strain, pressure = read_summary(path)
        np.testing.assert_allclose(strain, [0.2])
        np.testing.assert_allclose(pressure, [-0.5])

    def test_equal_weight_average_and_sample_standard_deviation(self):
        first = self.write("first.dat", "1.1 -2\n1.2 -4\n")
        second = self.write("second.dat", "1.1 -4\n1.2 -8\n")
        strain, mean, std = average_summaries([first, second])
        np.testing.assert_allclose(strain, [0.1, 0.2])
        np.testing.assert_allclose(mean, [-3, -6])
        np.testing.assert_allclose(std, [np.sqrt(2), np.sqrt(8)])

    def test_single_simulation_is_unchanged_without_a_variation_band(self):
        path = self.write("one.dat", "1.1 -2\n1.2 -4\n")
        _, mean, std = average_summaries([path])
        np.testing.assert_allclose(mean, [-2, -4])
        np.testing.assert_allclose(std, [0, 0])

    def test_rejects_different_strain_grids_or_point_counts(self):
        first = self.write("first.dat", "1.1 -2\n1.2 -4\n")
        for rows in ["1.1 -4\n1.3 -8\n", "1.1 -4\n"]:
            with self.subTest(rows=rows):
                second = self.write("second.dat", rows)
                with self.assertRaisesRegex(ValueError, "strain no coinciden"):
                    average_summaries([first, second])

    def test_rejects_invalid_summaries(self):
        for rows in ["", "1.0 0\n", "1.1 nan\n", "inf -1\n", "1.1 -1 3\n",
                     "1.1 -1\n1.1 -2\n", "1.2 -2\n1.1 -1\n"]:
            with self.subTest(rows=rows):
                path = self.write("invalid.dat", rows)
                with self.assertRaises(ValueError):
                    read_summary(path)
        with self.assertRaises(ValueError):
            average_summaries([])

    def test_cli_keeps_groups_separate_and_exports_their_means(self):
        from unittest.mock import patch

        self.write("network_1024_4.300_0000_summary.dat", "1.1 -2\n1.2 -4\n")
        self.write("network_1024_4.300_0001_summary.dat", "1.1 -4\n1.2 -8\n")
        self.write("network_4096_4.300_0000_summary.dat", "1.1 -10\n1.2 -20\n")
        destination = self.directory / "plots" / "average.png"
        argv = ["plot_resumen.py", "--sizes", "1024", "4096", "--z", "4.3",
                "--output-dir", str(self.directory), "--save", str(destination), "--dpi", "60"]
        with patch("sys.argv", argv):
            main()
        self.assertTrue(destination.is_file())
        small = np.loadtxt(destination.parent / "network_1024_4.300_mean.dat")
        large = np.loadtxt(destination.parent / "network_4096_4.300_mean.dat")
        np.testing.assert_allclose(small[:, 1], [-3, -6])
        np.testing.assert_allclose(large[:, 1], [-10, -20])
        np.testing.assert_allclose(small[:, 3], [2, 2])
        np.testing.assert_allclose(large[:, 3], [1, 1])
        self.assertEqual(small.shape, (2, 5))
        np.testing.assert_allclose(small[:, 4], [3.3, 7.2])
        np.testing.assert_allclose(large[:, 4], [11, 24])
        header = (destination.parent / "network_1024_4.300_mean.dat").read_text().splitlines()[0]
        self.assertEqual(header, "# strain pressure_mean pressure_std simulations force_mean")

    def test_default_destination_is_project_output_processed_even_with_custom_input(self):
        from unittest.mock import patch

        self.write("network_4096_4.300_0000_summary.dat", "1.1 -2\n1.2 -4\n")
        project = self.directory / "project"
        argv = ["plot_resumen.py", "--sizes", "4096", "--z", "4.3",
                "--output-dir", str(self.directory), "--dpi", "60"]
        with patch("sys.argv", argv), patch("plot_resumen.__file__", str(project / "plot_resumen.py")):
            main()
        destination = project / "output" / "processed"
        self.assertEqual(
            {path.name for path in destination.iterdir()},
            {"force_vs_strain.png", "network_4096_4.300_mean.dat"},
        )
        self.assertGreater((destination / "force_vs_strain.png").stat().st_size, 0)
        table = np.loadtxt(destination / "network_4096_4.300_mean.dat")
        np.testing.assert_allclose(table[:, 0], [0.1, 0.2])
        np.testing.assert_allclose(table[:, 1], [-2, -4])
        np.testing.assert_allclose(table[:, 4], [2.2, 4.8])
        self.assertFalse((self.directory / "plots").exists())

    def test_plot_uses_force_mean_and_scales_its_variation_band(self):
        from unittest.mock import patch
        from matplotlib.figure import Figure

        mean = np.array([-2.0, -4.0])
        groups = [(4096, 4.3, [Path("one"), Path("two")],
                   np.array([0.1, 0.2]), mean, np.array([0.2, 0.4]))]
        with patch.object(Figure, "savefig", autospec=True) as save:
            plot_averages(groups, self.directory / "plot.png")
            axes = save.call_args.args[0].axes[0]
            self.assertEqual(axes.get_xscale(), "log")
            self.assertEqual(axes.get_yscale(), "log")
            self.assertEqual(axes.get_ylabel(), "Mean tensile force")
            np.testing.assert_allclose(axes.lines[0].get_ydata(), [2.2, 4.8])
            vertices = axes.collections[0].get_paths()[0].vertices
            self.assertAlmostEqual(vertices[:, 1].min(), 1.98)
            self.assertAlmostEqual(vertices[:, 1].max(), 5.28)
        np.testing.assert_allclose(mean, [-2, -4])

    def test_nonpositive_values_leave_gaps_instead_of_using_absolute_values(self):
        from io import StringIO
        from unittest.mock import patch
        from matplotlib.figure import Figure

        groups = [(4096, 4.3, [Path("one")], np.array([0.1, 0.2, 0.3, 0.4]),
                   np.array([-2.0, 0.0, 1.0, -4.0]), np.zeros(4))]
        with patch.object(Figure, "savefig", autospec=True) as save, patch("sys.stderr", new_callable=StringIO) as warnings:
            plot_averages(groups, self.directory / "plot.png")
            axes = save.call_args.args[0].axes[0]
            y = axes.lines[0].get_ydata()
            np.testing.assert_array_equal(np.ma.getmaskarray(y), [False, True, True, False])
            np.testing.assert_allclose(y.compressed(), [2.2, 5.6])
            self.assertIn("2 puntos", warnings.getvalue())

    def test_band_crossing_zero_is_omitted_and_all_invalid_data_is_rejected(self):
        from io import StringIO
        from unittest.mock import patch
        from matplotlib.figure import Figure

        groups = [(4096, 4.3, [Path("one"), Path("two")], np.array([0.1, 0.2, 0.3]),
                   np.array([-0.1, -1.0, -2.0]), np.array([0.2, 0.1, 0.2]))]
        with patch.object(Figure, "savefig", autospec=True) as save, patch("sys.stderr", new_callable=StringIO) as warnings:
            plot_averages(groups, self.directory / "plot.png")
            axes = save.call_args.args[0].axes[0]
            for polygon in axes.collections[0].get_paths():
                self.assertTrue(np.all(polygon.vertices[:, 1] > 0))
                self.assertTrue(np.all(polygon.vertices[:, 0] >= 0.2))
            self.assertIn("banda", warnings.getvalue())
        invalid = [(4096, 4.3, [Path("one")], np.array([0.1]), np.array([1.0]), np.zeros(1))]
        with patch("sys.stderr", new_callable=StringIO), self.assertRaisesRegex(ValueError, "No hay valores"):
            plot_averages(invalid, self.directory / "invalid.png")


if __name__ == "__main__":
    unittest.main()

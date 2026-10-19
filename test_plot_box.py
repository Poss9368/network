"""Pruebas de los cortes periódicos y lectura de snapshots."""

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from plot_box import (
    find_link_intersections,
    find_steps,
    network_geometry,
    plot_step,
    read_step,
    split_periodic_contact,
)


class PeriodicContactTests(unittest.TestCase):
    def test_inside_box(self):
        segments = split_periodic_contact((0.2, 0.3), (0.4, 0.3), 1.0)
        np.testing.assert_allclose(segments, [[(0.2, 0.3), (0.4, 0.3)]])

    def test_crossing_x(self):
        segments = split_periodic_contact((0.95, 0.4), (0.05, 0.6), 1.0)
        np.testing.assert_allclose(segments, [[(0.95, 0.4), (1, 0.5)], [(0, 0.5), (0.05, 0.6)]])

    def test_crossing_y_in_box_of_length_two(self):
        segments = split_periodic_contact((0.8, 1.9), (1.2, 0.1), 2.0)
        np.testing.assert_allclose(segments, [[(0.8, 1.9), (1, 2)], [(1, 0), (1.2, 0.1)]])

    def test_two_separate_crossings(self):
        segments = split_periodic_contact((0.9, 0.8), (0.1, 0.1), 1.0)
        self.assertEqual(len(segments), 3)
        self.assertAlmostEqual(sum(math.dist(*segment) for segment in segments), math.hypot(0.2, 0.3))
        self.assertTrue(all(0 <= x <= 1 for segment in segments for point in segment for x in point))

    def test_corner_is_a_single_crossing(self):
        segments = split_periodic_contact((0.9, 0.9), (0.1, 0.1), 1.0)
        self.assertEqual(len(segments), 2)
        np.testing.assert_allclose(segments, [[(0.9, 0.9), (1, 1)], [(0, 0), (0.1, 0.1)]])

    def test_particle_exactly_on_border(self):
        segments = split_periodic_contact((0, 0.4), (0.9, 0.4), 1.0)
        np.testing.assert_allclose(segments, [[(1, 0.4), (0.9, 0.4)]])

    def test_reads_box_and_unique_contacts_and_sorts_steps_numerically(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            data = "2.0 -1.1\n1.9 1.0 1 -1\n0.1 1.0 0 -1\n"
            for step in [10, 2, 1]:
                (directory / f"sample_step_{step}.dat").write_text(data)
            (directory / "sample_resumen.dat").write_text("# ignored\n")
            name, steps = find_steps("sample.dat", directory)
            self.assertEqual(name, "sample")
            self.assertEqual([step for step, _ in steps], [1, 2, 10])
            length, positions, contacts = read_step(steps[0][1])
            self.assertEqual(length, 2.0)
            np.testing.assert_allclose(positions, [(1.9, 1), (0.1, 1)])
            self.assertEqual(contacts, [(0, 1)])


class LinkIntersectionTests(unittest.TestCase):
    def test_two_independent_links_cross(self):
        positions = [(0.3, 0.3), (0.7, 0.7), (0.3, 0.7), (0.7, 0.3)]
        _, crossings = network_geometry(positions, [(0, 1), (2, 3)], 1.0)
        np.testing.assert_allclose(crossings, [(0.5, 0.5)])

    def test_shared_particle_is_not_a_crossing(self):
        positions = [(0.5, 0.5), (0.3, 0.3), (0.7, 0.3)]
        _, crossings = network_geometry(positions, [(0, 1), (0, 2)], 1.0)
        self.assertEqual(crossings, [])

    def test_three_links_at_the_same_point_have_one_marker(self):
        positions = [(0.3, 0.5), (0.7, 0.5), (0.5, 0.3), (0.5, 0.7), (0.3, 0.3), (0.7, 0.7)]
        _, crossings = network_geometry(positions, [(0, 1), (2, 3), (4, 5)], 1.0)
        np.testing.assert_allclose(crossings, [(0.5, 0.5)])

    def test_periodic_link_alone_does_not_create_markers(self):
        segments, crossings = network_geometry([(0.95, 0.4), (0.05, 0.6)], [(0, 1)], 1.0)
        self.assertEqual(len(segments), 2)
        self.assertEqual(crossings, [])

    def test_intersection_with_part_of_a_periodic_link(self):
        positions = [(0.95, 0.4), (0.05, 0.6), (0.98, 0.42), (0.98, 0.52)]
        _, crossings = network_geometry(positions, [(0, 1), (2, 3)], 1.0)
        np.testing.assert_allclose(crossings, [(0.98, 0.46)])

    def test_true_crossing_at_boundary_has_both_visible_images(self):
        positions = [(0.95, 0.45), (0.05, 0.55), (0.95, 0.55), (0.05, 0.45)]
        _, crossings = network_geometry(positions, [(0, 1), (2, 3)], 1.0)
        np.testing.assert_allclose(sorted(crossings), [(0, 0.5), (1, 0.5)], atol=1e-12)

    def test_parallel_and_collinear_overlaps_do_not_create_spurious_points(self):
        segments = [((0.2, 0.2), (0.8, 0.2)), ((0.2, 0.3), (0.8, 0.3)), ((0.4, 0.2), (0.6, 0.2))]
        crossings = find_link_intersections(segments, [0, 1, 2], [(0, 1), (2, 3), (4, 5)], 1.0)
        self.assertEqual(crossings, [])

    def test_grid_does_not_miss_crossings_on_cell_boundaries(self):
        segments = [((0.88, 0.88), (0.92, 0.92)), ((0.88, 0.92), (0.92, 0.88))]
        segments.extend([((0.1, 0.1), (0.2, 0.1))] * 98)
        contacts = [(2 * i, 2 * i + 1) for i in range(len(segments))]
        crossings = find_link_intersections(segments, list(range(len(segments))), contacts, 1.0)
        np.testing.assert_allclose(crossings, [(0.9, 0.9)])

    def test_render_is_unit_sized_centered_and_has_no_axes_or_title(self):
        from unittest.mock import patch
        from matplotlib.figure import Figure

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample.dat"
            path.write_text("2.0 -1.1\n0.2 0.2 1\n0.4 0.2 0\n")
            with patch.object(Figure, "savefig", autospec=True) as save:
                plot_step(path, Path(temporary) / "sample.png", 100)
                self.assertEqual(save.call_count, 1)
                figure = save.call_args.args[0]
                axes = figure.axes[0]
                self.assertFalse(axes.axison)
                self.assertEqual(axes.get_title(), "")
                self.assertIsNone(figure._suptitle)
                self.assertIsNone(axes.get_legend())
                np.testing.assert_allclose(axes.get_xlim(), [-0.54, 0.54])
                np.testing.assert_allclose(axes.get_ylim(), [-0.54, 0.54])
                box = axes.patches[0]
                self.assertEqual(box.get_xy(), (-0.5, -0.5))
                self.assertEqual(box.get_width(), 1)
                self.assertEqual(box.get_height(), 1)
                np.testing.assert_allclose(axes.collections[0].get_segments(), [[(-0.4, -0.4), (-0.3, -0.4)]])

    def test_normalization_removes_affine_motion_but_keeps_relative_motion(self):
        from unittest.mock import patch
        from matplotlib.figure import Figure

        snapshots = [
            "1.0 -1.1\n0.2 0.2 1\n0.4 0.2 0\n",
            "2.0 -1.1\n0.4 0.4 1\n0.8 0.4 0\n",
            "2.0 -1.1\n0.4 0.4 1\n1.0 0.4 0\n",
        ]
        plotted = []
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample.dat"
            for snapshot in snapshots:
                path.write_text(snapshot)
                with patch.object(Figure, "savefig", autospec=True) as save:
                    plot_step(path, path.with_suffix(".png"), 100)
                    axes = save.call_args.args[0].axes[0]
                    plotted.append(axes.collections[1].get_offsets().copy())
        np.testing.assert_allclose(plotted[0], plotted[1])
        np.testing.assert_allclose(plotted[2], [(-0.3, -0.3), (0.0, -0.3)])
        self.assertNotAlmostEqual(plotted[2][1, 0], plotted[1][1, 0])

    def test_crossing_markers_are_normalized_and_centered(self):
        from unittest.mock import patch
        from matplotlib.figure import Figure

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample.dat"
            path.write_text("2.0 -1.1\n0.6 0.4 1\n1.4 1.2 0\n0.6 1.2 3\n1.4 0.4 2\n")
            with patch.object(Figure, "savefig", autospec=True) as save:
                plot_step(path, path.with_suffix(".png"), 100)
                axes = save.call_args.args[0].axes[0]
                np.testing.assert_allclose(axes.collections[2].get_offsets(), [(0.0, -0.1)], atol=1e-12)


if __name__ == "__main__":
    unittest.main()

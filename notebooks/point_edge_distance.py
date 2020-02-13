"""Test the point-edge distance computation."""

import unittest
import numpy


def point_point_distance(p, q):
    """Compute the distance from point p to point q."""
    return numpy.linalg.norm(p - q)


def point_edge_distance(p, v, w):
    """Compute the distance from point p to edge (v, w)."""
    v2p = p - v
    v2w = w - v
    edge_length_sq = v2w.T @ v2w
    if(edge_length_sq < 1e-10):
        return 0, point_point_distance(p, v)
    else:
        t = min(max((v2w.T @ v2p) / edge_length_sq, 0), 1)
        return t, point_point_distance(p, v + t * v2w)


class TestPointEdgeDistance(unittest.TestCase):
    """Test the point to edge distance computation."""

    def test_upper(self):
        """Test the point to edge distance computation."""
        for i in range(1000):
            p = numpy.random.random(2) * 100 - 50
            v = numpy.random.random(2) * 100 - 50
            w = numpy.random.random(2) * 100 - 50
            t, dist = point_edge_distance(p, v, w)

            more_t = min(t + numpy.random.random(), 1)
            dist_more_t = point_point_distance(p, v + more_t * (w - v))
            less_t = max(t - numpy.random.random(), 0)
            dist_less_t = point_point_distance(p, v + less_t * (w - v))

            self.assertTrue(dist <= dist_less_t)
            self.assertTrue(dist <= dist_more_t)


if __name__ == '__main__':
    unittest.main()

import argparse
from collections import Counter
from itertools import product
import json


OEIS_SHELLS = [
    1, 10, 33, 73, 128, 199, 285, 388, 506, 640, 789, 955, 1136,
    1333, 1545, 1774, 2018, 2278, 2553, 2845, 3152, 3475, 3813,
    4168, 4538, 4924, 5325, 5743, 6176, 6625, 7089, 7570, 8066,
    8578, 9105, 9649, 10208, 10783, 11373, 11980, 12602, 13240, 13893,
]
OEIS_BALLS = [
    1, 11, 44, 117, 245, 444, 729, 1117, 1623, 2263, 3052, 4007,
    5143, 6476, 8021, 9795, 11813, 14091, 16644, 19489, 22641,
    26116, 29929, 34097, 38635, 43559, 48884, 54627, 60803, 67428,
    74517, 82087, 90153, 98731, 107836, 117485, 127693, 138476, 149849,
]
OEIS_TCD_SHELLS = [
    1, 10, 33, 72, 126, 196, 281, 382, 498, 630, 777, 940, 1118,
    1312, 1521, 1746, 1986, 2242, 2513, 2800, 3102, 3420, 3753,
    4102, 4466, 4846, 5241, 5652, 6078, 6520, 6977, 7450, 7938,
    8442, 8961, 9496, 10046, 10612, 11193, 11790, 12402, 13030, 13673, 14332,
]
OEIS_TCD_BALLS = [
    1, 11, 44, 116, 242, 438, 719, 1101, 1599, 2229, 3006, 3946,
    5064, 6376, 7897, 9643, 11629, 13871, 16384, 19184, 22286,
    25706, 29459, 33561, 38027, 42873, 48114, 53766, 59844, 66364,
    73341, 80791, 88729, 97171, 106132, 115628, 125674, 136286, 147479, 159269,
]


def slab_counts(level, tiling='tea'):
    if tiling == 'tcd':
        return ((level + 1) // 2, 0) if level >= 0 else (0, (-level) // 2)
    if level >= 0:
        return (level + 3) // 4, (level + 1) // 4
    return (-level + 2) // 4, (-level) // 4


def residual(first, second, plus, minus):
    return max(0, first - plus, second - plus, first + second - plus,
               -minus - first, -minus - second, -minus - first - second)


def metric(vertex, tiling='tea'):
    first, second, level = vertex
    plus, minus = slab_counts(level, tiling)
    return abs(level) + residual(first, second, plus, minus)


def neighbors(vertex, tiling='tea'):
    first, second, level = vertex
    for delta_first, delta_second in [(1, 0), (-1, 0), (0, 1), (0, -1),
                                      (1, -1), (-1, 1)]:
        yield first + delta_first, second + delta_second, level
    yield first, second, level + 1
    yield first, second, level - 1
    if level % 2 == 0:
        sign = 1 if tiling == 'tcd' or level % 4 == 0 else -1
        yield first + sign, second, level + 1
        yield first, second + sign, level + 1
    if (level - 1) % 2 == 0:
        sign = -1 if tiling == 'tcd' or (level - 1) % 4 == 0 else 1
        yield first + sign, second, level - 1
        yield first, second + sign, level - 1


def shell_formula(radius, tiling='tea'):
    if radius == 0:
        return 1
    if tiling == 'tcd':
        numerator = 31 * radius**2 + [8, 9][radius % 2]
        assert numerator % 4 == 0
        return numerator // 4
    numerator = 63 * radius**2 + [16, 17, 12, 17][radius % 4]
    assert numerator % 8 == 0
    return numerator // 8


def ball_formula(radius, tiling='tea'):
    if tiling == 'tcd':
        numerator = 62 * radius**3 + 93 * radius**2 + 82 * radius + [24, 27][radius % 2]
        assert numerator % 24 == 0
        return numerator // 24
    numerator = (84 * radius**3 + 126 * radius**2 + 104 * radius
                 + [32, 38, 24, 30][radius % 4])
    assert numerator % 32 == 0
    return numerator // 32


def hexagon_count(plus, minus):
    numerator = plus**2 + 4 * plus * minus + minus**2 + 3 * plus + 3 * minus + 2
    assert numerator % 2 == 0
    return numerator // 2


def layer_shell(radius, tiling='tea'):
    total = 0
    for level in range(-radius, radius + 1):
        plus, minus = slab_counts(level, tiling)
        remaining = radius - abs(level)
        total += hexagon_count(plus, minus) if remaining == 0 else 3 * (plus + minus) + 6 * remaining
    return total


def convolve(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return result


def graph_checks(radius, tiling='tea'):
    distances = {(0, 0, 0): 0}
    frontier = [(0, 0, 0)]
    shells = [1]
    for distance in range(1, radius + 1):
        next_frontier = []
        for vertex in frontier:
            for adjacent in neighbors(vertex, tiling):
                if adjacent not in distances:
                    distances[adjacent] = distance
                    assert metric(adjacent, tiling) == distance, (adjacent, distance)
                    next_frontier.append(adjacent)
        frontier = next_frontier
        shells.append(len(frontier))
        assert shells[-1] == shell_formula(distance, tiling)
        assert len(distances) == ball_formula(distance, tiling)
    expected = OEIS_SHELLS if tiling == 'tea' else OEIS_TCD_SHELLS
    assert shells == expected[:radius + 1]
    assert Counter(distances.values()) == Counter(dict(enumerate(shells)))
    return len(distances)


def orientation_checks(radius=8):
    boundaries = list(range(-radius, radius, 2))
    graphs = 0
    checked_vertices = 0
    for choices in product([-1, 1], repeat=len(boundaries)):
        signs = dict(zip(boundaries, choices))

        def parameters(level):
            crossed = range(0, level) if level >= 0 else range(level, 0)
            effective = [signs[boundary] * (1 if level >= 0 else -1)
                         for boundary in crossed if boundary % 2 == 0]
            return effective.count(1), effective.count(-1)

        def adjacent(vertex):
            first, second, level = vertex
            for delta_first, delta_second in [(1, 0), (-1, 0), (0, 1), (0, -1),
                                              (1, -1), (-1, 1)]:
                yield first + delta_first, second + delta_second, level
            for direction in [-1, 1]:
                boundary = level if direction == 1 else level - 1
                yield first, second, level + direction
                if boundary % 2 == 0:
                    sign = direction * signs[boundary]
                    yield first + sign, second, level + direction
                    yield first, second + sign, level + direction

        distances = {(0, 0, 0): 0}
        frontier = [(0, 0, 0)]
        for distance in range(1, radius + 1):
            next_frontier = []
            for vertex in frontier:
                for neighbor in adjacent(vertex):
                    if neighbor not in distances:
                        distances[neighbor] = distance
                        plus, minus = parameters(neighbor[2])
                        assert abs(neighbor[2]) + residual(
                            neighbor[0], neighbor[1], plus, minus) == distance
                        next_frontier.append(neighbor)
            frontier = next_frontier
            plus, minus = parameters(distance)
            down_plus, down_minus = parameters(-distance)
            correction = plus * minus + down_plus * down_minus
            expected = shell_formula(distance, 'tcd') + correction
            assert len(frontier) == expected
            assert shell_formula(distance, 'tcd') <= expected <= shell_formula(distance)
            assert (expected == shell_formula(distance, 'tcd')) == (correction == 0)
            assert (expected == shell_formula(distance)) == (
                abs(plus - minus) <= 1 and abs(down_plus - down_minus) <= 1)
            assert ball_formula(distance, 'tcd') <= len(distances) <= ball_formula(distance)
        graphs += 1
        checked_vertices += len(distances)
    return {'boundary_levels': boundaries, 'sign_assignments': graphs,
            'bfs_radius': radius, 'total_vertices': checked_vertices,
            'coverage': 'Distances, exact orientation correction, sharp sphere and ball bounds, sphere equality conditions.'}


def main():
    parser = argparse.ArgumentParser(description='Exact finite checks for the tea graph formulas.')
    parser.add_argument('--radius', type=int, default=32)
    arguments = parser.parse_args()
    assert 1 <= arguments.radius <= 42
    vertices = graph_checks(arguments.radius)
    tcd_vertices = graph_checks(arguments.radius, 'tcd')
    orientations = orientation_checks()
    gap_sum = 0
    for radius in range(1001):
        gap = radius**2 // 8
        assert shell_formula(radius) - shell_formula(radius, 'tcd') == gap
        gap_sum += gap
        assert ball_formula(radius) - ball_formula(radius, 'tcd') == gap_sum
    for tiling in ['tea', 'tcd']:
        for radius in range(1001):
            assert layer_shell(radius, tiling) == shell_formula(radius, tiling)
    for radius, expected in enumerate(OEIS_SHELLS):
        assert shell_formula(radius) == expected
    for radius, expected in enumerate(OEIS_BALLS):
        assert ball_formula(radius) == expected
    for radius, expected in enumerate(OEIS_TCD_SHELLS):
        assert shell_formula(radius, 'tcd') == expected
    for radius, expected in enumerate(OEIS_TCD_BALLS):
        assert ball_formula(radius, 'tcd') == expected
    for plus in range(9):
        for minus in range(9):
            points = sum(-minus <= first + second <= plus
                         for first in range(-minus, plus + 1)
                         for second in range(-minus, plus + 1))
            assert points == hexagon_count(plus, minus)
            for distance in range(1, 7):
                assert hexagon_count(plus + distance, minus + distance) - hexagon_count(
                    plus + distance - 1, minus + distance - 1) == 3 * (plus + minus) + 6 * distance
    for level in range(-12, 13):
        for first in range(-8, 9):
            for second in range(-8, 9):
                vertex = (first, second, level)
                adjacent_vertices = list(neighbors(vertex))
                assert len(set(adjacent_vertices)) == 10
                for adjacent in adjacent_vertices:
                    assert vertex in set(neighbors(adjacent))
                    assert abs(metric(vertex) - metric(adjacent)) <= 1
                reflected = (-first, -second, 1 - level)
                translated = (-first, -second, level + 2)
                assert {(-adj[0], -adj[1], 1 - adj[2]) for adj in adjacent_vertices} == set(neighbors(reflected))
                assert {(-adj[0], -adj[1], adj[2] + 2) for adj in adjacent_vertices} == set(neighbors(translated))
    shell_denominator = [1, -2, 1, 0, -1, 2, -1]
    ball_denominator = convolve(shell_denominator, [1, -1])
    expected_numerator = [1, 8, 14, 17, 14, 8, 1]
    shells = [shell_formula(radius) for radius in range(1001)]
    balls = [ball_formula(radius) for radius in range(1001)]
    for series, denominator in [(shells, shell_denominator), (balls, ball_denominator)]:
        product = convolve(series, denominator)
        assert product[:7] == expected_numerator
        assert not any(product[7:1001])
    assert sum(expected_numerator) == 63
    assert sum(value * (-1)**index for index, value in enumerate(expected_numerator)) == -3
    remainder = [0, 0]
    for index, value in enumerate(expected_numerator):
        remainder[index % 2] += value * (-1)**(index // 2)
    assert remainder == [0, -1]
    tcd_numerator = [1, 8, 13, 8, 1]
    tcd_denominator = [1, -2, 0, 2, -1]
    for denominator, formula in [(tcd_denominator, shell_formula),
                                 (convolve(tcd_denominator, [1, -1]), ball_formula)]:
        series = [formula(radius, 'tcd') for radius in range(1001)]
        product = convolve(series, denominator)
        assert product[:5] == tcd_numerator
        assert not any(product[5:1001])
    assert sum(tcd_numerator) == 31
    assert sum(value * (-1)**index for index, value in enumerate(tcd_numerator)) == -1
    print(json.dumps({
        'status': 'pass',
        'coverage': 'Finite integer checks of tea, tcd and arbitrary slab-orientation graphs, distances, sphere and ball formulas.',
        'arbitrary_orientations': orientations,
        'tea_tcd_gap_radii': [0, 1000],
        'bfs_radius': arguments.radius,
        'bfs_vertices': vertices,
        'tcd_bfs_vertices': tcd_vertices,
        'distance_formula': f'Every vertex of graph distance at most {arguments.radius}.',
        'oeis_a299285_terms': len(OEIS_SHELLS),
        'oeis_a299286_terms': len(OEIS_BALLS),
        'oeis_a299287_terms': len(OEIS_TCD_SHELLS),
        'oeis_a299288_terms': len(OEIS_TCD_BALLS),
        'layer_sum_radii': [0, 1000],
        'hexagon_parameters': {'plus': [0, 8], 'minus': [0, 8], 'shell_depth': [1, 6]},
        'edge_parameter_box': {'first': [-8, 8], 'second': [-8, 8], 'level': [-12, 12]},
        'generating_function_coefficients': [0, 1000],
        'numerator_at_one': 63,
        'numerator_at_minus_one': -3,
        'numerator_mod_x_squared_plus_one': [0, -1],
    }, indent=2))


if __name__ == '__main__':
    main()

<?php

declare(strict_types=1);

namespace CheapRuler;

final class CheapRuler
{
    // values that define WGS84 ellipsoid model of the Earth
    private const RE = 6378.137; // equatorial radius
    private const FE = 1 / 298.257223563; // flattening
    private const E2 = self::FE * (2 - self::FE);
    private const RAD = M_PI / 180;

    public const KILOMETERS = 1;
    public const MILES = 2;
    public const NAUTICALMILES = 3;
    public const METERS = 4;
    public const METRES = 5;
    public const YARDS = 6;
    public const FEET = 7;
    public const INCHES = 8;

    private const FACTORS = [
        self::KILOMETERS =>  1,
        self::MILES => 1000 / 1609.344,
        self::NAUTICALMILES => 1000 / 1852,
        self::METERS => 1000,
        self::METRES => 1000,
        self::YARDS => 1000 / 0.9144,
        self::FEET => 1000 / 0.3048,
        self::INCHES => 1000 / 0.0254,
    ];

    private float $kx;

    private float $ky;

    /**
     * Creates a ruler object from tile coordinates (y and z).
     */
    public static function fromTile(float $y, float $z, int $units): self
    {
        $n = M_PI * (1 - 2 * ($y + 0.5) / (2 ** $z));
        $lat = atan(0.5 * (exp($n) - exp(-$n))) / self::RAD;
        return new self($lat, $units);
    }

    /**
     * Creates a ruler instance for very fast approximations to common geodesic measurements around a certain latitude.
     */
    public function __construct(float $lat, int $units)
    {
        $factor = self::FACTORS[$units] ?? null;
        if ($factor === null) {
            throw new InvalidArgumentException('Invalid units passed.');
        }

        // curvature formulas from https://en.wikipedia.org/wiki/Earth_radius#Meridional
        $m = self::RAD * self::RE * $factor;
        $coslat = cos($lat * self::RAD);
        $w2 = 1 / (1 - self::E2 * (1 - $coslat * $coslat));
        $w = sqrt($w2);

        // multipliers for converting longitude and latitude degrees into distance
        $this->kx = $m * $w * $coslat; // based on normal radius of curvature
        $this->ky = $m * $w * $w2 * (1 - self::E2); // based on meridonal radius of curvature
    }

    /**
     * Given two points of the form [longitude, latitude], returns the distance.
     * @param float[] $a a point [longitude, latitude]
     * @param float[] $b b point [longitude, latitude]
     * @return float distance
     */
    public function distance(array $a, array $b): float
    {
        $dx = self::wrap($a[0] - $b[0]) * $this->kx;
        $dy = ($a[1] - $b[1]) * $this->ky;
        return sqrt($dx * $dx + $dy * $dy);
    }

    /**
     * @param float[] $a a point [longitude, latitude]
     * @param float[] $b b point [longitude, latitude]
     * @return float
     */
    public function bearing(array $a, array $b): float
    {
        $dx = self::wrap($b[0] - $a[0]) * $this->kx;
        $dy = ($b[1] - $a[1]) * $this->ky;
        return atan2($dx, $dy) / self::RAD;
    }

    /**
     * Returns a new point given distance and bearing from the starting point.
     * @param float[] $p p point [longitude, latitude]
     * @param float $dist distance
     * @param float $bearing bearing
     * @return float[] point [longitude, latitude]
     */
    public function destination(array $p, float $dist, float $bearing): array
    {
        $a = $bearing * self::RAD;
        return $this->offset($p, sin($a) * $dist, cos($a) * $dist);
    }

    /**
     * Returns a new point given easting and northing offsets (in ruler units) from the starting point.
     * @param array $p p point [longitude, latitude]
     * @param float $dx easting
     * @param float $dy northing
     * @return float[] point [longitude, latitude]
     */
    public function offset(array $p, float $dx, float $dy): array
    {
        return [$p[0] + $dx / $this->kx, $p[1] + $dy / $this->ky];
    }

    /**
     * @param array<float[]> $points array of points [longitude, latitude]
     * @return float total line distance
     */
    public function lineDistance(array $points): float
    {
        $total = 0.0;
        $length = count($points);
        for ($i = 0; $i < $length - 1; $i++) {
            $total += $this->distance($points[$i], $points[$i + 1]);
        }
        return $total;
    }

    /**
     * Given a polygon (an array of rings, where each ring is an array of points), returns the area.
     * @param array<array<float[]>> $polygon polygon
     * @return float value in the specified units (square kilometers by default)
     */
    public function area(array $polygon): float
    {
        $sum = 0.0;
        foreach ($polygon as $i => $ring) {
            $len = count($ring);
            $k = $len -1;
            for ($j = 0; $j < $len; $k = $j++) {
                $sum += self::wrap($ring[$j][0] - $ring[$k][0]) * ($ring[$j][1] + $ring[$k][1]) * ($i ? -1 : 1);
            }
        }
        return (abs($sum) / 2) * $this->kx * $this->ky;
    }

    /**
     * Returns the point at a specified distance along the line
     * @param array<float[]> $line line
     * @param float $dist $distance
     * @return array
     */
    public function along(array $line, float $dist): array
    {
        $sum = 0;

        if ($dist <= 0) {
            return $line[0];
        }

        $len = count($line) -1;
        for ($i = 0; $i < $len - 1; $i++) {
            $p0 = $line[$i];
            $p1 = $line[$i + 1];
            $d = $this->distance($p0, $p1);
            $sum += $d;
            if ($sum > $dist) {
                return self::interpolate($p0, $p1, ($dist - ($sum - $d)) / $d);
            }
        }

        return $line[$len];
    }

    /**
     * Returns the distance from a point `p` to a line segment `a` to `b`
     * @param float[] $p p point [longitude, latitude]
     * @param float[] $a p1 segment point 1 [longitude, latitude]
     * @param float[] $b p2 segment point 2 [longitude, latitude]
     * @return float
     */
    public function pointToSegmentDistance(array $p, array $a, array $b): float
    {
        [$x, $y] = $a;
        $dx = self::wrap($b[0] - $x) * $this->kx;
        $dy = ($b[1] - $y) * $this->ky;

        if ($dx !== 0 || $dy !== 0) {
            $t = (self::wrap($p[0] - $x) * $this->kx * $dx + ($p[1] - $y) * $this->ky * $dy) / ($dx * $dx + $dy * $dy);

            if ($t > 1) {
                $x = $b[0];
                $y = $b[1];
            } elseif ($t > 0) {
                $x += ($dx / $this->kx) * $t;
                $y += ($dy / $this->ky) * $t;
            }
        }

        $dx = self::wrap($p[0] - $x) * $this->kx;
        $dy = ($p[1] - $y) * $this->ky;

        return sqrt($dx * $dx + $dy * $dy);
    }

    /**
     * Returns a part of the given line between the start and the stop points (or their closest points on the line).
     * @param float[] $start start point [longitude, latitude]
     * @param float[] $stop stop point [longitude, latitude]
     * @param array<float[]> $line line
     * @return array<float[]> line part of a line
     */
    public function lineSlice(array $start, array $stop, array $line): array
    {
        $p1 = $this->pointOnLine($line, $start);
        $p2 = $this->pointOnLine($line, $stop);

        if ($p1['index'] > $p2['index'] || ($p1['index'] === $p2['index'] && $p1['t'] > $p2['t'])) {
            $tmp = $p1;
            $p1 = $p2;
            $p2 = $tmp;
        }

        $slice = [$p1['point']];

        $l = $p1['index'] + 1;
        $r = $p2['index'];

        if (!self::equals($line[$l], $slice[0]) && $l <= $r) {
            $slice[] = $line[$l];
        }

        for ($i = $l + 1; $i <= $r; $i++) {
            $slice[] = $line[$i];
        }

        if (!self::equals($line[$r], $p2['point'])) {
            $slice[] = $p2['point'];
        }

        return $slice;
    }

    /**
     * @param float $start start distance
     * @param float $stop stop distance
     * @param array<float[]> $line line
     * @return array<float[]> line part of a line
     */
    public function lineSliceAlong(float $start, float $stop, array $line): array
    {
        $sum = 0;
        $slice = [];
        $length = count($line);

        for ($i = 0; $i < $length - 1; $i++) {
            $p0 = $line[$i];
            $p1 = $line[$i + 1];
            $d = $this->distance($p0, $p1);

            $sum += $d;

            if ($sum > $start && count($slice) === 0) {
                $slice[] = self::interpolate($p0, $p1, ($start - ($sum - $d)) / $d);
            }

            if ($sum >= $stop) {
                $slice[] = self::interpolate($p0, $p1, ($stop - ($sum - $d)) / $d);
                return $slice;
            }

            if ($sum > $start) {
                $slice[] = $p1;
            }
        }
        return $slice;
    }

    /**
     * Returns an object of the form {point, index, t}, where point is closest point on the line
     * from the given point, index is the start index of the segment with the closest point,
     * and t is a parameter from 0 to 1 that indicates where the closest point is on that segment.
     * @param array<float[]> $line line
     * @param float[] $p p point [longitude, latitude]
     * @return array Array of point, index, t
     */
    public function pointOnLine(array $line, array $p): array
    {
        $minDist = PHP_FLOAT_MAX; // Infinity;
        // $minX, $minY, $minI, $minT;
        $length = count($line);

        for ($i = 0; $i < $length - 1; $i++) {
            $x = $line[$i][0];
            $y = $line[$i][1];
            $dx = self::wrap($line[$i + 1][0] - $x) * $this->kx;
            $dy = ($line[$i + 1][1] - $y) * $this->ky;
            $t = 0;

            if ($dx !== 0 || $dy !== 0) {
                $t = (self::wrap($p[0] - $x) * $this->kx * $dx + ($p[1] - $y) * $this->ky * $dy) / ($dx * $dx + $dy * $dy);
                if ($t > 1) {
                    $x = $line[$i + 1][0];
                    $y = $line[$i + 1][1];
                } elseif ($t > 0) {
                    $x += ($dx / $this->kx) * $t;
                    $y += ($dy / $this->ky) * $t;
                }
            }

            $dx = self::wrap($p[0] - $x) * $this->kx;
            $dy = ($p[1] - $y) * $this->ky;

            $sqDist = $dx * $dx + $dy * $dy;
            if ($sqDist < $minDist) {
                $minDist = $sqDist;
                $minX = $x;
                $minY = $y;
                $minI = $i;
                $minT = $t;
            }
        }

        return ['point' => [$minX, $minY], 'index' => $minI, 't' => max(0, min(1, $minT))];
    }

    /**
     * @param float[] $p p point [longitude, latitude]
     * @param float $buffer buffer
     * @return float[] box object ([w, s, e, n])
     */
    public function bufferPoint(array $p, float $buffer): array
    {
        $v = $buffer / $this->ky;
        $h = $buffer / $this->kx;
        return [$p[0] - $h, $p[1] - $v, $p[0] + $h, $p[1] + $v];
    }

    /**
     * @param array<float[]> $bbox box object ([w, s, e, n])
     * @param float $buffer buffer
     * @return float[] box object ([w, s, e, n])
     */
    public function bufferBBox(array $bbox, float $buffer): array
    {
        $v = $buffer / $this->ky;
        $h = $buffer / $this->kx;

        return [
            $bbox[0] - $h,
            $bbox[1] - $v,
            $bbox[2] + $h,
            $bbox[3] + $v,
        ];
    }

    /**
     * @param float[] $p p point [longitude, latitude]
     * @param float[] $bbox box object ([w, s, e, n])
     * @return bool
     */
    public function insideBBox(array $p, array $bbox): bool
    {
        return self::wrap($p[0] - $bbox[0]) >= 0 && self::wrap($p[0] - $bbox[2]) <= 0 && $p[1] >= $bbox[1] && $p[1] <= $bbox[3];
    }

    /**
     * @param float[] $a a point [longitude, latitude]
     * @param float[] $b b point [longitude, latitude]
     * @return bool
     */
    private static function equals(array $a, array $b): bool
    {
        return $a[0] === $b[0] && $a[1] === $b[1];
    }

    private static function interpolate(array $a, array $b, float $t): array
    {
        $dx = self::wrap($b[0] - $a[0]);
        $dy = $b[1] - $a[1];
        return [$a[0] + $dx * $t, $a[1] + $dy * $t];
    }

    /**
     * Normalize a degree value into [-180..180] range
     * @param float $deg
     * @return float
     */
    private static function wrap(float $deg): float
    {
        while ($deg < -180) {
            $deg += 360;
        }
        while ($deg > 180) {
            $deg -= 360;
        }
        return $deg;
    }
}

using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using Windows.Devices.Geolocation;
using Windows.Foundation;
using Windows.UI.Xaml.Media;

namespace IlanCurves.CurveFunctions
{
    public static class CurveFunctions
    {
        /// <summary>
        /// For a given set of X control points, calculate the value t
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static IEnumerable<double> FindTForXCubic(this List<Point> points, double x)
        {
            return FindTForControlPointsCubic(points[0].X, points[1].X, points[2].X, points[3].X, x);
        }

        /// <summary>
        /// For a given set of X control points, calculate the value t
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static IEnumerable<double> FindTForYCubic(this List<Point> points, double y)
        {
            return FindTForControlPointsCubic(points[0].Y, points[1].Y, points[2].Y, points[3].Y, y);
        }

        /// <summary>
        /// Given cubic bezier control points on a plane, find all values of t for a given intersection.
        /// </summary>
        /// <param name="p0">Control point 0</param>
        /// <param name="p1">Control point 1</param>
        /// <param name="p2">Control point 2</param>
        /// <param name="p3">Control point 3</param>
        /// <param name="intersection">A list of all the t values</param>
        /// <returns></returns>
        public static IEnumerable<double> FindTForControlPointsCubic(double p0, double p1, double p2, double p3, double intersection)
        {
            var roots =
                FindRoots.Polynomial(
                    new Polynomial(
                        new double[]
                        {
                            p0 - intersection,
                            3.0 * (p1 - p0),
                            3.0 * (p0 - 2.0 * p1 + p2),
                            -p0 + 3.0 * p1 - 3.0 * p2 + p3
                        })
                    );

            return roots.Where(r => r.IsReal() && r.Real >= 0 && r.Real <= 1).Select(r => r.Real);
        }

        /// <summary>
        /// Finds all value Y values for a given X intersect coordinate within a curve
        /// </summary>
        /// <param name="points">The control points along the curve, must be n*3+1 points in length</param>
        /// <param name="x">The value defining the x intersect</param>
        /// <returns>The resulting y values or an empty IENumerable</returns>
        public static IEnumerable<double> FindYForX(this PointCollection points, double x)
        {
            List<double> results = new List<double>();

            for (int i = 0; i < points.Count - 3; i += 3)
            {
                var current = points.Skip(i).Take(4).ToList();
                var ts = current.FindTForXCubic(x);

                // TODO : find the right linq function for this
                foreach (var t in ts)
                {
                    var y = BezierSolveForT(current[0].Y, current[1].Y, current[2].Y, current[3].Y, t);
                    results.Add(y);
                }
            }

            return results;
        }


        /// <summary>
        /// Finds all value X values for a given Y intersect coordinate within a curve
        /// </summary>
        /// <param name="points">The control points along the curve, must be n*3+1 points in length</param>
        /// <param name="x">The value defining the x intersect</param>
        /// <returns>The resulting x values or an empty IENumerable</returns>
        public static IEnumerable<double> FindXForY(this PointCollection points, double y)
        {
            List<double> results = new List<double>();

            for (int i = 0; i < points.Count - 3; i += 3)
            {
                var current = points.Skip(i).Take(4).ToList();
                var ts = current.FindTForYCubic(y);

                // TODO : find the right linq function for this
                foreach (var t in ts)
                {
                    var x = BezierSolveForT(current[0].X, current[1].X, current[2].X, current[3].X, t);
                    results.Add(x);
                }
            }

            return results;
        }


        /// <summary>
        /// Calculates a point for a given t along a cubic bezier curve given the provided control points
        /// </summary>
        /// <param name="p0">Start control point</param>
        /// <param name="p1">Starting-inner control point</param>
        /// <param name="p2">Ending-inner control point</param>
        /// <param name="p3">End control point</param>
        /// <param name="t">t, must be between 0 and 1</param>
        /// <returns>The corresponding point along path.</returns>
        /// <seealso cref="https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Cubic_B%C3%A9zier_curves"/>
        public static double BezierSolveForT(double p0, double p1, double p2, double p3, double t)
        {
            return
                p0 * Math.Pow(1.0 - t, 3) +
                p1 * 3 * Math.Pow(1.0 - t, 2) * t +
                p2 * 3 * (1.0 - t) * Math.Pow(t, 2) +
                p3 * Math.Pow(t, 3)
                ;
        }

        /// <summary>
        /// Given a series of points, perform a bezier spline interpolation.
        /// </summary>
        /// <param name="points">The points to start at, end at and to pass through a an ordered list</param>
        /// <returns>A series of control points for drawing bezier curves</returns>
        public static PointCollection SplineInterpolate(this IEnumerable<Point> points)
        {
            var cpx = ComputeControlPoints(points.Select(p => p.X).ToArray());
            var cpy = ComputeControlPoints(points.Select(p => p.Y).ToArray());

            var result = new PointCollection { points.First() };

            int i = 0;
            foreach (var p in points.Skip(1))
            {
                result.Add(new Point(cpx[0][i], cpy[0][i]));
                result.Add(new Point(cpx[1][i], cpy[1][i]));
                i++;

                result.Add(p);
            }

            return result;
        }

        /// <summary>
        /// Performs control point estimation for cubic bezier spline interpolation through a series of points
        /// </summary>
        /// <param name="K">An array of input values to interpolate</param>
        /// <returns>A two dimension array containing first and second control points</returns>
        /// 
        /// This function is a direct port of the code presented on https://www.particleincell.com/2012/bezier-splines/
        /// TODO: Consider using matrices as provided by C#
        private static double[][] ComputeControlPoints(double[] K)
        {
            var n = K.Length - 1;

            // rhs vector
            var a = new double[K.Length];
            var b = new double[K.Length];
            var c = new double[K.Length];
            var rightHandSide = new double[K.Length];

            // left most segment
            a[0] = 0;
            b[0] = 2;
            c[0] = 1;
            rightHandSide[0] = K[0] + 2 * K[1];

            // internal segments
            for (var i = 1; i < n - 1; i++)
            {
                a[i] = 1;
                b[i] = 4;
                c[i] = 1;
                rightHandSide[i] = 4 * K[i] + 2 * K[i + 1];
            }

            // right segment
            a[n - 1] = 2;
            b[n - 1] = 7;
            c[n - 1] = 0;
            rightHandSide[n - 1] = 8 * K[n - 1] + K[n];

            // solves Ax=b with the Thomas algorithm (from Wikipedia)
            for (var i = 1; i < n; i++)
            {
                var m = a[i] / b[i - 1];
                b[i] -= m * c[i - 1];
                rightHandSide[i] = rightHandSide[i] - m * rightHandSide[i - 1];
            }

            var p1 = new double[K.Length];
            var p2 = new double[K.Length];

            p1[n - 1] = rightHandSide[n - 1] / b[n - 1];
            for (var i = n - 2; i >= 0; --i)
                p1[i] = (rightHandSide[i] - c[i] * p1[i + 1]) / b[i];

            /*we have p1, now compute p2*/
            for (var i = 0; i < n - 1; i++)
                p2[i] = 2 * K[i + 1] - p1[i + 1];

            p2[n - 1] = 0.5 * (K[n] + p1[n - 1]);

            return new double[][] { p1, p2 };
        }

        /// <summary>
        /// Calculates the bounding rectangle of a multisegment Bezier curve
        /// </summary>
        /// <param name="points">The control points of the curve</param>
        /// <returns>The bounding rectangle or null for invalid data</returns>
        public static Rect BezierExtents(PointCollection points)
        {
            if (points.Count == 0)
                return Rect.Empty;

            if ((points.Count - 1) % 3 != 0)
                //throw new ArgumentOutOfRangeException("points", "Cubic bezier curves must contain 1 + n*3 points");
                return Rect.Empty;

            var result = BezierExtentsSingle(points[0], points[1], points[2], points[3]);
            for (var i = 4; i < points.Count; i += 3)
                result.Union(BezierExtentsSingle(points[i - 1], points[i], points[i + 1], points[i + 2]));

            return result;
        }

        /// <summary>
        /// Calculates the bounding rectangle of the Bezier curve
        /// </summary>
        /// <param name="points">A list of 4 control points</param>
        /// <returns>The rectangle describing the boundaries of the given bezier curve</returns>
        public static Rect BezierExtentsSingle(Point p0, Point p1, Point p2, Point p3)
        {
            var xBounds = BezierExtentsSingle(p0.X, p1.X, p2.X, p3.X);
            var yBounds = BezierExtentsSingle(p0.Y, p1.Y, p2.Y, p3.Y);

            return new Rect
            {
                X = xBounds.Item1,
                Y = yBounds.Item1,
                Width = xBounds.Item2 - xBounds.Item1,
                Height = yBounds.Item2 - yBounds.Item1
            };
        }

        /// <summary>
        /// Finds the extents of the cubic bezier function for a given set of control points on a single axis
        /// </summary>
        /// <param name="p0">First control point</param>
        /// <param name="p1">Second control point</param>
        /// <param name="p2">Third control point</param>
        /// <param name="p3">Fourth control point</param>
        /// <returns>A tuple containing the ordered lower and upper bounds of the curve</returns>
        private static Tuple<double, double> BezierExtentsSingle(double p0, double p1, double p2, double p3)
        {
            // Start by preparing the results as a simple bounding box of the two ends of the line, this should handle linear conditions
            var lowerBound = p0;
            var upperBound = p0;

            if (p3 < lowerBound)
                lowerBound = p3;

            if (p3 > upperBound)
                upperBound = p3;

            // Coefficients for first derivative of the cubic bezier function
            var a = 3 * p3 - 9 * p2 + 9 * p1 - 3 * p0;
            var b = 6 * p0 - 12 * p1 + 6 * p2;
            var c = 3 * p1 - 3 * p0;

            // Quadratic formula to solve for first derivative
            var radicand = b * b - 4 * a * c;

            // Only solve if the radicand is zero or positive.
            if (radicand >= 0)
            {
                var root = Math.Sqrt(radicand);
                var t1 = (-b + root) / (2 * a);
                var t2 = (-b - root) / (2 * a);

                if (t1 > 0 && t1 < 1)
                {
                    var result1 = BezierSolveForT(p0, p1, p2, p3, t1);
                    if (result1 < lowerBound)
                        lowerBound = result1;

                    if (result1 > upperBound)
                        upperBound = result1;
                }

                if (t2 > 0 && t2 < 1)
                {
                    var result2 = BezierSolveForT(p0, p1, p2, p3, t2);

                    if (result2 < lowerBound)
                        lowerBound = result2;

                    if (result2 > upperBound)
                        upperBound = result2;
                }
            }

            // Return the upper and lower bounds
            return new Tuple<double, double>(lowerBound, upperBound);
        }
    }
}
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
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
            var roots =
                FindRoots.Polynomial(
                    new Polynomial(
                        new double[]
                        {
                            points[0].X - x,
                            3.0 * (points[1].X - points[0].X),
                            3.0 * (points[0].X - 2.0 * points[1].X + points[2].X),
                            -points[0].X + 3.0 * points[1].X - 3.0 * points[2].X + points[3].X
                        }
                        )
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
        /// TODO: Eliminate lists and use arrays instead.
        /// TODO: Consider using matrices as provided by C#
        private static double[][] ComputeControlPoints(double[] K)
        {
            var n = K.Length - 1;

            // rhs vector
            var a = new double[K.Length];
            var b = new double[K.Length];
            var c = new double[K.Length];
            var r = new double[K.Length];

            // left most segment
            a[0] = 0;
            b[0] = 2;
            c[0] = 1;
            r[0] = K[0] + 2 * K[1];

            // internal segments
            for (var i = 1; i < n - 1; i++)
            {
                a[i] = 1;
                b[i] = 4;
                c[i] = 1;
                r[i] = 4 * K[i] + 2 * K[i + 1];
            }

            // right segment
            a[n - 1] = 2;
            b[n - 1] = 7;
            c[n - 1] = 0;
            r[n - 1] = 8 * K[n - 1] + K[n];

            // solves Ax=b with the Thomas algorithm (from Wikipedia)
            for (var i = 1; i < n; i++)
            {
                var m = a[i] / b[i - 1];
                b[i] = b[i] - m * c[i - 1];
                r[i] = r[i] - m * r[i - 1];
            }

            var p1 = new double[K.Length];
            var p2 = new double[K.Length];

            p1[n - 1] = r[n - 1] / b[n - 1];
            for (var i = n - 2; i >= 0; --i)
                p1[i] = (r[i] - c[i] * p1[i + 1]) / b[i];

            /*we have p1, now compute p2*/
            for (var i = 0; i < n - 1; i++)
                p2[i] = 2 * K[i + 1] - p1[i + 1];

            p2[n - 1] = 0.5 * (K[n] + p1[n - 1]);

            return new double[][] { p1, p2 };
        }
    }
}
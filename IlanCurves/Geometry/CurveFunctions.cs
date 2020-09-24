using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using Windows.Foundation;
using Windows.UI.Xaml.Media;

namespace IlanCurves.Geometry
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
                MathNet.Numerics.FindRoots.Polynomial(
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

        public static Point NearestPointOnCurve(Point reflectFrom, Point[] points)
        {
            var P = new Vector2 { X = reflectFrom.X, Y = reflectFrom.Y };
            var V = points.Select(p => new Vector2 { X = p.X, Y = p.Y }).ToArray();

            var curveResults = new List<Vector2>();
            for (var i = 3; i < points.Length; i += 3)
                curveResults.Add(NearestPointOnCurve(P, V.Skip(i - 3).Take(4).ToArray()));

            var distances = curveResults.Select(x => (x - P).SquaredLength).ToList();
            var shortestDistance = double.MaxValue;
            Vector2 result = P;
            for(var i=0; i<distances.Count(); i++)
            {
                if(distances[i] < shortestDistance)
                {
                    result = curveResults[i];
                    shortestDistance = distances[i];
                }
            }

            return new Point { X = result.X, Y = result.Y };
        }

        // Cubic Bezier curve
        private static int Degree = 3;
        // Degree of eqn to find roots of
        private static int WDegree = 5;
        // Maximum depth for recursion
        private static int MaxDepth = 64;
        // Flatness control value 
        private static double Epsilon { get; } = Math.Pow(2, -MaxDepth);

        /*
         *  NearestPointOnCurve :
         *  	Compute the parameter value of the point on a Bezier
         *		curve segment closest to some arbtitrary, user-input point.
         *		Return the point on the curve at that parameter value.
         *
         */
        public static Vector2 NearestPointOnCurve(Vector2 P, Vector2[] V)
        {
            //  Convert problem to 5th-degree Bezier form
            var w = ConvertToBezierForm(P, V);

            // Find all possible roots of 5th-degree equation 
            var t_candidate = FindRoots(w, WDegree, 0);

            /* Compare distances of P to all candidates, and to t=0, and t=1 */

            double new_dist;

            // Check distance to beginning of curve, where t = 0	
            var v = P - V[0];
            var dist = v.SquaredLength;

            //Parameter value of closest pt
            var t = 0.0;

            /* Find distances for candidate points	*/
            for (var i = 0; i < t_candidate.Count(); i++)
            {
                var p = Curve(V, Degree, t_candidate[i]).Point;

                //p = Curve(V, Degree, t_candidate[i], null, null);
                new_dist = (P - p).SquaredLength;
                if (new_dist < dist)
                {
                    dist = new_dist;
                    t = t_candidate[i];
                }
            }

            // Finally, look at distance to end point, where t = 1.0
            new_dist = (P - V[Degree]).SquaredLength;
            if (new_dist < dist)
            {
                //dist = new_dist;
                t = 1.0;
            }

            // Return the point on the curve at parameter value t 
            //System.Diagnostics.Debug.WriteLine("t : " + t.ToString("N4.12"));
            return Curve(V, Degree, t).Point;
        }


        // Precalculate the choose component of the Berenstein Bezier polynomial
        // for all values of the equation for a cubic Bezier curve
        private static readonly double[,] BerensteinChooseCubic = new double[,]
        {	
            // Precomputed "z" for cubics
	        {1.0, 0.6, 0.3, 0.1},
            {0.4, 0.6, 0.6, 0.4},
            { 0.1, 0.3, 0.6, 1.0},
        };


        /// ConvertToBezierForm :
        ///		Given a point and a Bezier curve, generate a 5th-degree
        ///		Bezier-format equation whose solution finds the point on the
        ///     curve nearest the user-defined point.
        ///         - Point2 P;      The point to find t for	
        ///         - Point2* V;	 The control points		
        static Vector2[] ConvertToBezierForm(Vector2 P, Vector2[] V)
        {
            // Determine the c's -- these are vectors created by subtracting
            // point P from each of the control points				
            var c = new Vector2[Degree + 1];      // V(i)'s - P
            for (var i = 0; i <= Degree; i++)
                c[i] = V[i] - P;

            // Determine the d's -- these are vectors created by subtracting
            // each control point from the next					
            var d = new Vector2[Degree];    // V(i+1) - V(i)
            for (var i = 0; i <= Degree - 1; i++)
                d[i] = (V[i + 1] - V[i]).Scale(3.0);

            // Create the c,d table -- this is a table of dot products of the 
            // c's and d's							
            var cdTable = new double[3, 4]; // Dot product of c, d
            for (var row = 0; row <= Degree - 1; row++)
            {
                for (var column = 0; column <= Degree; column++)
                    cdTable[row, column] = d[row].DotProduct(c[column]);
            }

            // Now, apply the z's to the dot products, on the skew diagonal
            // Also, set up the x-values, making these "points"		
            var w = new Vector2[WDegree + 1];
            for (var i = 0; i <= WDegree; i++)
            {
                w[i] = new Vector2
                {
                    X = Convert.ToDouble(i) / Convert.ToDouble(WDegree),
                    Y = 0
                };
            }

            var n = Degree;
            var m = Degree - 1;
            for (var k = 0; k <= n + m; k++)
            {
                var lb = Math.Max(0, k - m);
                var ub = Math.Min(k, n);
                for (var i = lb; i <= ub; i++)
                {
                    var j = k - i;
                    w[i + j].Y += cdTable[j, i] * BerensteinChooseCubic[j, i];
                }
            }

            return (w);
        }

        // Ugly return value for Curve function
        internal class CurveResult
        {
            public Vector2[] LeftCurve = null;
            public Vector2[] RightRight = null;
            public Vector2 Point = null;
        }

        /// Bezier : 
        ///	Evaluate a Bezier curve at a particular parameter value
        ///     Fill in control points for resulting sub-curves if "Left" and
        ///	"Right" are non-null.
        ///
        ///     - int degree;      Degree of bezier curve	
        ///     - Point2* V;       Control pts			
        ///     - double t;        Parameter value		
        private static CurveResult Curve(Vector2[] V, int degree, double t)
        {
            const int WDegree = 5;
            Vector2[,] Vtemp = new Vector2[5 + 1, WDegree + 1];

            // Copy control points	
            for (var j = 0; j <= degree; j++)
                Vtemp[0, j] = V[j];

            // Triangle computation	
            for (var i = 1; i <= degree; i++)
            {
                for (var j = 0; j <= degree - i; j++)
                {
                    Vtemp[i, j] = new Vector2
                    {
                        X = (1.0 - t) * Vtemp[i - 1, j].X + t * Vtemp[i - 1, j + 1].X,
                        Y = (1.0 - t) * Vtemp[i - 1, j].Y + t * Vtemp[i - 1, j + 1].Y
                    };
                }
            }

            var result = new CurveResult
            {
                Point = Vtemp[degree, 0],
                LeftCurve = new Vector2[WDegree + 1],
                RightRight = new Vector2[WDegree + 1]
            };

            for (var i = 0; i <= degree; i++)
            {
                result.LeftCurve[i] = Vtemp[i, 0];
                result.RightRight[i] = Vtemp[degree - i, i];
            }

            return result;
        }

        static int CrossingCount(Vector2[] V, int degree)
        {
            int result = 0;    // Number of zero-crossings	

            //  Sign of coefficients	
            var sign = Math.Sign(V[0].Y);
            var old_sign = sign;

            for (var i = 1; i <= degree; i++)
            {
                sign = Math.Sign(V[i].Y);
                if (sign != old_sign)
                    result++;

                old_sign = sign;
            }

            return result;
        }

        /// ControlPolygonFlatEnough :
        ///   Check if the control polygon of a Bezier curve is flat enough
        ///    for recursive subdivision to bottom out.
        /// 
        /// Corrections by James Walker, jw@jwwalker.com, as follows:
        /// 
        /// There seem to be errors in the ControlPolygonFlatEnough function in the
        /// Graphics Gems book and the repository (NearestPoint.c). This function
        /// is briefly described on p. 413 of the text, and appears on pages 793-794.
        /// I see two main problems with it.
        /// 
        /// The idea is to find an upper bound for the error of approximating the x
        /// intercept of the Bezier curve by the x intercept of the line through the
        /// first and last control points. It is claimed on p. 413 that this error is
        /// bounded by half of the difference between the intercepts of the bounding
        /// box. I don't see why that should be true. The line joining the first and
        /// last control points can be on one side of the bounding box, and the actual
        /// curve can be near the opposite side, so the bound should be the difference
        /// of the bounding box intercepts, not half of it.
        /// 
        /// Second, we come to the implementation. The values distance[i] computed in
        /// the first loop are not actual distances, but squares of distances. I
        /// realize that minimizing or maximizing the squares is equivalent to
        /// minimizing or maximizing the distances.  But when the code claims that
        /// one of the sides of the bounding box has equation
        /// a * x + b * y + c + max_distance_above, where max_distance_above is one of
        /// those squared distances, that makes no sense to me.
        /// 
        /// I have appended my version of the function. If you apply my code to the
        /// cubic Bezier curve used to test NearestPoint.c,
        /// 
        /// static Point2 bezCurve[4] = {    /  A cubic Bezier curve    /
        ///   { 0.0, 0.0 },
        ///   { 1.0, 2.0 },
        ///   { 3.0, 3.0 },
        ///   { 4.0, 2.0 },
        /// };
        /// 
        /// my code computes left_intercept = -3.0 and right_intercept = 0.0, which you
        /// can verify by sketching a graph. The original code computes
        /// left_intercept = 0.0 and right_intercept = 0.9
        static bool ControlPolygonFlatEnough(Vector2[] V, int degree)
        {
            // Derive the implicit equation for line connecting first 
            //  and last control points 
            var a = V[0].Y - V[degree].Y;
            var b = V[degree].X - V[0].X;
            var c = V[0].X * V[degree].Y - V[degree].X * V[0].Y;

            var max_distance_above = 0.0;
            var max_distance_below = 0.0;

            for (var i = 1; i < degree; i++)
            {
                var value = a * V[i].X + b * V[i].Y + c;

                if (value > max_distance_above)
                    max_distance_above = value;
                else if (value < max_distance_below)
                    max_distance_below = value;
            }

            //  Implicit equation for zero line 
            const double a1 = 0.0;
            const double b1 = 1.0;
            const double c1 = 0.0;

            //  Implicit equation for "above" line 
            var a2 = a;
            var b2 = b;
            var c2 = c - max_distance_above;

            var det = a1 * b2 - a2 * b1;
            var dInv = 1.0 / det;

            var intercept_1 = (b1 * c2 - b2 * c1) * dInv;

            //  Implicit equation for "below" line 
            a2 = a;
            b2 = b;
            c2 = c - max_distance_below;

            det = a1 * b2 - a2 * b1;
            dInv = 1.0 / det;

            var intercept_2 = (b1 * c2 - b2 * c1) * dInv;

            // Compute intercepts of bounding box   
            var left_intercept = Math.Min(intercept_1, intercept_2);
            var right_intercept = Math.Max(intercept_1, intercept_2);

            //Precision of root
            var error = right_intercept - left_intercept;

            return error < Epsilon;
        }

        /// ComputeXIntercept :
        /// Compute intersection of chord from first control point to last
        ///  	with 0-axis.
        /// 
        /// NOTE: "T" and "Y" do not have to be computed, and there are many useless
        ///   * operations in the following(e.g. "0.0 - 0.0").
        /// Point2* V;          - Control points	
        /// int degree; 		- Degree of curve	
        static double ComputeXIntercept(Vector2[] V, int degree)
        {
            var XLK = 1.0 - 0.0;
            var YLK = 0.0 - 0.0;
            var XNM = V[degree].X - V[0].X;
            var YNM = V[degree].Y - V[0].Y;
            var XMK = V[0].X - 0.0;
            var YMK = V[0].Y - 0.0;

            var det = XNM * YLK - YNM * XLK;
            var detInv = 1.0 / det;

            var S = (XNM * YMK - YNM * XMK) * detInv;
            //  T = (XLK*YMK - YLK*XMK) * detInv; 

            var X = 0.0 + XLK * S;
            //  Y = 0.0 + YLK * S; 

            return X;
        }

        /// <summary>
        /// Given a 5th-degree equation in Bernstein-Bezier form, find
        /// all of the roots in the interval [0, 1]. 
        /// </summary>
        /// <param name="w">The control points representing the curve</param>
        /// <param name="degree">The degree of the polynomial</param>
        /// <param name="depth">The depth of the recurssion</param>
        /// <returns>The roots found</returns>
        private static double[] FindRoots(Vector2[] w, int degree, int depth)
        {

            int crossings = CrossingCount(w, degree);

            if (crossings == 0) // No solutions here	
                return new double[0];

            if (crossings == 1) // Unique solution
            {
                // Stop recursion when the tree is deep enough	
                // if deep enough, return 1 solution at midpoint 	
                if (depth >= MaxDepth)
                {
                    return new double[] {
                        (w[0].X + w[WDegree].X) / 2.0
                    };
                }

                if (ControlPolygonFlatEnough(w, degree))
                {
                    return new double[] {
                        ComputeXIntercept(w, degree)
                    };
                }
            }

            // Otherwise, solve recursively after	
            // subdividing control polygon		
            var curve = Curve(w, degree, 0.5);

            var leftT = FindRoots(curve.LeftCurve, degree, depth + 1);
            var rightT = FindRoots(curve.RightRight, degree, depth + 1);

            // Solutions from kids
            var result = new double[leftT.Count() + rightT.Count()];
            leftT.CopyTo(result, 0);
            rightT.CopyTo(result, leftT.Count());

            return result;
        }

        public static Point DarrensAlgorithm(double t, PointCollection points)
        {
            return new Point {
                X = DarrensAlgorithm(t, points[0].X, points[1].X, points[2].X, points[3].X),
                Y = DarrensAlgorithm(t, points[0].Y, points[1].Y, points[2].Y, points[3].Y)
            };
        }

        /// <summary>
        /// Solves the bezier curve using the same formula but written in a way that may be faster.
        /// </summary>
        /// <param name="t"></param>
        /// <param name="p0"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <returns></returns>
        private static double DarrensAlgorithm(double t, double p0, double p1, double p2, double p3)
        {
            //$$A = P_0 + t(P_1 - P_0)$$
            //$$B = P_1 + t(P_2 - P_1)$$
            //$$C = P_2 + t(P_3 - P_1)$$
            //$$D = A + t(B - A)$$
            //$$E = B + t(C - B)$$
            //$$F = D + t(E - D)$$

            var f =  t * ((3 * (p1 - p0)) + t * ((3 * (p2 - 2 * p1 + p0)) + t * (p3 - 3 * p2 + 3 * p1 - p0))) + p0;

            return f;
        }
    }
}
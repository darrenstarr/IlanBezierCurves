using System;
using Windows.Foundation;

namespace IlanCurves.Geometry
{
    public static class Clipping
    {
        public static Tuple<double, double, double, double> ClipLineToRect(Rect r, double x1, double y1, double x2, double y2)
        {
            if (x2 == x1)
            {
                if (y1 < r.Top) y1 = r.Top;
                if (y2 < r.Top) y2 = r.Top;
                if (y1 > r.Bottom) y1 = r.Bottom;
                if (y2 > r.Bottom) y2 = r.Bottom;
            }
            else
            {
                if (x1 > x2)
                {
                    var t = x1;
                    x1 = x2;
                    x2 = t;

                    t = y1;
                    y1 = y2;
                    y2 = t;
                }

                var m = (y2 - y1) / (x2 - x1);
                var b = -(m * x1 - y1);

                if (x1 < r.Left)
                {
                    x1 = r.Left;
                    y1 = m * r.Left + b;
                }

                if (x2 > r.Right)
                {
                    x2 = r.Right;
                    y2 = m * r.Right + b;
                }

                if (y1 > y2)
                {
                    var t = x1;
                    x1 = x2;
                    x2 = t;

                    t = y1;
                    y1 = y2;
                    y2 = t;

                    m = (y2 - y1) / (x2 - x1);
                    b = -(m * x1 - y1);
                }

                if (y1 < r.Top)
                {
                    y1 = r.Top;
                    x1 = (y1 - b) / m;
                }

                if (y2 > r.Bottom)
                {
                    y2 = r.Bottom;
                    x2 = (y2 - b) / m;
                }
            }

            return new Tuple<double, double, double, double>(x1, y1, x2, y2);
        }
    }
}

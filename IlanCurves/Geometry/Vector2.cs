using System;

namespace IlanCurves.Geometry
{
    public class Vector2
    {
        public double X { get; set; } = 0.0;
        public double Y { get; set; } = 0.0;

        public double SquaredLength
        {
            get { return (X * X) + (Y * Y); }
        }

        public double DistanceTo(Vector2 other)
        {
            return Math.Sqrt(Math.Pow(X - other.X, 2) + Math.Pow(Y - other.Y, 2));
        }

        public static Vector2 operator +(Vector2 left, Vector2 right)
        {
            return new Vector2 { X = left.X + right.X, Y = left.Y + right.Y };
        }

        public static Vector2 operator -(Vector2 left, Vector2 right)
        {
            return new Vector2 { X = left.X - right.X, Y = left.Y - right.Y };
        }

        public static Vector2 operator *(Vector2 left, Vector2 right)
        {
            return new Vector2 { X = left.X * right.X, Y = left.Y * right.Y };
        }

        public Vector2 Scale(double scale)
        {
            return new Vector2 { X = X * scale, Y = Y * scale };
        }

        public double DotProduct(Vector2 against)
        {
            var product = this * against;
            //System.Diagnostics.Debug.WriteLine($"{X}, {Y} times {against.X}, {against.Y} = {product}");
            return product.X + product.Y;
        }

        public override string ToString()
        {
            return "(" + X.ToString() + ", " + Y.ToString() + ")";
        }
    }
}
using IlanCurves.CurveFunctions;
using System;
using System.Collections.Generic;
using System.Linq;
using Windows.ApplicationModel;
using Windows.Foundation;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Shapes;

// The User Control item template is documented at https://go.microsoft.com/fwlink/?LinkId=234236

namespace IlanCurves
{
    public sealed partial class BezierView : UserControl
    {
        private Rect rawControlPointsBounds;
        private bool showControlPoints = true;
        private bool showGrid = true;
        private bool showCurve = true;
        private bool invertY = true;
        private bool mouseTracking = true;
        private bool showMouseXIntercepts = true;
        private bool showMouseYIntercepts = true;
        private bool showLinear = false;

        private Point graphLowerBounds = new Point(0, 0);
        private Point graphExtents = new Point(1, 1);

        private Point scale = new Point(1, 1);
        private readonly Point graphPaddingPercent = new Point(5, 5);
        private Point graphOffset = new Point();
        private Rect graphRect = new Rect();
        private readonly PointCollection scaledControlPoints = new PointCollection();

        private Point TrackingMousePosition;
        private Line MouseTrackingLineX;
        private Line MouseTrackingLineY;


        private List<Ellipse> MouseXInterceptBubbles = new List<Ellipse>();
        private List<Ellipse> MouseYInterceptBubbles = new List<Ellipse>();

        private readonly Windows.UI.Color[] BubbleColors = new Windows.UI.Color[]
        {
                    Windows.UI.Colors.Yellow,
                    Windows.UI.Colors.Orange,
                    Windows.UI.Colors.HotPink
        };

        public static readonly DependencyProperty ControlPointsProperty =
            DependencyProperty.Register("ControlPoints", typeof(PointCollection), typeof(BezierView), new PropertyMetadata(0));

        public PointCollection ControlPoints
        {
            get => (PointCollection)GetValue(ControlPointsProperty);
            set
            {
                SetValue(ControlPointsProperty, value);
                ScaleControlPoints();
                PaintCanvas();
            }
        }

        public bool ShowControlPoints
        {
            get => showControlPoints;
            set
            {
                showControlPoints = value;
                PaintCanvas();
            }
        }

        public bool ShowGrid
        {
            get => showGrid;
            set
            {
                showGrid = value;
                PaintCanvas();
            }
        }

        public bool ShowCurve
        {
            get => showCurve;
            set
            {
                showCurve = value;
                PaintCanvas();
            }
        }

        public bool ShowLinear
        {
            get => showLinear;
            set
            {
                showLinear = value;
                PaintCanvas();
            }
        }

        public bool InvertY
        {
            get => invertY;
            set
            {
                invertY = value;
                PaintCanvas();
            }
        }

        public bool MouseTracking
        {
            get => mouseTracking;
            set
            {
                mouseTracking = value;
                MouseTrackingLineX = null;
                MouseTrackingLineY = null;
                TrackingMousePosition = new Point();
                PaintCanvas();
            }
        }

        public bool ShowMouseXIntercepts
        {
            get => showMouseXIntercepts;
            set
            {
                showMouseXIntercepts = value;
                PaintCanvas();
            }
        }

        public bool ShowMouseYIntercepts
        {
            get => showMouseYIntercepts;
            set
            {
                showMouseYIntercepts = value;
                PaintCanvas();
            }
        }

        /// <summary>
        /// Size of the big bubbles on the UI
        /// </summary>
        public double BubbleSize { get; set; } = 9.0;

        /// <summary>
        /// Size of the little bubbles on the UI
        /// </summary>
        public double LittleBubbleSize { get; set; } = 5;


        public BezierView()
        {
            this.InitializeComponent();

            ControlPoints = new PointCollection();

            if (DesignMode.DesignModeEnabled)
            {
                var points = new PointCollection
                {
                    new Point(10, 10),
                    new Point(15, 5),
                    new Point(12, 3)
                };
                ControlPoints = points.SplineInterpolate();
            }
        }

        private double ScaleX(double x)
        {
            return (x - graphLowerBounds.X) * scale.X + graphOffset.X;
        }

        private double UnscaleX(double x)
        {
            return (x - graphOffset.X) / scale.X - graphLowerBounds.X;
        }

        private double ScaleY(double y)
        {
            if (invertY)
                return canvas.ActualHeight - (y - graphLowerBounds.Y) * scale.Y - graphOffset.Y;

            return (y - graphLowerBounds.Y) * scale.Y + graphOffset.Y;
        }

        private double UnscaleY(double y)
        {
            if (invertY)
                return (canvas.ActualHeight - y - graphOffset.Y) / scale.Y + graphLowerBounds.Y;
            
            return (y - graphOffset.Y) / scale.Y - graphLowerBounds.Y;
        }

        private Point ScaledPoint(Point p)
        {
            return new Point
            {
                X = ScaleX(p.X),
                Y = ScaleY(p.Y)
            };
        }

        private void ScaleControlPoints()
        {
            if (canvas.ActualWidth == 0 || canvas.ActualHeight == 0)
                return;

            rawControlPointsBounds = CurveFunctions.CurveFunctions.BezierExtents(ControlPoints);

            graphLowerBounds.X = rawControlPointsBounds.X < 0 ? Math.Floor(rawControlPointsBounds.X) : 0;
            graphLowerBounds.Y = rawControlPointsBounds.Y < 0 ? Math.Floor(rawControlPointsBounds.Y) : 0;
            graphExtents.X = Math.Ceiling(rawControlPointsBounds.Right - graphLowerBounds.X);
            graphExtents.Y = Math.Ceiling(rawControlPointsBounds.Bottom - graphLowerBounds.Y);

            scaledControlPoints.Clear();

            graphOffset.X = canvas.ActualWidth / 100 * graphPaddingPercent.X;
            graphOffset.Y = canvas.ActualHeight / 100 * graphPaddingPercent.Y;

            var graphWidth = canvas.ActualWidth - (graphOffset.X * 2);
            var graphHeight = canvas.ActualHeight - (graphOffset.Y * 2);

            graphRect.X = graphOffset.X;
            graphRect.Y = graphOffset.Y;
            graphRect.Width = graphWidth;
            graphRect.Height = graphHeight;

            scale.X = graphWidth / graphExtents.X;
            scale.Y = graphHeight / graphExtents.Y;

            foreach (var p in ControlPoints)
                scaledControlPoints.Add(ScaledPoint(p));

            PaintCanvas();
        }

        private void PaintCanvas()
        {
            canvas.Children.Clear();
            MouseXInterceptBubbles.Clear();
            MouseYInterceptBubbles.Clear();

            // Paint a border to the grid
            canvas.Children.Add(
                new Rectangle
                {
                    Width = graphRect.Width,
                    Height = graphRect.Height,
                    Margin = new Thickness(graphRect.X, graphRect.Y, 0, 0),
                    StrokeThickness = 1,
                    Stroke = new SolidColorBrush(Windows.UI.Colors.LightBlue),
                });

            if (ShowGrid)
                PaintGrid();

            if (MouseTracking)
                PaintMouseTracking();

            if (ShowLinear)
                PaintLinear();

            if (ShowCurve)
                PaintCurve();

            if (ShowControlPoints)
                PaintControlPoints();

            PaintEndPoints();
        }

        private void PaintGrid()
        {
            for (var i = Math.Ceiling(graphLowerBounds.X); i < rawControlPointsBounds.Right; i += 1.0)
            {
                var v = ScaleX(i);
                canvas.Children.Add(
                    new Line
                    {
                        X1 = v,
                        Y1 = graphRect.Top,
                        X2 = v,
                        Y2 = graphRect.Bottom,
                        StrokeThickness = (i == 0) ? 2 : 1,
                        Stroke = new SolidColorBrush(i == 0 ? Windows.UI.Colors.Blue : Windows.UI.Colors.LightBlue),
                    });
            }

            for (var i = Math.Ceiling(graphLowerBounds.Y); i < rawControlPointsBounds.Bottom; i += 1.0)
            {
                var v = ScaleY(i);
                canvas.Children.Add(
                    new Line
                    {
                        X1 = graphRect.Left,
                        Y1 = v,
                        X2 = graphRect.Right,
                        Y2 = v,
                        StrokeThickness = (i == 0) ? 2 : 1,
                        Stroke = new SolidColorBrush(i == 0 ? Windows.UI.Colors.Blue : Windows.UI.Colors.LightBlue),
                    });
            }
        }

        private void PaintCurve()
        {
            var pathSegmentCollection1 = new PathSegmentCollection();
            for (var i = 1; i < scaledControlPoints.Count; i += 3)
            {
                pathSegmentCollection1.Add(new BezierSegment
                {
                    Point1 = scaledControlPoints[i],
                    Point2 = scaledControlPoints[i + 1],
                    Point3 = scaledControlPoints[i + 2]
                });
            }

            canvas.Children.Add(new Path
            {
                Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                StrokeThickness = 1,
                Data = new GeometryGroup
                {
                    Children = {
                            new PathGeometry
                            {
                                Figures = new PathFigureCollection
                                {
                                    new PathFigure
                                    {
                                        StartPoint = scaledControlPoints[0],
                                        Segments = pathSegmentCollection1
                                    }
                                }
                            }
                        }
                }
            });
        }

        private void PaintMouseTracking()
        {
            MouseTrackingLineX = new Line
            {
                X1 = TrackingMousePosition.X,
                Y1 = graphRect.Top,
                X2 = TrackingMousePosition.X,
                Y2 = graphRect.Bottom,
                StrokeThickness = 1,
                Stroke = new SolidColorBrush(Windows.UI.Colors.LightCoral)
            };

            canvas.Children.Add(MouseTrackingLineX);

            MouseTrackingLineY = new Line
            {
                X1 = graphRect.Left,
                Y1 = TrackingMousePosition.Y,
                X2 = graphRect.Right,
                Y2 = TrackingMousePosition.Y,
                StrokeThickness = 1,
                Stroke = new SolidColorBrush(Windows.UI.Colors.LightCoral)
            };

            canvas.Children.Add(MouseTrackingLineY);
        }

        private void PaintControlPoints()
        {
            for (var t = 1; t < scaledControlPoints.Count; t++)
            {
                if ((t % 3) != 0)
                {
                    var p = scaledControlPoints[t];
                    var anchor = (t % 3) == 1 ? scaledControlPoints[t - 1] : scaledControlPoints[t + 1];

                    canvas.Children.Add(new Line
                    {
                        X1 = p.X,
                        Y1 = p.Y,
                        X2 = anchor.X,
                        Y2 = anchor.Y,
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                    });

                    canvas.Children.Add(
                        new Ellipse
                        {
                            Width = LittleBubbleSize,
                            Height = LittleBubbleSize,
                            Margin = new Thickness(p.X - LittleBubbleSize / 2, p.Y - LittleBubbleSize / 2, 0, 0),
                            StrokeThickness = 1,
                            Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                            Fill = new SolidColorBrush(BubbleColors[t % 3]),
                        });
                }
            }
        }

        private void PaintLinear()
        {
            var points = new PointCollection();
            for (var t = 0; t < scaledControlPoints.Count; t+=3)
                points.Add(scaledControlPoints[t]);

            canvas.Children.Add(
                new Polyline
                {
                    StrokeThickness = 1,
                    Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                    Points = points
                });
        }

        private void PaintEndPoints()
        {
            for (var t = 0; t < scaledControlPoints.Count; t++)
            {
                if ((t % 3) == 0)
                {
                    var p = scaledControlPoints[t];

                    canvas.Children.Add(
                        new Ellipse
                        {
                            Width = BubbleSize,
                            Height = BubbleSize,
                            Margin = new Thickness(p.X - BubbleSize / 2, p.Y - BubbleSize / 2, 0, 0),
                            StrokeThickness = 1,
                            Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                            Fill = new SolidColorBrush(Windows.UI.Colors.Yellow),
                        });
                }
            }
        }

        private void Canvas_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            ScaleControlPoints();
        }


        private void UpdateMouseTrackingLines()
        {
            if (MouseTrackingLineX != null)
            {
                MouseTrackingLineX.X1 = TrackingMousePosition.X;
                MouseTrackingLineX.X2 = TrackingMousePosition.X;
            }

            if (MouseTrackingLineY != null)
            {
                MouseTrackingLineY.Y1 = TrackingMousePosition.Y;
                MouseTrackingLineY.Y2 = TrackingMousePosition.Y;
            }
        }

        private void UpdateXIntercepts()
        {
            var currentX = UnscaleX(TrackingMousePosition.X);

            var xIntercepts = ControlPoints.FindYForX(currentX);

            while (xIntercepts.Count() > MouseXInterceptBubbles.Count)
            {
                var trackingBubble = new Ellipse
                {
                    Width = BubbleSize,
                    Height = BubbleSize,
                    StrokeThickness = 1,
                    Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                    Fill = new SolidColorBrush(Windows.UI.Colors.Green),
                };
                MouseXInterceptBubbles.Add(trackingBubble);
                canvas.Children.Add(trackingBubble);
            }

            foreach (var m in MouseXInterceptBubbles)
                m.Visibility = Visibility.Collapsed;

            var index = 0;
            foreach (var xIntercept in xIntercepts)
            {
                var scaledXIntercept = ScaleY(xIntercept);

                MouseXInterceptBubbles[index].Visibility = Visibility.Visible;
                MouseXInterceptBubbles[index].Margin =
                    new Thickness(
                        TrackingMousePosition.X - BubbleSize / 2,
                        scaledXIntercept - BubbleSize / 2,
                        0,
                        0);

                index++;
            }
        }

        private void UpdateYIntercepts()
        {
            var currentY = UnscaleY(TrackingMousePosition.Y);

            var yIntercepts = ControlPoints.FindXForY(currentY);

            while (yIntercepts.Count() > MouseYInterceptBubbles.Count)
            {
                var trackingBubble = new Ellipse
                {
                    Width = BubbleSize,
                    Height = BubbleSize,
                    StrokeThickness = 1,
                    Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                    Fill = new SolidColorBrush(Windows.UI.Colors.Red),
                };
                MouseYInterceptBubbles.Add(trackingBubble);
                canvas.Children.Add(trackingBubble);
            }

            foreach (var m in MouseYInterceptBubbles)
                m.Visibility = Visibility.Collapsed;

            var index = 0;
            foreach (var yIntercept in yIntercepts)
            {
                var scaledYIntercept = ScaleX(yIntercept);

                MouseYInterceptBubbles[index].Visibility = Visibility.Visible;
                MouseYInterceptBubbles[index].Margin =
                    new Thickness(
                        scaledYIntercept - BubbleSize / 2,
                        TrackingMousePosition.Y - BubbleSize / 2,
                        0,
                        0);

                index++;
            }
        }

        private void Canvas_PointerMoved(object sender, Windows.UI.Xaml.Input.PointerRoutedEventArgs e)
        {
            TrackingMousePosition = e.GetCurrentPoint(canvas).Position;

            if (MouseTracking)
                UpdateMouseTrackingLines();

            if (ShowMouseXIntercepts)
                UpdateXIntercepts();

            if (ShowMouseYIntercepts)
                UpdateYIntercepts();
        }
    }
}

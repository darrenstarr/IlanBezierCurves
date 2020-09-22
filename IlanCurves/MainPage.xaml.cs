using IlanCurves.CurveFunctions;
using Microsoft.Toolkit.Uwp.UI.Controls;
using System;
using System.Collections.Generic;
using System.Linq;
using Windows.Foundation;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Shapes;

// The Blank Page item template is documented at https://go.microsoft.com/fwlink/?LinkId=402352&clcid=0x409

namespace IlanCurves
{
    /// <summary>
    /// An empty page that can be used on its own or navigated to within a Frame.
    /// </summary>
    public sealed partial class MainPage : Page
    {
        /// <summary>
        /// Size of the big bubbles on the UI
        /// </summary>
        public double BubbleSize { get; set; } = 9.0;

        /// <summary>
        /// Size of the little bubbles on the UI
        /// </summary>
        public double LittleBubbleSize { get; set; } = 5;

        private readonly Windows.UI.Color[] BubbleColors = new Windows.UI.Color[]
        {
            Windows.UI.Colors.Yellow,
            Windows.UI.Colors.Orange,
            Windows.UI.Colors.HotPink
        };

        private Line trackingLine;
        private Ellipse trackingBubble;

        public MainPage()
        {
            InitializeComponent();
            PaintCanvas();
        }

        private void AddSampleClicked(object sender, RoutedEventArgs e)
        {
            var viewModel = DataContext as MainViewModel;

            var newSample = new GraphSample
            {
                DayNumber = viewModel.NewSampleDay,
                Value = viewModel.NewSampleValue
            };

            viewModel.AddSample(newSample);

            viewModel.NewSampleDay = 0;
            viewModel.NewSampleValue = 0.0;

            PaintCanvas();
        }

        private void canvas_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            PaintCanvas();
        }

        private void dataGrid_CellEditEnded(object sender, DataGridCellEditEndedEventArgs e)
        {
            PaintCanvas();
        }

        private IEnumerable<Point> SamplesPoints
        {
            get
            {
                return (DataContext as MainViewModel).Samples
                    .Select(x =>
                        new Point { X = x.DayNumber, Y = x.Value }
                    );
            }
        }

        private PointCollection ToPointCollection(IEnumerable<Point> points)
        {
            // TODO: There must be a better way than this
            var result = new PointCollection();
            foreach (var p in points)
                result.Add(p);

            return result;
        }

        private void PaintCanvas()
        {
            if (canvas == null || canvas.Children == null)
                return; 

            var samples = SamplesPoints;

            var minX = samples.Select(x => x.X).Min();
            var maxX = samples.Select(x => x.X).Max();
            var minY = samples.Select(x => x.Y).Min();
            var maxY = samples.Select(x => x.Y).Max();

            var horizontalExtent = maxX;
            var verticalExtent = maxY;

            trackingLine = null;
            trackingBubble = null;
            canvas.Children.Clear();
            
            // Scale all points by these values
            var scaleX = canvas.ActualWidth / horizontalExtent * 0.9;
            var scaleY = canvas.ActualHeight / verticalExtent * 0.9;

            // Translate all points by the values
            var translateX = (canvas.ActualWidth - (horizontalExtent * scaleX)) / 2.0;
            var translateY = (canvas.ActualHeight - (verticalExtent * scaleY)) / 2.0;

            // Paint a border to the grid
            canvas.Children.Add(new Rectangle
            {
                Width = horizontalExtent * scaleX,
                Height = verticalExtent * scaleY,
                Margin = new Thickness(translateX, translateY, 0, 0),
                StrokeThickness = 1,
                Stroke = new SolidColorBrush(Windows.UI.Colors.LightBlue),
            });

            // Paint the vertical grid lines
            for (double i = Math.Min(minX, 0.0); i <= maxX; i += 1.0)
            {
                canvas.Children.Add(
                    new Line
                    {
                        X1 = translateX + i * scaleX,
                        Y1 = canvas.ActualHeight - translateY,
                        X2 = translateX + i * scaleX,
                        Y2 = canvas.ActualHeight - translateY - verticalExtent * scaleY,
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(Windows.UI.Colors.LightBlue)
                    });
            }

            // Paint the horizontal grid lines
            for (double i = Math.Min(minY, 0.0); i <= maxY; i += 1.0)
            {
                canvas.Children.Add(
                    new Line
                    {
                        X1 = translateX,
                        Y1 = translateY + i * scaleY,
                        X2 = translateX + horizontalExtent * scaleX,
                        Y2 = translateY + i * scaleY,
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(Windows.UI.Colors.LightBlue)
                    });
            }

            // Scale the samples to line up with the grid
            var scaledSamples = samples.Select(p => new Point(p.X * scaleX + translateX, canvas.ActualHeight - p.Y * scaleY - translateY));

            if ((bool)ShowLinear.IsChecked)
            {
                // Draw a polyline through all the scaled points
                canvas.Children.Add(new Polyline
                {
                    Points = ToPointCollection(scaledSamples),
                    StrokeThickness = 1,
                    Stroke = new SolidColorBrush(Windows.UI.Colors.Black)
                });
            }

            var ScaledCurve = scaledSamples.SplineInterpolate();            

            // Draw bubbles for each scaled point as well as connectors for control points to anchors
            if ((bool)ShowControlPoints.IsChecked)
            {
                for (var t=1; t<ScaledCurve.Count; t++)
                { 
                    if ((t % 3) != 0)
                    {
                        var p = ScaledCurve[t];
                        var anchor = (t % 3) == 1 ? ScaledCurve[t - 1] : ScaledCurve[t + 1];

                        canvas.Children.Add(new Line
                            {
                                X1 = p.X,
                                Y1 = p.Y,
                                X2 = anchor.X,
                                Y2 = anchor.Y,
                                StrokeThickness = 1,
                                Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                            }
                            );

                        canvas.Children.Add(
                            new Ellipse
                            {
                                Width = LittleBubbleSize,
                                Height = LittleBubbleSize,
                                Margin = new Thickness(p.X - LittleBubbleSize / 2, p.Y - LittleBubbleSize / 2, 0, 0),
                                StrokeThickness = 1,
                                Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                                Fill = new SolidColorBrush(BubbleColors[t % 3]),
                            }
                            );
                    }
                }
            }

            // Draw bubbles for each scaled point
            foreach (var p in scaledSamples)
            {
                canvas.Children.Add(
                    new Ellipse
                    {
                        Width = BubbleSize,
                        Height = BubbleSize,
                        Margin = new Thickness(p.X - BubbleSize / 2, p.Y - BubbleSize / 2, 0, 0),
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                        Fill = new SolidColorBrush(Windows.UI.Colors.Yellow),
                    }
                    );
            }

            // Draw the bezier curve if it's selected
            if ((bool)ShowCurve.IsChecked)
            {
                var pathSegmentCollection1 = new PathSegmentCollection();
                for (var i = 1; i < ScaledCurve.Count; i += 3)
                {
                    pathSegmentCollection1.Add(new BezierSegment
                    {
                        Point1 = ScaledCurve[i],
                        Point2 = ScaledCurve[i + 1],
                        Point3 = ScaledCurve[i + 2]
                    });
                }

                var path1 = new Windows.UI.Xaml.Shapes.Path
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
                                        StartPoint = ScaledCurve[0],
                                        Segments = pathSegmentCollection1
                                    }
                                }
                            }
                        }
                    }
                };

                canvas.Children.Add(path1);
            }

            if ((bool)ShowTestPoints.IsChecked)
            {
                var curve = ToPointCollection(samples.SplineInterpolate());
                for (var n = minX; n <= maxX; n += 0.19)
                {
                    var ys = curve.FindYForX(n);
                    if (ys.Count() == 0)
                        continue;

                    var y = canvas.ActualHeight - (ys.First() * scaleY + translateY);
                    var x = n * scaleX + translateX;

                    canvas.Children.Add(
                        new Ellipse
                        {
                            Width = LittleBubbleSize,
                            Height = LittleBubbleSize,
                            Margin = new Thickness(x - LittleBubbleSize / 2, y - LittleBubbleSize / 2, 0, 0),
                            StrokeThickness = 1,
                            Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                            Fill = new SolidColorBrush(Windows.UI.Colors.Green)
                        }
                        );
                }
            }
        }

        private void CheckBox_Checked(object sender, RoutedEventArgs e)
        {
            PaintCanvas();
        }


        private void canvas_PointerMoved(object sender, PointerRoutedEventArgs e)
        {
            // this is super slow. there must be a better way
            var pt = e.GetCurrentPoint(canvas);

            var samples = SamplesPoints;

            var minX = samples.Select(x => x.X).Min();
            var maxX = samples.Select(x => x.X).Max();
            var minY = samples.Select(x => x.Y).Min();
            var maxY = samples.Select(x => x.Y).Max();

            var horizontalExtent = maxX;
            var verticalExtent = maxY;

            // Scale all points by these values
            var scaleX = canvas.ActualWidth / horizontalExtent * 0.9;
            var scaleY = canvas.ActualHeight / verticalExtent * 0.9;

            // Translate all points by the values
            var translateX = (canvas.ActualWidth - (horizontalExtent * scaleX)) / 2.0;
            var translateY = (canvas.ActualHeight - (verticalExtent * scaleY)) / 2.0;


            // Draw a line showing the X position of the mouse.
            if (trackingLine == null)
            {
                trackingLine = new Line
                {
                    X1 = pt.Position.X,
                    Y1 = 0,
                    X2 = pt.Position.X,
                    Y2 = canvas.ActualHeight,
                    StrokeThickness = 0.75,
                    Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                };
                canvas.Children.Add(trackingLine);
            }
            else
            {
                trackingLine.X1 = pt.Position.X;
                trackingLine.X2 = pt.Position.X;
            }

            var xPos = (pt.Position.X - translateX) / scaleX;

            // Draw a bubble on the Y intersect of the lowest value. This should be extended to show all Y values a the intersection
            var ys = ToPointCollection(samples.SplineInterpolate()).FindYForX(xPos);
            if (ys.Count() > 0)
            {
                var minYValue = ys.Min();
                var p = new Point(xPos * scaleX + translateX, canvas.ActualHeight - minYValue * scaleY - translateY);

                if (trackingBubble == null)
                {
                    trackingBubble = new Ellipse
                    {
                        Width = BubbleSize,
                        Height = BubbleSize,
                        Margin = new Thickness(p.X - BubbleSize / 2, p.Y - BubbleSize / 2, 0, 0),
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(Windows.UI.Colors.Black),
                        Fill = new SolidColorBrush(Windows.UI.Colors.Green),
                    };
                    canvas.Children.Add(trackingBubble);
                }
                else
                {
                    trackingBubble.Margin = new Thickness(p.X - BubbleSize / 2, p.Y - BubbleSize / 2, 0, 0);
                }
            }
        }

        private void dataGrid_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            var grid = sender as DataGrid;

            RemoveSample.IsEnabled = grid.SelectedItems.Count != 0;
        }

        private void RemoveItemClicked(object sender, RoutedEventArgs e)
        {
            var viewModel = DataContext as MainViewModel;

            var toDelete = new List<GraphSample>();
            foreach (var item in dataGrid.SelectedItems)
                toDelete.Add(item as GraphSample);

            foreach (var item in toDelete)
                viewModel.Samples.Remove(item);

            PaintCanvas();
        }
    }
}

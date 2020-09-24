using IlanCurves.Geometry;
using Microsoft.Toolkit.Uwp.UI.Controls;
using System.Collections.Generic;
using Windows.Foundation;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Media;

// The Blank Page item template is documented at https://go.microsoft.com/fwlink/?LinkId=402352&clcid=0x409

namespace IlanCurves
{
    /// <summary>
    /// An empty page that can be used on its own or navigated to within a Frame.
    /// </summary>
    public sealed partial class MainPage : Page
    {
        public MainPage()
        {
            InitializeComponent();
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

            UpdateBezierView();
        }

        void UpdateBezierView()
        {
            if (bezierView == null)
                return;

            var result = new PointCollection();
            foreach (var sample in (DataContext as MainViewModel).Samples)
                result.Add(
                    new Point
                    {
                        X = sample.DayNumber,
                        Y = sample.Value
                    });

            bezierView.ControlPoints = result.SplineInterpolate();
        }

        private void DataGrid_CellEditEnded(object sender, DataGridCellEditEndedEventArgs e)
        {
            UpdateBezierView();
        }

        private void DataGrid_SelectionChanged(object sender, SelectionChangedEventArgs e)
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

            UpdateBezierView();
        }

        private void ShowCurve_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowCurve = (bool)ShowCurve.IsChecked;
        }

        private void ShowLinear_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowLinear = (bool)ShowLinear.IsChecked;
        }

        private void ShowControlPoints_Checked(object sender, RoutedEventArgs e)
        {
            if(bezierView != null)
                bezierView.ShowControlPoints = (bool)ShowControlPoints.IsChecked;
        }

        private void ShowXIntersects_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowMouseXIntercepts = (bool)ShowXIntersect.IsChecked;
        }

        private void ShowYIntersects_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowMouseYIntercepts = (bool)ShowYIntersect.IsChecked;
        }

        private void ShowGrid_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowGrid = (bool)ShowGrid.IsChecked;
        }

        private void ShowMouseTrackingLines_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.MouseTracking = (bool)ShowMouseTrackingLines.IsChecked;
        }

        private void ShowShortestPath_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowShortestPath = (bool)ShowShortestPath.IsChecked;
        }

        private void ShowPointsAlongCurve_Checked(object sender, RoutedEventArgs e)
        {
            if (bezierView != null)
                bezierView.ShowPointsAlongCurve = (bool)ShowPointsAlongCurve.IsChecked;
        }
    }
}

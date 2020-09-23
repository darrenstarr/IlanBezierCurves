using IlanCurves.CurveFunctions;
using System;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Media;

namespace IlanCurves
{
    class GraphSampleToCurveConverter : IValueConverter
    {
        public object Convert(object value, Type targetType, object parameter, string language)
        {
            var result = new PointCollection();
            foreach (var sample in value as GraphSamples)
                result.Add(
                    new Windows.Foundation.Point
                    {
                        X = sample.DayNumber,
                        Y = sample.Value
                    });

            return result.SplineInterpolate();
        }

        public object ConvertBack(object value, Type targetType, object parameter, string language)
        {
            throw new NotImplementedException();
        }
    }
}

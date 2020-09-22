using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace IlanCurves
{
    public class GraphSample
    {
        public int DayNumber { get; set;  }
        public double Value { get; set; }
    }

    public class GraphSamples : ObservableCollection<GraphSample>, IEnumerable<GraphSample>
    {
    }
}

using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Windows.ApplicationModel;

namespace IlanCurves
{
    public class MainViewModel
    {
        public ObservableCollection<GraphSample> Samples { get; set; }

        public MainViewModel()
        {
            Samples = new ObservableCollection<GraphSample>();

            ///if(DesignMode.DesignModeEnabled)
            {
                Samples.Add(new GraphSample { DayNumber = 1, Value = 10 });
                Samples.Add(new GraphSample { DayNumber = 4, Value = 18 });
                Samples.Add(new GraphSample { DayNumber = 7, Value = 22 });
                Samples.Add(new GraphSample { DayNumber = 10, Value = 31 });
                Samples.Add(new GraphSample { DayNumber = 20, Value = 68 });
            }
        }
    }
}

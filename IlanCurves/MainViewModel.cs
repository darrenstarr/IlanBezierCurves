using System.Collections.ObjectModel;
using System.ComponentModel;

namespace IlanCurves
{
    public class MainViewModel : INotifyPropertyChanged
    {
        private int newSampleDay = 0;
        private double newSampleValue = 0;

        public GraphSamples Samples { get; set; }

        public int NewSampleDay
        {
            get => newSampleDay;
            set
            {
                newSampleDay = value;
                PropertyChanged?.Invoke(this, new PropertyChangedEventArgs("NewSampleDay"));
            }
        }

        public double NewSampleValue
        {
            get => newSampleValue;
            set
            {
                newSampleValue = value;
                PropertyChanged?.Invoke(this, new PropertyChangedEventArgs("NewSampleValue"));
            }
        }

        public void AddSample(GraphSample sample)
        {
            int i = 0;
            for (; i < Samples.Count; i++)
            {
                if (Samples[i].DayNumber > sample.DayNumber)
                    break;
            }
            Samples.Insert(i, sample);
        }

        public MainViewModel()
        {
            Samples = new GraphSamples();

            ///if(DesignMode.DesignModeEnabled)
            {
                Samples.Add(new GraphSample { DayNumber = 1, Value = 10 });
                Samples.Add(new GraphSample { DayNumber = 4, Value = 18 });
                Samples.Add(new GraphSample { DayNumber = 7, Value = 22 });
                Samples.Add(new GraphSample { DayNumber = 10, Value = 31 });
                Samples.Add(new GraphSample { DayNumber = 7, Value = 25 });
                Samples.Add(new GraphSample { DayNumber = 11, Value = 20 });
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;        
    }
}

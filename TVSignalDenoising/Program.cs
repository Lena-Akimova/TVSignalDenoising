
using System.IO;

namespace TVSignalDenoising
{
    class Program
    {
        static void Main(string[] args)
        {
            var location = System.Reflection.Assembly.GetExecutingAssembly().Location;
            var path = Path.GetDirectoryName(location);
            var p = path.Remove(path.IndexOf("bin"));

            var signal = FileIO.ReadSignal($"{p}Resources\\signal.txt");

            double[] x = new double[signal.Length];
            for(int i=0; i<signal.Length; i++)
                x[i] = signal[i] + 2;

            Regularisator reg = new Regularisator(2, signal);

            var x1 = reg.Projection(x, 20);
            reg.MinModelValue(x, x1);

        }
        
    }
}



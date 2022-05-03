
namespace TVSignalDenoising
{
    class Program
    {
        static void Main(string[] args)
        {
            var signal = FileIO.ReadSignal();

            double[] x = new double[signal.Length];
            for(int i=0; i<signal.Length; i++)
                x[i] = signal[i] + 2;

            Regularisator reg = new Regularisator(2, signal);

            var x1 = reg.Projection(x, 20);
            reg.MinModelValue(x, x1);

        }
        
    }
}



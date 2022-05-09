
using System;
using System.IO;

namespace TVSignalDenoising
{
    class Program
    {
        static void Main(string[] args)
        {
            //1.
            var ex1 = new Example05x2()
            {
                BoxUp = new double[] { 15 },
                BoxLow = new double[] { 10 },
                N = 1
            };
            var x01 = new double[] { 12 };
            var eps = 0.0001;
            var teta = 1-1/(2+Math.Sqrt(2));
            
            LevelSetOptimizator lsEx1 = new LevelSetOptimizator(ex1, -1, x01, eps, teta);

            var argmin1 = lsEx1.Minimize();
            var minVal1= ex1.GetValueAt(argmin1);
            Console.WriteLine($"\nОтвет: \neps = {eps} \nargmin = {alglib.ap.format(argmin1, 3)} \nminVal = {alglib.ap.format(minVal1.Value, 3)}  ");


            /*
            //2.
            var signal = FileIO.ReadSignal("");
            var noisedSignal = NoiseSignal(signal);
            var ex2 = new ExampleTV()
            {
                NoisedSignal = noisedSignal,
                N = signal.Length,
                Lambda = 2
            };

            var x02 = new double[signal.Length];
            LevelSetOptimizator lsEx2 = new LevelSetOptimizator(ex2, -1, x02, eps, teta);
            var argmin2 = lsEx2.Minimize();
            var minVal2= ex2.GetValueAt(argmin2);
            Console.WriteLine($"\nОтвет: \neps = {eps} \nargmin = {alglib.ap.format(argmin2, 3)} \nminVal = {alglib.ap.format(minVal2.Value, 3)}  ");
            */

        }

        private static double[] NoiseSignal(double[] signal, float k= 0.3f, float prob=0.03f)
        {
            var noised = new double[signal.Length];
            Random r = new Random();
            for(var i=0; i<signal.Length; i++)
            {
                var rand = (r.NextDouble() * 2 - 1) * k;
                if (r.NextDouble() < prob)
                    noised[i] = signal[i] + rand * 7;
                else
                    noised[i] = signal[i] + rand;
            }
            return noised;
        }

    }
}



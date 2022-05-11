
using System;
using System.IO;
using System.Linq;

namespace TVSignalDenoising
{
    class Program
    {
        static void Main(string[] args)
        {
            Example1();
            Example2();
        }





        static void Example1()
        {
            var ex1 = new Example05x2()
            {
                BoxUp = new double[] { 15 },
                BoxLow = new double[] { 10 },
                N = 1,
                Scale = 1
            };
            var x01 = new double[] { 12 };
            var eps = 0.0001;
            var teta = 0.7;// 1-1/(2+Math.Sqrt(2));

            LevelSetOptimizator lsEx1 = new LevelSetOptimizator(ex1, -1, x01, eps, teta);

            var argmin1 = lsEx1.Minimize();
            var minVal1 = ex1.GetValueAt(argmin1);
            Console.WriteLine($"\nОтвет: \neps = {eps} \nargmin = {alglib.ap.format(argmin1, 3)} \nminVal = {alglib.ap.format(minVal1.Value, 3)}  ");

        }


        static void Example2()
        {
            var eps = 0.000000001;
            var teta = 0.7;// 1-1/(2+Math.Sqrt(2));

            var ex2 = new ExampleTV() { Lambda = 1, Scale = 0.000000000000001 };
            var bUp = new double[ex2.N];
            var bLow = new double[ex2.N];
            var x0 = ex2.NoisedSignal;
            var mean = ex2.NoisedSignal.Sum() / ex2.N;

            x0 = new double[ex2.N];

            for (int i = 0; i < bUp.Length; i++)
            {
                bUp[i] = double.PositiveInfinity;
                bLow[i] = -double.PositiveInfinity;
                x0[i] = mean;
            }

            ex2.BoxUp = bUp;
            ex2.BoxLow = bLow;

            LevelSetOptimizator lsEx2 = new LevelSetOptimizator(ex2, -1, x0, eps, teta);
            ex2.DenoisedSignal = lsEx2.Minimize();
            var minTVVal = ex2.GetValueAt(ex2.DenoisedSignal);
            Console.WriteLine($"\nОтвет: \neps = {eps} \nargmin = {alglib.ap.format(ex2.DenoisedSignal, 3)} \n minTVVal = {alglib.ap.format(minTVVal.Value, 3)}  ");

            Display1DSignalTVRegResults(ex2);

        }

        static void Display1DSignalTVRegResults(ExampleTV ex)
        {

            var xPoints = new double[ex.N];
            for(int i=0; i<ex.N; i++)
                xPoints[i] = i;
            dislin.metafl("xwin");
            dislin.disini();
            dislin.titlin("1D SIGNAL TV-REGULARISATION", 1);
            dislin.titlin($"RMSE origin-noised {ExampleTV.RMSE(ex.NoisedSignal, ex.Signal)}", 2);
            dislin.titlin($"RMSE origin-tvRegulaised {ExampleTV.RMSE(ex.DenoisedSignal, ex.Signal)}", 3);

            dislin.graf(0, ex.N, 0.0, 20, -10.0, 10.0, -10.0, 1);
            dislin.title();

            //шумный
            dislin.color("red");
            dislin.curve(xPoints, ex.NoisedSignal, ex.N);
            //оригинальный
            dislin.color("green");
            dislin.curve(xPoints, ex.Signal, ex.N);
            //фильтрованный
            dislin.color("blue");
            dislin.curve(xPoints, ex.DenoisedSignal, ex.N);
            dislin.disfin();
        }





    }
}



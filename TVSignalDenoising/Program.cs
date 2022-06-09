using System;
using System.Diagnostics;
using System.Linq;


namespace TVSignalDenoising
{
    class Program
    {
        static void Main(string[] args)
        {
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

            var repeating = 100;
            var eps = 0.000001;
            var teta = 0.7;// 1-1/(2+Math.Sqrt(2));
            var Fkl = new (double[], int)[repeating];
            var Fkl1 = new (double[], int)[repeating];
            var Fks = new (double[], int)[repeating];
            var time1 = new double[repeating];
            var time2 = new double[repeating];
            var time3 = new double[repeating];
            var steps1 = new int[repeating];
            var steps2 = new int[repeating];
            var steps3 = new int[repeating];
            var min1 = new double[repeating];
            var min2 = new double[repeating];
            var min3 = new double[repeating];

            for (int e = 0; e < repeating; e++)
            {
                var x01 = GetRandomPoint(ex1.BoxLow, ex1.BoxUp, ex1.N);
                x01 =new double[] { 15 };
                Console.WriteLine($"x0: {alglib.ap.format(x01, 3)} ");

                //ПМУ с отсек. плоскостями
                LevelSetOptimizator lsEx1 = new LevelSetOptimizator(ex1, -1, x01, eps, teta, true);
                lsEx1.RemoveCuttingPlanes = true;

                Stopwatch stopWatch = Stopwatch.StartNew();
                var res1 = lsEx1.Minimize();
                stopWatch.Stop();

                var argmin1 = res1.Item1[res1.Item3];
                Fkl[e] = (res1.Item2, res1.Item3);
                var minVal1 = ex1.GetValueAt(argmin1);
                time1[e] = stopWatch.ElapsedMilliseconds / 1000.0;
                steps1[e] = res1.Item4;
                min1[e] = minVal1.Value;

                Console.WriteLine($"\nПМУ с отсек. плоскостями:\nОтвет: \neps = {eps} \nШагов = {res1.Item4} " +
                    $"\nargmin = {alglib.ap.format(argmin1, 3)} " +
                    $"\nminVal = {alglib.ap.format(minVal1.Value, 3)} " +
                    $"\nЗатраченное время = {time1[e]} сек ");
                PrintBorder();

                //ПМУ 
                lsEx1.RemoveCuttingPlanes = false;

                stopWatch.Restart();
                var res3 = lsEx1.Minimize();
                stopWatch.Stop();

                var argmin3 = res3.Item1[res3.Item3];
                Fkl1[e] = (res3.Item2, res3.Item4);
                var minVal3 = ex1.GetValueAt(argmin3);
                time3[e] = stopWatch.ElapsedMilliseconds / 1000.0;
                steps3[e] = res3.Item4;
                min3[e] = minVal3.Value;

                Console.WriteLine($"\nПМУ:\nОтвет: \neps = {eps} \nШагов = {res3.Item4} " +
                    $"\nargmin = {alglib.ap.format(argmin3, 3)} " +
                    $"\nminVal = {alglib.ap.format(minVal3.Value, 3)} " +
                    $"\nЗатраченное время = {stopWatch.ElapsedMilliseconds / 1000.0} сек ");
                PrintBorder();

                //Субградиентный
                var k = 100;
                var h = 0.01;
                SubgradientDescentOptimizator sd = new SubgradientDescentOptimizator(ex1, k, x01, h, eps);

                stopWatch.Restart();
                var res2 = sd.Minimize();
                stopWatch.Stop();
                time2[e] = stopWatch.ElapsedMilliseconds / 1000.0;

                Fks[e] = (res2.Item2, res2.Item4);
                var argmin2 = res2.Item1[res2.Item3];
                var minVal2 = ex1.GetValueAt(argmin2);
                steps2[e] = res2.Item4;
                min2[e] = minVal2.Value;

                Console.WriteLine($"\nСубградиентный:\nОтвет: \neps = {eps} \nШагов = {res2.Item4}  " +
                    $"\nargmin = {alglib.ap.format(argmin2, 3)} " +
                    $"\nminVal = {alglib.ap.format(minVal2.Value, 3)} " +
                    $"\nЗатраченное время = {stopWatch.ElapsedMilliseconds/1000.0} сек ");
                PrintBorder();
            }

            PrintBorder();
            Console.WriteLine($"Средние результаты за {repeating} итераций");
            PrintBorder();
            Console.WriteLine($"\nПМУ с отсек. плоскостями:" +
                    $"\nОтвет: " +
                    $"\neps = {eps} " +
                    $"\nСредний минимум = {min1.Sum() / repeating}" +
                    $"\nСреднее количество шагов = {steps1.Sum()/(double)repeating} " +
                    $"\nСреднее затраченное время = {time1.Sum() / (double)repeating} сек ");

            Console.WriteLine($"\nПМУ :" +
                $"\nОтвет: " +
                $"\neps = {eps} " +
                    $"\nСредний минимум = {min3.Sum()/repeating}" +
                    $"\nСреднее количество шагов = {steps3.Sum() / (double)repeating} " +
                    $"\nСреднее затраченное время = {time3.Sum() / (double)repeating} сек ");

            Console.WriteLine($"\nСубградиентный:" +
                $"\nОтвет: " +
                $"\neps = {eps} " +
                    $"\nСредний минимум = {min2.Sum() / repeating}" +
                    $"\nСреднее количество шагов = {steps2.Sum() / (double)repeating}  " +
                    $"\nСреднее затраченное время = {time2.Sum()/ (double)repeating} сек ");
            PrintBorder();

            DisplayFuncDeskending(Fkl, "Proximal Level-Set Method with cutting planes", ex1.GetValueAt(ex1.BoxLow).Value, ex1.GetValueAt(ex1.BoxUp).Value, 1,100);
            DisplayFuncDeskending(Fkl1, "Proximal Level-Set Method", ex1.GetValueAt(ex1.BoxLow).Value, ex1.GetValueAt(ex1.BoxUp).Value, 10,100);
            DisplayFuncDeskending(Fks, "Subgradient Method", ex1.GetValueAt(ex1.BoxLow).Value, ex1.GetValueAt(ex1.BoxUp).Value, 100,100 );


        }


        static void Example2()
        {
            var eps = 0.0000001;
            var teta = 0.7;// 1-1/(2+Math.Sqrt(2));
            var repeating = 1;
            var ex2 = new ExampleTV() { Lambda = 1, Scale = 0.00000000000000000001 };
            var bUp = new double[ex2.N];
            var bLow = new double[ex2.N];
            var time1 = new double[repeating];
            var time2 = new double[repeating];
            var time3 = new double[repeating];
            var steps1 = new int[repeating];
            var steps2 = new int[repeating];
            var steps3 = new int[repeating];
            var min1 = new double[repeating];
            var min2 = new double[repeating];
            var min3 = new double[repeating];
            var Fks = new (double[], int)[repeating];
            var Fkl = new (double[], int)[repeating];
            var Fkl1 = new (double[], int)[repeating];
            var SNR = 1.5;

            for (int i = 0; i < bUp.Length; i++)
            {
                bUp[i] = double.PositiveInfinity;
                bLow[i] = -double.PositiveInfinity;
            }

            ex2.BoxUp = bUp;
            ex2.BoxLow = bLow;


            for (int e = 0; e < repeating; e++) {
                PrintBorder();
                Console.WriteLine($"* ЭКСПЕРИМЕНТ {e+1}");
                PrintBorder();

                var x0 = ExampleTV.NoiseSignal(ex2.Signal, SNR);
                ex2.NoisedSignal = x0;

                //ПМУ с отсек. плоскостями
                LevelSetOptimizator lsEx2 = new LevelSetOptimizator(ex2, -1, x0, eps, teta, true);

                Stopwatch stopWatch = Stopwatch.StartNew();
                var res1 = lsEx2.Minimize();
                stopWatch.Stop();

                Fkl[e] = (res1.Item2, res1.Item3);
                ex2.DenoisedSignal = res1.Item1[res1.Item3];
                var minTVVal = ex2.GetValueAt(ex2.DenoisedSignal);
                time1[e] = stopWatch.ElapsedMilliseconds / 1000.0;
                steps1[e] = res1.Item4;
                min1[e] = minTVVal.Value;

                Console.WriteLine($"\nПМУ с отсек. плоскостями:" +
                    $"\nОтвет: " +
                    $"\nEps = {eps} " +
                    $"\nШагов = {res1.Item4} " +
                    $"\nArgmin = {alglib.ap.format(ex2.DenoisedSignal, 3)} " +
                    $"\nMin = {alglib.ap.format(minTVVal.Value, 3)} " +
                    $"\nЗатраченное время = {stopWatch.ElapsedMilliseconds / 1000.0} сек ");

                if (e == repeating - 1)
                    Display1DSignalTVRegResults(ex2, "Proximal level-set method with cutting planes removing");

                PrintBorder();

                //ПМУ 

                lsEx2.RemoveCuttingPlanes = false;

                stopWatch.Restart();
                var res3 = lsEx2.Minimize();
                stopWatch.Stop();

                var argmin3 = res3.Item1[res3.Item3];
                Fkl1[e] = (res3.Item2, res3.Item4);
                ex2.DenoisedSignal = res3.Item1[res3.Item3];
                var minVal3 = ex2.GetValueAt(argmin3);
                time3[e] = stopWatch.ElapsedMilliseconds / 1000.0;
                steps3[e] = res3.Item4;
                min3[e] = minVal3.Value;

                Console.WriteLine($"\nПМУ:" +
                    $"\nОтвет: " +
                    $"\nEps = {eps} " +
                    $"\nШагов = {res3.Item4} " +
                    $"\nArgmin = {alglib.ap.format(argmin3, 3)} " +
                    $"\nMin = {alglib.ap.format(minVal3.Value, 3)} " +
                    $"\nЗатраченное время = {stopWatch.ElapsedMilliseconds / 1000.0} сек ");

                if (e == repeating - 1)
                    Display1DSignalTVRegResults(ex2, "Proximal level-set method");
                PrintBorder();


                //Субградиентный
                var k = 100;
                var h = 0.01;
                SubgradientDescentOptimizator sd = new SubgradientDescentOptimizator(ex2, k, x0, h, eps);

                stopWatch.Restart();
                var res2 = sd.Minimize();
                stopWatch.Stop();

                Fks[e] = (res2.Item2, res2.Item4);
                var argmin2 = res2.Item1[res2.Item3];
                ex2.DenoisedSignal = argmin2;
                var minVal2 = ex2.GetValueAt(argmin2);
                time2[e] = stopWatch.ElapsedMilliseconds / 1000.0;
                steps2[e] = res2.Item4;
                min2[e] = minVal2.Value;

                Console.WriteLine($"\nСубградиентный:" +
                    $"\nОтвет: " +
                    $"\nEps = {eps} " +
                    $"\nШагов = {res2.Item4}  " +
                    $"\nArgmin = {alglib.ap.format(argmin2, 3)} " +
                    $"\nMin = {alglib.ap.format(minVal2.Value, 3)} " +
                    $"\nЗатраченное время = {stopWatch.ElapsedMilliseconds / 1000.0} сек ");

                if (e == repeating - 1)
                    Display1DSignalTVRegResults(ex2, "Subgradient method");
                PrintBorder();
            }

            PrintBorder();
            Console.WriteLine($"* Средние результаты за {repeating} итераций ");
            PrintBorder();

            Console.WriteLine($"\nПМУ с отсек. плоскостями:" +
                $"\nОтвет: " +
                $"\nEps = {eps} " +
                    $"\nСредний минимум = {min1.Sum() / repeating}" +
                    $"\nСреднее количество шагов = {steps1.Sum() / (double)repeating} " +
                    $"\nСреднее затраченное время = {time1.Sum() / (double)repeating} сек ");

            Console.WriteLine($"\nПМУ :" +
                $"\nОтвет: " +
                $"\nEps = {eps} " +
                    $"\nСредний минимум = {min3.Sum() / repeating}" +
                    $"\nСреднее количество шагов = {steps3.Sum() / (double)repeating} " +
                    $"\nСреднее затраченное время = {time3.Sum() / (double)repeating} сек ");

            Console.WriteLine($"\nСубградиентный:" +
                $"\nОтвет: " +
                $"\nEps = {eps} " +
                    $"\nСредний минимум = {min2.Sum() / repeating}" +
                    $"\nСреднее количество шагов = {steps2.Sum() / (double)repeating}  " +
                    $"\nСреднее затраченное время = {time2.Sum() / (double)repeating} сек ");


            DisplayFuncDeskending(Fkl, "Proximal Level-Set Method with cutting planes", 0, 100, 3, 30);
            DisplayFuncDeskending(Fkl1, "Proximal Level-Set Method", 0, 200, 1, 30);
            DisplayFuncDeskending(Fks, "Subgradient Method", 0, 100, 1000, 30);


        }

        /// <summary>
        /// Случайая тока их области определения
        /// </summary>
        static double[] GetRandomPoint(double[] boxLow, double[] boxUp, int n)
        {
            Random r = new Random();
            double[] x = new double[n];
            for(int i=0;i<n; i++)
            {
                x[i] = boxLow[i]+(boxUp[i]-boxLow[i])*r.NextDouble();
            }
            return x;
        }


        /// <summary>
        /// Средние убывающие значения функции на каждом шаге
        /// </summary>
        /// <param name="fk"></param>
        /// <returns></returns>
        static double[] MeanFk((double[], int)[] fk, int len)
        {
            var mean = new double[len];

            for(int i=0; i< len; i++)
            {
                var sum = 0.0;
                for(int j=0; j< fk.Length; j++)
                {
                    sum+=fk[j].Item1[i];
                }
                mean[i] = sum / fk.Length;
            }
            return mean;
        }

        /// <summary>
        /// Стандартное отклонение убывающих значений функции на каждом шаге
        /// </summary>
        /// <returns></returns>
        static double[] StdFk((double[], int)[] fk, int len)
        {
            var mean = MeanFk(fk, len);
            var n = mean.Length;
            var std = new double[n];

            for (int i = 0; i < n; i++)
            {
                var sum = 0.0;
                for (int j = 0; j < fk.Length; j++)
                {
                    sum += Math.Pow(fk[j].Item1[i]-mean[i], 2);
                }
                std[i] = Math.Sqrt(sum / fk.Length);
            }
            return std;
        }

        static void DisplayPoints( double[][] Xk,  string title, double xstep)
        {
            int k = Xk.Length;
            var xPoints = new double[k];
            var yPoints = new double[k];
            for (int s = 0; s < 8; s++)
                xPoints[s] = Xk[s][0];

            for (int s = 0; s < 8; s++)
                yPoints[s] = Xk[s][1];

            dislin.metafl("xwin");
            dislin.disini();
            dislin.titlin(title, 1);
            var xstart = xPoints.Min()-10;
            var xend = xPoints.Max() + 10;

            var ystart = yPoints.Min() - 10;
            var yend = yPoints.Max() + 10;

            dislin.graf(xstart, xend, xstart, xstep, ystart, yend, ystart, 1000);
            dislin.title();

            dislin.color("white");
            dislin.curve(xPoints, yPoints, k);

            dislin.disfin();
        }


        /// <summary>
        /// Вывод средних значений цф по итерациям
        /// </summary>
        /// <param name="Fk"></param>
        /// <param name="title"></param>
        /// <param name="yl"></param>
        /// <param name="yu"></param>
        /// <param name="xstep"></param>
        /// <param name="ystep"></param>
        static void DisplayFuncDeskending((double[], int)[] Fk, string title, double yl, double yu, double xstep,double ystep)
        {
            var averShod = Fk.Select(x => x.Item2).Average();
            Fk = Fk.Where(x => x.Item2 >= averShod).ToArray();
            var resIndex = Fk.Min(x => x.Item2);
            var xPoints = new double[resIndex+1];
            for (int s = 0; s < resIndex + 1; s++)
                xPoints[s] = s;

            dislin.metafl("xwin");
            dislin.scrmod("REVERSE");
            dislin.disini();
            dislin.name("F(x)", "Y");
            dislin.name("Iterations", "X");
            dislin.titlin(title, 1);

            dislin.graf(0, resIndex+1, 0.0, xstep, yl, yu, yl, ystep);
            dislin.title();

            var len = Math.Min(Fk.Length, resIndex + 1);
            var meanF = MeanFk(Fk, len);
            dislin.color("red");
            dislin.curve(xPoints, meanF, len);

            var std = StdFk(Fk, len);
            var upMean = meanF.Zip(std).Select(ms => ms.First + ms.Second).ToArray();
            var lowMean = meanF.Zip(std).Select(ms => ms.First - ms.Second).ToArray();

            dislin.color("blue");
            dislin.curve(xPoints, upMean, resIndex + 1);
            dislin.curve(xPoints, lowMean, resIndex + 1);

            dislin.disfin();
        }


        /// <summary>
        /// Представление сигналов
        /// </summary>
        /// <param name="ex"></param>
        static void Display1DSignalTVRegResults(ExampleTV ex, string title)
        {

            var xPoints = new double[ex.N];
            for(int i=0; i<ex.N; i++)
                xPoints[i] = i;
            dislin.metafl("xwin");
            dislin.scrmod("REVERSE");
            dislin.disini();
            dislin.titlin(title, 1);
            dislin.titlin($"RMSE origin-noised {ExampleTV.RMSE(ex.NoisedSignal, ex.Signal)}", 2);
            dislin.titlin($"RMSE origin-tvRegulaised {ExampleTV.RMSE(ex.DenoisedSignal, ex.Signal)}", 3);

            dislin.titlin($"TV origin = {ExampleTV.TV(ex.Signal)} TV noised = {ExampleTV.TV(ex.NoisedSignal)} TV tvRegulaised = {ExampleTV.TV(ex.DenoisedSignal) } ", 4); 

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

        static void PrintBorder(int len=100)
        {
            string s = "";
            for (int i = 0; i < len; i++)
                s += "*";
            Console.WriteLine(s);
        }
    }
}



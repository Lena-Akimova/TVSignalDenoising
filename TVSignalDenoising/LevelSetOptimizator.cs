using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TVSignalDenoising
{
    /// <summary>
    /// Оптимизатор проксимальным методом уровней
    /// </summary>
    public class LevelSetOptimizator: IOptimizator
    {
        private Example example;
        private double alfa; //начальное мин. значение модели функции
        private double[] x0; //первая точка модели функции
        private double eps; //требуемая точность
        private double teta; //коэффициент уровня
        public bool CuttingPlanes { get; set; }


        public LevelSetOptimizator(Example ex, double _alfa, double[] _x0, double _eps, double _teta, bool _cuttingPlanes)
        {
            example = ex;
            alfa = _alfa;
            x0 = _x0;
            eps = _eps;
            teta = _teta;
            CuttingPlanes = _cuttingPlanes;
        }

        public (double[][], double[], int, int) Minimize()
        {
            if (CuttingPlanes)
                return MinimizeWithCuttingPlanes();
            else
                return MinimizeWithoutCuttingPlanes();
        }

        /// <summary>
        /// Минимизация функции для <see cref="Example">примера</see>
        /// </summary>
        /// <returns></returns>
        public (double[][], double[], int, int) MinimizeWithCuttingPlanes()
        {

            var n = example.N;

            //var steps = FindMaxStepsCount();
            //Console.WriteLine($"Число шагов ПМУ не превышает {steps}\n"); //слишком большое

            //var c = (int)-Math.Log(eps, 10);
            var steps = 15 * example.N;
            var k = 0;
            var X = new double[steps][];
            var F = new double[steps];
            var Xk = new double[steps][];
            var Fk = new double[steps];
            var Ak = new double[steps];
            var Bk = new double[steps];
            var Sk = new double[steps][];
            var Lk = new double[steps];
            var i = 0;
            X[k] = x0;
            var A = alfa;
            var B = example.GetValueAt(X[k]).Value; //начальное рекордное значение модели ф-ии
            F[k] = B;
            var counter = 0;
            //Console.WriteLine($"x0={alglib.ap.format(X[k], 3)} \nf0= {B} ");

            try
            {
                //Пока ошибка модели больше требуемой точности
                while (B - A > eps)
                {
                    //Console.WriteLine($"k = {k}:");
                    Xk = new double[steps][];
                    Fk = new double[steps];
                    Ak = new double[steps];
                    Bk = new double[steps];
                    Sk = new double[steps][];
                    Lk = new double[steps];

                    i = 0;
                    Xk[i] = X[k];
                    Ak[i] = A;
                    Bk[i] = B;

                    //Уровень - между мин. и рекордным значениями
                    Lk[i] = (1 - teta) * Ak[i] + teta * Bk[i];
                    //Console.WriteLine($"L{k}{i} = {Lk[i]}");

                    //Касательная в последней построенной точке модели ф-ии
                    var S = example.GetSubGradAt(X[k]);
                    Sk[i] = S;
                    Fk[i] = example.GetValueAt(Xk[i]).Value;
                    var sum = Fk[i] + Sk[i].Zip(Xk[i]).Sum(sx => -sx.First * sx.Second);

                    double[,] linCon = new double[1, n];
                    for (int ii = 0; ii < n; ii++)
                        linCon[0, ii] = Sk[i][ii];

                    Ak[i + 1] = Lk[i] - 1;// чтобы зайти в цикл


                    //Если мин значение новой модели меньше уровня старой модели
                    while (Ak[i + 1] < Lk[i])
                    {
                        counter++;
                        //Console.WriteLine($"i = {i}:");
                        //3. Задаем множество уровней U
                        //строим проекцию на него - ищем точку, которая ближе всего к последней xk,
                        //но значение касательной к xk в ней не больше уровня
                        Xk[i + 1] = Projection(X[k], linCon, example.BoxUp, example.BoxLow, new double[] { Lk[i] - sum });

                        //Console.WriteLine($"Проекция точки {alglib.ap.format(X[k], 3)} \nна множество U{k},{i}: \n{alglib.ap.format(Xk[i + 1], 3)} ");
                        //Новая модель функции с большим на 1 кол-вом точек

                        //4. Проводим касательную в этой новой точке
                        Sk[i + 1] = example.GetSubGradAt(Xk[i + 1]);
                        Fk[i + 1] = example.GetValueAt(Xk[i + 1]).Value;

                        //Пересчитаем мин значение модели функции
                        var au1 = Fk[i] + Sk[i].Zip(Xk[i]).Sum(sx => -sx.First * sx.Second);
                        var au2 = Fk[i + 1] + Sk[i + 1].Zip(Xk[i + 1]).Sum(sx => -sx.First * sx.Second);
                        var au = new double[] { -au1, -au2 };

                        double[,] linKoeff = new double[2, n + 1]; //коэффициенты 2ух лин. ограничений 
                        for (int ii = 0; ii < 2; ii++)
                        {
                            for (int j = 0; j < n; j++)
                            {
                                linKoeff[ii, j] = Sk[ii][j];
                            }
                            linKoeff[ii, n] = -1; //для t
                        }

                        //Mаксимальное из минимальных значений двух функций касательных
                        var minMV = MinModelValue(linKoeff, au, example.BoxUp, example.BoxLow);
                        //Console.WriteLine($"Мин. значение модели ^f{k}*: {minMV}\n ");

                        Ak[i + 1] = Math.Max(Ak[i], minMV);
                        Fk[i + 1] = example.GetValueAt(Xk[i + 1]).Value;
                        //Пересчитываем рекордное значение модели
                        Bk[i + 1] = Math.Min(Bk[i], Fk[i + 1]);
                        //Console.WriteLine($"Pекордное значение модели f{k}*: {Bk[i + 1]}\n ");

                        if (Ak[i + 1] < Lk[i])
                        {
                            //опускаем уровень
                            Lk[i + 1] = (1 - teta) * A + teta * Bk[i + 1];
                            //Console.WriteLine($"L{k}{i+1} = {Lk[i+1]}");
                            i++;
                            //проектируем точку с новым l
                        }
                        else break;
                    }


                    //5
                    //Если мин знач новой модели ф-ии >= уровню
                    A = Ak[i + 1];
                    B = Bk[i + 1];

                    //текущая точка минимума выбирается из построенных 
                    X[k + 1] = Xk[MinIndex(Fk, i + 1)];
                    F[k + 1] = example.GetValueAt(X[k + 1]).Value;
                    k++;
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                var x1 = new double[2*X.Length][];
                Array.Copy(X, 0, x1, 0, k+1);
                Array.Copy(Xk, 1, x1, k+1, i);

                var f1 = new double[2*X.Length];
                Array.Copy(F, 0, f1, 0, k+1);
                Array.Copy(Fk, 1, f1, k+1, i);

                return (x1, f1, k+i, counter);
            }
            return (X, F, k, counter);
        }



        public (double[][], double[], int,int) MinimizeWithoutCuttingPlanes()
        {
            var n = example.N;
            //var c = (int)-Math.Log(eps, 10);
            var steps = 15* example.N;
            var k = 0;
            var X = new double[steps][];
            var F = new double[steps];
            var A = new double[steps];
            var B = new double[steps];
            var S = new double[steps][];
            var L = new double[steps];
            X[k] = x0;
            A[k] = alfa;
            B[k] = example.GetValueAt(X[k]).Value; //начальное рекордное значение модели ф-ии
            F[k] = B[k];
            try
            {
                //Пока ошибка модели больше требуемой точности
                while (B[k] - A[k] > eps)
                {
                    //Console.WriteLine($"k = {k}:");

                    //Уровень - между мин. и рекордным значениями
                    L[k] = (1 - teta) * A[k] + teta * B[k];
                    //Console.WriteLine($"L{k} = {L[k]}");

                    //Касательная в последней построенной точке модели ф-ии
                    S[k] = example.GetSubGradAt(X[k]);
                    F[k] = example.GetValueAt(X[k]).Value;

                    double[,] linCon = new double[1, n];
                    double[] l = new double[1];
                    for (int i = k; i <= k; i++)
                    {
                        var sum = F[i] + S[i].Zip(X[k]).Sum(sx => -sx.First * sx.Second);
                        for (int ii = 0; ii < n; ii++)
                            linCon[i-k, ii] = S[i][ii];
                        l[i-k] = L[k] - sum;
                    }
 
                    //строим проекцию на множество уровней  - ищем точку, которая ближе всего к последней xk,
                    //но значение всех касательных к xi-м в ней не больше уровня
                    X[k + 1] = Projection(X[k], linCon, example.BoxUp, example.BoxLow, l);

                    //Console.WriteLine($"Проекция точки {alglib.ap.format(X[k], 3)} \nна множество U{k}: \n{alglib.ap.format(X[k + 1], 3)} ");
                    //Новая модель функции с большим на 1 кол-вом точек

                    //4. Проводим касательную в этой новой точке
                    S[k + 1] = example.GetSubGradAt(X[k+1]);
                    F[k + 1] = example.GetValueAt(X[k + 1]).Value;

                    //Пересчитаем мин значение модели функции
                    var au = new double[2];
                    double[,] linKoeff = new double[2, n + 1]; //коэффициенты лин. ограничений 

                    for (int ii = k; ii < k+2; ii++)
                    {
                        au[ii-k] = -F[ii] - S[ii].Zip(X[ii]).Sum(sx => -sx.First * sx.Second);
                        for (int j = 0; j < n; j++)
                        {
                            linKoeff[ii-k, j] = S[ii][j];
                        }
                        linKoeff[ii-k, n] = -1; //для t
                    }

                    //Mаксимальное из минимальных значений функций касательных
                    var minMV = MinModelValue(linKoeff, au, example.BoxUp, example.BoxLow);
                    //Console.WriteLine($"Мин. значение модели ^f{k}*: {minMV}\n ");

                    A[k + 1] = Math.Max(A[k], minMV);
                    F[k + 1] = example.GetValueAt(X[k + 1]).Value;
                    //Пересчитываем рекордное значение модели
                    B[k + 1] = Math.Min(B[k], F[k + 1]);
                    //Console.WriteLine($"Pекордное значение модели f{k}*: {B[k + 1]}\n ");
                    k++;
                }
        }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            return (X, F, k, k+1);
        }

        /// <summary>
        /// Индекс минимального рекордного значения функции
        /// </summary>
        /// <param name="fk">значения функции в построенных точках</param>
        /// <param name="end">индекс последнего рекорд. значения</param>
        /// <returns></returns>
        private int MinIndex(double[] fk, int end)
        {
            var min = double.PositiveInfinity;
            var index = 0;
            for(int i=0; i <= end; i++)
            {
                if(fk[i]<min)
                {
                    min = fk[i];
                    index = i;
                }
            }
            return index;
        }



        /// <summary>
        /// Проекция точки на множество, заданное лин. неравенством.
        /// Задача условной квадратичной минимизации
        /// ||x-x0||^2 -> min = Sum((xi-x0i)^2)-> min
        /// f(x0)+<s,x-x0> <= l
        /// xi>=0
        /// </summary>
        /// <param name="x0">проектируемая точка</param>
        /// <param name="l">уровень</param>
        /// <returns></returns>
        public double[] Projection(double[] x0, double[,] linCon,double[] boxUp, double[] boxLow, double[] l)
        {
            int n = example.N;

            double[,] A = new double[n, n]; //гессиан
            double[] B = new double[n]; // вектор коэфф-в b
            double[] aUp = l; //ограничение лин. фнкц. сверху
            double[] aLow = new double[l.Length]; //снизу
            double[] scale = new double[n]; //масштаб

            for (int i = 0; i < l.Length; i++)
                aLow[i] = -double.PositiveInfinity;

            for (int i = 0; i < n; i++)
            {
                A[i, i] = 2.0; //в цф нет произведений разноименных переменных
                B[i] = -2 * x0[i];
                scale[i] = example.Scale; //один масштаб у всех иксов=
            }

            double[] x;
            alglib.minqpstate state;
            alglib.minqpreport rep;

            alglib.minqpcreate(n, out state);

            alglib.minqpsetquadraticterm(state, A);
            alglib.minqpsetlinearterm(state, B);

            alglib.minqpsetlc2dense(state, linCon, aLow, aUp);
            alglib.minqpsetbc(state, boxLow, boxUp);

            alglib.minqpsetscale(state, scale);

            alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
            alglib.minqpoptimize(state);
            alglib.minqpresults(state, out x, out rep);
            return x;
        }


        /// <summary>
        /// Оптимальное значение модели (симплексом)
        /// t-> min
        /// f(x1)+<s1, x-x1> <=t
        /// f(x2)+<s2, x-x2> <=t
        /// t >= 0
        /// xi=[box]
        /// </summary>
        /// <param name="x1">Итерационная точка</param>
        /// <param name="x2">Ее проекция на U</param>
        /// <returns></returns>
        public double MinModelValue(double[,] A, double[] au, double[] boxU, double[] boxL)
        {
            var n = example.N;
            double[] al = new double[] { -Double.PositiveInfinity, -Double.PositiveInfinity };

            double[] c = new double[n + 1]; //цф: 0*x1+..0*xn + t -> min
            c[n] = 1;

            double[] scale = new double[n + 1];
            double[] boxLow = new double[n + 1];
            double[] boxUp = new double[n + 1];
            for (int i = 0; i < n + 1; i++)
                scale[i] = example.Scale; //один масштаб у всех иксов
            for (int i = 0; i < n; i++)
            {
                boxUp[i] = boxU[i]; //иксы сверху
                boxLow[i] = boxL[i]; //иксы снизу
            }
            boxUp[n] = Double.PositiveInfinity; //t
            boxLow[n] = -Double.PositiveInfinity; //t

            double[] x;
            alglib.minlpstate state;
            alglib.minlpreport rep;

            alglib.minlpcreate(n + 1, out state);

            alglib.minlpsetcost(state, c);
            alglib.minlpsetbc(state, boxLow, boxUp);
            alglib.minlpsetlc2dense(state, A, al, au, 2);

            alglib.minlpsetscale(state, scale);

            //revised dual simplex method
            alglib.minlpsetalgodss(state, 0.0);
            alglib.minlpoptimize(state);
            alglib.minlpresults(state, out x, out rep);

            return x[n];
        }


        /// <summary>
        /// Диаметр области определения
        /// </summary>
        /// <returns></returns>
        public double FindDiametr(double[] boxLow, double[] boxUp)
        {
            var n = 2*example.N;
            double[,] A = new double[n,n]; //гессиан
            double[] B = new double[n]; // вектор коэфф-в b
            double[] scale = new double[n]; //масштаб
            double[] boxU = new double[n];
            double[] boxL = new double[n];

            for (int i = 0; i < n; i++)
            {
                A[i, i] = -2.0;
                if (i % 2 == 0)
                {
                    A[i, i + 1] = 2.0;
                    A[i + 1, i] = 2.0;
                }
                scale[i] = 1; //один масштаб у всех иксов
                boxU[i] = i % 2 != 0 ? double.PositiveInfinity : boxUp[i / 2];
                boxL[i] = i % 2 == 0 ? -double.PositiveInfinity : boxLow[i / 2];

            }

            double[] x;
            alglib.minqpstate state;
            alglib.minqpreport rep;
            alglib.minqpcreate(n, out state);

            alglib.minqpsetquadraticterm(state, A);
            alglib.minqpsetlinearterm(state, B);
            alglib.minqpsetscale(state, scale);


            alglib.minqpsetbc(state, boxL, boxU);
            alglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
            alglib.minqpoptimize(state);
            alglib.minqpresults(state, out x, out rep);

            var zn = 0.0;
            for(int i=0; i<n-1; i++)
                zn += i % 2 == 0 ? (x[i] - x[i + 1]) * (x[i] - x[i + 1]) : 0;
           
           // Console.WriteLine($"Диаметр множества =  {zn} ");
            return Math.Sqrt(zn);
        }

        /// <summary>
        /// Величина, ограничивающая сверху количество шагов ПМУ
        /// </summary>
        /// <returns></returns>
        public long FindMaxStepsCount()
        {
            var d = FindDiametr(example.BoxLow, example.BoxUp);
            var L = example.FindL();
            var ch = d * d * L * L;
            var zn = eps * eps * (1 - teta) * (1 - teta) * (1 - (1 - teta) * (1 - teta));
           
            //var ch1 = 4 * L * L * d * d;
            //var zn1 = eps * eps;
            //Console.WriteLine($"Число шагов ПМУ не превышает {(long)Math.Truncate(ch1 / zn1) + 1}\n");

            return (long)Math.Truncate(ch / zn) + 1;
        }

    }
}

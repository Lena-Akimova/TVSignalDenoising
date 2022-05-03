using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TVSignalDenoising
{
    public class Regularisator
    {
        public double Lambda { get; set; }
        public double[] Signal { get; set; }


        public Regularisator(double lambda, double[] signal)
        {
            Lambda = lambda;
            Signal = new double[signal.Length];
            signal.CopyTo(Signal, 0);
        }



        /// <summary>
        /// Значение функции Sum(si-xi)^2 + lambda*Sum(|xi_1 - xi|) в точке 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public double GetValueAt(double[] x)
        {
            var sumKvadr = Signal.Zip(x).Sum(sx => Math.Pow(sx.First - sx.Second, 2));
            var sumMod = 0.0;
            for (int i = 0; i < x.Length - 1; i++)
                sumMod += Math.Abs(x[i + 1] - x[i]);

            return sumKvadr + Lambda * sumMod;
        }

        /// <summary>
        /// Субградиент в точке
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public double[] GetSubGradAt(double[] x)
        {
            //var s = signal.Zip(x).Select(sx => 2 * sx.Second - 2 * sx.First + lambda*2);

            var s = new double[x.Length];
            for (int i = 0; i < s.Length; i++)
                s[i] = i + 1;
            return s;
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
        public double[] Projection(double[] x0, double l)
        {
            int n = x0.Length;

            var s = GetSubGradAt(x0);
            var f = GetValueAt(x0);
            var sum = f + s.Zip(x0).Sum(sx => -sx.First * sx.Second);

            double[,] linCon = new double[1, n]; //одно основное лин. ограничение
            for (int i = 0; i < n; i++)
                linCon[0, i] = s[i];

            double[,] A = new double[n, n]; //гессиан
            double[] B = new double[n]; // вектор коэфф-в b
            double[] aUp = new double[] { l - sum }; //ограничение лин. фнкц. сверху
            double[] aLow = new double[] { -Double.PositiveInfinity }; //снизу
            double[] boxUp = new double[n]; //ограничения на значение икса
            double[] boxLow = new double[n];
            double[] scale = new double[n]; //масштаб

            for (int i = 0; i < n; i++)
            {
                A[i, i] = 2.0; //в цф нет произведений разноименных переменных
                B[i] = -2 * x0[i];
                scale[i] = 1; //один масштаб у всех иксов
                boxUp[i] = Double.PositiveInfinity; //иксы неотрицательны
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
            Console.WriteLine($"Проекция точки {alglib.ap.format(x0, 3)} \nна множество U: \n{alglib.ap.format(x, 3)} ");

            return x;

        }


        /// <summary>
        /// Оптимальное значение модели (симплексом)
        /// t-> min
        /// f(x1)+<s1, x-x1> <=t
        /// f(x2)+<s2, x-x2> <=t
        /// xi,t >= 0
        /// </summary>
        /// <param name="x1">Итерационная точка</param>
        /// <param name="x2">Ее проекция на U</param>
        /// <returns></returns>
        public double MinModelValue(double[] x1, double[] x2)
        {
            var n = x1.Length;
            var s1 = GetSubGradAt(x1);
            var s2 = GetSubGradAt(x2);
            var s = new double[2][];
            s[0] = s1;
            s[1] = s2;

            var f1 = GetValueAt(x1);
            var f2 = GetValueAt(x2);

            //Ограничения лин. функц. сверху
            var au1 = f1 + s1.Zip(x1).Sum(sx => -sx.First * sx.Second);
            var au2 = f2 + s2.Zip(x2).Sum(sx => -sx.First * sx.Second);


            double[,] A = new double[2, n + 1]; //коэффициенты лин. ограничений
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = s[i][j];
                }
                A[i, n] = -1; //для t
            }

            double[] al = new double[] { -Double.PositiveInfinity, -Double.PositiveInfinity };
            double[] au = new double[] { -au1, -au2 };

            double[] c = new double[n + 1]; //цф: 0*x1+..0*xn + t -> min
            c[n] = 1;

            double[] scale = new double[n + 1];
            double[] boxLow = new double[n + 1];
            double[] boxUp = new double[n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                scale[i] = 1; //один масштаб у всех иксов
                boxUp[i] = Double.PositiveInfinity; //иксы неотрицательны
            }

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
            Console.WriteLine($"Мин. значение модели ^fk*: {x[n]}\n {alglib.ap.format(x, 3)}");

            return x[n];
        }

    }
}

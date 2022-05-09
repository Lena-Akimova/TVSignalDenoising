using System;
using System.Collections.Generic;
using System.Text;

namespace TVSignalDenoising
{
    interface IOptimizator
    {
        double[] Minimize();
    }
}

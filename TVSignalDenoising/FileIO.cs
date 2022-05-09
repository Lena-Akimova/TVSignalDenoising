using System;
using System.Collections.Generic;
using System.IO;

namespace TVSignalDenoising
{
    public class FileIO
    {
        public static double[] ReadSignal(string path = "Resources/signal.txt")
        {
            if(string.IsNullOrEmpty(path))
            {
                var location = System.Reflection.Assembly.GetExecutingAssembly().Location;
                var locName = Path.GetDirectoryName(location);
                var p = locName.Remove(locName.IndexOf("bin"));
                path = $"{p}Resources\\signal.txt";
            }

            using (StreamReader reader = new StreamReader(path))
            {
                var res = new List<double>();
                string? line;
                while ((line = reader.ReadLine()) != null)
                {
                    double s = 0;
                    if(double.TryParse(line, out s))
                        res.Add(s);
                }
                return res.ToArray();
            }
        }
    }
}

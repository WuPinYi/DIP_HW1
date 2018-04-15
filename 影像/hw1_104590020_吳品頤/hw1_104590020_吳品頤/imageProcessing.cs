using System;

using System.Runtime.InteropServices;
using Emgu.CV;
using Emgu.CV.Structure;
using System.Drawing;
using System.Drawing.Imaging;
using System.Drawing.Drawing2D;

namespace hw1_104590020_吳品頤
{
    public static class imageProcessing
    {
        public static Image<Bgr, Byte> ConvertToGray(Image<Bgr, Byte> source)  //灰階 function //Gray = R*0.299 + G*0.587 + B*0.114
        {
            Image<Bgr, Byte> result = new Image<Bgr, Byte>(source.Width, source.Height);
            int rows = source.Height;
            int cols = source.Width;
            byte r, g, b;

            for (int y = 0; y < rows; y++)
            {
                for (int x = 0; x < cols; x++)
                {
                    b = (byte)source.Data[y, x, 0];  //B
                    g = (byte)source.Data[y, x, 1];  //G
                    r = (byte)source.Data[y, x, 2];  //R
                    byte grayColor = (byte)(b * 0.114 + g * 0.587 + r * 0.299);
                    result.Data[y, x, 0] = grayColor;
                    result.Data[y, x, 1] = grayColor;
                    result.Data[y, x, 2] = grayColor;
                }
            }

            return result;
        }

        public static Image<Bgr, Byte> ConvertToMirror(Image<Bgr, Byte> source)
        {
            Image<Bgr, Byte> result = new Image<Bgr, Byte>(source.Width, source.Height);
            int rows = source.Height;
            int cols = source.Width;
            byte r, g, b;

            for (int y = 0; y < rows; y++)
            {
                for (int x = 0; x < cols; x++)
                {
                    b = (byte)source.Data[y, x, 0];  //B
                    g = (byte)source.Data[y, x, 1];  //G
                    r = (byte)source.Data[y, x, 2];  //R
                    result.Data[y, (cols - x - 1), 0] = b;
                    result.Data[y, (cols - x - 1), 1] = g;
                    result.Data[y, (cols - x - 1), 2] = r;
                }
            }

            return result;
        }

        public static Image<Bgr, Byte> Rotating(Image<Bgr, Byte> source, double theta)
        {
            Image<Bgr, Byte> result = new Image<Bgr, Byte>(source.Width, source.Height);
            double r = Math.Sqrt(Math.Pow((double)source.Width / 2d, 2d) + Math.Pow((double)source.Height / 2d, 2d));
            double OriginalAngle = Math.Acos((source.Width / 2d) / r) / Math.PI * 180d;
            double minW = 0d, maxW = 0d, minH = 0d, maxH = 0d;
            double[] drawPoint = new double[4];
            drawPoint[0] = (-OriginalAngle + theta) * Math.PI / 180d;
            drawPoint[1] = (OriginalAngle + theta) * Math.PI / 180d;
            drawPoint[2] = (180f - OriginalAngle + theta) * Math.PI / 180d;
            drawPoint[3] = (180f + OriginalAngle + theta) * Math.PI / 180d;

            foreach (double point in drawPoint)
            {
                double x = r * Math.Cos(point);
                double y = r * Math.Sin(point);

                if (x < minW)
                    minW = x;

                if (x > maxW)
                    maxW = x;

                if (y < minH)
                    minH = y;

                if (y > maxH)
                    maxH = y;
            }

            Size newSize = new Size((int)(maxW - minW), (int)(maxH - minH));
            ///Bitmap result = new Bitmap(newSize.Width, newSize.Height);
            PointF centerPoint = new PointF((float)result.Width / 2f, (float)result.Height / 2f);
            Graphics g = Graphics.FromImage(result.Bitmap);
            g.CompositingQuality = CompositingQuality.HighQuality;
            g.SmoothingMode = SmoothingMode.HighQuality;
            g.InterpolationMode = InterpolationMode.HighQualityBicubic;
            g.TranslateTransform(centerPoint.X, centerPoint.Y);
            g.RotateTransform((float)theta);
            g.TranslateTransform(-centerPoint.X, -centerPoint.Y);
            g.DrawImage(source.Bitmap, (float)(result.Width - source.Width) / 2f, (float)(result.Height - source.Height) / 2f, source.Width, source.Height);
            g.Dispose();
            return result;
        }

        public static Image<Bgr, Byte> imageBlending(Image<Bgr, Byte> source, Image<Bgr, Byte> source2, double Threshold)
        {
            Image<Bgr, Byte> result = new Image<Bgr, Byte>(source.Width, source.Height);
            int rows = source.Height;
            int cols = source.Width;
            int rows2 = source2.Height;
            int cols2 = source2.Width;
            int r, g, b;
            int r2, g2, b2;

            for (int y = 0, i = 0; y < rows; y++, i++)
            {
                for (int x = 0, j = 0; x < cols; x++, j++)
                {
                    b = (int)source.Data[y, x, 0];  //B
                    g = (int)source.Data[y, x, 1];  //G
                    r = (int)source.Data[y, x, 2];  //R
                    b2 = (int)source2.Data[i, j, 0];  //B2
                    g2 = (int)source2.Data[i, j, 1];  //G2
                    r2 = (int)source2.Data[i, j, 2];  //R2
                    result.Data[y, x, 0] = (byte)(b * Threshold + b2 * (1 - Threshold));
                    result.Data[y, x, 1] = (byte)(g * Threshold + g2 * (1 - Threshold));
                    result.Data[y, x, 2] = (byte)(r * Threshold + r2 * (1 - Threshold));
                }
            }

            return result;
        }

        public static Image<Bgr, Byte> ConvertToOtsu(Image<Bgr, Byte> source)
        {
            Image<Bgr, Byte> result = new Image<Bgr, Byte>(source.Width, source.Height);
            int rows = source.Height;
            int cols = source.Width;
            int PixNum = rows * cols;
            int GrayNum = 256;
            int[] PixArr = new int[PixNum];
            int[] GrayArr = new int[256];
            int r, g, b;

            for (int y = 0; y < rows; y++)
            {
                for (int x = 0; x < cols; x++)
                {
                    b = (int)source.Data[y, x, 0];  //B
                    g = (int)source.Data[y, x, 1];  //G
                    r = (int)source.Data[y, x, 2];  //R
                    int grayColor = (int)(b * 0.114 + g * 0.587 + r * 0.299);
                    PixArr[y * cols + x] = grayColor;
                    GrayArr[PixArr[y * cols + x]]++;
                }
            }

            double u0, u1;
            double maxVal = 0;
            int endval = 0;
            double w0, w1;

            for (int i = 0; i < GrayNum; i++)
            {
                w1 = w0 = u0 = u1 = 0;

                for (int j = 0; j < GrayNum; j++)
                {
                    if (j <= i)
                    {
                        w0 += GrayArr[j] / (PixNum * 1.0);
                        u0 += GrayArr[j] / (PixNum * 1.0) * j;
                    }
                    else
                    {
                        w1 += GrayArr[j] / (PixNum * 1.0);
                        u1 += GrayArr[j] / (PixNum * 1.0) * j;
                    }
                }

                u0 = u0 / w0;
                u1 = u1 / w1;
                double vala = w1 * w0 * Math.Pow(u0 - u1, 2);

                if (maxVal < vala)
                {
                    maxVal = vala;
                    endval = i;
                }
            }

            int val = endval;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    b = (int)source.Data[i, j, 0];  //B
                    g = (int)source.Data[i, j, 1];  //G
                    r = (int)source.Data[i, j, 2];  //R
                    int grayColor = (int)(b * 0.114 + g * 0.587 + r * 0.299);

                    if (grayColor > val)
                    {
                        result.Data[i, j, 0] = 255;
                        result.Data[i, j, 1] = 255;
                        result.Data[i, j, 2] = 255;
                    }
                    else
                    {
                        result.Data[i, j, 0] = 0;
                        result.Data[i, j, 1] = 0;
                        result.Data[i, j, 2] = 0;
                    }
                }
            }

            return result;
        }
        //

        public static Image<Bgr, Byte> HistogramEqualization(Image<Bgr, Byte> source)
        {
            Image<Bgr, Byte> result = new Image<Bgr, Byte>(source.Width, source.Height);
            int rows = source.Height;
            int cols = source.Width;
            int[,] pdf = new int[3, 256];
            int[,] cdf = new int[3, 256];

            // 算pdf
            for (int i = 0; i < 256; i++)   // i: RGB的0~255
            {
                for (int y = 0; y < rows; y++)
                {
                    for (int x = 0; x < cols; x++)
                    {
                        for (int j = 0; j < 3; j++) // j: 三個通道
                        {
                            if (source.Data[y, x, j] == i) pdf[j, i]++;
                        }
                    }
                }
            }

            // 算cdf
            for (int j = 0; j < 3; j++)
            {
                int temp = 0;

                for (int i = 0; i < 256; i++)
                {
                    temp += pdf[j, i];
                    cdf[j, i] = temp;
                }
            }

            // 結果
            for (int y = 0; y < rows; y++)
            {
                for (int x = 0; x < cols; x++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        byte original = source.Data[y, x, j];
                        int min = 0;

                        // 求CDF Min
                        for (int i = 0; i < 255; i++)
                            if (cdf[j, i] > min) min = cdf[j, i];

                        for (int i = 0; i < 255; i++)
                        {
                            if (pdf[j, i] != 0 && cdf[j, i] < min)
                                min = cdf[j, i];
                        }

                        result.Data[y, x, j] = H(cdf[j, original], min, rows, cols);
                    }
                }
            }

            return result;
        }
        private static byte H(int cdf, int cdfMin, int m, int n)  // m: 長, n: 寬
        {
            double num = (cdf - cdfMin) * 255 / ((m * n) - cdfMin) ;
            return (byte)Math.Round(num);
        }
    }
}

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using Emgu.CV;
using Emgu.CV.Util;
using Emgu.CV.Structure;
//http://www.emgu.com/wiki/index.php/Setting_up_EMGU_C_Sharp

namespace hw1_104590020_吳品頤
{
    public partial class Form1 : Form
    {
        Image<Bgr, Byte> _sourceImage = null;
        Image<Bgr, Byte> _sourceImage2 = null;
        Image<Bgr, Byte> _resultImage;
        Button nowClick = null;

        public Form1()
        {
            InitializeComponent();
        }

        private string LoadImageFile()
        {
            string fileName = "";
            OpenFileDialog dialog = new OpenFileDialog();
            DirectoryInfo dir = new DirectoryInfo(System.Windows.Forms.Application.StartupPath);
            dialog.Title = "Open an Image File";
            dialog.RestoreDirectory = true;
            dialog.InitialDirectory = dir.Parent.Parent.FullName;
            dialog.Filter = "PNG (*.png)|*.png|JPEG (*.jpeg;*.jpg)|*.jpg;*.jpeg|Bmp (*.bmp)|*.bmp";
            if (dialog.ShowDialog() == System.Windows.Forms.DialogResult.OK && dialog.FileName != null)
            {
                fileName = dialog.FileName;
            }
            return fileName;
        }

        private void _loadSourceImageButton_Click(object sender, EventArgs e)
        {
            string fileName = LoadImageFile();
            if (fileName != "")
            {
                _sourceImage = new Image<Bgr, Byte>(fileName);
                _sourcePictureBox.Image = _sourceImage.Bitmap;
            }
        }

        private void saveResultImageButton_Click(object sender, EventArgs e)
        {
            if (_resultImage == null) return;
            SaveFileDialog dialog = new SaveFileDialog();
            dialog.Filter = "PNG (*.png)|*.png|JPEG (*.jpeg;*.jpg)|*.jpg;*.jpeg|Bmp (*.bmp)|*.bmp";
            dialog.Title = "Save an Image File";
            if (dialog.ShowDialog() == System.Windows.Forms.DialogResult.OK && dialog.FileName != null)
            {
                string fileName = dialog.FileName;
                _resultImage.Save(fileName);
            }
        }

        /// //////////////////////////// gray ///////////////////////////////////////////
        private void _grayButton_Click(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;
            ifShow_Threshold_trackBar_Scroll(false);
            _resultImage = imageProcessing.ConvertToGray(_sourceImage);
            _resultPictureBox.Image = _resultImage.Bitmap;

            nowClick = _grayButton;
        }

        /// //////////////////////////// mirror ///////////////////////////////////////////
        private void _mirrorButton_Click(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;
            ifShow_Threshold_trackBar_Scroll(false);
            _resultImage = imageProcessing.ConvertToMirror(_sourceImage);
            _resultPictureBox.Image = _resultImage.Bitmap;

            nowClick = _mirrorButton;
        }

        /// //////////////////////////// Rotating ///////////////////////////////////////////
        private void _RotatingButton_Click(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;

            ifShow_Threshold_trackBar_Scroll(true);
            if (nowClick != _RotatingButton)
            {
                nowClick = _RotatingButton;
                setThreshold_trackBar(0, 36000, 7595);
            }
            _Threshold_Label.Text = (_Threshold_trackBar.Value / 100.0).ToString();

            _resultImage = imageProcessing.Rotating(_sourceImage, _Threshold_trackBar.Value / 100.0);
            _resultPictureBox.Image = _resultImage.Bitmap;
        }

        /// //////////////////////////// imageBlending ///////////////////////////////////////////
        private void _blendingButton_Click(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;

            ifShow_Threshold_trackBar_Scroll(true);
            if (nowClick != _blendingButton)
            {
                string fileName2 = LoadImageFile();
                if (fileName2 != "")
                {
                    _sourceImage2 = new Image<Bgr, Byte>(fileName2);
                }
                if (_sourceImage2 == null) return;

                nowClick = _blendingButton;
                setThreshold_trackBar(0, 100, 50);
            }
            _Threshold_Label.Text = (_Threshold_trackBar.Value / 100.0).ToString();

            _resultImage = imageProcessing.imageBlending(_sourceImage, _sourceImage2, _Threshold_trackBar.Value / 100.0);
            _resultPictureBox.Image = _resultImage.Bitmap;
        }

        /// //////////////////////////// ConvertToOtsu ///////////////////////////////////////////
        private void _OtsuButton_Click(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;
            ifShow_Threshold_trackBar_Scroll(false);
            _resultImage = imageProcessing.ConvertToOtsu(_sourceImage);
            _resultPictureBox.Image = _resultImage.Bitmap;

            nowClick = _OtsuButton;
        }

        /// //////////////////////////// HistogramEqualization ///////////////////////////////////////////

        private void _HistogramEqualizationButton_Click(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;
            ifShow_Threshold_trackBar_Scroll(false);
            _resultImage = imageProcessing.HistogramEqualization(_sourceImage);
            _resultPictureBox.Image = _resultImage.Bitmap;

            nowClick = _HistogramEqualizationButton;
        }

        /// //////////////////////////// _Threshold_trackBar_Scroll ///////////////////////////////////////////
        private void _Threshold_trackBar_Scroll(object sender, EventArgs e)
        {
            if (_sourceImage == null) return;
            nowClick.PerformClick();
        }
        void setThreshold_trackBar(int min, int max, int initialValue)
        {
            _Threshold_trackBar.Minimum = min;
            _Threshold_trackBar.Maximum = max;
            _Threshold_trackBar.Value = initialValue;
            _Threshold_min_Label.Text = (_Threshold_trackBar.Minimum / 100.0).ToString();
            _Threshold_max_Label.Text = (_Threshold_trackBar.Maximum / 100.0).ToString();
            _Threshold_Label.Text = (_Threshold_trackBar.Value / 100.0).ToString();
        }
        void ifShow_Threshold_trackBar_Scroll(bool s)
        {
            if (s)
            {
                _Threshold_trackBar.Show();
                _Threshold_min_Label.Show();
                _Threshold_max_Label.Show();
                _Threshold_Label.Show();
            }
            else
            {
                _Threshold_trackBar.Hide();
                _Threshold_min_Label.Hide();
                _Threshold_max_Label.Hide();
                _Threshold_Label.Hide();
            }
        }

    }
}


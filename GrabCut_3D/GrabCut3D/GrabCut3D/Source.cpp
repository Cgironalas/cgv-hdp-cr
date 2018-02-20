#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <vector>

using namespace cv;
using namespace std;

static bool go, comeBack1, comeBack2;

static Mat imageToShow, tempMat, tempMask;

static string sourceDir,	//Directory of the source images
maskDir,		//Directory of pre-existing masks or where they will be saved.
fileExtension,
tempDir;		//Temporal string for various usages.

static string *imgFiles,	// names of all the images from sourceDir
*maskFiles;	// names of all the masks from maskDir

static int index, currentImage, folderSize, indexBegin, middle, indexEnd, firstMask, lastMask, counter, rows, columns;

static unsigned char *matPixel;
static unsigned char *maskPixel;


static void help() {
	cout << "\nSelect a rectangular area around the object you want to segment\n" <<
		"\nHot keys: \n"
		"\tESC - quit the program\n"
		"\tr - restore the original image\n"
		"\tn - next iteration\n"
		"\n"
		"\tleft mouse button - set rectangle\n"
		"\n"
		"\tCTRL+left mouse button - set GC_BGD pixels\n"
		"\tSHIFT+left mouse button - set GC_FGD pixels\n"
		"\n"
		"\tCTRL+right mouse button - set GC_PR_BGD pixels\n"
		"\tSHIFT+right mouse button - set GC_PR_FGD pixels\n"
		"\n"
		"\tj - previous image\n"
		"\tk - next image\n"
		"\n"
		"\th - first image with mask\n"
		"\tl - last image with mask\n"
		"\n" << endl;
}

const Scalar RED = Scalar(0, 0, 255);
const Scalar PINK = Scalar(230, 130, 255);
const Scalar BLUE = Scalar(255, 0, 0);
const Scalar LIGHTBLUE = Scalar(255, 255, 160);
const Scalar GREEN = Scalar(0, 255, 0);

const int BGD_KEY = EVENT_FLAG_CTRLKEY;
const int FGD_KEY = EVENT_FLAG_SHIFTKEY;

static void getBinMask(const Mat& comMask, Mat& binMask) {
	if (comMask.empty() || comMask.type() != CV_8UC1)
		CV_Error(Error::StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)");
	if (binMask.empty() || binMask.rows != comMask.rows || binMask.cols != comMask.cols)
		binMask.create(comMask.size(), CV_8UC1);
	binMask = comMask & 1;
}

class GCApplication {
public:
	enum {
		NOT_SET = 0,
		IN_PROCESS = 1,
		SET = 2
	};
	static const int radius = 2;
	static const int thickness = -1;

	bool isInitialized;
	Mat mask;

	void reset();
	void setImageAndWinName(const Mat& _image, const string& _winName);
	void showImage() const;
	void mouseClick(int event, int x, int y, int flags, void* param);
	int nextIter();
	int getIterCount() const { return iterCount; }
	Mat getMask() const { return mask; }

private:
	void setRectInMask();
	void setLblsInMask(int flags, Point p, bool isPr);

	const string* winName;
	const Mat* image;
	Mat bgdModel, fgdModel;
	uchar rectState, lblsState, prLblsState;

	Rect rect;
	vector<Point> fgdPxls, bgdPxls, prFgdPxls, prBgdPxls;
	int iterCount;
};

//Set all the elements of the class as null or equivalent (mat with 0s, NOT_SET, 0, etc.)
void GCApplication::reset() {
	if (!mask.empty())
		mask.setTo(Scalar::all(GC_BGD));
	bgdPxls.clear();
	fgdPxls.clear();
	prBgdPxls.clear();
	prFgdPxls.clear();

	isInitialized = false;
	rectState = NOT_SET;
	lblsState = NOT_SET;
	prLblsState = NOT_SET;
	iterCount = 0;
}

//Kind of like a constructor, just sets the image and winName elements of the class
void GCApplication::setImageAndWinName(const Mat& _image, const string& _winName) {
	if (_image.empty() || _winName.empty())
		return;
	image = &_image;
	winName = &_winName;
	mask.create(image->size(), CV_8UC1);
	reset();
}

//
void GCApplication::showImage() const {
	if (image->empty() || winName->empty())
		return;

	Mat res;
	Mat binMask, colorMask;
	int r, g, b;
	Scalar intensityBin;
	Vec3b intensityRGB;
	if (!isInitialized)
		image->copyTo(res);
	else {
		image->copyTo(res);
		getBinMask(mask, binMask);
		
		cout << "Image size " << image->cols << " - " << image->rows;
		for (int x = 0; x < image->cols; x++) {
			for (int y = 0; y < image->rows; y++) {
				intensityBin = binMask.at<uchar>(y, x);
				intensityRGB = image->at<Vec3b>(y, x);
				b = intensityRGB.val[0];
				g = intensityRGB.val[1];
				r = intensityRGB.val[2];

				if (intensityBin.val[0] == 0) {
					r = (int)std::max(0.0, r - 50.0);
					g = (int)std::max(0.0, g - 50.0);
					b = (int)std::max(0.0, b - 50.0);
					res.at<Vec3b>(y, x) = Vec3b(r, g, b);
				}
			}
		}
	}

	vector<Point>::const_iterator it;
	for (it = bgdPxls.begin(); it != bgdPxls.end(); ++it)
		circle(res, *it, radius, BLUE, thickness);
	for (it = fgdPxls.begin(); it != fgdPxls.end(); ++it)
		circle(res, *it, radius, RED, thickness);
	for (it = prBgdPxls.begin(); it != prBgdPxls.end(); ++it)
		circle(res, *it, radius, LIGHTBLUE, thickness);
	for (it = prFgdPxls.begin(); it != prFgdPxls.end(); ++it)
		circle(res, *it, radius, PINK, thickness);

	if (rectState == IN_PROCESS || rectState == SET)
		rectangle(res, Point(rect.x, rect.y), Point(rect.x + rect.width, rect.y + rect.height), GREEN, 2);

	imshow(*winName, res);
}

//
void GCApplication::setRectInMask() {
	CV_Assert(!mask.empty());
	mask.setTo(GC_BGD);
	rect.x = max(0, rect.x);
	rect.y = max(0, rect.y);
	rect.width = min(rect.width, image->cols - rect.x);
	rect.height = min(rect.height, image->rows - rect.y);
	(mask(rect)).setTo(Scalar(GC_PR_FGD));
}

//
void GCApplication::setLblsInMask(int flags, Point p, bool isPr) {
	vector<Point> *bpxls, *fpxls;
	uchar bvalue, fvalue;
	if (!isPr) {
		bpxls = &bgdPxls;
		fpxls = &fgdPxls;
		bvalue = GC_BGD;
		fvalue = GC_FGD;
	}
	else {
		bpxls = &prBgdPxls;
		fpxls = &prFgdPxls;
		bvalue = GC_PR_BGD;
		fvalue = GC_PR_FGD;
	}
	if (flags & BGD_KEY) {
		bpxls->push_back(p);
		circle(mask, p, radius, bvalue, thickness);
	}
	if (flags & FGD_KEY) {
		fpxls->push_back(p);
		circle(mask, p, radius, fvalue, thickness);
	}
}

//Listener for mouse events
void GCApplication::mouseClick(int event, int x, int y, int flags, void*) {
	// TODO add bad args check
	switch (event) {
	case EVENT_LBUTTONDOWN: {/*set rect or GC_BGD(GC_FGD) labels*/
		bool isb = (flags & BGD_KEY) != 0,
			isf = (flags & FGD_KEY) != 0;
		if (rectState == NOT_SET && !isb && !isf) {
			rectState = IN_PROCESS;
			rect = Rect(x, y, 1, 1);
		}
		if ((isb || isf) && rectState == SET)
			lblsState = IN_PROCESS;
		break;
	}

	case EVENT_RBUTTONDOWN: {/*set GC_PR_BGD(GC_PR_FGD) labels*/
		bool isb = (flags & BGD_KEY) != 0,
			isf = (flags & FGD_KEY) != 0;
		if ((isb || isf) && rectState == SET)
			prLblsState = IN_PROCESS;
		break;
	}

	case EVENT_LBUTTONUP:
		if (rectState == IN_PROCESS) {
			rect = Rect(Point(rect.x, rect.y), Point(x, y));
			rectState = SET;
			setRectInMask();
			CV_Assert(bgdPxls.empty() && fgdPxls.empty() && prBgdPxls.empty() && prFgdPxls.empty());
			showImage();
		}
		if (lblsState == IN_PROCESS) {
			setLblsInMask(flags, Point(x, y), false);
			lblsState = SET;
			showImage();
		}
		break;
	case EVENT_RBUTTONUP:
		if (prLblsState == IN_PROCESS) {
			setLblsInMask(flags, Point(x, y), true);
			prLblsState = SET;
			showImage();
		}
		break;

	case EVENT_MOUSEMOVE:
		if (rectState == IN_PROCESS) {
			rect = Rect(Point(rect.x, rect.y), Point(x, y));
			CV_Assert(bgdPxls.empty() && fgdPxls.empty() && prBgdPxls.empty() && prFgdPxls.empty());
			showImage();
		}
		else if (lblsState == IN_PROCESS) {
			setLblsInMask(flags, Point(x, y), false);
			showImage();
		}
		else if (prLblsState == IN_PROCESS) {
			setLblsInMask(flags, Point(x, y), true);
			showImage();
		}
		break;
	}
}

//
int GCApplication::nextIter() {
	if (isInitialized)
		grabCut(*image, mask, rect, bgdModel, fgdModel, 1);
	else {
		if (rectState != SET)
			return iterCount;

		if (lblsState == SET || prLblsState == SET)
			grabCut(*image, mask, rect, bgdModel, fgdModel, 1, GC_INIT_WITH_MASK);
		else
			grabCut(*image, mask, rect, bgdModel, fgdModel, 1, GC_INIT_WITH_RECT);

		isInitialized = true;
	}
	iterCount++;

	bgdPxls.clear();
	fgdPxls.clear();
	prBgdPxls.clear();
	prFgdPxls.clear();
	return iterCount;
}

GCApplication gcapp;



//Passthrough function for the GCApplication mouseClick
static void on_mouse(int event, int x, int y, int flags, void* param) {
	gcapp.mouseClick(event, x, y, flags, param);
}



void readConfigurationFile(string configFile) {
	ifstream file;
	file.open(configFile);

	// Verify that the configuration file could be opened
	if (!file) {
		cerr << "Unable to read configuration file: " + configFile << endl;
		exit(-7);
	}

	int i = 0;
	string curLine;
	while (getline(file, curLine) && i < 5) {
		// Read source directory input
		if (i == 0) {
			sourceDir = curLine.substr(curLine.find(" = ") + 3);
			if (sourceDir.back() != '/') {
				sourceDir.push_back('/');
			}
		}
		// Read directory for binary masks
		else if (i == 1) {
			maskDir = curLine.substr(curLine.find(" = ") + 3);
			if (maskDir.back() != '/') {
				maskDir.push_back('/');
			}
		}
		// Read file extensions of the source files
		else if (i == 2) {
			fileExtension = curLine.substr(curLine.find(" = ") + 3);
		}
		// Read the beginning index of the source files
		else if (i == 3) {
			indexBegin = stoi(curLine.substr(curLine.find(" = ") + 3));
		}
		// Read the finishing index of the source files
		else if (i == 4) {
			indexEnd = stoi(curLine.substr(curLine.find(" = ") + 3));
		}

		i++;
	}

	// Verify that the configuration file could read all parameters expected 
	if (i < 5) {
		cerr << "The configuration file doesn't have enough parameters!" << endl;
		exit(-8);
	}
}

void printConfigurationValues() {
	cout << endl;
	cout << "Directory:\t" + sourceDir << endl;
	cout << "Mask dir:\t" + maskDir << endl;

	cout << "File Extension:\t" + fileExtension << endl;
	cout << "Index Begin:\t" + to_string(indexBegin) << endl;
	cout << "Index End:\t" + to_string(indexEnd) << endl;
}

//Read the image based on the current index
Mat readCurrentImage() {
	Mat image = imread(sourceDir + imgFiles[index], IMREAD_COLOR);
	if (image.empty()) {
		cout << endl << "Couldn't read image from file: " << imgFiles[index] << endl;
		return Mat();
	}

	return image;
}

// Attempts to read a mask for the image based on the current index
// If read correctly the mask will be transformed to GrabCut mask and used for the current image
// If mask is empty then it is reported and the variable go will be used
void readCurrentMask() {
	cout << maskDir + maskFiles[index];
	tempMat = imread(maskDir + maskFiles[index], IMREAD_GRAYSCALE);
	if (!tempMat.empty()) {
		cout << " - mask was read" << endl;
		// Threshold converts all the 255s in the binary mask to 3s, which GrabCut reads as foreground
		cv::threshold(tempMat, gcapp.mask, 250, 3, cv::THRESH_BINARY);
		gcapp.isInitialized = true;
		go = false;
	}
	else {
		cout << " - Mask was empty" << endl;
		go = true;
	}
	tempMat = Mat();
}

// Save the binary mask for the image based on the current index
string saveCurrentMask() {
	//Copy the grabCut mask
	Mat tempMask = gcapp.getMask().clone();
	//Create a dummy of same size and type as the GrabCut mask
	Mat binMask(tempMask.size(), tempMask.type());

	string maskFile = "rgb" + to_string(currentImage) + ".bmp";

	cv::threshold(tempMask, binMask, 2, 255, cv::THRESH_BINARY);

	imwrite(maskDir + maskFile, binMask);

	cout << "Saved mask for image: " + to_string(currentImage) << endl << endl;


	return maskFile;
}


// Read current image and mask, call showImage from gcapp
void showCurrentImage() {
	cout << endl << "Current Image: " + to_string(currentImage) << endl;
	cout << "Index: " + to_string(index) << endl;

	imageToShow = readCurrentImage();
	readCurrentMask();
	gcapp.showImage();
}

// Sets the current image to the one in the center of the directory and
// the index to the previously calculated "middle" value
void loadMiddleImage() {
	currentImage = indexBegin + middle;
	index = middle;
}

// Load and show first (middle) image
void resetImage() {
	if (index != middle) {
		loadMiddleImage();
		showCurrentImage();
	}
}

void loadNextImage() {
	if (currentImage == indexEnd) {
		cout << "Finished right side" << endl;
		loadMiddleImage();
	}
	else {
		currentImage++;
		index++;
	}
}

void loadPreviousImage() {
	if (currentImage == indexBegin) {
		cout << "Finished left side" << endl;
		loadMiddleImage();
	}
	else {
		currentImage--;
		index--;
	}
}


// Check if the mask of a given index exists
bool maskExists(int i) {
	//cout << "Mask to check: " + maskDir + maskFiles[i] << endl;
	Mat temp = imread(maskDir + maskFiles[i], IMREAD_GRAYSCALE);
	if (!temp.empty()) {
		go = false;
		return true;
	}
	else {
		go = true;
		return false;
	}
}

std::string hexStr(unsigned char *data, int len)
{
	std::stringstream ss;
	ss << std::hex;
	for (int i = 0; i < len; ++i)
		ss << std::setw(2) << std::setfill('0') << (int)data[i];
	return ss.str();
}


void generateVox() {
	cout << "Vox is being generated." << endl;
	fstream dataVox = fstream("data.bin", std::ios::out | std::ios::binary);
	fstream maskVox = fstream("mask.bin", std::ios::out | std::ios::binary);
	Mat tempMat;
	Mat tempMask;
	int current, maskCounter, dataCounter;
	maskCounter = 0;
	dataCounter = 0;
	for (current = folderSize - 1; current >= 0; current--) {
		dataCounter++;
		cout << "Index: " + to_string(current) << endl;
		tempMask = imread(maskDir + maskFiles[current], IMREAD_GRAYSCALE);
		//cout << "Read mask" << endl;
		tempMat = imread(sourceDir + imgFiles[current], IMREAD_COLOR);
		//cout << "Read source" << endl;
		if (!tempMask.empty()) {
			maskCounter++;
			for (int i = 0; i < tempMat.rows; i++) {
				for (int j = 0; j < tempMat.cols; j++) {
					matPixel = tempMat.ptr(i, j);
					maskPixel = tempMask.ptr(i, j);

					if (dataVox.is_open()) {

						ostringstream r, g, b;
						b << hex << std::setw(2) << std::setfill('0') << (int)maskPixel[0];
						g << hex << std::setw(2) << std::setfill('0') << (int)maskPixel[1];
						r << hex << std::setw(2) << std::setfill('0') << (int)maskPixel[2];

						if (b.str() == "00" && g.str() == "00" && r.str() == "00") {
							dataVox.write((char *)&maskPixel[2], 1);
							dataVox.write((char *)&maskPixel[1], 1);
							dataVox.write((char *)&maskPixel[0], 1);
						}
						else {
							dataVox.write((char *)&matPixel[2], 1);
							dataVox.write((char *)&matPixel[1], 1);
							dataVox.write((char *)&matPixel[0], 1);
						}
					}

					if (maskVox.is_open()) {
						maskVox.write((char *)&maskPixel[2], 1);
						maskVox.write((char *)&maskPixel[1], 1);
						maskVox.write((char *)&maskPixel[0], 1);
					}
				}
			}
			//cout << "Columns: " + to_string(tempMat.cols) + " - Rows: " + to_string(tempMat.rows) + "." << endl;
		}
		else {
			for (int i = 0; i < tempMat.rows; i++) {
				for (int j = 0; j < tempMat.cols; j++) {
					matPixel = tempMat.ptr(i, j);

					if (dataVox.is_open()) {
						dataVox.write((char *)&matPixel[2], 1);
						dataVox.write((char *)&matPixel[1], 1);
						dataVox.write((char *)&matPixel[0], 1);
					}
				}
			}

		}
	}
	dataVox.close();
	maskVox.close();
	cout << "Finished generating vox!" << endl;
	cout << "masks read: " << to_string(maskCounter) << " - Images read: " << to_string(dataCounter) << endl;
}


int main(int argc, char** argv) {

	//Input Error, no cofiguration file
	if (argc <= 1) {
		cout << "No configuration file in input!!!" << endl;
		exit(-1);
	}

	readConfigurationFile(argv[1]);
	printConfigurationValues();

	//Print program usage instructions
	help();


	folderSize = indexEnd - indexBegin + 1;
	cout << "Folder size: " + to_string(folderSize) << endl;


	imgFiles = new string[folderSize];
	maskFiles = new string[folderSize];

	for (index = 0, currentImage = indexBegin; index < folderSize; index++, currentImage++) {

		imgFiles[index] = "rgb" + to_string(currentImage) + "." + fileExtension;

		maskFiles[index] = "rgb" + to_string(currentImage) + ".bmp";
	}

	middle = (indexEnd - indexBegin) / 2;
	cout << "Middle: " + to_string(middle) << endl << endl << endl;

	loadMiddleImage();

	cout << "Current Image: " + to_string(currentImage) << endl;
	cout << "Index: " + to_string(index) + "\n" << endl;


	const string winName = "GrabCut";
	const string imageWin = "source";

	namedWindow(winName, WINDOW_AUTOSIZE);
	//Creates a window with name 'image', possible flags:
	//WINDOW_AUTOSIZE (get size of image, cannot be resized)
	//WINDOW_OPENGL (opengl compatible window)

	//namedWindow(imageWin, WINDOW_AUTOSIZE);

	setMouseCallback(winName, on_mouse, 0);
	//on_mouse is an in between function that they (in openCV in general) to manage mouse interaction
	//They just send the values to their own defined function gcapp.mouseClick
	//setMouseCallback(imageWin, on_mouse, 0);

	imageToShow = readCurrentImage();

	if (imageToShow.empty()) {
		cout << "First image empty!!!" << endl;
		return -2;
	}

	gcapp.setImageAndWinName(imageToShow, winName);

	readCurrentMask();
	gcapp.showImage();

	for (;;) {
		char c = (char)waitKey(0);
		switch (c) {
			// Close GrabCut Window (exits program)
		case '\x1b': {
			cout << "Exiting ..." << endl;
			goto exit_main;
		}

					// Run GrabCut through every image on both sides automatically
		case 'c': {
			cout << "Pressed c" << endl;
			comeBack1 = true;

			resetImage();
			loadPreviousImage();
			showCurrentImage();

			cout << "Current index: " + to_string(index) << endl;
			while (index >= indexBegin - 1 && index != middle) {
				if (!maskExists(index)) {
					goto nextIter;
				}
			come_back1:
				showCurrentImage();
				loadPreviousImage();
				cout << "Current index: " + to_string(index) << endl;
			}
			comeBack1 = false;

			showCurrentImage();

			comeBack2 = true;
			loadNextImage();
			showCurrentImage();

			cout << "Current index: " + to_string(index) << endl;
			while (index <= indexEnd - 1 && index != middle) {
				if (!maskExists(index)) {
					goto nextIter;
				}
			come_back2:
				showCurrentImage();
				loadNextImage();
				cout << "Current index: " + to_string(index) << endl;
			}
			comeBack2 = false;

			resetImage();
			break;
		}

					// Show first image (middle of the dataset)
		case 'f': {
			resetImage();
			break;
		}

					// Go to lowest image with a mask
		case 'h': {
			resetImage();
			go = false;

			cout << "Index: " << std::to_string(index) << " - IndexBegin: " << std::to_string(indexBegin) << endl;
			for (; index >= 0; index--, currentImage--) {
				if (!maskExists(index)) {
					cout << "Breaks at image: " + to_string(currentImage) << endl;
					go = true;
					break;
				}
			}

			if (!go) {
				cout << "All images on the left side have masks!" << endl;
			}

			loadNextImage();
			showCurrentImage();
			go = false;

			break;
		}

					// Go to the next lower image
		case 'j': {
			loadPreviousImage();
			showCurrentImage();

			if (go) {
				go = false;
				goto nextIter;
			}

			break;
		}

					// Go to the next higher image
		case 'k': {
			loadNextImage();
			showCurrentImage();

			if (go) {
				go = false;
				goto nextIter;
			}

			break;
		}

					// Go to highest image with a mask
		case 'l': {
			resetImage();
			go = false;

			for (; index < folderSize; index++, currentImage++) {
				if (!maskExists(index)) {
					cout << "Breaks at image: " + to_string(currentImage) << endl;
					go = true;
					break;
				}
			}

			if (!go) {
				cout << "All images on the right side have masks!" << endl;
			}

			loadPreviousImage();
			showCurrentImage();
			go = false;

			break;
		}

			  nextIter:
				  // Run an iteration of GrabCut in the current image
		case 'n': {
			int iterCount = gcapp.getIterCount();
			cout << "<" << iterCount << "... ";
			int newIterCount = gcapp.nextIter();
			if (newIterCount > iterCount) {
				gcapp.showImage();
				cout << iterCount << ">" << endl;
				maskFiles[index] = saveCurrentMask();
			}
			else
				cout << "rect must be determined>" << endl;

			if (comeBack1) {
				goto come_back1;
			}
			if (comeBack2) {
				goto come_back2;
			}
			break;
		}

				  // Reset current image GrabCut mask
		case 'r': {
			cout << "Mask was reset for image: " + to_string(currentImage) << endl
				<< "The saved mask will not be changed until a new iteration is ran." << endl;
			gcapp.reset();
			gcapp.showImage();
			break;
		}

				  // Generate vox file
		case 'v': {
			generateVox();
		}
		}
	}
exit_main:
	destroyWindow("image");
	return 0;
}